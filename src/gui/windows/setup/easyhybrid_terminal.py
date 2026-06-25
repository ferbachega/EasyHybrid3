#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  easyhybrid_terminal.py
#
#  Copyright 2022-2025 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
# ============================================================================
#  QUICK GUIDE TO THE TERMINAL  (DSL syntax:  command arg=value arg2=value2)
# ----------------------------------------------------------------------------
#  Discovery
#     help                         list all commands
#     list                         list loaded objects (with indices)
#
#  Representations  (rep = lines|sticks|spheres|dash|...)
#     show rep=sticks              applies to the ACTIVE SELECTION
#     show rep=sticks obj=1        targets only object 1 (selection unchanged)
#     show rep=spheres obj=0 chain=A
#     show rep=sticks  obj=0 resn=HIS name=CA
#     hide rep=lines   obj=0 resi=10-25
#
#  Active selection (persistent)  -- afterwards, commands without obj= act on it
#     select obj=0                 whole object
#     select obj=0 chain=A,B       chains A and B
#     select obj=0 resi=10-30      residue range
#     select obj=0 chain=A resn=HIS name=CA
#     deselect                     clear the selection
#
#  Filters (combinable, AND logic), accepted by show/hide/select:
#     chain = A   | A,B                 (chain name)
#     resi  = 45  | 10-20 | 1-5,8,12    (residue index; ranges and lists)
#     resn  = HIS                       (residue name)
#     name  = CA                        (atom name)
#
#  Trajectory (e.g. 2000 frames)
#     frame                        show the current frame
#     frame n=1000                 jump to frame 1000
#     next  /  prev                step forward / backward one frame
#
#  Camera
#     center                       center on the active selection
#     center obj=0                 center of mass of object 0
#     center obj=0 chain=A resn=HIS
#     zoom dir=in                  zoom in (in) or out (out)
#     zoom dir=out steps=10        several steps at once
#
#  Scene / file
#     axes show=true | false
#     load file=/path/system.pdb
#     load file="/with space/sys.xyz"   (use quotes for paths with spaces)
#
#  Keys: Up/Down arrows browse history; Tab completes the command;
#        double-Tab lists the options matching the current prefix.
#
#  API in scripts (same logic, without the DSL string):
#     cmd.show(rep="sticks", obj=1, chain="A")
#     cmd.select(obj=0, resi="10-30")
#     for i in range(0, 2000, 100): cmd.frame(n=i)
# ============================================================================
import gi
import sys
import io
import time
import shlex
import inspect
import traceback
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, GLib

import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


# ============================================================================
#  Command — unified command API for EasyHybrid
# ----------------------------------------------------------------------------
#  This class is the SINGLE SOURCE of commands. It replaces the three formerly
#  fragmented pieces (the terminal's stub Command, the CommandLine in eSession
#  and the raw eval/exec). The execution model is the requested DSL:
#
#        command arg1=value1 arg2=value2
#
#  Dual access path (hence "both" for the API):
#    - In the terminal (text):  cmd.run("show rep=sticks")
#    - In Python scripts:        cmd.show(rep="sticks")
#  Both routes land in the SAME cmd_* method, so there is no duplicated logic.
#
#  Each public command is a named method cmd_<name>. The cmd_ prefix serves
#  three purposes: (1) the parser knows what is dispatchable, (2) autocomplete
#  enumerates commands by introspection (fixing the old empty command_list),
#  (3) internal helper methods don't accidentally become commands.
# ============================================================================
class Command:
    def __init__(self, console, vm_session=None):
        self.console = console
        self.vm_session = vm_session

    # ----------------------------------------------------------------- parser
    def _parse(self, cmd_text):
        """
        Converts 'command arg=val arg2=val2' into (name, kwargs).

        Improvements over the old parser:
          - uses shlex to honor quotes ("name with space");
          - does NOT use str.replace(func, '') (which corrupted args that
            contained the command name);
          - tolerates empty input, extra spaces and '=' joined or separated;
          - coerces types: true/false -> bool, ints, floats, else string.
        """
        cmd_text = (cmd_text or "").strip()
        if not cmd_text:
            return None, {}
        try:
            tokens = shlex.split(cmd_text)
        except ValueError:
            # unbalanced quotes etc.: fall back to a simple split
            tokens = cmd_text.split()
        name = tokens[0]
        kwargs = {}
        for tok in tokens[1:]:
            if "=" in tok:
                k, v = tok.split("=", 1)
                kwargs[k.strip()] = self._coerce(v.strip())
            else:
                # loose positional argument -> accumulate in _args
                kwargs.setdefault("_args", []).append(self._coerce(tok))
        return name, kwargs

    @staticmethod
    def _coerce(value):
        """ Converts a DSL string to its most likely Python type. """
        low = value.lower()
        if low in ("true", "yes", "on"):
            return True
        if low in ("false", "no", "off"):
            return False
        if low in ("none", "null"):
            return None
        for cast in (int, float):
            try:
                return cast(value)
            except (ValueError, TypeError):
                pass
        return value

    # ------------------------------------------------------------- dispatcher
    def run(self, cmd_text):
        """
        TERMINAL entry point (DSL text). Parses, finds the cmd_<name> method
        and runs it with the kwargs. Returns a log string (or None) that the
        terminal prints. Errors become a readable log, never crashing the
        window.
        """
        name, kwargs = self._parse(cmd_text)
        if name is None:
            return None
        method = getattr(self, "cmd_" + name, None)
        if method is None or not callable(method):
            return "Command '{}' not recognized. Use 'help'.".format(name)
        try:
            return method(**kwargs)
        except TypeError as te:
            # typically a wrong argument: show the signature to help
            sig = inspect.signature(method)
            return "Uso: {} {}\n  ({})".format(name, self._sig_hint(sig), te)
        except Exception as exc:
            return "Erro em '{}': {}".format(name, exc)

    @staticmethod
    def _sig_hint(sig):
        parts = []
        for pname, p in sig.parameters.items():
            if pname == "_args":
                continue
            if p.default is inspect.Parameter.empty:
                parts.append(pname)
            else:
                parts.append("{}={!r}".format(pname, p.default))
        return " ".join(parts)

    # ---------------------------------------------------- introspection for TAB
    def command_names(self):
        """ Lists available command names (without the cmd_ prefix).
            Used by autocomplete and by 'help'. """
        return sorted(
            attr[len("cmd_"):]
            for attr in dir(self)
            if attr.startswith("cmd_") and callable(getattr(self, attr))
        )

    # ---------------------------------------------------- target helpers
    def _resolve_object(self, obj):
        """ Takes an index (int) and returns the matching vm_object, or None
            if it does not exist. Also accepts None (no target). """
        if obj is None or self.vm_session is None:
            return None
        try:
            return self.vm_session.vm_objects_dic.get(int(obj))
        except (ValueError, TypeError):
            return None

    def _selection_for_atoms(self, vm_object, atoms):
        """ Builds a temporary selection (VMSele) containing `vm_object` and
            the given `atoms` set, WITHOUT touching the user's active
            selection.

            Used by show/hide when they receive inline filters (obj=, chain=,
            ...), so the action is one-off and does not change selection state.
        """
        active = self.vm_session.selections[self.vm_session.current_selection]
        temp = active.__class__(self.vm_session)          # same VMSele type
        temp.selected_objects = {vm_object}
        temp.selected_atoms = set(atoms)
        temp.selected_atom_ids = set(
            getattr(a, "atom_id", getattr(a, "index", None)) for a in atoms)
        return temp

    @staticmethod
    def _match_field(value, spec):
        """ True if `value` matches `spec`, where spec may be:
              - single value:      "A"        -> value == "A"
              - comma list:         "A,B,C"    -> value in {A,B,C}
              - numeric range:      "10-20"    -> 10 <= int(value) <= 20
              - list of ranges:     "1-5,8,12" -> combines the above
            Text comparison is case-sensitive (chain/residue names usually
            are); ranges require value to be an integer.

            Examples:
              _match_field("A",  "A,B")   -> True
              _match_field("15", "10-20") -> True
              _match_field("CA", "CA")    -> True
        """
        value = str(value)
        for part in str(spec).split(","):
            part = part.strip()
            if "-" in part and part.replace("-", "").isdigit():
                lo, hi = part.split("-", 1)
                try:
                    if int(lo) <= int(value) <= int(hi):
                        return True
                except ValueError:
                    pass
            elif part == value:
                return True
        return False

    def _filter_atoms(self, vm_object, chain=None, resi=None, resn=None,
                      name=None):
        """ Returns the list of atoms of `vm_object` that satisfy ALL given
            filters (AND logic). A None filter is ignored.

            The compared attributes exist in VisMol's atom model:
              chain -> atom.chain.name     resn -> atom.residue.name
              resi  -> atom.residue.index  name -> atom.name

            NOTE: the session's selecting_by_* methods do NOT filter by name --
            they expand from a clicked atom. That is why we filter here
            directly.

            Accepted spec examples (via _match_field):
              chain="A"        chain="A,B"
              resi="45"        resi="10-20"      resi="1-5,8,12"
              resn="HIS"       name="CA"
        """
        def keep(atom):
            if chain is not None:
                cn = getattr(getattr(atom, "chain", None), "name", None)
                if cn is None or not self._match_field(cn, chain):
                    return False
            if resn is not None:
                rn = getattr(getattr(atom, "residue", None), "name", None)
                if rn is None or not self._match_field(rn, resn):
                    return False
            if resi is not None:
                ri = getattr(getattr(atom, "residue", None), "index", None)
                if ri is None or not self._match_field(ri, resi):
                    return False
            if name is not None:
                an = getattr(atom, "name", None)
                if an is None or not self._match_field(an, name):
                    return False
            return True
        return [a for a in vm_object.atoms.values() if keep(a)]

    # ========================================================================
    #  REAL COMMANDS  (each is also directly callable via the Python API)
    #  Anchored on vm_session's public API, not on stubs.
    # ========================================================================
    def cmd_help(self, **_):
        """ Lists available commands and the syntax of each one. """
        lines = ["Available commands:"]
        for name in self.command_names():
            method = getattr(self, "cmd_" + name)
            doc = (method.__doc__ or "").strip().split("\n")[0]
            lines.append("  {:<12} {}".format(name, doc))
        return "\n".join(lines)

    def cmd_list(self, **_):
        """ Lists the molecular objects loaded in the session. """
        if self.vm_session is None:
            return "Session unavailable."
        objs = self.vm_session.vm_objects_dic
        if not objs:
            return "No object loaded."
        return "\n".join("  [{}] {}".format(i, getattr(o, "name", "?"))
                         for i, o in objs.items())

    def cmd_show(self, rep="lines", obj=None, chain=None, resi=None,
                 resn=None, name=None, **_):
        """ Shows a representation.

            With no filters, acts on the ACTIVE SELECTION. With obj= (and
            optionally chain/resi/resn/name) it targets only those atoms,
            WITHOUT changing the active selection -- the action is one-off.

            rep  : lines | sticks | spheres | dash | ...
            obj  : object index (see 'list')
            chain: chain name        e.g. A   or A,B
            resi : residue index     e.g. 45  or 10-20  or 1-5,8
            resn : residue name      e.g. HIS
            name : atom name         e.g. CA

            Examples:
              show rep=sticks
              show rep=sticks obj=1
              show rep=spheres obj=0 chain=A
              show rep=sticks  obj=0 resn=HIS name=CA
              show rep=lines   obj=0 resi=10-25
        """
        return self._show_or_hide(True, rep, obj, chain, resi, resn, name)

    def cmd_hide(self, rep="lines", obj=None, chain=None, resi=None,
                 resn=None, name=None, **_):
        """ Hides a representation. Same filters as 'show'.

            Examples:
              hide rep=lines
              hide rep=lines   obj=0
              hide rep=spheres obj=0 chain=B
              hide rep=sticks  obj=0 resn=HOH
        """
        return self._show_or_hide(False, rep, obj, chain, resi, resn, name)

    def _show_or_hide(self, show, rep, obj, chain, resi, resn, name):
        """ Shared implementation of show/hide. Decides the target:
              - no filter              -> active selection (selection=None);
              - obj= [+ fine filters]  -> temporary selection with target only.
        """
        if self.vm_session is None:
            return "Session unavailable."
        verb = "Showing" if show else "Hiding"

        # Case 1: no explicit target -> use the user's active selection.
        if obj is None:
            self.vm_session.show_or_hide(rep_type=rep, selection=None, show=show)
            return "{} '{}' on the active selection.".format(verb, rep)

        # Case 2: explicit target by object (and maybe fine filters).
        target = self._resolve_object(obj)
        if target is None:
            return "Object {} not found. Use 'list'.".format(obj)

        has_fine = any(f is not None for f in (chain, resi, resn, name))
        if has_fine:
            atoms = self._filter_atoms(target, chain, resi, resn, name)
            if not atoms:
                return "No atom matches the filters in obj={}.".format(obj)
            selection = self._selection_for_atoms(target, atoms)
            scope = "obj={} ({} filtered atoms)".format(obj, len(atoms))
        else:
            selection = self._selection_for_atoms(target, target.atoms.values())
            scope = "obj={} (whole)".format(obj)

        self.vm_session.show_or_hide(rep_type=rep, selection=selection, show=show)
        return "{} '{}' on {}.".format(verb, rep, scope)

    def cmd_select(self, obj=None, chain=None, resi=None, resn=None,
                   name=None, **_):
        """ Sets the ACTIVE SELECTION (persistent). Subsequent commands
            without obj= act on it until a new 'select' or 'deselect'.

            Combinable filters (AND), same as show/hide:
              obj  : object index (required)
              chain: A   or A,B
              resi : 45  or 10-20  or 1-5,8
              resn : HIS
              name : CA

            Examples:
              select obj=0
              select obj=0 chain=A
              select obj=0 chain=A resn=HIS
              select obj=0 resi=10-30
              select obj=0 chain=A resn=HIS name=CA
        """
        if self.vm_session is None:
            return "Session unavailable."
        target = self._resolve_object(obj)
        if target is None:
            return "Usage: select obj=N [chain= resi= resn= name=]  ('list')"

        atoms = self._filter_atoms(target, chain, resi, resn, name)
        if not atoms:
            return "No atom matches the filters in obj={}.".format(obj)

        # Write into the active (persistent) selection.
        active = self.vm_session.selections[self.vm_session.current_selection]
        active.selected_objects = {target}
        active.selected_atoms = set(atoms)
        active.selected_atom_ids = set(
            getattr(a, "atom_id", getattr(a, "index", None)) for a in atoms)
        # Sync each object atom's .selected flag.
        for a in target.atoms.values():
            a.selected = False
        for a in atoms:
            a.selected = True

        filters = ", ".join(
            "{}={}".format(k, v) for k, v in
            (("chain", chain), ("resi", resi), ("resn", resn), ("name", name))
            if v is not None) or "whole object"
        return "Selected obj={} ({}): {} atoms.".format(obj, filters, len(atoms))

    def cmd_deselect(self, **_):
        """ Clears the active selection. """
        if self.vm_session is None:
            return "Session unavailable."
        active = self.vm_session.selections[self.vm_session.current_selection]
        active.selected_objects = set()
        active.selected_atoms = set()
        active.selected_atom_ids = set()
        return "Selection cleared."

    def cmd_frame(self, n=None, **_):
        """ Without n: shows the current frame. With n=int: jumps to it. """
        if self.vm_session is None:
            return "Session unavailable."
        if n is None:
            return "Current frame: {}".format(self.vm_session.get_frame())
        self.vm_session.set_frame(frame=int(n))
        return "Frame -> {}".format(int(n))

    def cmd_next(self, **_):
        """ Steps one frame forward in the trajectory. """
        if self.vm_session is None:
            return "Session unavailable."
        self.vm_session.forward_frame()
        return "Frame -> {}".format(self.vm_session.get_frame())

    def cmd_prev(self, **_):
        """ Steps one frame backward in the trajectory. """
        if self.vm_session is None:
            return "Session unavailable."
        self.vm_session.reverse_frame()
        return "Frame -> {}".format(self.vm_session.get_frame())

    def cmd_axes(self, show=True, **_):
        """ Shows/hides the axes (show=true|false). """
        if self.vm_session is None:
            return "Session unavailable."
        if show:
            self.vm_session.show_axes()
            return "Axes visible."
        self.vm_session.hide_axes()
        return "Axes hidden."

    def cmd_center(self, obj=None, chain=None, resi=None, resn=None,
                   name=None, **_):
        """ Centers the camera on a target (animates the translation to it).

            Without obj=: centers on the ACTIVE SELECTION (atom centroid).
            With obj= and no fine filters: object center of mass.
            With obj= and filters: centroid of the filtered atoms.

            obj=N  chain=A  resi=45|10-20  resn=HIS  name=CA

            Examples:
              center                       (on the active selection)
              center obj=0                 (center of mass of object 0)
              center obj=0 chain=A         (centroid of chain A)
              center obj=0 resn=HIS name=CA
        """
        if self.vm_session is None:
            return "Session unavailable."
        glcore = self.vm_session.vm_glcore

        # Case A: no object -> center on the active selection.
        if obj is None:
            active = self.vm_session.selections[self.vm_session.current_selection]
            atoms = list(active.selected_atoms)
            if not atoms:
                return "Nothing selected. Use 'select' or pass obj=N."
            vm_object = next(iter(active.selected_objects), None)
            target = self._centroid(atoms, vm_object)
            glcore.center_on_coordinates(vm_object, target)
            return "Centered on the active selection ({} atoms).".format(len(atoms))

        # Case B: explicit object.
        vm_object = self._resolve_object(obj)
        if vm_object is None:
            return "Object {} not found. Use 'list'.".format(obj)

        has_fine = any(f is not None for f in (chain, resi, resn, name))
        if has_fine:
            atoms = self._filter_atoms(vm_object, chain, resi, resn, name)
            if not atoms:
                return "No atom matches the filters in obj={}.".format(obj)
            target = self._centroid(atoms, vm_object)
            glcore.center_on_coordinates(vm_object, target)
            return "Centered on obj={} ({} atoms).".format(obj, len(atoms))
        else:
            # whole object -> center of mass (same pattern as autocenter)
            glcore.center_on_coordinates(vm_object, vm_object.mass_center)
            return "Centered on the center of mass of obj={}.".format(obj)

    def _centroid(self, atoms, vm_object):
        """ Centroid (mean of coordinates) of a list of atoms at the current
            frame. Uses atom.coords(frame) when available; falls back to
            atom.coords if it is an attribute. Returns np.array float32 [x,y,z]. """
        import numpy as np
        frame = 0
        try:
            frame = self.vm_session.get_frame()
        except Exception:
            pass
        pts = []
        for a in atoms:
            c = getattr(a, "coords", None)
            if callable(c):
                pts.append(np.asarray(c(frame), dtype=np.float32))
            elif c is not None:
                pts.append(np.asarray(c, dtype=np.float32))
        if not pts:
            # fallback: object center of mass
            return vm_object.mass_center
        return np.mean(np.vstack(pts), axis=0).astype(np.float32)

    def cmd_zoom(self, dir="in", steps=5, **_):
        """ Zooms the camera in or out (equivalent to turning the mouse wheel).

            dir  : in | out          (in zooms in, out zooms out)
            steps: how many scroll 'clicks' to apply (default 5)

            Examples:
              zoom dir=in
              zoom dir=out steps=10
        """
        if self.vm_session is None:
            return "Session unavailable."
        glcore = self.vm_session.vm_glcore
        direction = 1 if str(dir).lower() in ("in", "+", "1", "up") else -1
        try:
            n = max(1, int(steps))
        except (ValueError, TypeError):
            n = 5
        for _i in range(n):
            glcore.mouse_scroll(direction)
        return "Zoom {} x{}.".format("in" if direction == 1 else "out", n)

    def cmd_load(self, file=None, **_):
        """ Loads a molecule from file (file=/path/file). """
        if self.vm_session is None:
            return "Session unavailable."
        if file is None:
            return "Usage: load file=/path/to/file"
        self.vm_session.load_molecule(file)
        return "Loaded: {}".format(file)

    # --- atalhos de API direta (mesma logica, nomes "pythonicos") -----------
    #  Permitem cmd.show(rep="sticks") em scripts externos sem o prefixo cmd_.
    # --- direct API shortcuts (same logic, "pythonic" names) ----------------
    #  Allow use in external scripts without the cmd_ prefix or the DSL string:
    #      cmd.show(rep="sticks", obj=1, chain="A")
    #      cmd.select(obj=0, resi="10-30")
    def show(self, rep="lines", obj=None, chain=None, resi=None, resn=None, name=None):
        return self.cmd_show(rep=rep, obj=obj, chain=chain, resi=resi, resn=resn, name=name)
    def hide(self, rep="lines", obj=None, chain=None, resi=None, resn=None, name=None):
        return self.cmd_hide(rep=rep, obj=obj, chain=chain, resi=resi, resn=resn, name=name)
    def select(self, obj=None, chain=None, resi=None, resn=None, name=None):
        return self.cmd_select(obj=obj, chain=chain, resi=resi, resn=resn, name=name)
    def center(self, obj=None, chain=None, resi=None, resn=None, name=None):
        return self.cmd_center(obj=obj, chain=chain, resi=resi, resn=resn, name=name)
    def zoom(self, dir="in", steps=5):      return self.cmd_zoom(dir=dir, steps=steps)
    def frame(self, n=None):                return self.cmd_frame(n=n)
    def list(self):                         return self.cmd_list()


class TerminalWindow():
    """ EasyHybrid terminal window. """

    def open_window(self):
        """ """
        if self.visible == False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.main.home, 'src/gui/windows/setup/easyhybrid_terminal.glade'))
            self.builder.connect_signals(self)

            self.window = self.builder.get_object('window')
            self.window.set_default_size(650, 350)
            self.window.set_title('EasyHybrid Terminal')
            self.window.connect('destroy-event', self.close_window)
            self.window.set_keep_above(True)

            self.tag_table = self.textbuffer.get_tag_table()
            self.textbuffer.create_tag("green", foreground="green")
            self.textbuffer.create_tag("red", foreground="red")
            self.textbuffer.create_tag("blue", foreground="blue")

            self.locals = {}
            # Execution context: the unified Command
            self.cmd = Command(self, self.vm_session)
            self.locals['cmd'] = self.cmd
            # command_list is now actually populated (it used to be empty)
            self.command_list = self.cmd.command_names()

            self.entry_terminal = self.builder.get_object('entry_terminal')
            self.textview = self.builder.get_object('entry_text_buffer')
            self.textview.set_buffer(self.textbuffer)
            self.textview.set_wrap_mode(Gtk.WrapMode.WORD)
            self.textview.get_style_context().add_provider(self.css_provider, Gtk.STYLE_PROVIDER_PRIORITY_USER)
            self.entry_terminal.get_style_context().add_provider(self.css_provider, Gtk.STYLE_PROVIDER_PRIORITY_USER)

            self.window.show_all()
            self.visible = True
        else:
            self.window.present()

    def close_window(self, button, data=None):
        """ Function doc """
        self.window.destroy()
        self.visible = False
        self.main.cmd_terminal_button.set_active(False)

    def __init__(self, main=None):
        """ Class initialiser """
        self.main = main
        self.visible = False
        self.p_session = main.p_session
        self.vm_session = main.vm_session

        self.cmd_history = []
        self.cmd_history_counter = 0
        self.textbuffer = self.main.terminal_text_buffer
        self.command_list = []

        text = '''
        --------------------------------------------------------
               Welcome to the EasyHybrid Terminal.
        --------------------------------------------------------

                   Created by J.F.R Bachega


        '''
        end_iter = self.textbuffer.get_end_iter()
        self.textbuffer.insert(end_iter, text)

        self.last_tab_time = 0
        self.tab_timer_id = None

        self.css_provider = Gtk.CssProvider()
        self.css_provider.load_from_data(b"""
        textview {
            font-family: Monospace;
            font-size: 12pt;
        }

        entry {
            font-family: Monospace;
            font-size: 12pt;
        }
        """)

    def run_cmd(self, cmd):
        """ Runs a DSL command via Command.run and prints the log. """
        self.write_output(">" + cmd, "normal")
        log = self.cmd.run(cmd)
        if log is not None:
            self.write_output(log, "normal")

    def write_output(self, text: str, color: str = "normal"):
        end_iter = self.textbuffer.get_end_iter()
        tag = self.tag_table.lookup(color)
        if not tag:
            tag = self.textbuffer.create_tag(color, foreground=color)
        self.textbuffer.insert_with_tags(end_iter, text + "\n", tag)
        self.textview.scroll_to_iter(end_iter, 0.0, False, 0, 0)

    def on_entry_terminal(self, widget):
        """
        Enter in the input field. Now routes through the DSL (Command.run)
        instead of raw eval/exec -- safer and consistent with autocomplete.
        """
        command = self.entry_terminal.get_text()
        self.entry_terminal.set_text("")
        if not command.strip():
            return
        self.cmd_history.append(command)
        self.cmd_history_counter = 0
        self.write_output(">>> " + command, "normal")
        try:
            log = self.cmd.run(command)
            if log is not None:
                self.write_output(str(log), "normal")
        except Exception:
            self.write_output(traceback.format_exc().strip(), "red")

    # ------------------------------------------------------------ navigation
    def on_entry_terminal_backspace(self, widget):
        pass

    def on_entry_terminal_move_cursor(self, widget, data, data2, data3):
        pass

    def on_entry_terminal_change(self, widget):
        pass

    def update_window(self, system_names=True, coordinates=False, selections=True):
        pass

    def on_key_press_event(self, widget, event):
        """ Arrows: browse history. Tab: autocomplete. Double-Tab: list. """
        k_name = Gdk.keyval_name(event.keyval)
        size = -len(self.cmd_history)

        if k_name in ['Down', 'Up']:
            if k_name == 'Up':
                self.cmd_history_counter += -1
            if k_name == 'Down':
                self.cmd_history_counter += 1

            if self.cmd_history_counter >= 0:
                self.entry_terminal.set_text('')
                self.cmd_history_counter = 0
            elif self.cmd_history_counter < 0 and self.cmd_history_counter > size:
                self.entry_terminal.set_text(self.cmd_history[self.cmd_history_counter])
            else:
                self.cmd_history_counter = size
                self.entry_terminal.set_text(self.cmd_history[self.cmd_history_counter])
            # move the cursor to the end of the recovered text
            self.entry_terminal.set_position(-1)
            return True

        if event.keyval == Gdk.KEY_Tab:
            now = time.time()
            if now - self.last_tab_time < 0.2:
                # double-Tab: cancel the timer and LIST all options
                if self.tab_timer_id is not None:
                    GLib.source_remove(self.tab_timer_id)
                    self.tab_timer_id = None
                self.last_tab_time = 0
                self._list_completions()
            else:
                # single Tab (scheduled): try to complete the current prefix
                self.last_tab_time = now
                self.tab_timer_id = GLib.timeout_add(200, self._simple_Tab)
            # returning True prevents GTK from moving focus out of the entry
            return True
        return False

    # ----------------------------------------------------------- autocomplete
    def _current_word(self):
        """ The first typed word (the command name being built). """
        return self.entry_terminal.get_text().strip().split(" ")[0]

    def _matches(self, prefix):
        return [c for c in self.command_list if c.startswith(prefix)]

    def _simple_Tab(self):
        """
        Single Tab: completes the command name.
          - 1 match   -> complete and add a space;
          - several   -> complete up to the longest common prefix;
          - none       -> nothing.
        """
        self.tab_timer_id = None
        full = self.entry_terminal.get_text()
        # only complete while the user is still typing the NAME (no space)
        if " " in full.strip():
            return False
        prefix = full.strip()
        matches = self._matches(prefix)
        if not matches:
            return False
        if len(matches) == 1:
            self.entry_terminal.set_text(matches[0] + " ")
        else:
            common = os.path.commonprefix(matches)
            if len(common) > len(prefix):
                self.entry_terminal.set_text(common)
            else:
                self._list_completions()
        self.entry_terminal.set_position(-1)
        return False  # remove o timer

    def _list_completions(self):
        """ Shows in the buffer the options matching the current prefix. """
        prefix = self._current_word()
        matches = self._matches(prefix) if prefix else list(self.command_list)
        if matches:
            self.write_output("  ".join(matches), "blue")

    def _simple_Tab_legacy(self):
        # kept only as a reference; not used
        self.tab_timer_id = None
        return False
