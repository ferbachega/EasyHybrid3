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
#  Builder -- bonds (pick 2 atoms as pk1/pk2 with the measurement
#  picking tool -- same pk1..pk4 used for distances/angles/dihedrals --
#  or address atoms directly by obj=/atom1=/atom2=, no picking needed)
#     bond                         bonds pk1-pk2 (single bond)
#     bond order=2                 bonds pk1-pk2 as a double bond
#     bond order=3                 bonds pk1-pk2 as a triple bond
#     bond obj=0 atom1=3 atom2=7 order=2   bonds by atom_id directly
#     unbond                       removes the bond between pk1-pk2
#     unbond obj=0 atom1=3 atom2=7         removes by atom_id directly
#
#  Builder -- Dynamic Bonds (same 'bond'/'unbond', add frame=...):
#     bond order=2 frame=true      edits ONLY the current frame
#     bond order=2 frame=12        edits ONLY frame 12
#     bond order=2 frame=1:20      edits frames 1 through 20 (inclusive)
#     bond order=1 frame=all       edits every frame
#     unbond frame=true/12/1:20/all   same frame= forms, for removing
#     *** frame=... is REPRESENTATION-ONLY: it changes what is drawn for
#     the selected frame(s), NOT a real chemical bond -- it does NOT
#     touch the pDynamo system's topology/force field/force constants.
#     See 'help bond' for the full explanation.
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
#     cmd.bond(order=2)                       # bonds pk1-pk2, double bond
#     cmd.bond(obj=0, atom1=3, atom2=7, order=2)
#     cmd.bond(order=2, frame="1:20")         # Dynamic Bonds, frames 1-20
#     cmd.unbond()                            # unbonds pk1-pk2
#     for i in range(0, 2000, 100): cmd.frame(n=i)
# ============================================================================
import gi
import sys
import io
import time
import shlex
import inspect
import traceback
import numpy as np
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
# ============================================================================
#  Named colours for the 'color' command below. Atom colours in this
#  codebase are stored as RGB floats in the 0-1 range (see Atom._init_color()
#  in vismol/model/atom.py: "the returned value is in scale of 0 to 1") --
#  matching that convention here, rather than the more common 0-255 range,
#  to stay consistent with everything colours already do in this project.
#  This is a small, deliberately short list of common names (not a full
#  CSS/X11 palette) -- explicit r=/g=/b= is always available in 'color' for
#  anything not in here.
# ============================================================================
NAMED_COLORS = {
    "red":     (1.0, 0.0, 0.0),
    "green":   (0.0, 1.0, 0.0),
    "blue":    (0.0, 0.0, 1.0),
    "yellow":  (1.0, 1.0, 0.0),
    "cyan":    (0.0, 1.0, 1.0),
    "magenta": (1.0, 0.0, 1.0),
    "white":   (1.0, 1.0, 1.0),
    "black":   (0.0, 0.0, 0.0),
    "gray":    (0.5, 0.5, 0.5),
    "grey":    (0.5, 0.5, 0.5),
    "orange":  (1.0, 0.5, 0.0),
    "purple":  (0.5, 0.0, 0.5),
    "pink":    (1.0, 0.6, 0.7),
    "brown":   (0.6, 0.3, 0.1),
}


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
            return "Usage: {} {}\n  ({})".format(name, self._sig_hint(sig), te)
        except Exception as exc:
            return "Error in '{}': {}".format(name, exc)

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
                      name=None, symbol=None,
                      not_chain=None, not_resi=None, not_resn=None,
                      not_name=None, not_symbol=None):
        """ Returns the list of atoms of `vm_object` that satisfy ALL given
            INCLUDE filters (AND logic) and NONE of the given EXCLUDE
            filters (each not_* rules an atom out if it matches -- e.g.
            not_symbol=H removes every hydrogen, regardless of what the
            include filters selected). A None filter (include or exclude)
            is ignored.

            The compared attributes exist in VisMol's atom model:
              chain  -> atom.chain.name     resn -> atom.residue.name
              resi   -> atom.residue.index  name -> atom.name
              symbol -> atom.symbol (the chemical element, e.g. "H", "C",
                        "O" -- NOT the same as `name`, which can be things
                        like "HA1"/"CB"/"OXT" even for the same element)

            NOTE: the session's selecting_by_* methods do NOT filter by name --
            they expand from a clicked atom. That is why we filter here
            directly.

            Accepted spec examples (via _match_field), for BOTH the
            include and the not_* exclude versions of each field:
              chain="A"        chain="A,B"
              resi="45"        resi="10-20"      resi="1-5,8,12"
              resn="HIS"       name="CA"         symbol="H"
              symbol="H,D"     (hydrogen AND deuterium)

            Examples (see cmd_select's own docstring for the full
            command-line syntax):
              everything except hydrogens:        not_symbol="H"
              everything except waters:            not_resn="HOH"
              chain A except its hydrogens:        chain="A", not_symbol="H"
        """
        def keep(atom):
            # --- filtros de INCLUSAO (E) -- todos os dados precisam bater ---
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
            if symbol is not None:
                sy = getattr(atom, "symbol", None)
                if sy is None or not self._match_field(sy, symbol):
                    return False
            # --- filtros de EXCLUSAO -- bater em QUALQUER um ja tira o atomo ---
            if not_chain is not None:
                cn = getattr(getattr(atom, "chain", None), "name", None)
                if cn is not None and self._match_field(cn, not_chain):
                    return False
            if not_resn is not None:
                rn = getattr(getattr(atom, "residue", None), "name", None)
                if rn is not None and self._match_field(rn, not_resn):
                    return False
            if not_resi is not None:
                ri = getattr(getattr(atom, "residue", None), "index", None)
                if ri is not None and self._match_field(ri, not_resi):
                    return False
            if not_name is not None:
                an = getattr(atom, "name", None)
                if an is not None and self._match_field(an, not_name):
                    return False
            if not_symbol is not None:
                sy = getattr(atom, "symbol", None)
                if sy is not None and self._match_field(sy, not_symbol):
                    return False
            return True
        return [a for a in vm_object.atoms.values() if keep(a)]

    # ========================================================================
    #  REAL COMMANDS  (each is also directly callable via the Python API)
    #  Anchored on vm_session's public API, not on stubs.
    # ========================================================================
    def cmd_help(self, cmd=None, **_):
        """ Lists available commands. With cmd=<nome>, shows the full
        description and usage examples for that specific command.

        Examples:
          help
          help cmd=show
          help cmd=bond
        """
        if cmd is not None:
            method = getattr(self, "cmd_" + cmd, None)
            if method is None:
                return ("Command '{}' does not exist. Available commands: {}"
                        .format(cmd, ", ".join(self.command_names())))
            doc = inspect.getdoc(method) or "(no description)"
            return "{}\n{}\n{}".format(cmd, "-" * len(cmd), doc)

        lines = ["Available commands (use 'help cmd=<name>' for details and examples):"]
        for name in self.command_names():
            method = getattr(self, "cmd_" + name)
            doc = (method.__doc__ or "").strip().split("\n")[0]
            lines.append("  {:<12} {}".format(name, doc))
        return "\n".join(lines)

    def cmd_list(self, **_):
        """ Lists the molecular objects loaded in the session, with the
        index each one has (use that index as obj= in every other
        command that needs one -- show, hide, select, center, add,
        delete, bond, placemode...).

        Examples:
          list
        """
        if self.vm_session is None:
            return "Session unavailable."
        objs = self.vm_session.vm_objects_dic
        if not objs:
            return "No object loaded."
        return "\n".join("  [{}] {}".format(i, getattr(o, "name", "?"))
                         for i, o in objs.items())

    def cmd_show(self, rep="lines", obj=None, chain=None, resi=None,
                 resn=None, name=None, symbol=None,
                 not_chain=None, not_resi=None, not_resn=None,
                 not_name=None, not_symbol=None, **_):
        """ Shows a representation.

            With no filters, acts on the ACTIVE SELECTION. With obj= (and
            optionally chain/resi/resn/name/symbol/not_*) it targets only
            those atoms, WITHOUT changing the active selection -- the
            action is one-off.

            rep    : lines | sticks | spheres | dash | ...
            obj    : object index (see 'list')
            chain  : chain name        e.g. A   or A,B
            resi   : residue index     e.g. 45  or 10-20  or 1-5,8
            resn   : residue name      e.g. HIS
            name   : atom name         e.g. CA
            symbol : chemical element  e.g. H   or C,N,O
            not_*  : exclude filter, same fields (not_chain, not_resi,
                     not_resn, not_name, not_symbol) -- drops any atom
                     that matches, e.g. not_symbol=H to skip hydrogens

            Examples:
              show rep=sticks
              show rep=sticks obj=1
              show rep=spheres obj=0 chain=A
              show rep=sticks  obj=0 resn=HIS name=CA
              show rep=lines   obj=0 resi=10-25
              show rep=sticks  obj=0 not_symbol=H     (everything except hydrogens)
        """
        return self._show_or_hide(True, rep, obj, chain, resi, resn, name, symbol,
                                   not_chain, not_resi, not_resn, not_name, not_symbol)

    def cmd_hide(self, rep="lines", obj=None, chain=None, resi=None,
                 resn=None, name=None, symbol=None,
                 not_chain=None, not_resi=None, not_resn=None,
                 not_name=None, not_symbol=None, **_):
        """ Hides a representation. Same filters as 'show' (including
            symbol= and the not_* exclude filters).

            Examples:
              hide rep=lines
              hide rep=lines   obj=0
              hide rep=spheres obj=0 chain=B
              hide rep=sticks  obj=0 resn=HOH
              hide rep=sticks  obj=0 symbol=H          (hide only hydrogens)
        """
        return self._show_or_hide(False, rep, obj, chain, resi, resn, name, symbol,
                                   not_chain, not_resi, not_resn, not_name, not_symbol)

    def _show_or_hide(self, show, rep, obj, chain, resi, resn, name, symbol=None,
                       not_chain=None, not_resi=None, not_resn=None,
                       not_name=None, not_symbol=None):
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

        has_fine = any(f is not None for f in
                       (chain, resi, resn, name, symbol,
                        not_chain, not_resi, not_resn, not_name, not_symbol))
        if has_fine:
            atoms = self._filter_atoms(target, chain, resi, resn, name, symbol,
                                        not_chain, not_resi, not_resn, not_name, not_symbol)
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
                   name=None, symbol=None,
                   not_chain=None, not_resi=None, not_resn=None,
                   not_name=None, not_symbol=None, **_):
        """ Sets the ACTIVE SELECTION (persistent). Subsequent commands
            without obj= act on it until a new 'select' or 'deselect'.

            Combinable INCLUDE filters (AND), same as show/hide:
              obj    : object index (required)
              chain  : A   or A,B
              resi   : 45  or 10-20  or 1-5,8
              resn   : HIS
              name   : CA
              symbol : H   or C,N,O    (chemical element -- NOT the same
                       as name: an atom named "HA1" still has symbol "H")

            EXCLUDE filters (subtract from whatever the includes above
            matched -- an atom is dropped if it matches ANY not_* given):
              not_chain, not_resi, not_resn, not_name, not_symbol
              (same value syntax as their include counterparts)

            "Select everything except a group of atoms" (e.g. hydrogens)
            is exactly what not_symbol is for -- see the last 3 examples.

            Examples:
              select obj=0
              select obj=0 chain=A
              select obj=0 chain=A resn=HIS
              select obj=0 resi=10-30
              select obj=0 chain=A resn=HIS name=CA
              select obj=0 not_symbol=H                    (everything except hydrogens)
              select obj=0 chain=A not_symbol=H             (chain A, no hydrogens)
              select obj=0 not_resn=HOH not_symbol=H        (no waters, no hydrogens)
        """
        if self.vm_session is None:
            return "Session unavailable."
        target = self._resolve_object(obj)
        if target is None:
            return "Usage: select obj=N [chain= resi= resn= name= symbol= not_chain= not_resi= not_resn= not_name= not_symbol=]  ('list')"

        atoms = self._filter_atoms(target, chain, resi, resn, name, symbol,
                                    not_chain, not_resi, not_resn, not_name, not_symbol)
        if not atoms:
            return "No atom matches the filters in obj={}.".format(obj)

        # Write into the active (persistent) selection.
        active = self.vm_session.selections[self.vm_session.current_selection]
        # [EN] BUG FIX (user reported: selection worked logically -- show/
        # hide/center/bond all correctly acted on it -- but nothing
        # highlighted on screen). Root cause: this used to just set
        # active.selected_atoms / atom.selected directly, which is only
        # HALF of what the app's own click-based selection path does.
        # selection_function_viewing_set() (called by the normal picking
        # flow, see vismol_session.py's _selection_function_set()) always
        # finishes with
        # _build_selected_atoms_coords_and_selected_objects_from_selected_atoms(),
        # which populates vm_object.selected_atom_ids (a PER-OBJECT set,
        # different from active.selected_atom_ids) -- and that is what the
        # on-screen selection highlight actually reads from. Skipping it
        # meant the selection was logically correct (everything keyed off
        # active.selected_atoms already worked) but invisible.
        #
        # SECOND bug found testing this fix (not just reading code):
        # selection_function_viewing_set(None) -- the official "clear"
        # path -- does NOT reset atom.selected on the atoms that were
        # previously selected; it only clears active.selected_atoms and
        # every vm_object.selected_atom_ids. Confirmed live: selecting the
        # 4 hydrogens of a synthetic methane right after selecting just
        # the carbon left the CARBON's own atom.selected still True
        # (only vm_object.selected_atom_ids had correctly dropped it) --
        # a real risk of stale highlighting for whatever was selected
        # before. Fixed by capturing the OLD active.selected_atoms BEFORE
        # clearing, and explicitly setting .selected = False on each of
        # them.
        old_atoms = set(active.selected_atoms)  # copia -- selection_function_viewing_set(None) limpa o set original NO LUGAR
        active.selection_function_viewing_set(None)
        for a in old_atoms:
            a.selected = False
        active.selected_atoms = set(atoms)
        for a in atoms:
            a.selected = True
        active._build_selected_atoms_coords_and_selected_objects_from_selected_atoms()
        active.selected_atom_ids = set(
            getattr(a, "atom_id", getattr(a, "index", None)) for a in atoms)
        self.vm_session.vm_glcore.queue_draw()

        filters = ", ".join(
            "{}={}".format(k, v) for k, v in
            (("chain", chain), ("resi", resi), ("resn", resn), ("name", name),
             ("symbol", symbol), ("not_chain", not_chain), ("not_resi", not_resi),
             ("not_resn", not_resn), ("not_name", not_name), ("not_symbol", not_symbol))
            if v is not None) or "whole object"
        return "Selected obj={} ({}): {} atoms.".format(obj, filters, len(atoms))

    def cmd_deselect(self, **_):
        """ Clears the active (persistent) selection set by 'select'.
        Commands that act on the active selection when no obj= is given
        (show, hide, center) go back to having nothing to act on until
        a new 'select'.

        Examples:
          deselect
        """
        if self.vm_session is None:
            return "Session unavailable."
        active = self.vm_session.selections[self.vm_session.current_selection]
        # [EN] same bug/fix as cmd_select() above -- selection_function_
        # viewing_set(None) is the OFFICIAL clear path (also used by
        # deselecting via a normal click), and it correctly clears
        # vm_object.selected_atom_ids for every object, not just this
        # selection's own bookkeeping -- BUT it does NOT reset
        # atom.selected on the atoms that were previously selected
        # (confirmed live, not just by reading the code -- see the note
        # in cmd_select()). Captured here the same way: grab the OLD
        # selected_atoms before clearing, then explicitly set
        # .selected = False on each.
        old_atoms = set(active.selected_atoms)  # copia -- selection_function_viewing_set(None) limpa o set original NO LUGAR
        active.selection_function_viewing_set(None)
        for a in old_atoms:
            a.selected = False
        self.vm_session.vm_glcore.queue_draw()
        return "Selection cleared."

    def cmd_frame(self, n=None, **_):
        """ Without n=: shows the current frame number. With n=<int>:
        jumps directly to that frame of the trajectory (applies to
        whichever object(s) have a loaded trajectory).

        Examples:
          frame
          frame n=0
          frame n=50
        """
        if self.vm_session is None:
            return "Session unavailable."
        if n is None:
            return "Current frame: {}".format(self.vm_session.get_frame())
        self.vm_session.set_frame(frame=int(n))
        return "Frame -> {}".format(int(n))

    def cmd_next(self, **_):
        """ Steps one frame forward in the trajectory (same as clicking
        the ">" button on the trajectory player).

        Examples:
          next
        """
        if self.vm_session is None:
            return "Session unavailable."
        self.vm_session.forward_frame()
        return "Frame -> {}".format(self.vm_session.get_frame())

    def cmd_prev(self, **_):
        """ Steps one frame backward in the trajectory (same as clicking
        the "<" button on the trajectory player).

        Examples:
          prev
        """
        if self.vm_session is None:
            return "Session unavailable."
        self.vm_session.reverse_frame()
        return "Frame -> {}".format(self.vm_session.get_frame())

    def cmd_axes(self, show=True, **_):
        """ Shows or hides the XYZ axes gizmo in the 3D view.

        show : true | false   (default true)

        Examples:
          axes show=true
          axes show=false
        """
        if self.vm_session is None:
            return "Session unavailable."
        if show:
            self.vm_session.show_axes()
            return "Axes visible."
        self.vm_session.hide_axes()
        return "Axes hidden."

    def cmd_center(self, obj=None, chain=None, resi=None, resn=None,
                   name=None, symbol=None,
                   not_chain=None, not_resi=None, not_resn=None,
                   not_name=None, not_symbol=None, **_):
        """ Centers the camera on a target (animates the translation to it).

            Without obj=: centers on the ACTIVE SELECTION (atom centroid).
            With obj= and no fine filters: object center of mass.
            With obj= and filters: centroid of the filtered atoms.

            obj=N  chain=A  resi=45|10-20  resn=HIS  name=CA  symbol=C
            not_chain=  not_resi=  not_resn=  not_name=  not_symbol=

            Examples:
              center                       (on the active selection)
              center obj=0                 (center of mass of object 0)
              center obj=0 chain=A         (centroid of chain A)
              center obj=0 resn=HIS name=CA
              center obj=0 not_symbol=H    (centroid of everything except hydrogens)
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

        has_fine = any(f is not None for f in
                       (chain, resi, resn, name, symbol,
                        not_chain, not_resi, not_resn, not_name, not_symbol))
        if has_fine:
            atoms = self._filter_atoms(vm_object, chain, resi, resn, name, symbol,
                                        not_chain, not_resi, not_resn, not_name, not_symbol)
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

    def cmd_color(self, obj=None, color=None, r=None, g=None, b=None,
                  chain=None, resi=None, resn=None, name=None, symbol=None,
                  not_chain=None, not_resi=None, not_resn=None,
                  not_name=None, not_symbol=None, **_):
        """ Changes the colour of atoms, permanently (overrides the
        default periodic-table colour for those atoms specifically --
        not just a temporary highlight like 'select').

        Give the colour either by name (color=) or by explicit RGB
        (r=/g=/b=, each 0.0-1.0) -- not both.

        Two ways to pick WHICH atoms, same pattern as show/hide/center:
          - no obj=  -> colors the ACTIVE SELECTION (from 'select').
          - obj=N [+ chain=/resi=/resn=/name=/symbol=/not_*] -> colors
            exactly those atoms, one-off, without touching the active
            selection.

        obj    : object index (see 'list') -- omit to use the active selection
        color  : a name (see the list at the end of this help text)
        r,g,b  : explicit RGB, each 0.0 to 1.0 (use instead of color=)
        chain, resi, resn, name, symbol : same filters as 'select'
        not_chain, not_resi, not_resn, not_name, not_symbol : same
                 exclude filters as 'select'

        Examples:
          select obj=0 not_symbol=H
          color color=gray                              (colors the selection above)

          color obj=0 color=red                         (whole object red)
          color obj=0 symbol=H color=white               (all hydrogens white)
          color obj=0 not_symbol=H color=gray             (everything but H gray)
          color obj=0 chain=A resn=HIS color=yellow
          color obj=0 r=0.2 g=0.6 b=1.0                  (custom RGB, whole object)

        Available color names: (see NAMED_COLORS)
        """
        if self.vm_session is None:
            return "Session unavailable."

        if color is not None:
            rgb = NAMED_COLORS.get(str(color).lower())
            if rgb is None:
                return ("Unknown color name '{}'. Available: {}"
                        .format(color, ", ".join(sorted(NAMED_COLORS.keys()))))
        elif r is not None or g is not None or b is not None:
            rgb = (float(r or 0.0), float(g or 0.0), float(b or 0.0))
        else:
            return "Usage: give color=<name> OR r=/g=/b= (0.0-1.0 each)."

        # --- Decide WHICH atoms: active selection, or obj=+filters. ---
        if obj is None:
            active = self.vm_session.selections[self.vm_session.current_selection]
            atoms = list(active.selected_atoms)
            if not atoms:
                return "Nothing selected. Use 'select' or pass obj=N."
            source_desc = "the active selection"
        else:
            target = self._resolve_object(obj)
            if target is None:
                return "Object {} not found. Use 'list'.".format(obj)
            atoms = self._filter_atoms(target, chain, resi, resn, name, symbol,
                                        not_chain, not_resi, not_resn, not_name, not_symbol)
            if not atoms:
                return "No atom matches the filters in obj={}.".format(obj)
            source_desc = "obj={}".format(obj)

        for a in atoms:
            a.color = np.array(rgb, dtype=np.float32)

        # Rebuild the colours array read at draw time (same pattern used by
        # gui/windows/builder/atom_ops.py's add_atom()/remove_atom() -- see
        # VismolObject._generate_color_vectors(), which fills
        # self.colors[i] = atom.color for every atom), then force the
        # representations that actually carry a colour VBO to rebuild it.
        # Done per AFFECTED OBJECT, not just one: the active selection can
        # in principle span more than one vismol_object (selected_objects
        # is a set, not a single object), so grouping atoms by their own
        # .vm_object first avoids only recoloring whichever object
        # happened to be checked first.
        affected = {}
        for a in atoms:
            affected.setdefault(a.vm_object, []).append(a)
        for vm_object in affected:
            vm_object._generate_color_vectors(self.vm_session.atom_id_counter)
            vm_object.create_representation(rep_type="lines")
            vm_object.create_representation(rep_type="nonbonded")
        self.vm_session.vm_glcore.queue_draw()

        return "Colored {} atom(s) in {} -> RGB({:.2f}, {:.2f}, {:.2f}).".format(
            len(atoms), source_desc, rgb[0], rgb[1], rgb[2])

    # [EN] BUG CAUGHT BEFORE SHIPPING (while writing this method, not by
    # a later test): tried to build the "Available color names" part of
    # the docstring above by writing """...""".format(names=...) directly
    # -- but that turns the docstring from a bare string literal into a
    # Call expression, and Python ONLY recognises a function's __doc__
    # from a bare string literal as the first statement in its body. The
    # .format() call would have silently evaluated and discarded its
    # result, leaving cmd_color.__doc__ as None -- breaking
    # 'help cmd=color' (cmd_help() uses inspect.getdoc(), which returns
    # None/"(sem descricao)" for an undocumented method) without any
    # error or warning anywhere. Fixed by keeping the docstring a plain
    # literal (with a placeholder line instead of an f-string/.format
    # call) and appending the actual, dynamically-built colour list
    # here instead, right after the method is defined -- assigning
    # cmd_color.__doc__ directly still works from inside the class body,
    # since `cmd_color` is just a regular local name (the function object)
    # at this point in class construction.
    cmd_color.__doc__ = cmd_color.__doc__.replace(
        "Available color names: (see NAMED_COLORS)",
        "Available color names: " + ", ".join(sorted(NAMED_COLORS.keys())))

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

    def cmd_new(self, name="new_molecule", **_):
        """ Creates an empty, editable molecule object -- the starting
        point of the Builder (draw a molecule from scratch). Also
        creates a matching (initially empty) pDynamo System linked to
        it right away (see gui/windows/builder/empty_object.py's
        sync_pdynamo_system() for how, and its own module docstring for
        why this is the LEAST-verified part of the whole Builder feature
        set -- test this specific piece first if anything looks off).
        The new object appears at index len(list)-1 in 'list', nested
        under its own system row in the main treeview.

        name : the object's name (default "new_molecule")

        Examples:
          new
          new name=my_ligand
        """
        if self.vm_session is None:
            return "Session unavailable."
        from gui.windows.builder.empty_object import create_empty_vismol_object
        vobj = create_empty_vismol_object(self.vm_session, name=name)
        return "Created empty builder object: {} (index {}, 0 atoms)".format(vobj.name, vobj.index)

    def cmd_placemode(self, on=True, obj=None, symbol="C", **_):
        """ Toggles click-to-place-atom mode: while on, plain left
        clicks (no shift, no drag) on the 3D view add new atoms of
        `symbol` to `obj`, at the clicked position, instead of
        selecting whatever is already there. Starts in the "add" tool
        (see 'tool' to switch to "delete"). Shift-click still selects
        normally while this mode is on (needed to pick 2 atoms for
        'bond').

        on     : true | false           (default true)
        obj    : object index (see 'list') -- required when on=true
        symbol : element to place        (default "C")

        Examples:
          placemode on=true obj=0 symbol=C
          placemode on=true obj=0 symbol=O
          placemode on=false
        """
        if self.vm_session is None:
            return "Session unavailable."
        from gui.windows.builder.click_mode import enable_atom_placement_mode, disable_atom_placement_mode
        if not on or (isinstance(on, str) and on.lower() in ("false", "off", "0")):
            disable_atom_placement_mode(self.vm_session)
            return "Place-atom mode OFF."
        if obj is None:
            return "Usage: placemode on=true obj=<index> symbol=C"
        try:
            vobj = self.vm_session.vm_objects_dic[int(obj)]
        except (KeyError, ValueError):
            return "Object '{}' not found. Use 'list' to see the indices.".format(obj)
        enable_atom_placement_mode(self.vm_session, vobj, symbol=symbol)
        return ("Place-atom mode ON -- object '{}', element {}. "
                "LEFT-click on the 3D view (no dragging) to "
                "add atoms. Use 'placemode on=false' to turn it off.").format(vobj.name, symbol)

    def cmd_tool(self, name="add", **_):
        """ Switches which action a plain click performs while
        'placemode' is on -- same as pressing 'a'/'d' on the keyboard.
        Requires placemode to already be on (use 'placemode' first).

        name : add    -- plain click places a new atom (the default)
               delete -- plain click removes the clicked atom

        Examples:
          tool name=add
          tool name=delete
        """
        if self.vm_session is None:
            return "Session unavailable."
        if not getattr(self.vm_session, "builder_atom_mode", False):
            return "Builder is not on -- use 'placemode' first."
        from gui.windows.builder.click_mode import set_tool
        try:
            set_tool(self.vm_session, name)
        except ValueError as e:
            return str(e)
        return "Builder tool = '{}'.".format(name)

    def cmd_bond(self, obj=None, atom1=None, atom2=None, order=1, frame=None, **_):
        """ Adds a bond, or updates an existing one's order, at a chosen
        bond order. Two ways to use, tried in this order:
              bond obj=<indice> atom1=<id> atom2=<id> order=1
                                            -- bonds two atoms directly by
                                               atom_id, no picking needed
                                               (handy for terminal-only
                                               use).
              bond                          -- bonds the 2 atoms currently
                                                marked pk1/pk2 in the
                                                MEASUREMENT PICKING tool
                                                (the same pk1..pk4 used for
                                                on-screen distances/angles/
                                                dihedrals). This is the
                                                default when obj/atom1/
                                                atom2 aren't given.

        order : 1 (single, default), 2 (double) or 3 (triple) -- the
                bond's persisted order (see atom_ops.set_bond_order()'s
                own docstring).

        frame : NOT GIVEN (default) -> edits the object's STATIC
                topology (vismol_object.bonds), exactly as described
                above -- valid for every frame, this is what gets
                saved/exported.
                GIVEN -> edits the DYNAMIC BONDS representation instead
                (the per-frame connectivity used by the Dynamic Bonds
                display, typically the QC region of a QM/MM trajectory),
                for the frame(s) selected by the value:
                  frame=true            current frame only
                  frame=12              frame 12 only
                  frame=1:20            frames 1 through 20, INCLUSIVE
                  frame=all             every frame

                *** WARNING / READ BEFORE USING frame=...: this ONLY
                changes what is DRAWN on screen for the selected
                frame(s) -- it does NOT create a real chemical bond, and
                does NOT touch the linked pDynamo system's topology,
                force field, bond force constants, charges, or anything
                else used for an actual QM/MM calculation. It is meant
                for visually inspecting/correcting the Dynamic Bonds
                display (e.g. forcing a bond the automatic distance-based
                detection missed on one particular frame) -- NOT for
                defining the chemistry of the system. To create a real
                bond with proper force-field parameters, the pDynamo
                system's topology needs to be edited through other means.

        [EN] BUG FIX: previously used atom_ops.add_bond(), same as the
        Builder's own 'b' keyboard shortcut. Confirmed via live testing
        that this is WRONG here: add_bond() is built around the Builder
        canvas's own design (vismol_object.manual_bonds as the ONLY
        source of truth for connectivity -- see add_atom()'s docstring),
        which works for an empty object grown atom-by-atom, but silently
        DROPS every other bond when used on a normally-loaded structure
        (e.g. a PDB file), since those bonds were never registered in
        manual_bonds. Switched to atom_ops.set_bond_order(), which edits
        only the one pair involved (adds it if missing, or just updates
        its order if it already exists) and leaves every other bond
        exactly as it was, regardless of where it originally came from.
        Also no longer auto-adjusts hydrogens (that was a Builder-canvas
        convenience, not appropriate for editing bond order on an
        already-complete loaded structure) -- see set_bond_order()'s own
        docstring in atom_ops.py for the full explanation.

        Syncs the linked pDynamo system afterwards (static path only --
        frame=... is representation-only, see the warning above).

        Examples:
          bond
          bond order=2
          bond obj=0 atom1=3 atom2=7 order=2
          bond order=2 frame=true
          bond obj=0 atom1=3 atom2=7 order=2 frame=12
          bond order=2 frame=1:20
          bond order=1 frame=all
        """
        if self.vm_session is None:
            return "Session unavailable."
        try:
            order = int(order)
        except (TypeError, ValueError):
            return "order must be 1 (single), 2 (double) or 3 (triple)."
        if order not in (1, 2, 3):
            return "order must be 1 (single), 2 (double) or 3 (triple)."

        if obj is not None and atom1 is not None and atom2 is not None:
            try:
                vobj = self.vm_session.vm_objects_dic[int(obj)]
            except (KeyError, ValueError):
                return "Object '{}' not found. Use 'list' to see the indices.".format(obj)

            from gui.windows.builder.atom_ops import resolve_frame_arg
            try:
                frames = resolve_frame_arg(vobj, frame)
            except ValueError as e:
                return str(e)

            if frames is not None:
                from gui.windows.builder.atom_ops import set_dynamic_bond_order, push_undo_snapshot
                push_undo_snapshot(vobj)
                try:
                    n_created = set_dynamic_bond_order(vobj, int(atom1), int(atom2),
                                                        bond_order=order, frames=frames)
                except ValueError as e:
                    return str(e)
                frame_desc = "frame {}".format(frames[0]) if len(frames) == 1 else "{} frames".format(len(frames))
                return ("[Dynamic Bonds] Bond order between atom {} and atom {} set to {} for {} "
                        "({} new pair(s) added). NOTE: representation-only, does not change the "
                        "pDynamo system's real topology.".format(atom1, atom2, order, frame_desc, n_created))

            from gui.windows.builder.atom_ops import set_bond_order, push_undo_snapshot
            push_undo_snapshot(vobj)
            try:
                created = set_bond_order(vobj, int(atom1), int(atom2), bond_order=order)
            except ValueError as e:
                return str(e)
            from gui.windows.builder.empty_object import sync_pdynamo_system
            sync_pdynamo_system(vobj)
            return ("Bond created between atom {} and atom {} (order={}).".format(atom1, atom2, order) if created
                    else "Bond order between atom {} and atom {} set to {}.".format(atom1, atom2, order))

        from gui.windows.builder.click_mode import handle_bond_picking
        return handle_bond_picking(self.vm_session, bond_order=order, frame=frame)

    def cmd_unbond(self, obj=None, atom1=None, atom2=None, frame=None, **_):
        """ Removes a bond. Two ways to use, tried in this order:
              unbond obj=<indice> atom1=<id> atom2=<id>
                                            -- removes the bond between two
                                               atoms directly by atom_id,
                                               no picking needed.
              unbond                        -- removes the bond between
                                                the 2 atoms currently
                                                marked pk1/pk2 in the
                                                MEASUREMENT PICKING tool
                                                (same pk1..pk4 used for
                                                on-screen distances/
                                                angles/dihedrals). This is
                                                the default when obj/
                                                atom1/atom2 aren't given.

        frame : NOT GIVEN (default) -> edits the object's STATIC topology,
                exactly as described above.
                GIVEN -> removes the pair from the DYNAMIC BONDS
                representation instead, for the frame(s) selected by the
                value -- same accepted forms as 'bond' (true/N/'A:B'/'all').

        *** SAME WARNING AS 'bond': frame=... only changes what is drawn
        for the selected frame(s) -- it does NOT change the pDynamo
        system's real topology/force field. See 'bond's own docstring for
        the full explanation.

        [EN] BUG FIX: previously used atom_ops.remove_bond(). Same problem
        as 'bond' above, worse: remove_bond() refuses to touch any bond
        not explicitly recorded in manual_bonds AND resets
        vismol_object.index_bonds = None before rebuilding purely from
        manual_bonds -- confirmed via live testing that this drops every
        OTHER bond in a normally-loaded structure, not just fails to
        remove the intended one. Switched to atom_ops.unset_bond(), which
        works on any bond present in vismol_object.bonds regardless of
        its origin, and leaves everything else untouched.

        Syncs the linked pDynamo system afterwards, if a bond was removed
        (static path only -- frame=... is representation-only).

        Examples:
          unbond
          unbond obj=0 atom1=3 atom2=7
          unbond frame=true
          unbond obj=0 atom1=3 atom2=7 frame=12
          unbond frame=1:20
          unbond frame=all
        """
        if self.vm_session is None:
            return "Session unavailable."
        if obj is not None and atom1 is not None and atom2 is not None:
            try:
                vobj = self.vm_session.vm_objects_dic[int(obj)]
            except (KeyError, ValueError):
                return "Object '{}' not found. Use 'list' to see the indices.".format(obj)

            from gui.windows.builder.atom_ops import resolve_frame_arg
            try:
                frames = resolve_frame_arg(vobj, frame)
            except ValueError as e:
                return str(e)

            if frames is not None:
                from gui.windows.builder.atom_ops import unset_dynamic_bond, push_undo_snapshot
                push_undo_snapshot(vobj)
                n_removed = unset_dynamic_bond(vobj, int(atom1), int(atom2), frames=frames)
                frame_desc = "frame {}".format(frames[0]) if len(frames) == 1 else "{} frames".format(len(frames))
                return ("[Dynamic Bonds] Bond removed between atom {} and atom {} in {} of {} requested. "
                        "NOTE: representation-only, does not change the pDynamo system's real "
                        "topology.".format(atom1, atom2, n_removed, frame_desc))

            from gui.windows.builder.atom_ops import unset_bond, push_undo_snapshot
            push_undo_snapshot(vobj)
            removed = unset_bond(vobj, int(atom1), int(atom2))
            if removed:
                from gui.windows.builder.empty_object import sync_pdynamo_system
                sync_pdynamo_system(vobj)
            return ("Bond removed between atom {} and atom {}.".format(atom1, atom2) if removed
                    else "No bond existed between these 2 atoms (nothing removed).")

        from gui.windows.builder.click_mode import handle_unbond_picking
        return handle_unbond_picking(self.vm_session, frame=frame)

    def cmd_add(self, obj=None, symbol="C", x=0.0, y=0.0, z=0.0, **_):
        """ Adds one atom to a Builder object, at an explicit position
        (in Angstrom) -- the terminal-only equivalent of clicking on the
        3D view in 'placemode'. Adjusts this atom's hydrogens afterwards
        (it starts with zero bonds, so this just gives it its full
        complement -- see gui/windows/builder/atom_ops.
        adjust_hydrogens()) and syncs the linked pDynamo system.

        obj    : object index (see 'list')
        symbol : element symbol            (default "C")
        x, y, z: position in Angstrom       (default 0.0, 0.0, 0.0)

        Examples:
          add obj=0 symbol=C x=0.0 y=0.0 z=0.0
          add obj=0 symbol=O x=1.2 y=0.0 z=0.0
        """
        if self.vm_session is None:
            return "Session unavailable."
        if obj is None:
            return "Usage: add obj=<index> symbol=C x=0.0 y=0.0 z=0.0"
        try:
            vobj = self.vm_session.vm_objects_dic[int(obj)]
        except (KeyError, ValueError):
            return "Object '{}' not found. Use 'list' to see the indices.".format(obj)
        from gui.windows.builder.atom_ops import add_atom, adjust_hydrogens, push_undo_snapshot
        push_undo_snapshot(vobj)
        atom = add_atom(vobj, symbol=symbol, x=float(x), y=float(y), z=float(z))
        adjust_hydrogens(vobj, atom.atom_id)
        from gui.windows.builder.empty_object import sync_pdynamo_system
        sync_pdynamo_system(vobj)
        return "Atom added: {} #{} at ({}, {}, {}) -- object now has {} atom(s)".format(
            atom.symbol, atom.atom_id, x, y, z, len(vobj.atoms))

    def cmd_delete(self, obj=None, atom=None, **_):
        """ Removes one atom by atom_id from a Builder object --
        terminal-only equivalent of clicking on it in the "delete" tool
        ('tool name=delete'). Remaining atoms with a higher atom_id get
        renumbered down by one automatically (atom_id must stay a dense,
        contiguous 0..N-1 index). Adjusts the FORMER neighbours'
        hydrogens afterwards (unless the removed atom was itself a
        hydrogen -- see gui/windows/builder/atom_ops.adjust_hydrogens()'s
        call site in vismol_glcore.py for why: instantly replacing a
        deliberately-deleted H would make the deletion a no-op) and
        syncs the linked pDynamo system.

        obj  : object index (see 'list')
        atom : atom_id to remove (0-based)

        Examples:
          delete obj=0 atom=2
        """
        if self.vm_session is None:
            return "Session unavailable."
        if obj is None or atom is None:
            return "Usage: delete obj=<index> atom=<atom_id>"
        try:
            vobj = self.vm_session.vm_objects_dic[int(obj)]
        except (KeyError, ValueError):
            return "Object '{}' not found. Use 'list' to see the indices.".format(obj)
        from gui.windows.builder.atom_ops import remove_atom, adjust_hydrogens, push_undo_snapshot

        deleted_id = int(atom)
        if deleted_id not in vobj.atoms:
            return "Atom {} does not exist in this object.".format(deleted_id)

        deleted_symbol = vobj.atoms[deleted_id].symbol
        neighbor_objs = []
        for bond in vobj.bonds.values():
            if bond.atom_index_i != deleted_id and bond.atom_index_j != deleted_id:
                continue
            other_id = bond.atom_index_j if bond.atom_index_i == deleted_id else bond.atom_index_i
            neighbor = vobj.atoms[other_id]
            if neighbor.symbol != 'H':
                neighbor_objs.append(neighbor)

        push_undo_snapshot(vobj)
        try:
            remove_atom(vobj, deleted_id)
        except ValueError as e:
            return str(e)

        if deleted_symbol != 'H':
            for neighbor in neighbor_objs:
                adjust_hydrogens(vobj, neighbor.atom_id)

        from gui.windows.builder.empty_object import sync_pdynamo_system
        sync_pdynamo_system(vobj)
        return "Atom {} removed -- object now has {} atom(s).".format(atom, len(vobj.atoms))

    def cmd_load(self, file=None, **_):
        """ Loads a molecule from a file on disk (PDB, XYZ, and any
        other format the existing file loaders support), adding it as a
        new object in the session (see 'list' afterwards for its index).

        file : full path to the file

        Examples:
          load file=/home/user/structure.pdb
          load file=/home/user/molecule.xyz
        """
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
    def bond(self, obj=None, atom1=None, atom2=None, order=1, frame=None):
        return self.cmd_bond(obj=obj, atom1=atom1, atom2=atom2, order=order, frame=frame)
    def unbond(self, obj=None, atom1=None, atom2=None, frame=None):
        return self.cmd_unbond(obj=obj, atom1=atom1, atom2=atom2, frame=frame)
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
