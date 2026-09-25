#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Packmol component editor dialog
#
#  Copyright 2022-2026 Fernando Bachega
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
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import math
import os

from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import CoordinatesComboBox
from pdynamo.pDynamo2EasyHybrid.helpers import export_special_PDB
from pBabel import ImportSystem

# . Same order as the glade file's combo_constraint_type <items> list --
#   index i here must match that combo's row i (Gtk.ComboBoxText has no
#   per-row id lookup by string once built, only by index).
_CONSTRAINT_TYPES = ["inside_box", "outside_box", "inside_sphere", "outside_sphere", "fixed"]

# . Same order as the glade file's combo_source_type <items> list.
_SOURCE_TYPES = ["file", "vobject"]

# . Same order as the glade file's combo_number_mode <items> list.
_NUMBER_MODES = ["fixed", "density"]

# . Density-based counts only make sense for the two constraints that
#   enclose a well-defined, finite volume -- "outside" constraints are
#   open-ended relative to the packing region, and "fixed" places
#   exactly one copy by definition.
_VOLUME_CONSTRAINT_TYPES = ("inside_box", "inside_sphere")

_AVOGADRO_NUMBER = 6.02214076e23  # mol^-1
_ANGSTROM3_TO_CM3 = 1.0e-24


def _constraint_volume_A3(constraint_type, params):
    """ Volume enclosed by an inside_box/inside_sphere constraint, in
        cubic Angstrom -- None for any other constraint_type (no finite,
        well-defined volume to compute a density against).
    """
    if constraint_type == "inside_box":
        x1, y1, z1, x2, y2, z2 = params
        return abs(x2 - x1) * abs(y2 - y1) * abs(z2 - z1)
    if constraint_type == "inside_sphere":
        _cx, _cy, _cz, radius = params
        return (4.0 / 3.0) * math.pi * (radius ** 3)
    return None


def _molar_mass_g_per_mol(structure_path):
    """ Sums atomic weights (g/mol) straight from the structure file's
        own atoms via pDynamo3's ImportSystem() -- the same reader this
        whole app already relies on for PDB/XYZ/MOL2/... (rather than a
        bespoke element-guessing parser here), so this works uniformly
        whether structure_path is a user-picked file or one this dialog
        just exported from a loaded object. Returns None (component's
        caller shows the dialog) if the file can't be read at all.
    """
    try:
        system = ImportSystem(structure_path)
        return sum(atom.mass for atom in system.atoms)
    except Exception:
        return None


class AddPackmolComponentDialog:
    """ Modal dialog for adding or editing a single Packmol component
        (structure file + copy count + geometric constraint) -- one row
        of PreparePackmolWindow's component treeview. Same blocking
        `dialog.run()` + `self.ok`/result-attributes pattern as
        AddHarmonicRestraintDialog (dialogs/restraint_dialog.py).

        Pass `component=None` for "Add" (defaults, title/button say
        "Add"); pass an existing component dict (see
        packmol_runner.build_packmol_input()'s docstring for its shape)
        for "Edit" (fields prefilled, title/button say "Save").
    """

    def __init__(self, main=None, component=None):
        self.main = main
        self.home = main.home

        self.builder = Gtk.Builder()
        self.builder.add_from_file(
            os.path.join(self.home, 'src/gui/windows/setup/add_packmol_component_dialog.glade')
        )
        self.builder.connect_signals(self)

        self.dialog = self.builder.get_object('dialog')
        self.label_structure_path = self.builder.get_object('label_structure_path')
        self.spinbtn_component_number = self.builder.get_object('spinbtn_component_number')
        self.combo_constraint_type = self.builder.get_object('combo_constraint_type')
        self.grid_box_params = self.builder.get_object('grid_box_params')
        self.grid_sphere_params = self.builder.get_object('grid_sphere_params')
        self.grid_fixed_params = self.builder.get_object('grid_fixed_params')

        # . "Set number by: Density" -- see _NUMBER_MODES/_constraint_volume_A3()
        # /_molar_mass_g_per_mol() above.
        self.combo_number_mode      = self.builder.get_object('combo_number_mode')
        self.label_component_number = self.builder.get_object('label_component_number')
        self.label_density          = self.builder.get_object('label_density')
        self.spinbtn_density        = self.builder.get_object('spinbtn_density')
        self.label_computed_number  = self.builder.get_object('label_computed_number')

        self.spinbtn_box = [self.builder.get_object(name) for name in
                             ('spinbtn_box_x1', 'spinbtn_box_y1', 'spinbtn_box_z1',
                              'spinbtn_box_x2', 'spinbtn_box_y2', 'spinbtn_box_z2')]
        self.spinbtn_sphere = [self.builder.get_object(name) for name in
                                ('spinbtn_sphere_cx', 'spinbtn_sphere_cy',
                                 'spinbtn_sphere_cz', 'spinbtn_sphere_radius')]
        self.spinbtn_fixed = [self.builder.get_object(name) for name in
                               ('spinbtn_fixed_x', 'spinbtn_fixed_y', 'spinbtn_fixed_z',
                                'spinbtn_fixed_a', 'spinbtn_fixed_b', 'spinbtn_fixed_g')]

        self.structure_filechooser = FolderChooserButton(self.main, 'file', self.home)
        self.builder.get_object('box_structure_file').pack_start(
            self.structure_filechooser.btn, False, False, 0)
        self.structure_filechooser.btn.show_all()

        # . "Loaded object" source: a snapshot of an already-open vismol
        # object's CURRENT coordinates, exported to a PDB in the Packmol
        # job's own working folder the moment Add/Save is clicked (see
        # on_button_ok_clicked) -- from that point on it behaves exactly
        # like a file-based component (build_packmol_input() never knows
        # the difference), it just doesn't need the user to have saved a
        # file to disk first.
        self.combo_source_type   = self.builder.get_object('combo_source_type')
        self.label_structure_file = self.builder.get_object('label_structure_file')
        self.box_structure_file   = self.builder.get_object('box_structure_file')
        self.label_vobject        = self.builder.get_object('label_vobject')
        self.box_vobject_combo    = self.builder.get_object('box_vobject_combo')

        self.combobox_vobject = CoordinatesComboBox()
        e_id = self.main.p_session.active_id
        self.combobox_vobject.set_model(self.main.vobject_liststore_dict[e_id])
        size = len(self.main.vobject_liststore_dict[e_id])
        if size > 0:
            self.combobox_vobject.set_active(size - 1)
        self.box_vobject_combo.pack_start(self.combobox_vobject, True, True, 0)
        self.combobox_vobject.show_all()

        self.combo_constraint_type.connect('changed', self._on_constraint_type_changed)
        self.combo_constraint_type.connect('changed', self._on_density_inputs_changed)
        for spin in self.spinbtn_box + self.spinbtn_sphere:
            spin.connect('value-changed', self._on_density_inputs_changed)

        editing = component is not None
        self.dialog.set_title('Edit Packmol Component' if editing else 'Add Packmol Component')
        self.builder.get_object('button_ok').set_label('Save' if editing else 'Add')

        if editing:
            self._prefill(component)
        else:
            self._on_constraint_type_changed(self.combo_constraint_type)
            self._on_number_mode_changed(self.combo_number_mode)

        # . Editing always resolves to the plain "file" view: we don't
        # track whether an existing component's structure_path
        # originally came from a vobject export, so Edit just shows the
        # already-resolved path -- switching to "Loaded object" and
        # picking one re-exports and overwrites it on Save, same as Add.
        self._on_source_type_changed(self.combo_source_type)

        self.builder.get_object('button_cancel').connect('clicked', self.close_window)
        self.builder.get_object('button_ok').connect('clicked', self.on_button_ok_clicked)
        self.dialog.connect('destroy', self.close_window)

        self.component = None   # result dict, filled in on confirm
        self.ok = False

        self.dialog.run()

    def _prefill(self, component):
        structure_path = component.get('structure_path')
        if structure_path:
            self.structure_filechooser.set_folder(folder=structure_path)
            self.label_structure_path.set_text(structure_path)

        self.spinbtn_component_number.set_value(float(component.get('number', 1)))

        number_mode = component.get('number_mode', 'fixed')
        try:
            self.combo_number_mode.set_active(_NUMBER_MODES.index(number_mode))
        except ValueError:
            self.combo_number_mode.set_active(0)
        if 'density' in component:
            self.spinbtn_density.set_value(float(component['density']))

        constraint_type = component.get('constraint_type', 'inside_box')
        try:
            index = _CONSTRAINT_TYPES.index(constraint_type)
        except ValueError:
            index = 0
        self.combo_constraint_type.set_active(index)

        params = component.get('params', ())
        if constraint_type in ('inside_box', 'outside_box'):
            for spin, value in zip(self.spinbtn_box, params):
                spin.set_value(value)
        elif constraint_type in ('inside_sphere', 'outside_sphere'):
            for spin, value in zip(self.spinbtn_sphere, params):
                spin.set_value(value)
        elif constraint_type == 'fixed':
            for spin, value in zip(self.spinbtn_fixed, params):
                spin.set_value(value)

        self._on_constraint_type_changed(self.combo_constraint_type)
        self._on_number_mode_changed(self.combo_number_mode)

    def _on_constraint_type_changed(self, combo):
        constraint_type = _CONSTRAINT_TYPES[combo.get_active()]
        is_box = constraint_type in ('inside_box', 'outside_box')
        is_sphere = constraint_type in ('inside_sphere', 'outside_sphere')
        is_fixed = constraint_type == 'fixed'

        self.grid_box_params.set_visible(is_box)
        self.grid_sphere_params.set_visible(is_sphere)
        self.grid_fixed_params.set_visible(is_fixed)

        # . Packmol errors out if `number` > 1 is given for a fixed
        #   structure (see packmol_runner.py's CONSTRAINT_TYPES docstring)
        #   -- lock the spinbutton at 1 while "Fixed" is selected so the
        #   UI can't produce an invalid .inp.
        self.spinbtn_component_number.set_sensitive(not is_fixed)
        if is_fixed:
            self.spinbtn_component_number.set_value(1)

    def _on_source_type_changed(self, combo):
        source_type = _SOURCE_TYPES[combo.get_active()]
        is_file = source_type == 'file'

        self.label_structure_file.set_visible(is_file)
        self.box_structure_file.set_visible(is_file)
        self.label_vobject.set_visible(not is_file)
        self.box_vobject_combo.set_visible(not is_file)

    def _on_number_mode_changed(self, combo):
        number_mode = _NUMBER_MODES[combo.get_active()]
        is_density = number_mode == 'density'

        self.label_component_number.set_visible(not is_density)
        self.spinbtn_component_number.set_visible(not is_density)
        self.label_density.set_visible(is_density)
        self.spinbtn_density.set_visible(is_density)
        self.label_computed_number.set_visible(is_density)

        self._on_density_inputs_changed()

    def _on_density_inputs_changed(self, widget=None):
        """ Live preview of the density -> molecule-count conversion --
            recomputed on every constraint-type/box/sphere/density change
            (see the signal connections in __init__). Purely informative:
            on_button_ok_clicked() always redoes this calculation itself
            against the FINAL resolved structure_path rather than trust
            whatever this label last showed.
        """
        if _NUMBER_MODES[self.combo_number_mode.get_active()] != 'density':
            return
        count, error = self._compute_density_count()
        if error:
            self.label_computed_number.set_text('({})'.format(error))
        else:
            self.label_computed_number.set_text('≈ {:d} molecules'.format(count))

    def _current_molar_mass_g_per_mol_for_preview(self):
        """ Molar mass for the live preview -- unlike the final,
            authoritative calculation in on_button_ok_clicked() (which
            always has a real structure_path to read, since the vobject
            source is already exported to PDB by that point), the preview
            can run before any export has happened, so a "Loaded object"
            source is summed directly from the vobject's own atoms
            instead of requiring a throwaway export just to compute this.
        """
        source_type = _SOURCE_TYPES[self.combo_source_type.get_active()]
        if source_type == 'vobject':
            tree_iter = self.combobox_vobject.get_active_iter()
            if tree_iter is None:
                return None
            model = self.combobox_vobject.get_model()
            _name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            periodic_table = self.main.vm_session.periodic_table
            try:
                return sum(periodic_table.get_atomic_mass(symbol=atom.symbol) for atom in vobject.atoms.values())
            except Exception:
                return None
        else:
            structure_path = self.structure_filechooser.get_folder()
            if not structure_path or not os.path.isfile(structure_path):
                return None
            return _molar_mass_g_per_mol(structure_path)

    def _compute_density_count(self, structure_path=None):
        """ Returns (count, error_message) -- exactly one of the two is
            None. Pass structure_path once it's known (from
            on_button_ok_clicked, after the "Loaded object" source, if
            any, has already been exported) to compute the molar mass
            from that exact file instead of the live-preview shortcut in
            _current_molar_mass_g_per_mol_for_preview().
        """
        constraint_type = _CONSTRAINT_TYPES[self.combo_constraint_type.get_active()]
        if constraint_type not in _VOLUME_CONSTRAINT_TYPES:
            return None, 'Density only works with the "Inside box" or "Inside sphere" constraints.'

        if constraint_type == 'inside_box':
            params = tuple(spin.get_value() for spin in self.spinbtn_box)
        else:
            params = tuple(spin.get_value() for spin in self.spinbtn_sphere)
        volume_A3 = _constraint_volume_A3(constraint_type, params)
        if not volume_A3 or volume_A3 <= 0:
            return None, 'The constraint has zero volume.'

        if structure_path is not None:
            molar_mass = _molar_mass_g_per_mol(structure_path)
        else:
            molar_mass = self._current_molar_mass_g_per_mol_for_preview()
        if not molar_mass:
            return None, 'Choose a valid structure first.'

        density = self.spinbtn_density.get_value()
        volume_cm3 = volume_A3 * _ANGSTROM3_TO_CM3
        mass_g = density * volume_cm3
        moles = mass_g / molar_mass
        count = int(round(moles * _AVOGADRO_NUMBER))
        return max(count, 1), None

    def _export_selected_vobject_to_pdb(self):
        """ Writes the vobject currently selected in combobox_vobject to a
            PDB file inside the Packmol job's own working folder (falling
            back to home if no folder has been chosen there yet), and
            returns the path -- or None (having already shown an error
            dialog) if nothing valid is selected.
        """
        tree_iter = self.combobox_vobject.get_active_iter()
        if tree_iter is None:
            self.main.simple_dialog.info(msg='Please choose a loaded object.')
            return None
        model = self.combobox_vobject.get_model()
        name, vobject_id = model[tree_iter][:2]
        vobject = self.main.vm_session.vm_objects_dic[vobject_id]

        work_folder = self.main.prepare_packmol_window.work_folder_chooser.get_folder()
        if not work_folder or not os.path.isdir(work_folder):
            work_folder = self.home
        export_dir = os.path.join(work_folder, 'packmol_objects')
        os.makedirs(export_dir, exist_ok=True)

        safe_name = ''.join(c if (c.isalnum() or c in '-_') else '_' for c in name)
        output_path = os.path.join(export_dir, safe_name + '.pdb')

        # . frame=-1 (last frame) matches every other export_special_PDB()
        # caller in this codebase (prepare_namd_run.py, prepare_smd_run.py,
        # prepare_amber_system.py, ...).
        export_special_PDB(vobject=vobject, frame=-1, output=output_path)
        return output_path

    def on_button_ok_clicked(self, widget):
        source_type = _SOURCE_TYPES[self.combo_source_type.get_active()]
        if source_type == 'vobject':
            structure_path = self._export_selected_vobject_to_pdb()
            if structure_path is None:
                return
            self.label_structure_path.set_text(structure_path)
        else:
            structure_path = self.structure_filechooser.get_folder()
            if not structure_path or not os.path.isfile(structure_path):
                self.main.simple_dialog.info(msg='Please choose a valid structure file.')
                return

        constraint_type = _CONSTRAINT_TYPES[self.combo_constraint_type.get_active()]
        if constraint_type in ('inside_box', 'outside_box'):
            params = tuple(spin.get_value() for spin in self.spinbtn_box)
        elif constraint_type in ('inside_sphere', 'outside_sphere'):
            params = tuple(spin.get_value() for spin in self.spinbtn_sphere)
        else:
            params = tuple(spin.get_value() for spin in self.spinbtn_fixed)

        number_mode = _NUMBER_MODES[self.combo_number_mode.get_active()]
        density = None
        if number_mode == 'density':
            # . Recompute fresh against the now-final structure_path
            # (for a "Loaded object" source, this is the file just
            # exported above -- not the live vobject-atoms shortcut the
            # preview label uses) rather than trust anything already
            # showing in label_computed_number.
            count, error = self._compute_density_count(structure_path=structure_path)
            if error:
                self.main.simple_dialog.info(msg=error)
                return
            number = count
            density = self.spinbtn_density.get_value()
        else:
            # . Packmol errors out if `number` > 1 is given for a fixed
            #   structure -- the spinbutton is desensitized at 1 while
            #   "Fixed" is selected (see _on_constraint_type_changed()), but
            #   clamp here too rather than trust that as the only guard.
            number = 1 if constraint_type == 'fixed' else int(self.spinbtn_component_number.get_value())

        self.component = {
            'structure_path': structure_path,
            'number': number,
            'number_mode': number_mode,
            'constraint_type': constraint_type,
            'params': params,
        }
        if density is not None:
            self.component['density'] = density
        self.ok = True
        self.dialog.destroy()

    def close_window(self, widget, data=None):
        self.dialog.destroy()
