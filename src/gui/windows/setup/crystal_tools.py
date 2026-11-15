"""
crystal_tools
=============

System > Cell and Symmetry:

  Center in Cell...     pSimulation.CrystalCenterCoordinates - every molecule
                        (or atom) moved into the primary cell, optionally with a
                        selection (e.g. the protein) at the centre of the box.
                        The result is a new object of the system.
  Expand to P1 / Supercell...
                        pSimulation.CrystalExpandToP1 - a new system with P1
                        symmetry, replicated na x nb x nc times.
"""

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk


def _active(main):
    p = main.p_session
    system = p.psystem.get(p.active_id)
    if system is None:
        raise ValueError("There is no active system.")
    if system.symmetry is None or system.symmetryParameters is None:
        raise ValueError("'{}' has no periodic cell (System > Cell and Symmetry > Edit Cell, or Solvate System).".format(system.label))
    objs = [v for v in main.vm_session.vm_objects_dic.values() if getattr(v, "e_id", None) == system.e_id]
    vobject = objs[-1] if objs else None
    return system, vobject


def _sync_from_view(main, system, vobject):
    """The system takes the coordinates of the frame shown."""
    if vobject is None or len(vobject.atoms) != len(system.atoms):
        return
    frame = min(int(main.vm_session.get_frame() or 0), len(vobject.frames) - 1)
    X = vobject.frames[frame]
    for i in range(len(system.atoms)):
        for k in range(3):
            system.coordinates3[i, k] = float(X[i][k])


def _error(main, text):
    main.simple_dialog.error(msg=text)


def center_in_cell(main):
    try:
        system, vobject = _active(main)
    except ValueError as error:
        _error(main, str(error))
        return
    vm = main.vm_session
    sel = vm.selections[vm.current_selection]
    selected = sorted({a.index - 1 for a in sel.selected_atoms if vobject is not None and a.vm_object is vobject})

    dialog = Gtk.Dialog(title="Center in Cell", transient_for=main.window, modal=True)
    dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, "Center", Gtk.ResponseType.OK)
    box = dialog.get_content_area()
    box.set_spacing(6)
    box.set_border_width(10)
    combo = Gtk.ComboBoxText()
    for text in ("Whole molecules (bonded groups)", "Atoms"):
        combo.append_text(text)
    combo.set_active(0)
    check = Gtk.CheckButton(label="Put the selection ({} atoms) at the centre of the cell".format(len(selected)))
    check.set_active(bool(selected))
    check.set_sensitive(bool(selected))
    grid = Gtk.Grid(column_spacing=8, row_spacing=6)
    grid.attach(Gtk.Label(label="Wrap", xalign=0), 0, 0, 1, 1)
    grid.attach(combo, 1, 0, 1, 1)
    grid.attach(check, 0, 1, 2, 1)
    note = Gtk.Label(xalign=0)
    note.set_markup("<small>pDynamo CrystalCenterCoordinates, on the frame shown. The result is a new\n"
                    "object of '{}'; QC atoms are kept together.</small>".format(GLibEscape(system.label)))
    grid.attach(note, 0, 2, 2, 1)
    box.add(grid)
    dialog.show_all()
    answer = dialog.run()
    mode = "byMMIsolate" if combo.get_active() == 0 else "byAtom"
    to_center = check.get_active() and bool(selected)
    dialog.destroy()
    if answer != Gtk.ResponseType.OK:
        return
    try:
        from pCore import Selection
        from pSimulation import CrystalCenterCoordinates
        _sync_from_view(main, system, vobject)
        if mode == "byMMIsolate" and system.mmModel is None:
            mode = "byIsolate"
        coords, _translations = CrystalCenterCoordinates(system, log=None, mode=mode,
                                                         toCenter=Selection.FromIterable(selected) if to_center else None)
        if coords is None:
            raise ValueError("pDynamo returned no coordinates.")
        system.coordinates3 = coords
        main.p_session._add_vismol_object_to_easyhybrid_session(system, True, name="centered")
        main.vm_session.vm_glcore.queue_draw()
    except Exception as error:
        _error(main, "Could not center the coordinates:\n{}".format(error))


def expand_to_p1(main):
    try:
        system, vobject = _active(main)
    except ValueError as error:
        _error(main, str(error))
        return
    sp = system.symmetryParameters
    dialog = Gtk.Dialog(title="Expand to P1 / Supercell", transient_for=main.window, modal=True)
    dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, "Create System", Gtk.ResponseType.OK)
    box = dialog.get_content_area()
    box.set_spacing(6)
    box.set_border_width(10)
    grid = Gtk.Grid(column_spacing=8, row_spacing=6)
    spins = []
    for k, (axis, length) in enumerate((("a", sp.a), ("b", sp.b), ("c", sp.c))):
        s = Gtk.SpinButton.new_with_range(1, 20, 1)
        s.set_value(1)
        spins.append(s)
        grid.attach(Gtk.Label(label="cells along {} ({:.2f} A)".format(axis, length), xalign=0), 0, k, 1, 1)
        grid.attach(s, 1, k, 1, 1)
    label = Gtk.Label(xalign=0)
    grid.attach(label, 0, 3, 2, 1)

    def update(*_):
        n = int(spins[0].get_value() * spins[1].get_value() * spins[2].get_value())
        label.set_text("{} cells x {} atoms = {} atoms".format(n, len(system.atoms), n * len(system.atoms)))
    for s in spins:
        s.connect("value-changed", update)
    update()
    note = Gtk.Label(xalign=0)
    note.set_markup("<small>pDynamo CrystalExpandToP1: every symmetry image of the asymmetric unit, in\n"
                    "na x nb x nc cells. A new system; the force field and the NB model are kept.</small>")
    grid.attach(note, 0, 4, 2, 1)
    box.add(grid)
    dialog.show_all()
    answer = dialog.run()
    na, nb, nc = (int(s.get_value()) for s in spins)
    dialog.destroy()
    if answer != Gtk.ResponseType.OK:
        return
    try:
        from pCore import Clone
        from pSimulation import CrystalExpandToP1
        _sync_from_view(main, system, vobject)
        triclinic = len({na, nb, nc}) > 1
        new = CrystalExpandToP1(system, arange=range(na), brange=range(nb), crange=range(nc),
                                imposeTriclinic=triclinic, log=None)
        if new is None:
            raise ValueError("The system has no crystal symmetry.")
        if new is system:
            _error(main, "'{}' is already P1 and one cell was asked: nothing to do.".format(system.label))
            return
        if system.nbModel is not None and new.nbModel is None:
            new.DefineNBModel(Clone(system.nbModel))
        p = main.p_session
        name = "{}_P1_{}x{}x{}".format(system.label, na, nb, nc)
        new = p.add_new_system_to_psession(system=new, name=name, tag="P1")
        main.main_treeview.add_new_system_to_treeview(new)
        ff = getattr(new.mmModel, "forceField", "None")
        main.bottom_notebook.status_teeview_add_new_item(
            message="New System:  {} ({}) - Force Field:  {}".format(new.label, new.e_tag, ff), system=new)
        p._add_vismol_object_to_easyhybrid_session(new, True)
        main.vm_session.vm_glcore.queue_draw()
    except Exception as error:
        _error(main, "Could not expand the system:\n{}".format(error))


def GLibEscape(text):
    from gi.repository import GLib
    return GLib.markup_escape_text(str(text))
