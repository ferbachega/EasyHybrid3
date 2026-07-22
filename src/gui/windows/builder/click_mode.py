#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Molecule Builder -- click-to-place-atom edit mode
#
#  Description:
#      Third building block of the "Builder" tool. Lets the user toggle a
#      mode where a plain left click on the 3D view (a click, not a drag --
#      same "was this a drag or a click" distinction the existing picking/
#      selection code already makes) places a new atom at the clicked
#      position instead of selecting whatever is already there. If the
#      click instead lands ON an atom that already belongs to the target
#      object, that atom's ELEMENT is changed in place to the currently
#      selected symbol instead (see handle_click_to_place_atom() /
#      atom_ops.set_atom_element()) -- no duplicate atom stacked on top.
#
#      State (vm_session.builder_atom_mode / .builder_atom_symbol /
#      .builder_target_object) lives directly on vm_session, read
#      defensively via getattr(..., default) everywhere -- deliberately
#      NOT declared in VismolSession.__init__, to keep this feature's
#      footprint out of vismol_session.py entirely (see the hook added to
#      vismol_glcore.py's mouse_released() for the one unavoidable touch
#      point: the click handler itself has to live where clicks already
#      are handled).
#
#      3D position from a 2D click: world_pos_from_mouse() below does
#      PROPER perspective unprojection using the camera's real
#      view_matrix/projection parameters (field_of_view, aspect,
#      dist_cam_zrp) -- NOT the earlier approach (VismolGLCore._mouse_pos(),
#      an orthographic-style approximation always at a fixed z_near
#      depth), which was tried first, reusing an existing utility already
#      used for drag_pos_x/y/z, but which the user confirmed in practice
#      placed atoms away from the actual click position, as expected for
#      a scene rendered in perspective (see resize_window(),
#      mop.my_glPerspectivef) but unprojected as if orthographic. The
#      current math -- pixel -> NDC -> invert the perspective x/y scaling
#      at a chosen depth -> transform by inv(view_matrix) -- was verified
#      numerically offline (recovers an exact known world point, both for
#      a plain and for a rotated+translated synthetic camera) BEFORE
#      being wired into the real code, and again with the real
#      VismolGLCore/GLCamera classes (headless, Xvfb + a manual
#      glcore.initialize() call -- see the Builder session's test notes)
#      -- a click near the viewport centre with the default camera
#      (looking at the origin from (0,0,10), dist_cam_zrp=10) landed the
#      new atom almost exactly at (0,0,0), as expected.
#
#      Depth choice: FIRST tries to read the REAL depth buffer under the
#      cursor (see _read_depth_and_atom_at_pixel() below) -- the same technique
#      PyMOL/Avogadro use for "click to place/pick": a dedicated
#      colour+depth render pass identical in structure to the existing
#      _pick() atom-selection code (same clear, same
#      draw_background_sel_representation() loop, same BACKGROUND_ID
#      white-pixel test for "nothing rendered here"), just ALSO reading
#      GL_DEPTH_COMPONENT at the same pixel and converting it back to a
#      distance-from-camera via the exact inverse of the perspective
#      projection's z-mapping (verified numerically offline first, see
#      _read_depth_and_atom_at_pixel()'s docstring). If the clicked pixel is
#      background (nothing rendered there), falls back to
#      vm_glcore.dist_cam_zrp (the camera's pivot/focus distance) instead
#      -- so clicking ON existing atoms/bonds places the new atom glued
#      to that surface's real depth, and clicking on empty space still
#      gives a sensible default instead of failing.
#
#      NOT implemented yet: choosing the element interactively (symbol is
#      whatever enable_atom_placement_mode() was last called with -- a
#      periodic-table panel or similar is a separate, later step), visual
#      feedback for which mode is currently active, undo.
#
from util.debug import dprint
import numpy as np
from OpenGL import GL
import ctypes


def _current_frame_position ( vismol_object, atom_id ):
    """ [EN] Returns atom_id's position at the CURRENTLY DISPLAYED frame
    (vismol_object.vm_session.frame), NOT always frame 0.

    BUG FIXED (reported by the user: the hover ring/info-text drawn by
    draw_hover_highlight()/draw_hover_info_text() appeared "outside the
    atom" after navigating a real, multi-frame trajectory): several
    functions in this file used to read vismol_object.frames[0, atom_id]
    directly -- harmless for Builder objects specifically (which only
    ever have exactly one frame, so index 0 IS always "the" frame), but
    silently wrong the moment those same functions got reused for
    hovering/interacting with ANY object in the session, including
    real, multi-frame trajectories, where frame 0 is only correct until
    the user navigates away from it.

    Same clamping VismolGLCore._safe_frame_coords() itself uses (clamp
    to the last frame if vm_session.frame has somehow gone out of THIS
    object's own range -- e.g. two objects loaded with different
    trajectory lengths, current frame beyond the shorter one's count). """
    n_frames = vismol_object.frames.shape[0]
    frame_idx = vismol_object.vm_session.frame
    if frame_idx < 0:
        frame_idx = 0
    elif frame_idx >= n_frames:
        frame_idx = n_frames - 1
    return vismol_object.frames[frame_idx, atom_id]


def enable_atom_placement_mode ( vm_session, vismol_object, symbol = "C" ):
    """ Turns on Builder editing mode for vismol_object, starting in the
    "add" tool (plain left clicks on the 3D view place new atoms of
    `symbol` at the clicked position). Call disable_atom_placement_mode()
    to go back to normal click-to-select behaviour.

    [EN] Generalised (was a single on/off flag) to support the 'a'/'d'/
    'b' keyboard shortcuts: vm_session.builder_tool now tracks WHICH
    action a plain click performs while editing is on ('add' or
    'delete') -- see set_tool() below and the keyboard handlers added to
    VismolGTKWidget (_pressed_a / _pressed_d / _pressed_b). """
    vm_session.builder_atom_mode = True
    vm_session.builder_atom_symbol = symbol
    vm_session.builder_target_object = vismol_object
    vm_session.builder_tool = "add"


def disable_atom_placement_mode ( vm_session ):
    """ Restores normal click-to-select behaviour. """
    vm_session.builder_atom_mode = False
    vm_session.builder_target_object = None
    vm_session.builder_tool = "add"


def set_atom_placement_symbol ( vm_session, symbol ):
    """ Changes which element gets placed by the NEXT click, without
    turning the mode off/on again (e.g. switching from "C" to "O"
    mid-session while still in placement mode). """
    vm_session.builder_atom_symbol = symbol


def set_tool ( vm_session, tool ):
    """ [EN] Switches the active Builder tool -- 'add' (plain click
    places a new atom, the default) or 'delete' (plain click removes the
    clicked atom instead). Called by the 'a'/'d' keyboard shortcuts (see
    VismolGTKWidget._pressed_a / _pressed_d). Does nothing (no error) if
    Builder editing isn't currently on -- just sets the tool for whenever
    it next is. """
    if tool not in ( "add", "delete" ):
        raise ValueError ( "set_tool: tool deve ser 'add' ou 'delete', recebido: {}".format ( tool ) )
    vm_session.builder_tool = tool


def handle_bond_shortcut ( vm_session ):
    """ [EN] Called by the 'b' keyboard shortcut (VismolGTKWidget._pressed_b).
    NOT a persistent tool/mode -- a one-shot ACTION: looks at whatever is
    CURRENTLY selected (vm_session.selections[vm_session.current_selection]
    .selected_atoms -- the app's existing, already-accumulating multi-
    select mechanism; a plain click toggles an atom's membership in this
    set without clearing the rest, see
    VismolViewingSelection.selecting_by_atom() in vismol_selections.py),
    and if EXACTLY two atoms belonging to the SAME vismol_object are
    selected, adds a bond between them via atom_ops.add_bond(). Clears
    the selection afterwards on success, so the next two clicks start a
    fresh pair instead of accidentally including a third atom.

    Deliberately requires the SAME object for both atoms -- bonding
    across two different VismolObjects isn't something add_bond() (or
    the rest of this Builder) supports, and silently doing something
    unexpected there would be worse than just refusing with a clear
    message.

    Returns a short status string (meant to be logged/printed by the
    caller), rather than raising, for the common "not exactly 2 atoms
    selected yet" case -- that's an expected, frequent state while the
    user is still clicking atoms to select, not an error. """
    sel = vm_session.selections[vm_session.current_selection]
    selected = list ( sel.selected_atoms )

    if len ( selected ) != 2:
        return "Select exactly 2 atoms before pressing 'b' (currently selected: {}).".format ( len ( selected ) )

    atom_a, atom_b = selected
    if atom_a.vm_object is not atom_b.vm_object:
        return "The 2 selected atoms must belong to the same object."

    from gui.windows.builder.atom_ops import add_bond, adjust_hydrogens, push_undo_snapshot
    push_undo_snapshot ( atom_a.vm_object )
    created = add_bond ( atom_a.vm_object, atom_a.atom_id, atom_b.atom_id )
    adjust_hydrogens ( atom_a.vm_object, atom_a.atom_id )
    adjust_hydrogens ( atom_a.vm_object, atom_b.atom_id )

    from gui.windows.builder.empty_object import sync_pdynamo_system
    sync_pdynamo_system ( atom_a.vm_object )

    sel.selection_function_viewing_set ( None )   # limpa a selecao pro proximo par

    if created:
        return "Bond created between atom {} and atom {}.".format ( atom_a.atom_id, atom_b.atom_id )
    else:
        return "A bond between these 2 atoms already existed."


def handle_click_to_delete_atom ( vm_glcore ):
    """ [EN] Called from render() -- see the hook added there, right next
    to the existing "if self.picking: self._pick()" / new
    "builder_placing_atom" checks -- ONLY when builder_tool == 'delete'
    and a plain (non-shift) click just ran through NORMAL picking
    (reusing _pick()'s own atom identification instead of re-implementing
    it: simpler and safer than a second, parallel colour-picking pass
    just for this). If _pick() found an atom under the cursor
    (vm_glcore.atom_picked), removes it via atom_ops.remove_atom() and
    clears atom_picked (so the normal, non-Builder selection machinery
    downstream in mouse_released() doesn't also try to act on an atom
    that no longer exists). """
    atom = getattr ( vm_glcore, "atom_picked", None )
    if atom is None:
        return None

    from gui.windows.builder.atom_ops import remove_atom
    vismol_object = atom.vm_object
    atom_id = atom.atom_id
    remove_atom ( vismol_object, atom_id )
    vm_glcore.atom_picked = None
    return atom_id


def _read_depth_and_atom_at_pixel ( vm_glcore, mouse_x, mouse_y ):
    """ [EN] Reads the REAL depth buffer under the cursor via a dedicated
    render pass, mirroring VismolGLCore._pick() (same clear, same
    draw_background_sel_representation() loop over active objects/
    representations, same BACKGROUND_ID convention for "nothing here")
    -- just reading GL_DEPTH_COMPONENT in addition to (instead of) the
    RGBA colour ID. This is the standard technique used by PyMOL/
    Avogadro-style "click to place/pick" tools: place new geometry at
    the REAL depth of whatever is under the cursor, not a guessed plane.

    Returns (atom, distance):
      - atom     : the Atom under the cursor (decoded from the picking
                   colour ID via vm_session.atom_dic_id, same lookup
                   VismolGLCore._pick() itself uses), or None if the
                   pixel is background.
      - distance : positive distance-from-camera (float), or None if
                   the pixel is background -- caller should fall back
                   to a default depth (e.g. dist_cam_zrp) in that case.

    [EN] Extended (was depth-only, called _read_depth_and_atom_at_pixel) so the
    Builder's "add" tool can tell whether a click landed ON an existing
    atom (to replace its element in place -- see
    handle_click_to_place_atom() / atom_ops.set_atom_element()) using
    the SAME render pass already being done for the depth read, instead
    of a second, separate picking pass -- the colour ID was already
    being computed here just to test against BACKGROUND_ID and then
    thrown away; decoding it into an actual Atom is the same
    pickedID -> vm_session.atom_dic_id[pickedID] lookup _pick() does,
    just reusing this pass's own pixel read instead of triggering a new
    one.

    Depth-buffer-value -> distance-from-camera conversion verified
    numerically offline first (a known distance run through the forward
    projection math, then back through this inverse, recovered the
    original value to 5+ significant figures) before being wired in
    here. """
    BACKGROUND_ID = 16777215

    GL.glClearColor ( 1, 1, 1, 1 )
    GL.glClear ( GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT )
    vm_glcore.update_camera_ubo ( )

    for vm_object in vm_glcore.vm_session.vm_objects_dic.values ( ):
        if not vm_object.active:
            continue
        for rep in vm_object.representations.values ( ):
            if rep and rep.active:
                rep.draw_background_sel_representation ( )

    GL.glPixelStorei ( GL.GL_PACK_ALIGNMENT, 1 )

    x = int ( mouse_x )
    y = int ( vm_glcore.height - mouse_y )   # GTK (top-left) -> OpenGL (bottom-left)

    color_data = GL.glReadPixels ( x, y, 1, 1, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE )
    pickedID = color_data[0] + color_data[1] * 256 + color_data[2] * 256 * 256

    if pickedID == BACKGROUND_ID:
        #print ( "DEBUG click_mode: depth buffer at clicked pixel = background (nothing rendered there)" )
        return None, None

    atom = vm_glcore.vm_session.atom_dic_id.get ( pickedID )

    depth_raw = GL.glReadPixels ( x, y, 1, 1, GL.GL_DEPTH_COMPONENT, GL.GL_FLOAT )
    depth_buffer_value = float ( np.frombuffer ( depth_raw, dtype = np.float32 ) [0] )

    proj = vm_glcore.glcamera.projection_matrix
    p22, p32 = float ( proj[2,2] ), float ( proj[3,2] )
    ndc_z = 2.0 * depth_buffer_value - 1.0
    distance = p32 / ( ndc_z + p22 )

    dprint ( "DEBUG click_mode: depth buffer at clicked pixel = {:.5f}  ->  distance from camera = {:.3f}  ->  atom = {}".format (
            depth_buffer_value, distance, ( "#{} ({})".format ( atom.atom_id, atom.symbol ) if atom is not None else None ) ) )
    return atom, float ( distance )


def world_pos_from_mouse ( vm_glcore, mouse_x, mouse_y, depth = None ):
    """ [EN] PROPER perspective unprojection -- replaces the earlier
    approach of reusing VismolGLCore._mouse_pos() (an orthographic-style
    approximation at a fixed z_near depth), after the user reported
    atoms landing away from the actual click position.

    Converts a 2D pixel coordinate (mouse_x, mouse_y, GTK convention:
    origin top-left) into a 3D WORLD position, at a chosen distance
    ("depth") from the camera along the view direction. If `depth` isn't
    given explicitly, tries _read_depth_and_atom_at_pixel() first (the REAL depth
    under the cursor, PyMOL/Avogadro-style), falling back to
    vm_glcore.dist_cam_zrp (the camera's current pivot/focus distance)
    only if that pixel is background.

    Math (verified numerically, offline, against both a plain and a
    rotated+translated synthetic camera before touching the real code --
    recovers the exact original world point in both cases):
      1. pixel -> NDC:  ndc_x = 2*mouse_x/width  - 1
                        ndc_y = 1 - 2*mouse_y/height   (Y flip: GTK is
                        top-left origin, NDC/OpenGL is bottom-left)
      2. Invert the perspective projection's x/y scaling (see
         matrix_operations.pyx: my_glPerspectivef -- f = 1/tan(fovy in
         radians), matrix entries P[0,0]=f/aspect, P[1,1]=f) to get the
         VIEW-space x/y at the chosen depth:
           view_x = ndc_x * depth * aspect / f
           view_y = ndc_y * depth / f
           view_z = -depth   (camera looks down -Z in view space)
      3. Transform the view-space point back to WORLD space using the
         INVERSE of glcamera.view_matrix. This codebase stores the
         translation component in the LAST ROW of its 4x4 matrices (see
         matrix_operations.pyx: my_glTranslatef, get_xyz_coords) --  a
         row-vector convention, i.e. `view_row = world_row @ view_matrix`
         -- so the inverse relation is `world_row = view_row @
         inv(view_matrix)`, NOT `inv(view_matrix) @ view_row`.

    Prints the pixel/NDC/world coordinates every time it's called -- the
    user asked to see exactly what position is being computed, to check
    it against where they actually clicked. """
    width  = float ( vm_glcore.width )
    height = float ( vm_glcore.height )

    if depth is None:
        _atom_unused, depth = _read_depth_and_atom_at_pixel ( vm_glcore, mouse_x, mouse_y )
        if depth is None:
            depth = float ( abs ( vm_glcore.dist_cam_zrp ) )

    ndc_x = 2.0 * mouse_x / width  - 1.0
    ndc_y = 1.0 - 2.0 * mouse_y / height

    fovy   = float ( vm_glcore.glcamera.field_of_view )
    aspect = float ( vm_glcore.glcamera.viewport_aspect_ratio )
    f = 1.0 / np.tan ( fovy * np.pi / 180.0 )   # MESMA formula de my_glPerspectivef (matrix_operations.pyx)

    view_x = ndc_x * depth * aspect / f
    view_y = ndc_y * depth / f
    view_z = -depth

    view_point = np.array ( [ view_x, view_y, view_z, 1.0 ], dtype = np.float32 )
    inv_view   = np.linalg.inv ( vm_glcore.glcamera.view_matrix )
    world_point = view_point @ inv_view

    dprint ( "DEBUG click_mode: mouse=({:.1f}, {:.1f})  viewport=({:.0f}x{:.0f})  ndc=({:.3f}, {:.3f})".format (
            mouse_x, mouse_y, width, height, ndc_x, ndc_y ) )
    dprint ( "DEBUG click_mode: fovy={:.2f} aspect={:.3f} depth_used={:.3f}".format ( fovy, aspect, depth ) )
    dprint ( "DEBUG click_mode: view_point=({:.3f}, {:.3f}, {:.3f})".format ( view_x, view_y, view_z ) )
    dprint ( "DEBUG click_mode: world_point=({:.3f}, {:.3f}, {:.3f})  <- new atom position".format (
            world_point[0], world_point[1], world_point[2] ) )

    return float ( world_point[0] ), float ( world_point[1] ), float ( world_point[2] )


def handle_click_to_place_atom ( vm_glcore, mouse_x, mouse_y ):
    """ Called from VismolGLCore.mouse_released() -- see the hook added
    there -- ONLY when builder_atom_mode is active and the mouse event
    was a genuine click (not a drag, not a rotate/pan/zoom gesture).

    Two behaviours, decided by what's under the cursor (read ONCE via
    _read_depth_and_atom_at_pixel(), reused for whichever branch below
    is taken -- no second render pass):

      1. Click lands ON an atom that already belongs to THIS Builder's
         target object: that atom's element is changed in place to the
         currently selected symbol (atom_ops.set_atom_element()) --
         same atom_id, same position, no duplicate atom stacked on top.
         A no-op if the atom already has the selected symbol (nothing
         to replace). [EN] Deliberately NOT done for an atom belonging
         to a DIFFERENT vismol_object (not the current Builder target)
         -- falls through to case 2 instead, same reasoning as
         add_bond()'s existing same-object-only restriction elsewhere
         in the Builder: silently mutating a different object the user
         didn't ask to edit would be worse than the (harmless) fallback
         of placing a new atom at that position in the CORRECT object.

      2. Click lands on empty space (or on a different object): computes
         a 3D position from the 2D click (world_pos_from_mouse(), proper
         perspective unprojection -- see its own docstring), converts it
         into the TARGET OBJECT's own local/model space (see below --
         this step was missing at first, and caused new atoms to land in
         the wrong place as soon as the model had been rotated), and
         adds one atom there via atom_ops.add_atom().

    Returns the affected Atom (new or replaced), or None if the mode
    isn't actually usable right now (no target object set).

    [EN] model_mat correction (bug fixed after the user reported atoms
    landing correctly at first, but in the wrong place after rotating
    the model): this codebase implements "rotate the view" by rotating
    EVERY vismol_object's OWN model_mat (see VismolGLCore._rotate_view())
    -- the camera's view_matrix itself never changes. world_pos_from_mouse()
    only inverts view_matrix/projection_matrix, recovering a position in
    the space view_matrix operates in -- call it "world space". But the
    coordinates actually stored in vismol_object.frames (what add_atom()
    writes to) live in that object's OWN LOCAL space, one more transform
    away: rendered_position = local_row @ vismol_object.model_mat @
    view_matrix @ projection_matrix (row-vector convention, matching
    vertex_shader_surface's "view_mat * model_mat * vec4(vert_coord,1.0)"
    written the other way around in GLSL's column-vector notation).
    So converting a clicked "world space" point into the LOCAL coordinate
    that needs to be stored requires ONE more inverse step:
    local = world_row @ inverse(vismol_object.model_mat). Verified
    numerically offline first (a known local point, forward-projected
    through a ROTATED model_mat + fixed camera, then recovered exactly
    through this same inverse chain) before adding it here -- without
    this step, the recovered point was off by exactly the model_mat
    rotation (matching the reported symptom precisely: correct when
    model_mat is still identity/unrotated, wrong after rotating). """
    vm_session = vm_glcore.vm_session
    vismol_object = getattr ( vm_session, "builder_target_object", None )
    if vismol_object is None:
        return None

    symbol = getattr ( vm_session, "builder_atom_symbol", "C" )

    picked_atom, depth = _read_depth_and_atom_at_pixel ( vm_glcore, mouse_x, mouse_y )
    if depth is None:
        depth = float ( abs ( vm_glcore.dist_cam_zrp ) )

    if picked_atom is not None:
        dprint ( "DEBUG click_mode: an atom WAS clicked -- atom #{} ({}), object '{}'".format (
                picked_atom.atom_id, picked_atom.symbol, picked_atom.vm_object.name ) )
    
        dprint(picked_atom.atom_id, picked_atom.bonds )
    else:
        dprint ( "DEBUG click_mode: no atom was clicked (empty space)." )

    if picked_atom is not None and picked_atom.vm_object is vismol_object:
        if picked_atom.symbol == symbol:
            dprint ( "DEBUG click_mode: clicked atom #{} is already '{}' -- nothing to replace.".format (
                    picked_atom.atom_id, symbol ) )
            return picked_atom

        from gui.windows.builder.atom_ops import set_atom_element, push_undo_snapshot, adjust_hydrogens
        push_undo_snapshot ( vismol_object )
        atom = set_atom_element ( vismol_object, picked_atom.atom_id, symbol )
        # [EN] The new element's standard valence is (almost always)
        # different from the old one's -- adjust THIS atom's own
        # hydrogens to match (its neighbours' bond orders to it are
        # unchanged, so THEY don't need adjusting, only this atom does).
        adjust_hydrogens ( vismol_object, atom.atom_id )
        from gui.windows.builder.empty_object import sync_pdynamo_system
        sync_pdynamo_system ( vismol_object )
        dprint ( "DEBUG click_mode: replaced atom #{} -> '{}'".format ( atom.atom_id, symbol ) )
        return atom

    wx, wy, wz = world_pos_from_mouse ( vm_glcore, mouse_x, mouse_y, depth = depth )

    world_point = np.array ( [ wx, wy, wz, 1.0 ], dtype = np.float32 )
    inv_model = np.linalg.inv ( vismol_object.model_mat )
    local_point = world_point @ inv_model
    x, y, z = float ( local_point[0] ), float ( local_point[1] ), float ( local_point[2] )

    dprint ( "DEBUG click_mode: world_point=({:.3f}, {:.3f}, {:.3f})  -> local_point (after inv(model_mat))=({:.3f}, {:.3f}, {:.3f})".format (
            wx, wy, wz, x, y, z ) )

    from gui.windows.builder.atom_ops import add_atom, push_undo_snapshot
    push_undo_snapshot ( vismol_object )
    atom = add_atom ( vismol_object, symbol = symbol, x = x, y = y, z = z )

    # [EN] DESIGN CHANGE: this used to have its own ad-hoc hydrogenation
    # logic here (a `tmp` dict of fixed C/N/O direction templates),
    # ONLY ever applied to a brand-new atom placed in empty space (never
    # bonded to anything else). Replaced with a call to the new, general
    # atom_ops.adjust_hydrogens() (see its own docstring), which works the
    # same way here (atom has zero bonds yet -> needed_h == full target
    # valence -> same fixed templates, same result) but ALSO now runs
    # from every other place an atom's bonding can change (replace
    # element, bond-drag finish, cycle bond order, delete atom/bond --
    # see each of those call sites' own comment), instead of being a
    # one-off special case just for this one interaction.
    from gui.windows.builder.atom_ops import adjust_hydrogens
    adjust_hydrogens ( vismol_object, atom.atom_id )

    from gui.windows.builder.empty_object import sync_pdynamo_system
    sync_pdynamo_system ( vismol_object )

    return atom


# =====================================================================================
#   Click-and-drag to create a bonded atom
#   ------------------------------------------------------------------------------
#   Third click interaction for the "add" tool (alongside "click on empty
#   space -> new atom" and "click on an existing atom -> replace its
#   element", both above): press-and-HOLD on an existing atom, then drag
#   without releasing -- creates a new atom, already bonded to the one
#   pressed, that follows the cursor live while the button stays down.
#   Releasing the button drops the new atom at its current position and
#   finalises the bond. A plain click (press+release with no real drag)
#   on an atom still means "replace", exactly as before -- see the
#   mouse_pressed()/mouse_motion()/mouse_released() hooks in
#   vismol_glcore.py for how the two interactions are told apart.
#
#   State lives on vm_session (same getattr(..., default)-everywhere,
#   nothing declared in VismolSession.__init__ convention already used
#   for builder_atom_mode/builder_target_object above):
#     builder_bond_drag_active       : bool, is a drag currently happening
#     builder_bond_drag_origin_atom  : the Atom the drag started FROM
#     builder_bond_drag_new_atom     : the Atom being dragged
#     builder_bond_drag_object       : the vismol_object both belong to
#     builder_bond_drag_depth        : FIXED distance-from-camera the new
#                                       atom is kept at for the whole drag
#                                       (see start_bond_drag()) -- picked
#                                       once, at the origin atom's own
#                                       depth, rather than re-reading
#                                       whatever's under the cursor each
#                                       motion event (which would jitter
#                                       every time the cursor crosses over
#                                       a DIFFERENT atom/bond mid-drag).
# =====================================================================================

def start_bond_drag ( vm_glcore, origin_atom, depth ):
    """ Begins a click-and-drag-to-create-a-bonded-atom interaction.
    Called from mouse_motion() (see the hook added there) the FIRST time
    real mouse movement is detected after a press that landed on an atom
    belonging to the current Builder target object -- render() only
    records that atom as a CANDIDATE (see the hook added there, next to
    builder_placing_atom); the drag doesn't actually start until this
    function runs, precisely so that a plain click (press+release, no
    real movement) never reaches here at all and still falls through to
    the "replace element" interaction instead (handle_click_to_place_atom()).

    Creates the new atom immediately, AT THE SAME POSITION as
    origin_atom (distance zero -- it hasn't been dragged anywhere yet),
    bonded to origin_atom via add_atom()'s own bonded_to parameter --
    an EXPLICIT bond, which is now the ONLY kind that exists (distance-
    based auto-detection is off entirely for the Builder -- see atom_ops.
    add_atom()'s docstring). From here, update_bond_drag() repositions
    this same new atom every motion event, and finish_bond_drag()
    re-confirms the bond once the drag ends (see its own docstring for
    why that re-confirmation still matters even though the bond was
    already explicit from the very start).

    `depth` is the REAL distance-from-camera of the pressed pixel
    (already computed by the caller via _read_depth_and_atom_at_pixel()
    -- not re-read here) -- stored on vm_session so update_bond_drag()
    can keep the dragged atom at that SAME depth for the whole gesture
    (see the module-level note above for why a fixed depth, not a
    per-motion-event re-read, is what gives smooth/predictable
    dragging).

    Returns the new Atom. """
    vm_session    = vm_glcore.vm_session
    vismol_object = origin_atom.vm_object
    symbol        = getattr ( vm_session, "builder_atom_symbol", "C" )

    ox, oy, oz = [ float ( c ) for c in _current_frame_position ( vismol_object, origin_atom.atom_id ) ]

    from gui.windows.builder.atom_ops import add_atom, push_undo_snapshot
    push_undo_snapshot ( vismol_object )
    new_atom = add_atom ( vismol_object, symbol = symbol, x = ox, y = oy, z = oz,
                           bonded_to = origin_atom.atom_id )

    vm_session.builder_bond_drag_active      = True
    vm_session.builder_bond_drag_origin_atom = origin_atom
    vm_session.builder_bond_drag_new_atom    = new_atom
    vm_session.builder_bond_drag_object      = vismol_object
    vm_session.builder_bond_drag_depth       = depth

    dprint ( "DEBUG click_mode: bond-drag started -- new atom #{} ('{}') bonded to atom #{} ('{}')".format (
            new_atom.atom_id, symbol, origin_atom.atom_id, origin_atom.symbol ) )

    vm_glcore.updated_coords = True
    vm_glcore.queue_draw ( )
    return new_atom


def update_bond_drag ( vm_glcore, mouse_x, mouse_y ):
    """ Called from mouse_motion() on every motion event while
    vm_session.builder_bond_drag_active is True. Repositions the
    dragged atom (atom_ops.move_atom() -- cheap, no bond recompute) to
    the current cursor position, unprojected at the FIXED depth captured
    by start_bond_drag() (see module-level note above).

    Pure math -- world_pos_from_mouse() only touches the GPU when its
    `depth` argument is None, and here it never is -- safe to call
    directly from mouse_motion() (a plain GTK handler, NOT render()),
    unlike handle_click_to_place_atom() which needs a real depth-buffer
    read and must stay deferred to render(). """
    vm_session = vm_glcore.vm_session
    if not getattr ( vm_session, "builder_bond_drag_active", False ):
        return None

    vismol_object = vm_session.builder_bond_drag_object
    new_atom      = vm_session.builder_bond_drag_new_atom
    depth         = vm_session.builder_bond_drag_depth

    wx, wy, wz = world_pos_from_mouse ( vm_glcore, mouse_x, mouse_y, depth = depth )

    world_point = np.array ( [ wx, wy, wz, 1.0 ], dtype = np.float32 )
    inv_model   = np.linalg.inv ( vismol_object.model_mat )
    local_point = world_point @ inv_model
    x, y, z = float ( local_point[0] ), float ( local_point[1] ), float ( local_point[2] )

    from gui.windows.builder.atom_ops import move_atom
    move_atom ( vismol_object, new_atom.atom_id, x, y, z )
    return new_atom


def _find_bond_snap_target ( vismol_object, dragged_atom, exclude_ids, tolerance = 1.3 ):
    """ [EN] Looks for an existing atom close enough to `dragged_atom`'s
    CURRENT position to "snap" onto -- used by finish_bond_drag() (see
    its own docstring) to decide whether dropping the dragged atom near/
    onto another EXISTING atom should connect the two directly instead
    of keeping the dragged (temporary) atom.

    Deliberately a PURE 3D-DISTANCE check (position vs position), NOT a
    screen-space click/pick: a pixel/depth-based pick at the cursor's
    CURRENT position would almost always just hit the dragged atom
    ITSELF (it's rendered exactly at the cursor, being actively dragged
    there every motion event) -- so proximity in the object's own
    coordinate space is what actually detects "dropped onto atom X"
    here, not what pixel colour is on top.

    Threshold is the same kind of covalent-radius-sum heuristic the old
    (now-removed, see atom_ops.add_atom()'s docstring) distance-based
    auto-detector used, just computed here on demand for this one
    comparison rather than for the whole object: two atoms within
    (cov_rad_a + cov_rad_b) * tolerance of each other are considered
    "aimed at the same spot", not two unrelated atoms that merely ended
    up somewhat close.

    Returns the closest matching Atom (excluding anything in
    `exclude_ids`, normally {origin_atom.atom_id, dragged_atom.atom_id}),
    or None if nothing is within range. """
    position = _current_frame_position ( vismol_object, dragged_atom.atom_id )

    best_atom = None
    best_dist = None
    for candidate_id, candidate in vismol_object.atoms.items ( ):
        if candidate_id in exclude_ids:
            continue
        candidate_pos = _current_frame_position ( vismol_object, candidate_id )
        dist = float ( np.linalg.norm ( candidate_pos - position ) )
        threshold = ( dragged_atom.cov_rad + candidate.cov_rad ) * tolerance
        if dist <= threshold and ( best_dist is None or dist < best_dist ):
            best_atom = candidate
            best_dist = dist

    return best_atom


def finish_bond_drag ( vm_glcore ):
    """ Called from mouse_released() (checked FIRST, before anything
    else -- see the hook added there) once the button comes back up
    while a bond-drag is active. Finalises the dragged atom at its
    current (already up to date, from the last update_bond_drag() call)
    position.

    NEW: if that final position is close enough to an EXISTING atom
    (other than the one the drag started from) -- see
    _find_bond_snap_target() -- the temporary dragged atom is REMOVED,
    and the bond is made directly between the origin atom and that
    existing atom instead. This is what makes "drag from atom A, aim at
    already-existing atom B, release" connect A and B directly rather
    than leaving a new, redundant atom sitting on top of B.

    [EN] DESIGN CHANGE: this used to re-run distance-based bond
    detection from scratch here (vismol_object.find_bonded_and_
    nonbonded_atoms()) to pick up any OTHER bond the dropped atom might
    have landed close enough to, then re-add the origin<->new-atom bond
    explicitly afterwards since that recompute could silently drop it.
    Distance-based auto-detection is now turned OFF entirely for the
    Builder (see atom_ops.add_atom()'s docstring for why) -- the
    "landed close to another atom" case is now handled explicitly and
    deliberately above (_find_bond_snap_target()), not as a side effect
    of a general-purpose distance recompute; and the origin<->new-atom
    bond, when the dragged atom IS kept, was already explicit from the
    moment start_bond_drag() created it (via add_atom()'s bonded_to
    parameter), so it was never at risk of being dropped in the first
    place.

    Clears all builder_bond_drag_* state and turns vm_glcore.updated_coords
    back off (start_bond_drag()/update_bond_drag() need it on for the
    live-follow effect during the drag; leaving it on afterwards would
    just mean every representation's coordinates get needlessly re-
    uploaded on every future frame, forever -- see the note in
    atom_ops.move_atom()).

    Returns the finalised Atom (the ORIGINAL dragged atom, or the
    existing atom it snapped onto if the dragged one was removed), or
    None if no drag was active. """
    vm_session = vm_glcore.vm_session
    if not getattr ( vm_session, "builder_bond_drag_active", False ):
        return None

    vismol_object = vm_session.builder_bond_drag_object
    origin_atom   = vm_session.builder_bond_drag_origin_atom
    new_atom      = vm_session.builder_bond_drag_new_atom

    from gui.windows.builder.atom_ops import add_bond, remove_atom, adjust_hydrogens

    snap_target = _find_bond_snap_target (
        vismol_object, new_atom,
        exclude_ids = { origin_atom.atom_id, new_atom.atom_id }
    )

    if snap_target is not None:
        # [EN] Dropped onto an existing atom -- connect origin directly
        # to IT, and throw away the temporary dragged atom. Capture both
        # ids BEFORE removing new_atom: new_atom is always the MOST
        # RECENTLY created atom in this object (nothing else has been
        # added since start_bond_drag() created it, only moved via
        # move_atom()), so its atom_id is the HIGHEST one currently in
        # use -- removing it therefore can NEVER renumber origin_atom's
        # or snap_target's ids out from under us (remove_atom() only
        # shifts ids ABOVE the removed one down by one).
        origin_id = origin_atom.atom_id
        target_id = snap_target.atom_id
        remove_atom ( vismol_object, new_atom.atom_id )
        add_bond ( vismol_object, origin_id, target_id )

        dprint ( "DEBUG click_mode: bond-drag finished -- snapped onto existing atom #{} ('{}'), connected to atom #{} ('{}'); temporary dragged atom removed".format (
                target_id, snap_target.symbol, origin_id, origin_atom.symbol ) )

        finalised_atom = vismol_object.atoms[target_id]
    else:
        add_bond ( vismol_object, origin_atom.atom_id, new_atom.atom_id )

        dprint ( "DEBUG click_mode: bond-drag finished -- atom #{} ('{}') dropped, bonded to atom #{} ('{}')".format (
                new_atom.atom_id, new_atom.symbol, origin_atom.atom_id, origin_atom.symbol ) )

        finalised_atom = new_atom

    # [EN] Adjust BOTH atoms' hydrogens now that their bonding changed --
    # origin_atom just gained a bond (may now have too MANY H's to still
    # be correct), and finalised_atom either just got created (needs its
    # full set of H's) or -- in the snap case -- also just gained a bond
    # (same "too many H's" possibility). Always re-read .atom_id from the
    # live Atom OBJECT here (not a cached int): remove_atom() -- possibly
    # triggered by the FIRST adjust_hydrogens() call below, if it happens
    # to remove some hydrogens -- mutates every SURVIVING atom's
    # .atom_id IN PLACE, so reading it fresh off the object is always
    # correct regardless of what the first call did, even though a
    # lower-numbered atom may have been removed in between the two calls.
    adjust_hydrogens ( vismol_object, origin_atom.atom_id )
    adjust_hydrogens ( vismol_object, finalised_atom.atom_id )

    from gui.windows.builder.empty_object import sync_pdynamo_system
    sync_pdynamo_system ( vismol_object )

    vm_session.builder_bond_drag_active      = False
    vm_session.builder_bond_drag_origin_atom = None
    vm_session.builder_bond_drag_new_atom    = None
    vm_session.builder_bond_drag_object      = None
    vm_session.builder_bond_drag_depth       = None
    vm_session.builder_press_candidate_atom  = None
    vm_session.builder_press_candidate_depth = None

    vm_glcore.updated_coords = False
    vm_glcore.queue_draw ( )
    return finalised_atom


# =====================================================================================
#   Ctrl+click on a bond -- cycle its order (single -> double -> triple -> single)
#   ------------------------------------------------------------------------------
#   Fourth click interaction for the "add" tool (alongside "click on empty
#   space -> new atom", "click on an existing atom -> replace its
#   element", and "click+drag from an atom -> new bonded atom", all
#   above): Ctrl+click (no drag) on an existing BOND cycles that bond's
#   order. Deliberately scoped to ONLY the current Builder target object
#   (vm_session.builder_target_object) -- by design, only one object is
#   ever being edited at a time, so there's no reason for this to affect,
#   or even look at, any other object in the session.
#
#   DESIGN CHOICE: bond "picking" here is a plain 2D-projection distance
#   check (project each candidate bond's two endpoints to screen pixels,
#   see _project_local_point_to_pixel(), then find the closest bond's
#   on-screen line segment to the click), NOT a GPU colour-ID pick like
#   atom picking (vm_session.atom_dic_id / VismolGLCore._pick()). Building
#   a parallel colour-ID system for bonds would need a genuinely separate
#   VAO/VBO (each bond needs its OWN 2 vertices with a UNIFORM colour --
#   the existing "lines" representation shares vertex data with atoms,
#   colouring each line endpoint by its OWN atom's colour, which can't
#   also encode "this bond" since one atom can be an endpoint of several
#   different bonds at once) and a brand new GPU render pass, which isn't
#   feasible to get right without a live GL context to test against. The
#   projection-distance approach only needs plain Python/numpy -- the
#   exact same matrices world_pos_from_mouse() already uses, just applied
#   in the FORWARD direction (3D -> 2D) instead of backward (2D -> 3D).
# =====================================================================================

def _project_local_point_to_pixel ( vm_glcore, vismol_object, local_xyz ):
    """ Forward projection: a LOCAL 3D point (vismol_object's own
    coordinate space, i.e. straight out of vismol_object.frames) -> 2D
    pixel coordinates (GTK convention, origin top-left). The exact
    inverse chain of world_pos_from_mouse()'s own unprojection math (see
    that function's docstring for the verified formulas) plus the
    model_mat step documented in handle_click_to_place_atom()'s
    docstring, run FORWARDS instead of backwards.

    Returns (pixel_x, pixel_y, view_z). view_z is returned too so the
    caller can discard points BEHIND the camera (view_z >= 0, i.e. depth
    <= 0) instead of projecting them to a nonsensical pixel position --
    in that case pixel_x/pixel_y come back as None. """
    width  = float ( vm_glcore.width )
    height = float ( vm_glcore.height )

    local_point = np.array ( [ local_xyz[0], local_xyz[1], local_xyz[2], 1.0 ], dtype = np.float32 )
    world_point = local_point @ vismol_object.model_mat
    view_point  = world_point @ vm_glcore.glcamera.view_matrix

    view_z = float ( view_point[2] )
    depth = -view_z   # camera olha para -Z no espaco de view (ver docstring de world_pos_from_mouse)

    if depth <= 1e-6:
        return None, None, view_z   # atras da camera -- nao ha pixel sensato

    fovy   = float ( vm_glcore.glcamera.field_of_view )
    aspect = float ( vm_glcore.glcamera.viewport_aspect_ratio )
    f = 1.0 / np.tan ( fovy * np.pi / 180.0 )   # mesma formula de my_glPerspectivef (matrix_operations.pyx)

    ndc_x = float ( view_point[0] ) * f / ( aspect * depth )
    ndc_y = float ( view_point[1] ) * f / depth

    pixel_x = ( ndc_x + 1.0 ) * width  / 2.0
    pixel_y = ( 1.0 - ndc_y ) * height / 2.0

    return pixel_x, pixel_y, view_z


def _point_to_segment_distance_2d ( px, py, ax, ay, bx, by ):
    """ Standard 2D point-to-line-SEGMENT distance: perpendicular
    distance if the closest point falls within the segment, distance to
    the nearest endpoint otherwise. Used to measure how close a click
    is to a bond's projected on-screen line. """
    abx, aby = bx - ax, by - ay
    seg_len_sq = abx * abx + aby * aby
    if seg_len_sq < 1e-9:
        return float ( np.hypot ( px - ax, py - ay ) )   # bond endpoints coincide on-screen -- treat as a point

    t = ( ( px - ax ) * abx + ( py - ay ) * aby ) / seg_len_sq
    t = max ( 0.0, min ( 1.0, t ) )
    closest_x = ax + t * abx
    closest_y = ay + t * aby
    return float ( np.hypot ( px - closest_x, py - closest_y ) )


def find_bond_at_pixel ( vm_glcore, vismol_object, mouse_x, mouse_y, pixel_threshold = 10.0 ):
    """ Finds the bond of `vismol_object` whose on-screen line segment is
    closest to (mouse_x, mouse_y), within pixel_threshold screen pixels.
    See the module-level note above for why this is a plain 2D-projection
    distance check rather than a GPU colour-ID pick.

    Only ever called with vismol_object = vm_session.builder_target_object
    (see cycle_bond_order()'s caller in vismol_glcore.py) -- bonds of any
    OTHER object are never even considered, by construction, matching the
    "only one object is editable at a time" design.

    Returns the closest Bond within range, or None. """
    best_bond = None
    best_dist = None

    for bond in vismol_object.bonds.values ( ):
        pos_i = _current_frame_position ( vismol_object, bond.atom_index_i )
        pos_j = _current_frame_position ( vismol_object, bond.atom_index_j )

        px_i, py_i, _z_i = _project_local_point_to_pixel ( vm_glcore, vismol_object, pos_i )
        px_j, py_j, _z_j = _project_local_point_to_pixel ( vm_glcore, vismol_object, pos_j )

        if px_i is None or px_j is None:
            continue   # uma das pontas esta atras da camera -- ignora esse bond

        dist = _point_to_segment_distance_2d ( mouse_x, mouse_y, px_i, py_i, px_j, py_j )
        if dist <= pixel_threshold and ( best_dist is None or dist < best_dist ):
            best_bond = bond
            best_dist = dist

    return best_bond


def find_atom_at_pixel_2d ( vm_glcore, vismol_object, mouse_x, mouse_y, pixel_threshold = 12.0 ):
    """ [EN] Pure-CPU hover test: projects every atom of `vismol_object`
    to a screen pixel (_project_local_point_to_pixel() -- no GPU calls
    at all, same technique find_bond_at_pixel() already uses for bonds)
    and returns whichever one is closest to (mouse_x, mouse_y), within
    pixel_threshold pixels.

    Deliberately NOT a GPU colour-ID pick (glReadPixels forces the GPU
    to finish and sync with the CPU before returning -- a real,
    well-documented stall, not just "one more call" -- see the
    conversation this was designed in for the full reasoning): this
    needs to be cheap enough to run on EVERY mouse_motion event, even
    pure hovering with no button held (POINTER_MOTION_MASK is enabled on
    this widget -- see vismol_gtkwidget.py -- so motion events fire
    constantly while the cursor is over the view, not just while
    dragging). Projecting N atoms is a handful of matrix multiplies,
    still far cheaper than one synchronous GPU readback, regardless of
    whether the camera happens to be moving.

    Only checks THIS ONE object -- see find_atom_at_pixel_2d_any_object()
    below for the general, "hover works everywhere" version actually
    wired into mouse_motion() now; this single-object version is kept
    around since it's still what the Builder-specific bond/atom-drag
    features (start_bond_drag() and friends) conceptually only ever
    needed one object for.

    Returns the closest Atom within range, or None. """
    best_atom = None
    best_dist = None

    for atom_id, atom in vismol_object.atoms.items ( ):
        pos = _current_frame_position ( vismol_object, atom_id )
        px, py, _z = _project_local_point_to_pixel ( vm_glcore, vismol_object, pos )
        if px is None:
            continue   # atras da camera

        dist = float ( np.hypot ( mouse_x - px, mouse_y - py ) )
        if dist <= pixel_threshold and ( best_dist is None or dist < best_dist ):
            best_atom = atom
            best_dist = dist

    return best_atom


def find_atom_at_pixel_2d_any_object ( vm_glcore, mouse_x, mouse_y, pixel_threshold = 12.0 ):
    """ [EN] Same idea as find_atom_at_pixel_2d() above (pure-CPU,
    no-GPU-calls screen-space hover test), but searches across EVERY
    active VismolObject in the session (vm_session.vm_objects_dic),
    not just one -- this is what actually makes hovering work "a
    qualquer momento" (over any loaded molecule, any time), not only
    while editing in the Builder.

    PERFORMANCE NOTE: this is a plain Python for-loop over every atom of
    every active object, run on every single mouse_motion event
    (including pure hover, no button held). For a Builder-sized molecule
    (a handful to a few dozen atoms) this is trivial. For a large,
    normally-loaded system (a protein with thousands of atoms), this
    could start to add up, unlike the GPU colour-ID pick it replaces
    (which is O(1) per pixel regardless of atom count -- its cost comes
    entirely from the CPU/GPU sync, not from iterating atoms). If hover
    ever feels laggy on a large system, the fix is to VECTORISE this
    with numpy -- project every atom's position through the view/
    projection matrices in one batched matrix multiply (instead of
    _project_local_point_to_pixel()'s current one-atom-at-a-time Python
    loop) and use a single vectorised distance comparison -- rather than
    reintroducing a GPU readback. Not done here to keep this change
    small and match exactly what was asked; flagged for later if it
    turns out to matter in practice.

    Returns the closest Atom within range across ALL active objects, or
    None. """
    best_atom = None
    best_dist = None

    for vm_object in vm_glcore.vm_session.vm_objects_dic.values ( ):
        if not vm_object.active:
            continue
        for atom_id, atom in vm_object.atoms.items ( ):
            pos = _current_frame_position ( vm_object, atom_id )
            px, py, _z = _project_local_point_to_pixel ( vm_glcore, vm_object, pos )
            if px is None:
                continue   # atras da camera

            dist = float ( np.hypot ( mouse_x - px, mouse_y - py ) )
            if dist <= pixel_threshold and ( best_dist is None or dist < best_dist ):
                best_atom = atom
                best_dist = dist

    return best_atom


def cycle_bond_order ( vm_glcore, vismol_object, bond ):
    """ Cycles bond.bond_order: 1 -> 2 -> 3 -> 1 (single -> double ->
    triple -> single). Triggered by Ctrl+click on an existing bond (see
    the mouse_released() hook in vismol_glcore.py).

    Persists the new order in vismol_object.manual_bond_orders (keyed by
    the normalized (min,max) atom-id pair), which atom_ops.
    _reapply_manual_bonds() now feeds into vismol_object.
    _bonds_from_pair_of_indexes_list() as `external_orders` -- REQUIRED,
    not optional: bonds get recreated FROM SCRATCH (fresh Bond() objects)
    every time ANYTHING changes on this object (add_atom(), remove_atom(),
    add_bond()...), so without persisting it somewhere durable, cycling a
    bond's order here would get silently overwritten back to the default
    the very next time anything else is edited -- the exact same class of
    bug manual_bonds itself had to be fixed for (see
    _reapply_manual_bonds()'s own docstring).

    [EN] BUG FIX, found while wiring this up: vismol_object.
    _bonds_from_pair_of_indexes_list()'s external_orders parameter
    already existed but its actual assignment (bond.bond_order = ...) was
    commented out (a dead "pass" in its place) -- passing external_orders
    had literally no effect before. Fixed there directly (see that
    method's own updated comment) since there was no way to make this
    feature work otherwise. """
    pair = ( bond.atom_index_i, bond.atom_index_j )
    pair = ( min ( pair ), max ( pair ) )

    new_order = ( bond.bond_order % 3 ) + 1

    from gui.windows.builder.atom_ops import _reapply_manual_bonds, push_undo_snapshot, adjust_hydrogens
    push_undo_snapshot ( vismol_object )

    # capturados ANTES de qualquer ajuste de hidrogenio -- ver comentario
    # abaixo sobre reler .atom_id do objeto ao vivo em vez de reusar pair[0]/pair[1]
    atom_a_obj = vismol_object.atoms[pair[0]]
    atom_b_obj = vismol_object.atoms[pair[1]]

    if not hasattr ( vismol_object, "manual_bond_orders" ) or vismol_object.manual_bond_orders is None:
        vismol_object.manual_bond_orders = { }
    vismol_object.manual_bond_orders[pair] = new_order

    _reapply_manual_bonds ( vismol_object )

    vismol_object.create_representation ( rep_type = "lines" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    # [EN] The bond's order changed, so BOTH atoms' valence sums changed
    # (e.g. single -> double frees up one unit of valence on each side,
    # which could now need one FEWER hydrogen apiece). Same "re-read
    # .atom_id from the live object" reasoning as finish_bond_drag()'s
    # own hydrogen-adjustment call -- adjusting atom_a_obj first might
    # remove a lower-numbered hydrogen than atom_b_obj, which would shift
    # atom_b_obj's own id if we used a stale cached int instead.
    adjust_hydrogens ( vismol_object, atom_a_obj.atom_id )
    adjust_hydrogens ( vismol_object, atom_b_obj.atom_id )

    from gui.windows.builder.empty_object import sync_pdynamo_system
    sync_pdynamo_system ( vismol_object )

    dprint ( "DEBUG click_mode: bond order cycled -- atoms #{} <-> #{} now order {}".format (
            pair[0], pair[1], new_order ) )

    vm_glcore.queue_draw ( )
    return new_order


# =====================================================================================
#   Ctrl+drag to reposition an EXISTING atom
#   ------------------------------------------------------------------------------
#   Fifth click interaction for the "add" tool (alongside "click on empty
#   space -> new atom", "click on an existing atom -> replace its
#   element", "click+drag from an atom -> new bonded atom", and
#   "Ctrl+click a bond -> cycle its order"): Ctrl+press-and-HOLD on an
#   existing atom, then drag -- moves THAT SAME atom to follow the
#   cursor, live, with NO new atom and NO new bond created (unlike the
#   plain, non-Ctrl click-and-drag feature above). Releasing drops it at
#   its current position.
#
#   Deliberately reuses the exact same "lazy start on real motion" +
#   "fixed depth for the whole gesture" design as start_bond_drag()/
#   update_bond_drag()/finish_bond_drag() above (see those functions'
#   own docstrings for the reasoning -- it's identical here): a plain
#   Ctrl+CLICK (no real drag) on a BOND still means "cycle its order"
#   (cycle_bond_order() above) -- these two Ctrl interactions never
#   conflict because they target different things (an ATOM vs. a BOND's
#   on-screen LINE), and because the atom-drag only ever actually
#   STARTS on real mouse movement, exactly like the plain click-and-
#   drag-to-create-a-bond feature does for the non-Ctrl case.
# =====================================================================================

def start_atom_drag ( vm_glcore, atom, depth ):
    """ Begins a Ctrl+drag-to-reposition-an-existing-atom interaction.
    Called from mouse_motion() (see the hook added there) the FIRST time
    real mouse movement is detected after a Ctrl+press that landed on an
    atom belonging to the current Builder target object -- render() only
    records that atom as a CANDIDATE (mirroring start_bond_drag()'s own
    candidate mechanism, see the render() hook in vismol_glcore.py); the
    drag doesn't actually start until this runs, so a plain Ctrl+click
    (press+release, no real movement) never reaches here.

    Pushes an undo snapshot immediately (BEFORE any repositioning
    happens) -- same "snapshot once per logical gesture, not once per
    motion event" reasoning as start_bond_drag(). """
    vm_session    = vm_glcore.vm_session
    vismol_object = atom.vm_object

    from gui.windows.builder.atom_ops import push_undo_snapshot
    push_undo_snapshot ( vismol_object )

    vm_session.builder_ctrl_drag_active = True
    vm_session.builder_ctrl_drag_atom   = atom
    vm_session.builder_ctrl_drag_object = vismol_object
    vm_session.builder_ctrl_drag_depth  = depth

    dprint ( "DEBUG click_mode: atom-drag started -- moving atom #{} ('{}')".format (
            atom.atom_id, atom.symbol ) )

    vm_glcore.updated_coords = True
    vm_glcore.queue_draw ( )


def update_atom_drag ( vm_glcore, mouse_x, mouse_y ):
    """ Called from mouse_motion() on every motion event while
    vm_session.builder_ctrl_drag_active is True. Repositions the dragged
    atom (atom_ops.move_atom() -- cheap, no bond recompute, exactly like
    update_bond_drag() uses for the atom it creates) to the current
    cursor position, unprojected at the FIXED depth captured by
    start_atom_drag() -- same reasoning as update_bond_drag()'s own
    fixed-depth choice: keeps the atom moving in a plane parallel to the
    screen, rather than jumping depth if the cursor crosses over some
    other atom/bond mid-drag.

    Pure math -- safe to call directly from mouse_motion() (a plain GTK
    handler), not render() -- see update_bond_drag()'s own docstring for
    why (world_pos_from_mouse() only touches the GPU when its `depth`
    argument is None, and here it never is). """
    vm_session = vm_glcore.vm_session
    if not getattr ( vm_session, "builder_ctrl_drag_active", False ):
        return None

    vismol_object = vm_session.builder_ctrl_drag_object
    atom          = vm_session.builder_ctrl_drag_atom
    depth         = vm_session.builder_ctrl_drag_depth

    wx, wy, wz = world_pos_from_mouse ( vm_glcore, mouse_x, mouse_y, depth = depth )

    world_point = np.array ( [ wx, wy, wz, 1.0 ], dtype = np.float32 )
    inv_model   = np.linalg.inv ( vismol_object.model_mat )
    local_point = world_point @ inv_model
    x, y, z = float ( local_point[0] ), float ( local_point[1] ), float ( local_point[2] )

    from gui.windows.builder.atom_ops import move_atom
    move_atom ( vismol_object, atom.atom_id, x, y, z )
    return atom


def finish_atom_drag ( vm_glcore ):
    """ Called from mouse_released() (checked FIRST, before the existing
    Ctrl+click-a-bond handling -- see the hook added there) once the
    button comes back up while an atom-drag is active.

    Unlike finish_bond_drag(), there's no bonding to reconcile here at
    all -- moving an atom doesn't change which OTHER atoms it's bonded
    to, or their bond orders, so no add_bond()/adjust_hydrogens() calls
    are needed. Just syncs the linked pDynamo system (positions changed,
    even though bonds/topology didn't -- see empty_object.
    sync_pdynamo_system()'s own docstring: it rebuilds coordinates3 every
    time regardless) and clears the drag state.

    Returns the repositioned Atom, or None if no drag was active. """
    vm_session = vm_glcore.vm_session
    if not getattr ( vm_session, "builder_ctrl_drag_active", False ):
        return None

    vismol_object = vm_session.builder_ctrl_drag_object
    atom          = vm_session.builder_ctrl_drag_atom

    dprint ( "DEBUG click_mode: atom-drag finished -- atom #{} ('{}') repositioned".format (
            atom.atom_id, atom.symbol ) )

    vm_session.builder_ctrl_drag_active           = False
    vm_session.builder_ctrl_drag_atom             = None
    vm_session.builder_ctrl_drag_object           = None
    vm_session.builder_ctrl_drag_depth            = None
    vm_session.builder_ctrl_press_candidate_atom  = None
    vm_session.builder_ctrl_press_candidate_depth = None

    vm_glcore.updated_coords = False

    from gui.windows.builder.empty_object import sync_pdynamo_system
    sync_pdynamo_system ( vismol_object )

    vm_glcore.queue_draw ( )
    return atom


# =====================================================================================
#   Hover highlight -- a flat, camera-facing (billboard) yellow ring
#   ------------------------------------------------------------------------------
#   Drawn around whichever atom find_atom_at_pixel_2d_any_object() reports
#   as hovered (see the mouse_motion() hook in vismol_glcore.py) -- makes
#   the SAME atom that's about to be printed also visible, and doubles as
#   a sanity check that the pure-CPU hover test agrees with what the
#   normal GPU pick would select (which is what it's actually built out
#   of now -- see find_atom_at_pixel_2d_any_object()'s own docstring:
#   the hover test used to disagree with real click-picking because it
#   ignored occlusion; drawing the highlighted atom directly on top of
#   the real atom makes any remaining disagreement immediately obvious
#   -- if the ring ever visibly sits on the WRONG atom, that's the thing
#   to re-check).
#
#   "Camera-facing, doesn't rotate with the model" -- built directly in
#   WORLD space (not the hovered object's own local/model space): the
#   ring's plane is spanned by the camera's OWN right/up axes (extracted
#   from view_matrix -- this codebase implements "rotate the view" by
#   rotating every vismol_object's OWN model_mat, see _rotate_view(), so
#   the camera's view_matrix itself never changes -- meaning these two
#   axes are effectively FIXED for the whole session, not something that
#   needs recomputing defensively every frame, though it's cheap enough
#   to just do it every time regardless), not the object's local axes --
#   so the ring's ORIENTATION stays exactly the same regardless of how
#   the model has been rotated, only its POSITION follows the atom.
#   Verified numerically offline first (see the conversation this was
#   built in): projected a ring built this way through several different
#   camera orientations and confirmed the projected shape's variation
#   was IDENTICAL across all of them (proving the ring's screen-space
#   shape doesn't depend on model_mat rotation) before wiring this in.
#
#   [EN] FILLED + TRANSPARENT variant added afterwards (the user asked
#   for a filled, partially-transparent disk instead of just an
#   outline): the "lines_sel" shader used for the outline ring hardcodes
#   alpha=1.0 in its OWN fragment shader (sel_fragment_shader_lines --
#   confirmed by reading shaders/lines.py directly: `final_color =
#   vec4(frag_color, 1.0)`), so no amount of vertex-colour or
#   glBlendFunc tweaking on this end could ever make it transparent --
#   the shader itself throws the alpha away. Rather than risk modifying
#   an existing, shared shader (used by real bond-line rendering
#   elsewhere) to add alpha support it was never built for, a small,
#   dedicated shader pair (_HOVER_FILL_VERTEX_SHADER/
#   _HOVER_FILL_FRAGMENT_SHADER below) was written instead, taking a
#   genuine vec4 (RGBA) per-vertex colour straight through to
#   final_color, compiled once via vm_glcore.load_shaders() (which
#   already handles linking AND binding the shared camera UBO to the new
#   program automatically -- confirmed by reading load_shaders() itself,
#   no extra setup needed on this end) and cached on vm_glcore so it
#   only compiles once, not every time an atom is hovered.
# =====================================================================================

_HOVER_FILL_VERTEX_SHADER = """
#version 330
precision highp float;
precision highp int;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
uniform mat4 model_mat;

in vec3 vert_coord;
in vec4 vert_color;

out vec4 frag_color;

void main(){
    frag_color = vert_color;
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

_HOVER_FILL_FRAGMENT_SHADER = """
#version 330
precision highp float;
precision highp int;

in vec4 frag_color;
out vec4 final_color;

void main(){
    final_color = frag_color;
}
"""


def _get_hover_fill_program ( vm_glcore ):
    """ Lazily compiles (once) and caches the small dedicated shader used
    by the FILLED hover disk -- see the module-level note above for why
    this couldn't just reuse "lines_sel". """
    program = getattr ( vm_glcore, "builder_hover_fill_program", None )
    if program is None:
        program = vm_glcore.load_shaders ( _HOVER_FILL_VERTEX_SHADER, _HOVER_FILL_FRAGMENT_SHADER )
        vm_glcore.builder_hover_fill_program = program
    return program


def draw_hover_highlight ( vm_glcore, atom, n_segments = 24, filled = True,
                            color = ( 1.0, 1.0, 0.0 ), alpha = 0.10 ):
                            #color = ( 1.0, 1.0, 0.0 ), alpha = 0.15 ):
    """ Draws a flat, camera-facing highlight around `atom`'s CURRENT
    position -- a semi-transparent FILLED disk by default (filled=True),
    or a solid outline ring (filled=False, the original style, opaque --
    see the module-level note above for why that style can't support
    transparency). Called from render() (see the hook added there),
    every frame that vm_glcore.builder_hover_atom is set -- rebuilds a
    tiny, throwaway VAO/VBO pair each time (a couple dozen vertices at
    most) rather than caching one: the position needs to track the
    hovered atom live anyway (in case it moves while hovered), and a
    buffer this small costs nothing worth optimising away.

    color : RGB tuple, 0-1 range.
    alpha : 0 (fully transparent) - 1 (fully opaque). Only meaningful
            when filled=True -- the outline ring is always fully opaque
            (see above).

    Returns (world_center, up, radius) -- the same camera-facing "up"
    vector and world position this function already had to compute,
    handed back so draw_hover_info_text() (see below) can position the
    info line just below this same highlight without recomputing any of
    it itself. """
    vismol_object = atom.vm_object
    view_matrix   = vm_glcore.glcamera.view_matrix

    right = view_matrix[0:3, 0]
    up    = view_matrix[0:3, 1]
    right = right / np.linalg.norm ( right )
    up    = up    / np.linalg.norm ( up )

    local_pos = _current_frame_position ( vismol_object, atom.atom_id )
    local_pos_h = np.array ( [ local_pos[0], local_pos[1], local_pos[2], 1.0 ], dtype = np.float32 )
    world_center = ( local_pos_h @ vismol_object.model_mat )[:3]

    #radius = float ( getattr ( atom, "vdw_rad", None ) or 0.4 ) * 1.3
    radius = 0.47

    thetas = np.linspace ( 0, 2 * np.pi, n_segments, endpoint = False )
    ring_points = [ world_center + right * np.cos ( t ) * radius + up * np.sin ( t ) * radius for t in thetas ]

    #GL.glEnable ( GL.GL_DEPTH_TEST )
    GL.glDisable ( GL.GL_DEPTH_TEST )
    GL.glEnable ( GL.GL_BLEND )
    GL.glBlendFunc ( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA )

    if filled:
        # [EN] GL_TRIANGLE_FAN: centre vertex first, then every
        # perimeter point in order, then the FIRST perimeter point
        # again to close the fan -- standard technique for a filled
        # N-gon (here, N=n_segments, close enough to a circle).
        vertices = [ world_center ] + ring_points + [ ring_points[0] ]
        coords = np.array ( vertices, dtype = np.float32 )
        rgba = ( color[0], color[1], color[2], alpha )
        colors = np.tile ( np.array ( rgba, dtype = np.float32 ), ( len ( vertices ), 1 ) )

        program = _get_hover_fill_program ( vm_glcore )
        GL.glUseProgram ( program )
        vm_glcore.load_matrices ( program, np.identity ( 4, dtype = np.float32 ) )   # ja em world space -- model_mat = identidade

        vao = GL.glGenVertexArrays ( 1 )
        GL.glBindVertexArray ( vao )

        coord_vbo = GL.glGenBuffers ( 1 )
        GL.glBindBuffer ( GL.GL_ARRAY_BUFFER, coord_vbo )
        GL.glBufferData ( GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_DYNAMIC_DRAW )
        att_position = GL.glGetAttribLocation ( program, "vert_coord" )
        GL.glEnableVertexAttribArray ( att_position )
        GL.glVertexAttribPointer ( att_position, 3, GL.GL_FLOAT, GL.GL_FALSE, 3 * coords.itemsize, ctypes.c_void_p ( 0 ) )

        col_vbo = GL.glGenBuffers ( 1 )
        GL.glBindBuffer ( GL.GL_ARRAY_BUFFER, col_vbo )
        GL.glBufferData ( GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_DYNAMIC_DRAW )
        att_color = GL.glGetAttribLocation ( program, "vert_color" )
        GL.glEnableVertexAttribArray ( att_color )
        GL.glVertexAttribPointer ( att_color, 4, GL.GL_FLOAT, GL.GL_FALSE, 4 * colors.itemsize, ctypes.c_void_p ( 0 ) )

        GL.glDrawArrays ( GL.GL_TRIANGLE_FAN, 0, len ( vertices ) )

        GL.glBindVertexArray ( 0 )
        GL.glDeleteBuffers ( 2, [ coord_vbo, col_vbo ] )
        GL.glDeleteVertexArrays ( 1, [ vao ] )
        GL.glUseProgram ( 0 )

    else:
        # [EN] Original outline-ring style -- "lines_sel" shader, always
        # fully opaque (see the module-level note above). sel_vertex_
        # shader_lines' geometry shader consumes GL_LINES (pairs of
        # vertices, one segment per pair) -- NOT GL_LINE_LOOP -- so each
        # consecutive pair of ring points (wrapping back to the first)
        # needs to be laid out as its OWN pair in the vertex buffer.
        
        
        GL.glEnable ( GL.GL_LINE_SMOOTH )
        GL.glHint ( GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST )

        vertices = [ ]
        for i in range ( n_segments ):
            vertices.append ( ring_points[i] )
            vertices.append ( ring_points[ ( i + 1 ) % n_segments ] )

        coords = np.array ( vertices, dtype = np.float32 )
        colors = np.tile ( np.array ( color, dtype = np.float32 ), ( len ( vertices ), 1 ) )

        program = vm_glcore.shader_programs["lines_sel"]
        GL.glUseProgram ( program )
        GL.glLineWidth ( 2.0 )
        vm_glcore.load_matrices ( program, np.identity ( 4, dtype = np.float32 ) )

        vao = GL.glGenVertexArrays ( 1 )
        GL.glBindVertexArray ( vao )

        coord_vbo = GL.glGenBuffers ( 1 )
        GL.glBindBuffer ( GL.GL_ARRAY_BUFFER, coord_vbo )
        GL.glBufferData ( GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_DYNAMIC_DRAW )
        att_position = GL.glGetAttribLocation ( program, "vert_coord" )
        GL.glEnableVertexAttribArray ( att_position )
        GL.glVertexAttribPointer ( att_position, 3, GL.GL_FLOAT, GL.GL_FALSE, 3 * coords.itemsize, ctypes.c_void_p ( 0 ) )

        col_vbo = GL.glGenBuffers ( 1 )
        GL.glBindBuffer ( GL.GL_ARRAY_BUFFER, col_vbo )
        GL.glBufferData ( GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_DYNAMIC_DRAW )
        att_color = GL.glGetAttribLocation ( program, "vert_color" )
        GL.glEnableVertexAttribArray ( att_color )
        GL.glVertexAttribPointer ( att_color, 3, GL.GL_FLOAT, GL.GL_FALSE, 3 * colors.itemsize, ctypes.c_void_p ( 0 ) )

        GL.glDrawArrays ( GL.GL_LINES, 0, len ( vertices ) )

        GL.glBindVertexArray ( 0 )
        GL.glDeleteBuffers ( 2, [ coord_vbo, col_vbo ] )
        GL.glDeleteVertexArrays ( 1, [ vao ] )

        GL.glDisable ( GL.GL_LINE_SMOOTH )
        GL.glLineWidth ( 1 )
        GL.glUseProgram ( 0 )

    GL.glDisable ( GL.GL_BLEND )

    return world_center, up, radius


def draw_hover_info_text ( vm_glcore, atom, world_center, up, radius ):
    """ [EN] Draws a short info line (atom index, element, residue/chain,
    object name) just below the hover ring built by draw_hover_highlight()
    (see its own docstring -- world_center/up/radius are exactly what
    that function already computed, handed back so this doesn't need to
    redo any of it). Called from render() right after the ring.

    Reuses this codebase's OWN existing freetype text pipeline (the
    ACTIVE LabelRepresentation class, gui/libgl/representations.py --
    confirmed by checking which of the two same-named classes in that
    file is actually live code and which is dead, commented-out legacy:
    counted every triple-quote in the file to track the string-literal
    state and found the SECOND LabelRepresentation sits entirely inside
    an unclosed docstring block, so it's inert -- the FIRST one, and its
    same vm_font.* calls used here, are the real, working API).

    Text is BILLBOARD by construction, same as the ring, but for a
    different, simpler reason: read shaders/vm_freetype.py directly and
    confirmed its vertex shader applies ONLY view_mat (no model_mat) to
    each character's world-space anchor point, and its geometry shader
    expands each point into a quad by offsetting X/Y in VIEW SPACE (not
    world space) before applying proj_mat -- so every character quad
    faces the camera automatically, with no billboard math needed on
    this end beyond supplying a single WORLD-space position per
    character (same "verify by reading the actual shader source, not by
    assuming" approach already used for the ring's own camera-facing
    math).

    [EN] NOT copied verbatim from LabelRepresentation: that class
    transforms its anchor point through `vm_glcore.model_mat` before
    upload, which is a DIFFERENT matrix from any given atom's own
    `vismol_object.model_mat` (the one this whole Builder feature set
    has consistently used everywhere else, verified numerically more
    than once earlier in this project). Using vismol_object.model_mat
    here instead keeps this consistent with everything else already
    built, rather than reusing a call that looks likely to be a
    pre-existing, rarely-exercised bug in that class. """
    
    
    #this fuction is not been used -but could
    return False
    
    
    
    font = getattr ( vm_glcore, "builder_hover_font", None )
    if font is None:
        from vismol.libgl.vismol_font import VismolFont
        font = VismolFont ( color = [ 1.0, 1.0, 0.0, 1.0 ] )
        font.set_dimensions ( width = 0.12, height = 0.12 )
        font.make_freetype_font ( )
        font.make_freetype_texture ( vm_glcore.core_shader_programs["freetype"] )
        vm_glcore.builder_hover_font = font

    residue = getattr ( atom, "residue", None )
    chain   = getattr ( atom, "chain", None )
    resn    = residue.name  if residue is not None else "?"
    resi    = residue.index if residue is not None else "?"
    chain_name = chain.name if chain is not None else "?"

    text = "#{}/{}/{} - {}/{}".format (
    #text = "#{} {} {}{}/{} @ {}".format (
    #        atom.atom_id, atom.symbol, resn, resi, chain_name, atom.vm_object.name )
            atom.atom_id, atom.symbol, atom.name, resn, resi)

    anchor = world_center - up * ( radius + 0.35 )   # um pouco abaixo do anel

    xyz_pos   = [ ]
    uv_coords = [ ]
    chars     = 0

    GL.glBindTexture ( GL.GL_TEXTURE_2D, font.texture_id )
    for i, c in enumerate ( text ):
        chars += 1
        c_id = ord ( c )
        x = c_id % 16
        y = c_id // 16 - 2
        xyz_pos.append ( anchor[0] + i * font.char_width - ( len ( text ) * font.char_width ) / 2.0 )
        xyz_pos.append ( anchor[1] )
        xyz_pos.append ( anchor[2] )
        uv_coords.append ( x * font.text_u )
        uv_coords.append ( y * font.text_v )
        uv_coords.append ( ( x + 1 ) * font.text_u )
        uv_coords.append ( ( y + 1 ) * font.text_v )

    xyz_pos   = np.array ( xyz_pos, dtype = np.float32 )
    uv_coords = np.array ( uv_coords, dtype = np.float32 )

    GL.glBindBuffer ( GL.GL_ARRAY_BUFFER, font.coord_vbo )
    GL.glBufferData ( GL.GL_ARRAY_BUFFER, xyz_pos.itemsize * len ( xyz_pos ), xyz_pos, GL.GL_DYNAMIC_DRAW )
    GL.glBindBuffer ( GL.GL_ARRAY_BUFFER, font.text_vbo )
    GL.glBufferData ( GL.GL_ARRAY_BUFFER, uv_coords.itemsize * len ( uv_coords ), uv_coords, GL.GL_DYNAMIC_DRAW )
    GL.glBindBuffer ( GL.GL_ARRAY_BUFFER, 0 )

    GL.glDisable ( GL.GL_DEPTH_TEST )   # texto sempre legivel, mesmo atras de outra geometria -- mesma convencao ja usada por LabelRepresentation
    GL.glEnable ( GL.GL_BLEND )
    GL.glBlendFunc ( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA )
    GL.glUseProgram ( vm_glcore.core_shader_programs["freetype"] )

    font.load_matrices ( vm_glcore.core_shader_programs["freetype"],
                          vm_glcore.glcamera.view_matrix, vm_glcore.glcamera.projection_matrix )
    font.load_font_params ( vm_glcore.core_shader_programs["freetype"] )

    GL.glBindVertexArray ( font.vao )
    GL.glDrawArrays ( GL.GL_POINTS, 0, chars )
    GL.glDisable ( GL.GL_BLEND )
    GL.glEnable ( GL.GL_DEPTH_TEST )
    GL.glBindVertexArray ( 0 )
    GL.glUseProgram ( 0 )
