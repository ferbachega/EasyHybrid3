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
import numpy as np
from OpenGL import GL


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

    from gui.windows.builder.atom_ops import add_bond
    created = add_bond ( atom_a.vm_object, atom_a.atom_id, atom_b.atom_id )

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
        print ( "DEBUG click_mode: depth buffer at clicked pixel = background (nothing rendered there)" )
        return None, None

    atom = vm_glcore.vm_session.atom_dic_id.get ( pickedID )

    depth_raw = GL.glReadPixels ( x, y, 1, 1, GL.GL_DEPTH_COMPONENT, GL.GL_FLOAT )
    depth_buffer_value = float ( np.frombuffer ( depth_raw, dtype = np.float32 ) [0] )

    proj = vm_glcore.glcamera.projection_matrix
    p22, p32 = float ( proj[2,2] ), float ( proj[3,2] )
    ndc_z = 2.0 * depth_buffer_value - 1.0
    distance = p32 / ( ndc_z + p22 )

    print ( "DEBUG click_mode: depth buffer at clicked pixel = {:.5f}  ->  distance from camera = {:.3f}  ->  atom = {}".format (
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

    print ( "DEBUG click_mode: mouse=({:.1f}, {:.1f})  viewport=({:.0f}x{:.0f})  ndc=({:.3f}, {:.3f})".format (
            mouse_x, mouse_y, width, height, ndc_x, ndc_y ) )
    print ( "DEBUG click_mode: fovy={:.2f} aspect={:.3f} depth_used={:.3f}".format ( fovy, aspect, depth ) )
    print ( "DEBUG click_mode: view_point=({:.3f}, {:.3f}, {:.3f})".format ( view_x, view_y, view_z ) )
    print ( "DEBUG click_mode: world_point=({:.3f}, {:.3f}, {:.3f})  <- new atom position".format (
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
        print ( "DEBUG click_mode: an atom WAS clicked -- atom #{} ({}), object '{}'".format (
                picked_atom.atom_id, picked_atom.symbol, picked_atom.vm_object.name ) )
    
        print(picked_atom.atom_id, picked_atom.bonds )
    else:
        print ( "DEBUG click_mode: no atom was clicked (empty space)." )

    if picked_atom is not None and picked_atom.vm_object is vismol_object:
        if picked_atom.symbol == symbol:
            print ( "DEBUG click_mode: clicked atom #{} is already '{}' -- nothing to replace.".format (
                    picked_atom.atom_id, symbol ) )
            return picked_atom

        from gui.windows.builder.atom_ops import set_atom_element
        atom = set_atom_element ( vismol_object, picked_atom.atom_id, symbol )
        print ( "DEBUG click_mode: replaced atom #{} -> '{}'".format ( atom.atom_id, symbol ) )
        return atom

    wx, wy, wz = world_pos_from_mouse ( vm_glcore, mouse_x, mouse_y, depth = depth )

    world_point = np.array ( [ wx, wy, wz, 1.0 ], dtype = np.float32 )
    inv_model = np.linalg.inv ( vismol_object.model_mat )
    local_point = world_point @ inv_model
    x, y, z = float ( local_point[0] ), float ( local_point[1] ), float ( local_point[2] )

    print ( "DEBUG click_mode: world_point=({:.3f}, {:.3f}, {:.3f})  -> local_point (after inv(model_mat))=({:.3f}, {:.3f}, {:.3f})".format (
            wx, wy, wz, x, y, z ) )

    from gui.windows.builder.atom_ops import add_atom
    atom = add_atom ( vismol_object, symbol = symbol, x = x, y = y, z = z )
    
    
    
    
    tmp = {'C': [[-0.785298, 0.243518, -0.653254], [0.322015, -0.981331, -0.189814], [-0.334691, 0.073016, 0.992665], [0.798227, 0.665009, -0.149645]], 'N': [[-0.785298, 0.243518, -0.653254], [0.322015, -0.981331, -0.189814], [-0.334691, 0.073016, 0.992665]], 'O': [[-0.785298, 0.243518, -0.653254], [0.322015, -0.981331, -0.189814]]}
    
    #tmp = {'C' : [
    #             [-0.785298,  0.243518, -0.653254],
    #             [ 0.322015, -0.981331, -0.189814],
    #             [-0.334691,  0.073016,  0.992665],
    #             [ 0.798227,  0.665009, -0.149645]
    #             ],
    #       
    #       'N' : [
    #             [-0.785298,  0.243518, -0.653254],
    #             [ 0.322015, -0.981331, -0.189814],
    #             [-0.334691,  0.073016,  0.992665],
    #             ],
    #       
    #       'O' : [
    #             [-0.785298,  0.243518, -0.653254],
    #             [ 0.322015, -0.981331, -0.189814],
    #             ]
    #      }
    
    if picked_atom:
       
       pass
    
    
    
    # aqui é quando adicionamos um átomo novo  não ligado a uma subestrutura existente
    else:
        if symbol in tmp.keys():
            H_list = tmp[symbol]
            
            for xyz  in H_list:
                atm_tmp = add_atom ( vismol_object, symbol = "H", 
                                             x = xyz[0]+x, 
                                             y = xyz[1]+y, 
                                             z = xyz[2]+z 
                                             )
    
    
        
    for bond in atom.bonds:
        print(atom,  bond.atom_i.symbol, bond.atom_j.symbol)
    return atom
