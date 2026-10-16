# Making EasyHybrid render natively on macOS (GTK3 + Quartz)

## Summary

EasyHybrid's main window rendered **completely blank/white on macOS**
(native Quartz GTK3 backend) — no menu bar, no panels, nothing, even
though the process was running fine and not crashed. Root cause and fix
below. After these changes the app runs natively on macOS (Apple
Silicon, tested on an M4) with no X11/XQuartz dependency, and molecules
render correctly (verified visually with `examples/pkl/1UAO_OPLS.pkl`).

All changes are contained to `src/graphics_engine` (the Vismol
submodule) plus one small, backwards-compatible addition to
`src/gui/main/main_window.py`. Nothing about the on-screen behavior,
controls, or file formats changed for other platforms.

## Root cause

GTK3's Quartz (native macOS) backend has a long-standing, apparently
never-fixed bug: realizing a `Gtk.GLArea` corrupts Cairo/window
compositing for the **entire** containing `NSWindow`, leaving it
permanently blank. This reproduces with a minimal 10-line script (a
bare `Gtk.Window` containing nothing but a `Gtk.GLArea` — no
EasyHybrid/Vismol code involved at all — already renders fully blank,
including plain `Gtk.Label` widgets sharing the same window that have
nothing to do with the GLArea).

This is GTK3-specific. GTK4 has a real fix (`GSK_RENDERER=gl`), but GSK
doesn't exist in GTK3 and there is no equivalent. Running under
XQuartz/X11 instead of Quartz was also tried and rejected: GTK3's X11
GL path segfaults on first draw on macOS. Neither backend has a working
`GtkGLArea` on this platform.

## The fix: render OpenGL outside of GTK's own compositing

`src/graphics_engine/src/vismol/gui/vismol_gtkwidget.py`:

- `VismolGTKWidget` no longer subclasses `Gtk.GLArea`. It's now a
  `Gtk.DrawingArea`.
- A new `_OffscreenGLContext` class creates a **hidden GLFW window**
  (via the `glfw` Python package + native GLFW library) purely to own
  an OpenGL 3.3 core context and its default framebuffer. This context
  is never attached to any GTK/AppKit view, so it never touches GTK's
  Quartz compositing.
- On every `"draw"` signal (Cairo callback), the widget: makes the
  offscreen context current, resizes the hidden GLFW window to match
  the widget's allocation, runs the *exact same* `VismolGLCore.render()`
  call as before, reads the finished frame back with `glReadPixels`,
  and blits it into the Cairo context via `cairo.ImageSurface` +
  `cr.paint()`.
- `"realize"`/`"render"`/`"resize"` GLArea signals were replaced with
  `"draw"` and `"size-allocate"` (DrawingArea equivalents).
- `set_has_depth_buffer()` / `set_has_alpha()` compatibility shims were
  added on the widget (the offscreen GLFW window is always created with
  depth+alpha, so these are no-ops) since `VismolGLCore.initialize()`
  calls them on `self.parent_widget` expecting the old GtkGLArea API.
- `capture_screenshot()` (used by PNG export / the preview window) now
  calls `self._gl_ctx.make_current()` instead of the old
  `Gtk.GLArea.make_current()`/`get_error()`.

New Python dependency: **`glfw`** (pip) + the **`glfw`** native library
(available via conda-forge, `mamba install -c conda-forge glfw`). No
other platform needs this — it's only imported by this one file, and
only actually exercised on Quartz.

Everything else in Vismol's render loop (`vismol_glcore.py`,
`representations.py`, shaders, camera, picking, trajectory playback,
export, etc.) is untouched — the fix is purely about *how the finished
frame reaches the screen*, not how it's produced.

### External API surface (unaffected)

Only two things outside `vismol_gtkwidget.py` touch this widget
directly, and both still work unchanged: `glwidget.queue_draw()`
(`treeview_menu.py`, `treeview_menu_new.py`) and
`glwidget.open_export_preview()` (`main_window.py`).

## Follow-on bugs this exposed

Once the window actually started rendering, three **pre-existing**
OpenGL correctness bugs surfaced — invisible before only because
nothing ever got far enough to hit them. These are not
macOS/Quartz-specific; a strict/modern OpenGL 3.3 **core profile**
driver (which is what Apple's GL always was) exposes them. They would
likely also surface on Linux the day its GL driver stops being lenient
about them (e.g. under any implementation that enforces the spec
strictly, such as Mesa's "no error" build turned off, ANGLE, or a
newer NVIDIA/AMD driver in core mode).

### 1. `gl_FragColor` in a core-profile shader

`src/vismol/libgl/shaders/pick.py`:
`fragment_shader_picking_dots_safe` used the legacy `gl_FragColor`
built-in, which does not exist under GLSL 330 core (removed after
GLSL 1.30/1.40). Fixed by declaring `out vec4 frag_color;` and writing
to that instead, matching the convention already used by every other
shader in the codebase.

(There is a second, unrelated `gl_FragColor` use in
`src/vismol/libgl/shaders/surface.py`, inside a `vertex`/`fragment`
shader pair that is never imported/used anywhere — left as-is,
out of scope.)

### 2. `varying` used as a plain global inside a single shader stage

`src/vismol/libgl/shaders/lines.py` and
`src/vismol/libgl/shaders/sticks.py`: several geometry shaders declared
helper variables (`rot_mat`, `mid_point`, `p_00`..`p_26`, `point_a/b/c`)
with the old `varying` qualifier, used purely as ordinary global
variables shared between functions *within the same shader stage* —
not as an actual vertex→geometry or geometry→fragment interface. `varying`
was removed from GLSL 330 core entirely. Fixed by dropping the
qualifier (making them plain globals, exactly like the neighboring
`center0`/`center1`/`mid_c`/`eff_rad` declarations already in the same
shaders).

### 3. Vertex attribute arrays disabled after unbinding the VAO

Five call sites, in `glaxis.py` (x2), `selection_box.py`,
`shapes.py` (x2), and `vismol_font.py`, built a VAO like:

```c
GL.glBindVertexArray(0)              // unbind FIRST
GL.glDisableVertexAttribArray(...)   // then try to disable — invalid!
```

Under a strict core-profile driver, `glDisableVertexAttribArray` while
*no* VAO is bound (id 0 has no meaning in core profile) raises
`GL_INVALID_OPERATION`. PyOpenGL turns that into a hard Python
exception, which aborted every single draw call that touched these
objects (the coordinate-axis gizmo, the selection box, primitive
shapes, on-screen text). Fixed by reordering: disable the attribute
arrays *first*, while the real VAO is still bound, then unbind.

### 4. `glLineWidth()` beyond what the driver supports

Apple's GL only guarantees `glLineWidth(1.0)` under a core profile
(`GL_ALIASED_LINE_WIDTH_RANGE` is commonly just `[1.0, 1.0]`). The
codebase calls `glLineWidth()` with wider values (2–5) in ~34 places
across `glaxis.py`, `representations.py`, and `selection_box.py` for
the axis gizmo, selection box, and wire representations.
Rather than edit all 34 call sites, this is patched centrally in
`vismol_gtkwidget.py` (`_patch_gl_line_width()`, called once at
startup): every module does `from OpenGL import GL` and calls
`GL.glLineWidth(...)`, i.e. all of them look up the function by
attribute on the *same shared* `OpenGL.GL` module object at call time
rather than binding their own copy at import time — so patching that
one attribute once covers every call site. The patched version clamps
the requested width to whatever the driver actually reports.

## Files changed

Main repo (`EasyHybrid3`):
- `src/gui/main/main_window.py` — minor: `window.present()` +
  `GLib.idle_add()` initial-redraw nudge (harmless on all platforms,
  belt-and-suspenders for the very first frame).
- `easyhybrid.py` — `gi.disable_legacy_autoinit()` +
  explicit `Gtk.init([])` before importing `Gtk`/`Gdk` (avoids
  PyGObject's automatic `Gdk.init_check()` running before the Quartz
  backend is fully loaded).

Submodule (`src/graphics_engine`, Vismol):
- `src/vismol/gui/vismol_gtkwidget.py` — the main fix (GLArea →
  DrawingArea + offscreen GLFW context), described above.
- `src/vismol/libgl/shaders/pick.py`, `lines.py`, `sticks.py` — GLSL
  core-profile compliance fixes (#1, #2 above).
- `src/vismol/libgl/glaxis.py`, `selection_box.py`, `shapes.py`,
  `vismol_font.py` — VAO/attribute-array ordering fixes (#3 above).

## New dependency

```
mamba install -c conda-forge glfw   # native library
pip install glfw                     # Python ctypes bindings
```

Only used on the code path that creates `VismolGTKWidget`
(`vismol_gtkwidget.py`); no impact on non-GUI/headless usage of
pDynamo3/Vismol.

## Known pre-existing, unrelated issue (not fixed, out of scope)

The CLI convenience shortcut in `easyhybrid.py` (`easyhybrid3
somefile.pdb` / `.pkl` as the last argument) calls
`vm_session.load_molecule(filein)` directly, which — for both raw PDB
*and* `.pkl` pDynamo system files — never correctly registers the
object's pDynamo `e_id`, and fails inside `eSession.py:
_add_vismol_object` (`KeyError: None`) — silently, because the call
site wraps it in a bare `except: pass`. The real "File → Open" /
"File → New" GUI menu paths (`vm_session.load(...)`, which is what
actually gets exercised in normal use) work correctly and are what was
used to verify the rendering fix end-to-end. This bug is
platform-independent (same code path would fail identically on Linux)
and unrelated to the macOS rendering work above — flagging it here
only because it was noticed along the way.
