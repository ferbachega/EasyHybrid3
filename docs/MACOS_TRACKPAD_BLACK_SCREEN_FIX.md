# macOS trackpad scroll could freeze the viewer or hide the molecule

## Symptom

When opening a PKL structure and scrolling quickly on a Mac trackpad, the
Vismol viewer could turn black or stop responding. Some runs also printed
`RuntimeWarning: overflow encountered in cast` while converting fog values
to OpenGL `float32` uniforms. Window resizing could separately assert while
building the perspective matrix if the clipping planes were invalid.

## Cause

GTK reports trackpad motion as smooth scroll deltas. The widget converted the
entire accumulated delta into discrete scroll steps and called
`mouse_scroll()` once for every step in a synchronous loop. Each outward
clipping step widened the camera slab exponentially. A large delta or fast
series of events could therefore stall the GTK draw loop and make the depth
range too large for useful molecular rendering. Fog start/end values are
derived from the camera clipping planes, so invalid or extreme camera/fog
settings could also overflow during OpenGL's `float32` conversion.

## Fix

### Limit smooth-scroll work per event

In `src/graphics_engine/src/vismol/gui/vismol_gtkwidget.py`, cap each smooth
scroll event to four steps and discard excess movement from that event. This
keeps a single trackpad event from creating a long blocking loop.

```python
steps = int(total)
applied_steps = max(-4, min(steps, 4))
self._scroll_remainder = 0.0 if applied_steps != steps else total - steps
if applied_steps == 0:
    return

direction = 1 if applied_steps > 0 else -1
for _ in range(abs(applied_steps)):
    self.vm_glcore.mouse_scroll(direction)
```

### Bound the clipping slab

In `src/graphics_engine/src/vismol/libgl/vismol_glcore.py`, cap the clipping
slab width at 1,000 Å. If a non-finite or wider slab is encountered, restore
a finite slab around the current camera distance. Outward zoom also caps its
exponential width increase at the same limit.

```python
max_slab_width = 1000.0
slab_width = max(float(self.glcamera.z_far - self.glcamera.z_near), 0.05)
if not np.isfinite(slab_width):
    center = max(float(self.dist_cam_zrp), self.glcamera.min_zfar)
    self.glcamera.z_near = center - max_slab_width / 2.0
    self.glcamera.z_far = center + max_slab_width / 2.0
    slab_width = max_slab_width
elif slab_width > max_slab_width:
    center = 0.5 * (float(self.glcamera.z_near) + float(self.glcamera.z_far))
    self.glcamera.z_near = center - max_slab_width / 2.0
    self.glcamera.z_far = center + max_slab_width / 2.0
    slab_width = max_slab_width

if up:
    new_width = min(slab_width * np.exp(rate), max_slab_width)
```

### Keep projection and fog inputs valid

The same `vismol_glcore.py` resize path now ignores zero-sized allocations
and keeps the effective far plane beyond the near plane before constructing
the projection matrix. In `src/graphics_engine/src/vismol/libgl/glcamera.py`,
fog range scale and center offset are checked for non-finite values; invalid
settings fall back to defaults, and range scale is capped at 1,000.

## Verification

Opening `examples/pkl/1UAO_OPLS.pkl` and scrolling quickly with the trackpad
no longer makes the viewer go black or hang. `git diff --check` passes. The
manual GUI check was performed on macOS; the code changes are in the Vismol
submodule at `src/graphics_engine`.
