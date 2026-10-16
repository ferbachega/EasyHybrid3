# Adding a trackpad-friendly pan binding on macOS (Cmd+Right-drag)

## Summary

EasyHybrid's camera panning is bound to a middle-mouse-button drag
(`vismol_glcore.py`). MacBook trackpads have no middle button and no
default gesture that emulates one, so panning was effectively
unreachable for anyone driving the app from a trackpad — rotate
(left-drag) and zoom (right-drag, reachable via a two-finger
secondary-click) both worked fine, but translating the view did not.

This adds a second way to trigger the same pan action: **holding Cmd
while right-dragging** (i.e. Cmd + two-finger-tap-drag on a trackpad),
matching PyMOL's own `Cmd+Right` = translate convention on macOS, so
the muscle memory carries over. The existing plain-right-drag = zoom
and plain-middle-drag = pan bindings are untouched; this is additive
only.

Verified end-to-end on an Apple Silicon Mac (M4): launched the app,
confirmed via a temporary debug print that GTK reports the Cmd key as
keyval `Meta_L` on this platform, then confirmed Cmd+trackpad-drag
visually pans the loaded structure instead of zooming it.

## Change

`src/graphics_engine/src/vismol/libgl/vismol_glcore.py`:
- Added a `self.cmd` flag (alongside the existing `self.ctrl` /
  `self.shift` modifier flags), defaulting to `False`.
- `mouse_pressed()`: right-button-drag now resolves to **pan** instead
  of **zoom** when `self.cmd` is `True`. Plain right-drag (no Cmd) still
  zooms, unchanged.

```python
self.mouse_rotate = left   and not (middle or right)
self.mouse_zoom   = right  and not (middle or left) and not self.cmd
self.mouse_pan    = (middle and not (right or left)) \
                     or (right and self.cmd and not (middle or left))
```

`src/graphics_engine/src/vismol/gui/vismol_gtkwidget.py`:
- Added key-press/release handlers for `Meta_L`, `Meta_R`, `Super_L`,
  and `Super_R`, all toggling the same `self.vm_glcore.cmd` flag. All
  four keyval names are wired (rather than just `Meta_L`) because which
  one GTK reports for the Cmd key can vary by GTK version/backend/
  keymap; only `Meta_L` was actually observed on the tested setup, but
  binding all four is cheap insurance against a different GTK build
  reporting one of the others.

This mirrors the existing pattern already used for `Control_L` /
`Shift_L` in the same file (`self.vm_glcore.ctrl` / `self.vm_glcore.shift`).

## Files changed

Submodule (`src/graphics_engine`, Vismol):
- `src/vismol/libgl/vismol_glcore.py` — `self.cmd` flag + pan/zoom
  gating in `mouse_pressed()`.
- `src/vismol/gui/vismol_gtkwidget.py` — `Meta_L/R` and `Super_L/R`
  key handlers.

## Not changed

- Left-drag rotate, middle-drag pan, and plain right-drag zoom are all
  unchanged for every platform.
- No new dependencies.
- Scroll-wheel zoom and the existing Ctrl/Shift-gated builder and
  selection-box behavior are untouched.
