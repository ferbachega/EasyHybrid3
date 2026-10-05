#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  theme.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  [EN] 2026-10-02 user request: try EasyHybrid with the look of the
#  "Scientific Ribbon" GTK3 prototype (~/Downloads/gtk3_ribbon_dark.py --
#  dark navy palette, cyan accent), in a way that is trivially reversible.
#
#  Purely a CSS layer on top of the user's GTK theme: no .glade file and no
#  widget is changed, so removing the provider restores the original look
#  exactly. Selected by:
#    - environment variable EASYHYBRID_THEME (wins; handy to just try it:
#        EASYHYBRID_THEME=ribbon_dark python3 easyhybrid.py )
#    - gl_parameters["interface_theme"] in .config.json (persistent)
#  "default" (or anything unknown) = no extra CSS = the original EasyHybrid.
#  set_theme() can also switch live (it removes the previous provider).
# ============================================================================

#    ┌───────────────┬────────────────────────────────────────────────────────┐
#    │     Nome      │                         Visual                         │
#    ├───────────────┼────────────────────────────────────────────────────────┤
#    │ default       │ EasyHybrid original, sem CSS extra                     │
#    ├───────────────┼────────────────────────────────────────────────────────┤
#    │ ribbon_dark   │ azul-marinho com destaque ciano (o protótipo)          │
#    ├───────────────┼────────────────────────────────────────────────────────┤
#    │ graphite_dark │ cinza neutro com destaque laranja                      │
#    ├───────────────┼────────────────────────────────────────────────────────┤
#    │ nord          │ paleta Nord, cinza-azulado com destaque azul-gelo      │
#    ├───────────────┼────────────────────────────────────────────────────────┤
#    │ slate_light   │ claro, com destaque azul; bom para projetor ou figuras │
#    ├───────────────┼────────────────────────────────────────────────────────┤
#    │ high_contrast │ preto e branco com destaque amarelo; acessibilidade    │
#    └───────────────┴────────────────────────────────────────────────────────┘
#    
#    Como usar
#    - Escolher na combo muda na hora. Para manter na próxima sessão, use "Apply and Save".
#    - "Reset parameters" volta para default.
#    - A variável EASYHYBRID_THEME=nome continua funcionando para testar sem salvar nada.
#    
#    Dois problemas visuais corrigidos
#    1. Radio e checkbox: o seu tema do sistema (Mint-L-Dark) desenha esses botões com imagens. Por isso o CSS não conseguia mudar a cor: no tema claro os radios desmarcados ficavam pretos e o marcado ficava verde-oliva. Troquei por um "check" que segue a cor do tema e um ponto de radio feito só em CSS. Agora ficam com a cor de destaque de cada tema.
#    2. Ícones da barra principal nos temas escuros: o traço preto de ícones mistos (documento escuro com átomo vermelho) sumia no fundo escuro. Isso já acontecia no ribbon_dark. Agora só os pixels pretos/cinza-escuros clareiam e as cores ficam. A mudança é reversível: ao passar para um tema claro, os ícones originais voltam sem reiniciar.
#    
#    Para criar outro tema, basta uma paleta e uma linha no theme.py. As cores que você não definir vêm do ribbon_dark:
#    _MEU = dict(window="#...", panel="#...", accent="#...", selection="#...", ...)
#    _register("meu_tema", "Meu tema", _MEU, dark=True)   # dark=True: clareia ícones pretos
#    No topo do theme.py há a lista das cores que dá para mudar e o que cada uma colore.
#    
#    Arquivos alterados
#    - src/gui/theme.py
#    - src/gui/config.py: novo padrão interface_theme
#    - src/gui/windows/setup/setup_interface.glade e setup_interface.py: a combo nova, que só acrescenta uma linha à grade da aba Startup
#    - easyhybrid.py: adapta os ícones da janela principal
#    
#    O fundo da área 3D continua controlado pela cor de fundo do Viewer. Se quiser, posso fazer cada tema sugerir uma cor de fundo combinando.


import os

import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk, Gdk


DEFAULT_THEME = "default"

# ----------------------------------------------------------------------------
#  [EN] 2026-10-02: more themes. Every theme is a PALETTE (a dict with the
#  keys of _BASE_KEYS) poured into the single CSS template below, so a new
#  theme = a new palette + one line in _register(). A theme that needs more
#  than colours (radius, padding, fonts) can still pass its own CSS string.
#  Keys:
#    window/text/text_dim/caption   base background, text, dim text, frame titles
#    bar/ribbon/panel/workspace     menubar, toolbars+buttons, lists, notebook pages
#    border/border_soft             outlines
#    hover/hover_border/hover_text  hovered buttons / menu items
#    accent/accent_text             highlight colour (checked tab, checks, sliders)
#    tab_checked/tab_hover          notebook tabs
#    entry/entry_border/entry_text  text inputs
#    selection/on_selection         selected rows, checked buttons, and their text
#    status/focus                   statusbar, focus ring
#    check_mark/knob                glyph on a checked check/radio; switch+scale knob
#    danger_bg/danger_border/danger_text   destructive buttons
# ----------------------------------------------------------------------------

# palette of the prototype (gtk3_ribbon_dark.py)
_RIBBON_DARK = {
    "window":      "#0b1220",
    "text":        "#dce8f7",
    "text_dim":    "#91b2d4",
    "caption":     "#7fa7d0",
    "bar":         "#0d1728",
    "ribbon":      "#101d31",
    "panel":       "#101d30",
    "workspace":   "#080f1b",
    "border":      "#24415f",
    "border_soft": "#203653",
    "hover":       "#193453",
    "hover_border":"#2d5d8d",
    "hover_text":  "#ffffff",
    "accent":      "#16a3ff",
    "accent_text": "#62c4ff",
    "tab_checked": "#132b47",
    "tab_hover":   "#142944",
    "entry":       "#14263d",
    "entry_border":"#31577e",
    "entry_text":  "#e5f2ff",
    "selection":   "#174b7a",
    "on_selection":"#ffffff",
    "status":      "#0d192a",
    "focus":       "#19a7ff",
    "check_mark":  "#ffffff",
    "knob":        "#dce8f7",
    "danger_bg":   "#5a1f2a",
    "danger_border":"#8c3443",
    "danger_text": "#ffffff",
}
_BASE_KEYS = tuple ( _RIBBON_DARK )

# neutral dark grey, orange accent
_GRAPHITE_DARK = dict ( _RIBBON_DARK,
    window = "#1e1f22", text = "#e3e3e3", text_dim = "#a9abb0", caption = "#c8b08a",
    bar = "#18191b", ribbon = "#26282c", panel = "#222326", workspace = "#151618",
    border = "#3a3d42", border_soft = "#303236",
    hover = "#34373c", hover_border = "#6a6e76", hover_text = "#ffffff",
    accent = "#ff9f1c", accent_text = "#ffb84d", tab_checked = "#2d2a24", tab_hover = "#2b2d31",
    entry = "#2a2c30", entry_border = "#4a4e55", entry_text = "#f0f0f0",
    selection = "#7a4a0e", on_selection = "#ffffff", status = "#1a1b1d", focus = "#ff9f1c",
    check_mark = "#1e1f22", knob = "#e3e3e3" )

# light blue-grey, blue accent (projectors, screenshots for papers)
_SLATE_LIGHT = dict ( _RIBBON_DARK,
    window = "#eef1f5", text = "#1f2933", text_dim = "#52606d", caption = "#3e4c59",
    bar = "#e1e6ec", ribbon = "#f7f9fb", panel = "#ffffff", workspace = "#e4e8ee",
    border = "#c3ccd6", border_soft = "#d5dce4",
    hover = "#dde8f5", hover_border = "#8fb3dc", hover_text = "#1f2933",
    accent = "#2f6fdb", accent_text = "#1d5bbf", tab_checked = "#ffffff", tab_hover = "#e9eef4",
    entry = "#ffffff", entry_border = "#b4bfcc", entry_text = "#1f2933",
    selection = "#2f6fdb", on_selection = "#ffffff", status = "#e1e6ec", focus = "#2f6fdb",
    check_mark = "#ffffff", knob = "#ffffff",
    danger_bg = "#fde8ea", danger_border = "#d64550", danger_text = "#a3212c" )

# Nord (https://www.nordtheme.com) -- polar night + frost
_NORD = dict ( _RIBBON_DARK,
    window = "#2e3440", text = "#e5e9f0", text_dim = "#b4bccb", caption = "#88c0d0",
    bar = "#272c36", ribbon = "#3b4252", panel = "#323845", workspace = "#242933",
    border = "#4c566a", border_soft = "#434c5e",
    hover = "#434c5e", hover_border = "#5e81ac", hover_text = "#eceff4",
    accent = "#88c0d0", accent_text = "#8fbcbb", tab_checked = "#3b4252", tab_hover = "#384050",
    entry = "#3b4252", entry_border = "#4c566a", entry_text = "#eceff4",
    selection = "#5e81ac", on_selection = "#ffffff", status = "#2a2f3a", focus = "#88c0d0",
    check_mark = "#2e3440", knob = "#d8dee9",
    danger_bg = "#5c2f35", danger_border = "#bf616a", danger_text = "#ffffff" )

# black/white/yellow, for accessibility and very bright rooms
_HIGH_CONTRAST = dict ( _RIBBON_DARK,
    window = "#000000", text = "#ffffff", text_dim = "#e6e6e6", caption = "#ffd700",
    bar = "#000000", ribbon = "#000000", panel = "#000000", workspace = "#000000",
    border = "#ffffff", border_soft = "#9a9a9a",
    hover = "#333333", hover_border = "#ffd700", hover_text = "#ffffff",
    accent = "#ffd700", accent_text = "#ffd700", tab_checked = "#000000", tab_hover = "#262626",
    entry = "#000000", entry_border = "#ffffff", entry_text = "#ffffff",
    selection = "#ffd700", on_selection = "#000000", status = "#000000", focus = "#ffd700",
    check_mark = "#000000", knob = "#ffffff",
    danger_bg = "#000000", danger_border = "#ff4d4d", danger_text = "#ff8080" )

_CSS_TEMPLATE = """
/* ---- base ---- */
window, dialog, .background, messagedialog {{ background-color: {window}; color: {text}; }}
label {{ color: {text}; }}
label:disabled, *:disabled {{ color: alpha({text}, 0.45); }}
separator {{ background-color: {border}; min-width: 1px; min-height: 1px; }}
tooltip, tooltip.background {{ background-color: {ribbon}; color: {text}; border: 1px solid {hover_border}; }}
tooltip * {{ color: {text}; }}

/* ---- menus / bars ---- */
menubar {{ background-color: {bar}; border-bottom: 1px solid {border_soft}; }}
menubar > menuitem {{ color: {text_dim}; padding: 4px 10px; }}
menubar > menuitem:hover {{ background-color: {tab_checked}; color: {accent_text}; box-shadow: inset 0 -2px {accent}; }}
menu, .menu, .context-menu, popover {{ background-color: {ribbon}; color: {text}; border: 1px solid {border}; }}
menu menuitem {{ color: {text}; padding: 4px 10px; }}
menu menuitem:hover {{ background-color: {hover}; color: {hover_text}; }}
menu menuitem:hover label {{ color: {hover_text}; }}
menu menuitem:disabled label {{ color: alpha({text}, 0.4); }}
toolbar, .toolbar, .primary-toolbar {{ background-color: {ribbon}; border-bottom: 1px solid {border}; padding: 2px; }}
toolbar button, toolbutton button {{ background: transparent; border: 1px solid transparent; border-radius: 3px; }}
toolbar button:hover, toolbutton button:hover {{ background-color: {hover}; border-color: {hover_border}; }}
toolbar button:checked, toolbutton button:checked {{ background-color: {tab_checked}; border-color: {accent}; }}
headerbar, .titlebar {{ background: {bar}; color: {text}; border-bottom: 1px solid {border_soft}; }}
statusbar {{ background-color: {status}; color: {text_dim}; border-top: 1px solid {border_soft}; }}
statusbar label {{ color: {text_dim}; }}

/* ---- notebooks (tabs) ---- */
notebook > header {{ background-color: {bar}; border-color: {border_soft}; }}
notebook > header tab {{ color: {text_dim}; padding: 4px 12px; border: none; }}
notebook > header tab:hover {{ background-color: {tab_hover}; color: {text}; }}
notebook > header tab:checked {{ color: {accent_text}; background-color: {tab_checked}; box-shadow: inset 0 -3px {accent}; }}
notebook > header tab label {{ color: inherit; }}
notebook > stack {{ background-color: {workspace}; }}

/* ---- panels / frames / panes ---- */
frame > border {{ border: 1px solid {border}; }}
frame > label {{ color: {caption}; }}
paned > separator {{ background-color: {border}; }}
scrolledwindow {{ background-color: {panel}; }}
viewport {{ background-color: {panel}; }}
box.panel, .panel {{ background-color: {panel}; border: 1px solid {border}; }}

/* ---- tree / list views ---- */
treeview, treeview.view, iconview, list, .view {{ background-color: {panel}; color: {text}; }}
treeview.view:hover {{ background-color: alpha({hover}, 0.6); }}
treeview.view:selected, treeview:selected, list row:selected, .view:selected, iconview:selected {{
    background-color: {selection}; color: {on_selection}; }}
list row:selected label {{ color: {on_selection}; }}
treeview.view header button, treeview header button {{
    background: {bar}; color: {caption}; border: none; border-right: 1px solid {border_soft};
    border-bottom: 1px solid {border}; padding: 2px 6px; }}
treeview.view header button:hover {{ background: {hover}; color: {text}; }}
list row:hover {{ background-color: alpha({hover}, 0.6); }}

/* ---- buttons ---- */
button {{ background-image: none; background-color: {ribbon}; color: {text};
          border: 1px solid {border}; border-radius: 3px; box-shadow: none; text-shadow: none; }}
button:hover {{ background-color: {hover}; border-color: {hover_border}; }}
button:active, button:checked {{ background-color: {selection}; border-color: {accent}; color: {on_selection};
                                  box-shadow: inset 0 -2px {accent}; }}
button:active label, button:checked label {{ color: {on_selection}; }}
button:disabled {{ background-color: alpha({ribbon}, 0.5); border-color: alpha({border}, 0.5); }}
button:focus {{ outline-color: {focus}; }}
button.flat {{ background: transparent; border-color: transparent; }}
button.flat:hover {{ background-color: {hover}; border-color: {hover_border}; }}
button.suggested-action {{ background-color: {selection}; border-color: {accent}; color: {on_selection}; }}
button.suggested-action label {{ color: {on_selection}; }}
button.destructive-action {{ background-color: {danger_bg}; border-color: {danger_border}; color: {danger_text}; }}
button.destructive-action label {{ color: {danger_text}; }}

/* ---- inputs ---- */
entry, spinbutton, searchentry, spinbutton entry {{
    background-color: {entry}; color: {entry_text}; border: 1px solid {entry_border}; border-radius: 3px;
    box-shadow: none; caret-color: {accent_text}; }}
entry:focus, spinbutton:focus {{ border-color: {accent}; }}
entry selection, textview text selection {{ background-color: {selection}; color: {on_selection}; }}
spinbutton button {{ background: {entry}; border: none; border-left: 1px solid {entry_border}; border-radius: 0; color: {entry_text}; }}
combobox button, combobox box button {{ background-color: {entry}; border-color: {entry_border}; color: {entry_text}; }}
combobox button label, combobox box button label {{ color: {entry_text}; }}
textview, textview text {{ background-color: {panel}; color: {text}; }}
check, radio {{ background-image: none; background-color: {entry}; border: 1px solid {entry_border}; color: {check_mark};
                box-shadow: none; -gtk-icon-source: none; min-width: 14px; min-height: 14px; }}
check {{ border-radius: 3px; }}
radio {{ border-radius: 50%; }}
check:hover, radio:hover {{ border-color: {accent}; }}
/* the system theme draws these with images: replace them by a symbolic
   check (tinted with `color`) and a pure-CSS radio dot (no icon needed) */
check:checked {{ background-color: {accent}; border-color: {accent}; -gtk-icon-source: -gtk-icontheme("object-select-symbolic"); }}
check:indeterminate {{ background-color: {accent}; border-color: {accent}; -gtk-icon-source: -gtk-icontheme("list-remove-symbolic"); }}
radio:checked {{ border-color: {accent}; background-image: radial-gradient(circle, {check_mark} 0%, {check_mark} 28%, {accent} 34%, {accent} 100%); }}
switch {{ background-color: {entry}; border: 1px solid {entry_border}; }}
switch:checked {{ background-color: {accent}; }}
switch slider {{ background-color: {knob}; border: 1px solid {entry_border}; }}
scale trough, progressbar trough {{ background-color: {entry}; border: 1px solid {entry_border}; }}
scale highlight, progressbar progress {{ background-color: {accent}; }}
scale slider {{ background-color: {knob}; border: 1px solid {hover_border}; }}

/* ---- scrollbars ---- */
scrollbar {{ background-color: {workspace}; border: none; }}
scrollbar slider {{ background-color: {border}; border-radius: 6px; min-width: 6px; min-height: 6px; }}
scrollbar slider:hover {{ background-color: {hover_border}; }}
"""
_RIBBON_DARK_CSS = _CSS_TEMPLATE       # old name, kept for anything importing it

THEMES   = { }      # name -> ( label shown in Preferences, CSS )
PALETTES = { }      # name -> palette (None for a theme registered with raw CSS)
# themes with a dark background: monochrome dark icons are lightened for them
DARK_THEMES = set ( )


def _register ( name, label, palette = None, css = None, dark = False ):
    """ Adds a theme. Give a `palette` (rendered through _CSS_TEMPLATE) or a
        ready `css` string. Missing palette keys fall back to ribbon_dark. """
    if palette is not None:
        palette = dict ( _RIBBON_DARK, **palette )
        css = _CSS_TEMPLATE.format ( **palette )
    THEMES[name]   = ( label, css )
    PALETTES[name] = palette
    if dark:
        DARK_THEMES.add ( name )


_register ( "ribbon_dark",   "Scientific Ribbon (dark navy)",    _RIBBON_DARK,   dark = True )
_register ( "graphite_dark", "Graphite (dark, orange accent)",   _GRAPHITE_DARK, dark = True )
_register ( "nord",          "Nord (dark, frost accent)",        _NORD,          dark = True )
_register ( "slate_light",   "Slate (light, blue accent)",       _SLATE_LIGHT )
_register ( "high_contrast", "High contrast (black / yellow)",   _HIGH_CONTRAST, dark = True )


def theme_choices ( ):
    """ [(name, label), ...] for a selector, "default" first. """
    return [ ( DEFAULT_THEME, "EasyHybrid default (system GTK theme)" ) ] + \
           [ ( name, label ) for name, ( label, _css ) in THEMES.items ( ) ]


_current_provider = None
_current_name = DEFAULT_THEME
# widget trees handed to adapt_icons_for_theme(): re-adapted on a live switch
import weakref
_adapted_roots = weakref.WeakSet ( )


def theme_from_config ( gl_parameters = None ):
    """ Theme name to use: EASYHYBRID_THEME environment variable first, then
        gl_parameters["interface_theme"], else "default". """
    env = os.environ.get ( "EASYHYBRID_THEME" )
    if env:
        return env.strip ( )
    if gl_parameters:
        return str ( gl_parameters.get ( "interface_theme", DEFAULT_THEME ) or DEFAULT_THEME )
    return DEFAULT_THEME


def set_theme ( name = DEFAULT_THEME ):
    """ Applies theme `name` to the whole application (every window, current
        and future). "default" / unknown names remove any extra CSS, i.e. the
        original EasyHybrid look. Returns the name actually applied. Icons of
        already-open windows registered via adapt_icons_for_theme() follow. """
    global _current_provider, _current_name
    screen = Gdk.Screen.get_default ( )
    if screen is None:
        return _current_name
    if _current_provider is not None:
        Gtk.StyleContext.remove_provider_for_screen ( screen, _current_provider )
        _current_provider = None
    if name not in THEMES:
        _current_name = DEFAULT_THEME
    else:
        provider = Gtk.CssProvider ( )
        provider.load_from_data ( THEMES[name][1].encode ( "utf-8" ) )
        Gtk.StyleContext.add_provider_for_screen ( screen, provider, Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION )
        _current_provider = provider
        _current_name = name
    for root in list ( _adapted_roots ):
        adapt_icons_for_theme ( root )
    return _current_name


def current_theme ( ):
    return _current_name


def is_dark ( ):
    return _current_name in DARK_THEMES


def _hex_rgb ( color ):
    color = color.lstrip ( "#" )
    return [ int ( color[i:i + 2], 16 ) for i in ( 0, 2, 4 ) ]


def _lightened_pixbuf ( pixbuf, light_hex = "#dce8f7" ):
    """ A copy of `pixbuf` with its dark, monochrome strokes turned light
        (alpha kept), or None when the icon is coloured / not dark -- e.g. the
        Builder's ring-fragment icons are black line drawings that vanish on
        a dark navy button. Original pixbuf untouched. """
    try:
        import numpy as np
        from gi.repository import GdkPixbuf, GLib
        if not pixbuf.get_has_alpha ( ) or pixbuf.get_n_channels ( ) != 4:
            return None
        w, h, stride = pixbuf.get_width ( ), pixbuf.get_height ( ), pixbuf.get_rowstride ( )
        raw = np.frombuffer ( pixbuf.get_pixels ( ), dtype = np.uint8 )
        img = np.zeros ( ( h, w, 4 ), dtype = np.uint8 )
        for y in range ( h ):
            img[y] = raw[y * stride: y * stride + w * 4].reshape ( w, 4 )
        opaque = img[:, :, 3] > 40
        if opaque.sum ( ) < 10:
            return None
        rgb3   = img[:, :, :3].astype ( np.int32 )
        spread = rgb3.max ( axis = 2 ) - rgb3.min ( axis = 2 )
        # dark, unsaturated pixels = the "ink" of the drawing
        ink = opaque & ( rgb3.max ( axis = 2 ) < 110 ) & ( spread < 40 )
        # Only icons that are mostly ink: pure line drawings (Builder rings)
        # and line drawings with coloured details (main toolbar: dark
        # document + red/green atom). Colour icons are left alone, and in
        # mixed ones only the ink pixels change -- the colours stay.
        if ink.sum ( ) < 0.30 * opaque.sum ( ):
            return None
        out = img.copy ( )
        light = np.array ( _hex_rgb ( light_hex ), dtype = np.int32 )
        shade = 1.0 - rgb3.mean ( axis = 2 ) / 255.0     # 1 = black stroke
        for c in range ( 3 ):
            out[:, :, c] = np.where ( ink, ( light[c] * np.clip ( shade, 0.35, 1.0 ) ).astype ( np.uint8 ), img[:, :, c] )
        data = GLib.Bytes.new ( out.tobytes ( ) )
        return GdkPixbuf.Pixbuf.new_from_bytes ( data, GdkPixbuf.Colorspace.RGB, True, 8, w, h, w * 4 )
    except Exception:
        return None


def adapt_icons_for_theme ( widget ):
    """ Makes the GtkImages under `widget` match the current theme: with a
        dark theme, dark monochrome icons get a lightened copy (tinted with
        the theme's text colour, see _lightened_pixbuf()); with a light or
        the default theme, any icon lightened before is put back. Reversible:
        the original pixbuf is kept on the image. `widget` is remembered so
        a later set_theme() re-adapts it live. """
    if widget is None:
        return
    try:
        _adapted_roots.add ( widget )
    except TypeError:
        pass
    dark  = is_dark ( )
    palette = PALETTES.get ( _current_name ) or _RIBBON_DARK
    light_hex = palette["text"]
    stack = [ widget ]
    while stack:
        w = stack.pop ( )
        if isinstance ( w, Gtk.Image ):
            original = getattr ( w, "_eh_original_pixbuf", None )
            if original is not None:                     # undo the previous theme's tint
                w.set_from_pixbuf ( original )
                w._eh_original_pixbuf = None
            if dark and w.get_storage_type ( ) == Gtk.ImageType.PIXBUF:
                pixbuf = w.get_pixbuf ( )
                light  = _lightened_pixbuf ( pixbuf, light_hex )
                if light is not None:
                    w._eh_original_pixbuf = pixbuf
                    w.set_from_pixbuf ( light )
        if isinstance ( w, Gtk.Container ):
            stack.extend ( w.get_children ( ) )
