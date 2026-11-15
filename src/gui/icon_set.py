#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  icon_set.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  2026-11-01: ícones da barra principal à escolha do usuário.
#    "modern"  = os ícones da 3.2 (icons/eh_color/png/32, definidos no glade)
#    "classic" = os ícones da 3.1.3 (icons/classic, 24 px)
#  Só as três barras que já existiam na 3.1.3 (toolbar1, toolbar2_selections,
#  toolbar4_pdynamo_tools) têm versão antiga; a barra do Builder (toolbar_edit)
#  fica sempre com os ícones novos. Selecionado por gl_parameters["toolbar_icons"]
#  (Preferences > Startup, aplicado na hora) ou pela variável EASYHYBRID_ICONS.
#  Reversível: o pixbuf do glade fica guardado na própria GtkImage.
# ============================================================================

import os

import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk, GdkPixbuf

DEFAULT_ICON_SET = "classic"

# GtkImage id em MainWindow.glade -> arquivo em src/gui/icons/classic
CLASSIC_ICONS = {
    "image1"  : "icon_new_project.png",
    "image2"  : "icon_open.png",
    "image3"  : "icon_save.png",
    "image4"  : "icon_save_as.png",
    "image20" : "icon_terminal.png",
    "image21" : "icon_player_play.png",
    "image10" : "draw_structure.png",
    "image22" : "icon_selection_list.png",
    "image23" : "icon_go_to_atom.png",
    "image24" : "icon_pDynamo_selections2.2.png",
    "image25" : "icon_QCSetup1.3.png",
    "image14" : "icon_single_point1.2.png",
    "image15" : "icon_minimization_24x24.png",
    "image13" : "icon_molecular_dynamics_24x24.png",
    "image12" : "icon_normal_modes.png",
    "image5"  : "icon_scan_24x24.png",
    "image6"  : "icon_umbrella_sampling_24x24.png",
    "image7"  : "icon_NEB_24x24.png",
    "image8"  : "icon_save_image.png",
    "image9"  : "icon_process_manager.png",
}

_ICON_SETS = { "modern"  : "Modern (EasyHybrid 3.2)",
               "classic" : "Classic (EasyHybrid 3.1)" }


def icon_set_choices ( ):
    """ [(name, label), ...] for a selector, the default first. """
    return list ( _ICON_SETS.items ( ) )


_current_name = DEFAULT_ICON_SET


def current_icon_set ( ):
    return _current_name


def icon_set_from_config ( gl_parameters = None ):
    """ EASYHYBRID_ICONS environment variable first, then gl_parameters["toolbar_icons"]. """
    name = os.environ.get ( "EASYHYBRID_ICONS" )
    if not name and gl_parameters:
        name = gl_parameters.get ( "toolbar_icons" )
    name = str ( name or DEFAULT_ICON_SET ).strip ( )
    return name if name in _ICON_SETS else DEFAULT_ICON_SET


def apply_icon_set ( builder, home, name = DEFAULT_ICON_SET ):
    """ Puts icon set `name` on the main window toolbars of `builder` (the
        MainWindow.glade builder). The theme's tint of dark icons is redone
        afterwards (gui.theme.adapt_icons_for_theme). Returns the name applied. """
    global _current_name
    if name not in _ICON_SETS:
        name = DEFAULT_ICON_SET
    folder = os.path.join ( home, "src", "gui", "icons", "classic" )
    window = None
    for image_id, filename in CLASSIC_ICONS.items ( ):
        image = builder.get_object ( image_id )
        if image is None:
            continue
        # undo the dark theme's tint first, so the right pixbuf is kept
        tinted = getattr ( image, "_eh_original_pixbuf", None )
        if tinted is not None:
            image.set_from_pixbuf ( tinted )
            image._eh_original_pixbuf = None
        modern = getattr ( image, "_eh_modern_pixbuf", None )
        if modern is None:
            modern = image._eh_modern_pixbuf = image.get_pixbuf ( )
        pixbuf = modern
        if name == "classic":
            try:
                pixbuf = GdkPixbuf.Pixbuf.new_from_file ( os.path.join ( folder, filename ) )
            except Exception:
                pixbuf = modern
        if pixbuf is not None:
            image.set_from_pixbuf ( pixbuf )
        window = image.get_toplevel ( )
    if window is not None:
        try:
            from gui.theme import adapt_icons_for_theme
            adapt_icons_for_theme ( window )
        except Exception:
            pass
    _current_name = name
    return name
