#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: mini trajectory player under the main treeview
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-21, user's request): a small, simplified version of
#  the Trajectory Player at the bottom of the Objects treeview: first,
#  previous, play/pause, next, last (16 px icons), a compact frame slider and
#  "frame / total". The range is the largest frame count of the visible
#  objects (as main_treeview.refresh_trajectory_scalebar). Playback runs on
#  a GLib timer in the GTK thread (the full player moves widgets from a
#  worker thread), loops at the end, and a slow timer keeps the slider in
#  step with frame changes made elsewhere (full player, terminal, ...).
#  Mouse wheel on the slider: one frame per step.
#

import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk, Gdk, GLib

_CSS = b"""
.eh-mini-player button { padding: 1px 3px; min-height: 0; min-width: 0; }
.eh-mini-player label  { font-size: 85%; }
.eh-mini-player scale  { padding: 0 2px; min-height: 0; }
.eh-mini-player scale slider { min-width: 10px; min-height: 10px; margin: -4px; }
"""


class MiniPlayer ( Gtk.Box ):

    FPS = 20

    def __init__ ( self, main ):
        super ( ).__init__ ( orientation = Gtk.Orientation.HORIZONTAL, spacing = 1 )
        self.main = main
        self.vm_session = main.vm_session
        self._play_id = None
        self._syncing = False
        self.n_frames = 1
        self.set_border_width ( 2 )
        self.get_style_context ( ).add_class ( "eh-mini-player" )
        provider = Gtk.CssProvider ( ); provider.load_from_data ( _CSS )
        Gtk.StyleContext.add_provider_for_screen ( Gdk.Screen.get_default ( ), provider, Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION )

        def button ( icon, tip, callback ):
            b = Gtk.Button ( )
            b.set_image ( Gtk.Image.new_from_icon_name ( icon, Gtk.IconSize.MENU ) )
            b.set_relief ( Gtk.ReliefStyle.NONE )
            b.set_tooltip_text ( tip )
            b.set_can_focus ( False )
            b.connect ( "clicked", callback )
            self.pack_start ( b, False, False, 0 )
            return b
        button ( "media-skip-backward-symbolic", "First frame", lambda b: self.go ( 0 ) )
        button ( "media-seek-backward-symbolic", "Previous frame", lambda b: self.step ( -1 ) )
        self.button_play = button ( "media-playback-start-symbolic", "Play / pause (loops)", self.on_play )
        button ( "media-seek-forward-symbolic", "Next frame", lambda b: self.step ( +1 ) )
        button ( "media-skip-forward-symbolic", "Last frame", lambda b: self.go ( self.n_frames - 1 ) )

        self.adjustment = Gtk.Adjustment ( value = 0, lower = 0, upper = 0, step_increment = 1, page_increment = 10 )
        self.scale = Gtk.Scale ( orientation = Gtk.Orientation.HORIZONTAL, adjustment = self.adjustment )
        self.scale.set_draw_value ( False )
        self.scale.set_digits ( 0 )
        self.scale.set_hexpand ( True )
        self.scale.set_can_focus ( False )
        self.scale.set_tooltip_text ( "Frame (mouse wheel: one frame per step)" )
        self.scale.connect ( "value-changed", self.on_scale_changed )
        self.scale.connect ( "scroll-event", self.on_scroll )
        self.pack_start ( self.scale, True, True, 0 )

        self.label = Gtk.Label ( label = "1 / 1" )
        self.label.set_width_chars ( 9 )
        self.label.set_xalign ( 1.0 )
        self.pack_start ( self.label, False, False, 2 )

        GLib.timeout_add ( 400, self._sync )

    # ------------------------------------------------------------------
    def _count_frames ( self ):
        higher = 1
        for vobject in self.vm_session.vm_objects_dic.values ( ):
            if getattr ( vobject, "active", False ):
                higher = max ( higher, len ( vobject.frames ) )
        return higher

    def refresh ( self ):
        """ Range and label from the session (called when objects change). """
        self.n_frames = self._count_frames ( )
        frame = min ( int ( self.vm_session.get_frame ( ) ), self.n_frames - 1 )
        self._syncing = True
        self.adjustment.set_upper ( self.n_frames - 1 )
        self.adjustment.set_value ( frame )
        self._syncing = False
        self.label.set_text ( "{} / {}".format ( frame + 1, self.n_frames ) )
        self.set_sensitive ( self.n_frames > 1 )
        if self.n_frames <= 1 and self._play_id is not None:
            self.pause ( )

    def _sync ( self ):
        if self._play_id is None:
            try:
                self.refresh ( )
            except Exception:
                pass
        return True

    def go ( self, frame ):
        frame = max ( 0, min ( int ( frame ), self.n_frames - 1 ) )
        self.vm_session.set_frame ( frame )
        self._syncing = True
        self.adjustment.set_value ( frame )
        self._syncing = False
        self.label.set_text ( "{} / {}".format ( frame + 1, self.n_frames ) )

    def step ( self, delta ):
        self.n_frames = self._count_frames ( )
        self.go ( int ( self.vm_session.get_frame ( ) ) + delta )

    # ------------------------------------------------------------------
    def on_scale_changed ( self, scale ):
        if not self._syncing:
            self.go ( int ( round ( scale.get_value ( ) ) ) )

    def on_scroll ( self, widget, event ):
        direction = event.direction
        if direction == Gdk.ScrollDirection.SMOOTH:
            ok, dx, dy = event.get_scroll_deltas ( )
            delta = 1 if dy > 0 else -1 if dy < 0 else 0
        else:
            delta = 1 if direction in ( Gdk.ScrollDirection.DOWN, Gdk.ScrollDirection.RIGHT ) else -1
        if delta: self.step ( delta )
        return True

    def on_play ( self, button ):
        if self._play_id is None: self.play ( )
        else: self.pause ( )

    def play ( self ):
        self.n_frames = self._count_frames ( )
        if self.n_frames <= 1: return
        self.button_play.set_image ( Gtk.Image.new_from_icon_name ( "media-playback-pause-symbolic", Gtk.IconSize.MENU ) )
        self._play_id = GLib.timeout_add ( int ( 1000 / self.FPS ), self._tick )

    def pause ( self ):
        if self._play_id is not None:
            GLib.source_remove ( self._play_id )
            self._play_id = None
        self.button_play.set_image ( Gtk.Image.new_from_icon_name ( "media-playback-start-symbolic", Gtk.IconSize.MENU ) )

    def _tick ( self ):
        frame = int ( self.vm_session.get_frame ( ) ) + 1
        self.go ( 0 if frame >= self.n_frames else frame )
        return True
