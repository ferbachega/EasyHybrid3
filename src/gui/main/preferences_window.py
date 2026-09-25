#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2025 Fernando Bachega
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
#  Description:
#      Provides functions for selecting atoms and residues in pDynamo systems
#      to facilitate QC/MM partitioning and molecular simulations.
#
# [EN] This file used to open with ~75 lines importing nearly every
# window class in the app (simulation/analysis windows, pCore, numpy,
# ...) -- a copy-paste leftover from main_window.py, none of it actually
# used anywhere below: PreferencesWindow only ever needs Gtk (it builds
# its own small inline glade XML). Removed as dead weight/misleading
# coupling; if something here genuinely needs one of those again, import
# it explicitly at that point instead of restoring this block wholesale.
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

class PreferencesWindow:
    """ Class doc """
    
    def __init__ (self, main = None, e_id = None , v_id = None):
        """ Class initialiser """
        
        self.xml='''
<?xml version="1.0" encoding="UTF-8"?>
<!-- Generated with glade 3.22.2 -->
<interface>
  <requires lib="gtk+" version="3.20"/>
  <object class="GtkWindow" id="window">
    <property name="can_focus">False</property>
    <property name="default_width">300</property>
    <child type="titlebar">
      <placeholder/>
    </child>
    <child>
      <object class="GtkAlignment">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="top_padding">5</property>
        <property name="bottom_padding">5</property>
        <property name="left_padding">5</property>
        <property name="right_padding">5</property>
        <child>
          <object class="GtkBox">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="orientation">vertical</property>
            <property name="spacing">5</property>
            <child>
              <object class="GtkGrid">
                <property name="visible">True</property>
                <property name="can_focus">False</property>
                <property name="row_spacing">5</property>
                <property name="column_spacing">5</property>
                <child>
                  <object class="GtkEntry" id="entry_name">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="hexpand">True</property>
                    <property name="text" translatable="yes">new pDynamo system</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_name">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="halign">end</property>
                    <property name="label" translatable="yes">Name:</property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_tag">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="halign">end</property>
                    <property name="label" translatable="yes">Tag: </property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">1</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkEntry" id="entry_tag">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="max_length">15</property>
                    <property name="width_chars">10</property>
                    <property name="text" translatable="yes">molsys</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">1</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_color">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="label" translatable="yes">Color:</property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">2</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkColorButton" id="button_color">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                    <property name="rgba">rgb(138,226,52)</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">2</property>
                  </packing>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">True</property>
                <property name="position">0</property>
              </packing>
            </child>
            <child>
              <object class="GtkButtonBox" id="dialog-action_area2">
                <property name="visible">True</property>
                <property name="can_focus">False</property>
                <property name="layout_style">end</property>
                <child>
                  <object class="GtkButton" id="button_cancel">
                    <property name="label" translatable="yes">Cancel</property>
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                  </object>
                  <packing>
                    <property name="expand">False</property>
                    <property name="fill">False</property>
                    <property name="position">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkButton" id="button_apply">
                    <property name="label" translatable="yes">Apply</property>
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                  </object>
                  <packing>
                    <property name="expand">False</property>
                    <property name="fill">False</property>
                    <property name="position">1</property>
                  </packing>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">False</property>
                <property name="position">2</property>
              </packing>
            </child>
          </object>
        </child>
      </object>
    </child>
  </object>
</interface>
'''        
        
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        
        self.e_id  = e_id 
        self.v_id  = v_id 
         
        
        self.builder = Gtk.Builder()
        self.builder.add_from_string(self.xml)

        self.window = self.builder.get_object('window')

        self.entry_name   = self.builder.get_object('entry_name')
        self.entry_tag    = self.builder.get_object('entry_tag')
        self.button_color = self.builder.get_object('button_color')
        
        self.button_cancel = self.builder.get_object('button_cancel')
        self.button_apply  = self.builder.get_object('button_apply')
        
        self.button_apply.connect('clicked', self.on_button_apply )
        self.button_cancel.connect('clicked', self.on_button_cancel )
        # [EN] BUG FIX: this used to be commented out, so closing the
        # window via the window manager's own close button (instead of
        # clicking Apply/Cancel) never reset rename_window_visible --
        # left it stuck True forever, so the NEXT "Rename" click found
        # rename_window_visible still True and called set_names() on
        # this already-destroyed window instead of opening a fresh one.
        # Connecting to "destroy" itself (which fires no matter HOW the
        # window closes) instead of duplicating the reset in both button
        # handlers guarantees it's reset exactly once, every time.
        self.window.connect("destroy", self._on_window_destroyed)
        self.main.main_treeview.treeview_menu.rename_window_visible = True

        self.window.set_resizable(False)
        self.window.show_all()
        self.button_color.hide()
        self.builder.get_object('label_color').hide()

    def set_names (self, name, tag):
        """ Function doc """

        self.entry_name.set_text(name)
        self.entry_tag .set_text(tag)

    def _on_window_destroyed (self, widget):
        """ Function doc """
        self.main.main_treeview.treeview_menu.rename_window_visible = False

    def on_button_apply (self, button):
        """ Function doc """
        name = self.entry_name.get_text()
        tag  = self.entry_tag .get_text()

        # [EN] BUG FIX: rename() can reject the new name (e.g. another
        # vobject already has it) and used to just dprint() the reason
        # with no dialog -- this code never checked the return value, so
        # the window closed anyway, looking like a successful rename that
        # silently did nothing. Now: only apply the tag and close the
        # window if the name change actually went through; on failure,
        # rename() has already shown the user why (see main_window.py),
        # and the window stays open so they can pick a different name.
        if not self.main.rename(e_id = self.e_id, v_id = self.v_id, name = name):
            return

        self.main.rename_tag(e_id = self.e_id, tag = tag)
        self.window.destroy()

    def on_button_cancel (self, button):
        """ Function doc """
        self.window.destroy()
