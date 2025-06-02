#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sequence_plot.py
#  
#  Copyright 2023 Fernando <fernando@fernando-Inspiron-7537>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
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
#  
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, cairo


one_letter_res_dict = {'C':'CYS', 
                       'D':'ASP', 
                       'S':'SER', 
                       'Q':'GLN', 
                       'K':'LYS',
                       'I':'ILE', 
                       'P':'PRO', 
                       'T':'THR', 
                       'F':'PHE', 
                       'N':'ASN', 
                       'G':'GLY', 
                       'H':'HIS', 
                       'L':'LEU', 
                       'R':'ARG', 
                       'W':'TRP', 
                       'A':'ALA', 
                       'V':'VAL', 
                       'E':'GLU', 
                       'Y':'TYR', 
                       'M':'MET',
                       '-':'gap',}

three_letter_res_dict = {'CYS': 'C', 
                       'ASP': 'D', 
                       'SER': 'S', 
                       'GLN': 'Q', 
                       'LYS': 'K',
                       'ILE': 'I', 
                       'PRO': 'P', 
                       'THR': 'T', 
                       'PHE': 'F', 
                       'ASN': 'N', 
                       'GLY': 'G', 
                       'HIS': 'H', 
                       
                       # amber
                       "HID": "H",
                       "HIE": "H",
                       "HIP": "H",
                       "ASH": "D",
                       "GLH": "E",
                       "CYX": "C",
                       
                       # charmm
                       "HSD": "H", 
                       "HSE": "H", 
                       "HSP": "H", 
                       
                       
                       'LEU': 'L', 
                       'ARG': 'R', 
                       'TRP': 'W', 
                       'ALA': 'A', 
                       'VAL': 'V', 
                       'GLU': 'E', 
                       'TYR': 'Y', 
                       'MET': 'M'}

class SeqResidue:
    """ Class doc """
    def __init__ (self, rname = "UNK" , rcode = "X", rnum = 0, chain = "X", color = [1,1,1]):
        """ Class initialiser """
        

        self.rname  = rname
        self.rcode  = rcode
        self.rnum   = rnum 
        self.chain  = chain
        self.vres   = False #vismol residues object
        
        self.is_selected = False
        self.rtype_color = color 
        self.is_solvent  = False 
        
    def change_rtype_color (self):
        """ Function doc """
        pass
    
    def set_active ( self, is_selected = False):
        """ Function doc """
        self.is_selected = is_selected

    def _is_selected (self, ):
        """ Function doc """
        # Checks if **only one** item has selected = True
        self.is_selected = False
        
        if self.vres:
            
            for i in self.vres.atoms.values():
                if i.selected:
                    self.is_selected = True
                else:
                    pass
                
                #print (i, i.selected)
            #self.is_selected = sum(obj.selected for obj in self.vres.atoms.values()) == 1
            
            #print(self.vres.atoms.values())
        else:
            pass
        
        
class TextDrawingArea(Gtk.DrawingArea):

    def __init__(self, main):
        self.main = main
        super(TextDrawingArea, self).__init__()
        self.connect('draw', self.on_draw), 
        self.canvas_width , self.canvas_height = 2200, 2800
        self.set_size_request(self.canvas_width, self.canvas_height)
        
        
        self.font_size  = 20
        self.fonte_type = "DejaVu Sans Mono"
        
        self.char_width  = None
        self.char_height = None
        
        self.gap_width  = 0
        self.gap_height = 6
        
        self.grid_x = 0.0 
        self.grid_y = 0.0 
    
        self.cr = None
    
        self.gap_color    = (0.5, 0.5, 0.5)
        self.marker_color = (1  , 1  , 1)
        self.bg_color     = (0  , 0  , 0)
        #self.bg_color     = (1  , 1  , 1)
        #self.connect("motion-notify-event", self.on_motion)
        #self.connect("button_press_event" , self.button_press)
        #self.connect("clicked" , self.button_press)
        #sefl.on_motion = main.on_motion
    
    def on_motion (self, widget, event):
        """ Function doc """
    
    def get_char_size (self,cr):
        """ Function doc """
        
        cr.select_font_face(self.fonte_type)
        cr.set_font_size(self.font_size)
        
        text = "X"
        x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
        
        #print(x_bearing, y_bearing, 
        #          width, height   , 
        #      x_advance, y_advance)
    
        
        self.char_width  = width
        self.char_height = height
        
    def on_draw(self, widget, cr):
        #print (self.bg_color, self.marker_color)
        self.cr = cr
        #cr.set_source_rgb(0, 0, 0)  # Cor de fundo branca
        #---------------------------------------------------------------
        # Background color
        cr.set_source_rgb(self.bg_color[0], 
                          self.bg_color[1],
                          self.bg_color[2],)  
        cr.paint()
        #---------------------------------------------------------------

        self.get_char_size (cr)
        cr.set_source_rgb(0, 0, 1)  # Cor do texto preto

        i = (self.char_height + self.gap_height)
        
        marker_counter = 0 
        
        '''
        xmax = 0
        ymax = 0
        
        These are the canvas boundaries, which will be used 
        to ensure that the scrollbars function properly.
        '''
        xmax = 0
        ymax = 0
        
        for index, sequence in enumerate(self.main.sequences):
            #print(index, 'index')

            if self.main.mode == 0:
                '''       Draw the placeholders.          '''
                k = 0
                for marker in range(10, len(sequence), 10):
                    
                    k =   12*marker -12
                    
                    cr.set_source_rgb(self.marker_color[0],
                                      self.marker_color[1],
                                      self.marker_color[2],)
                    text =  str(marker)
                    x, y = k, i  # Text position
                    cr.move_to(x, y)
                    cr.show_text(text)
                    
                    
                    x, y = k, i*2
                    cr.move_to(x, y)
                    cr.show_text('|')
                    
                    if y > ymax:
                        ymax = y
                    if x > xmax:
                        xmax = x
                
                #print('xmax', xmax,'ymax', ymax)
                #print(marker, k, x, y)
                show_selection = True
                if self.main.main.vm_session.selections[self.main.main.vm_session.current_selection].active:
                    pass
                else:
                    show_selection = False
                
                
                '''   Draw the resiues - one letter code.     '''
                k = 0
                for residue in sequence:
                    
                    #if self.main.main.vm_session.selections[self.main.vm_session.current_selection].active:
                    #print(self.main.main.vm_session.selections[self.main.main.vm_session.current_selection].active)
                    #if show_selection:
                    
                    x, y = k, i*3 +i*index
                    
                    if show_selection:
                        residue._is_selected()
                        if residue.is_selected:
                            #- - - - - - selection box - - - - - - - -
                            cr.set_source_rgb(0, 0, 1)
                            cr.set_source_rgb(residue.rtype_color[0],residue.rtype_color[1],residue.rtype_color[2] )
                            cr.rectangle(
                                        x,
                                        y-16 ,
                                        12 ,
                                        18 )
                            
                            #cr.rectangle(
                            #             k,
                            #             i-12, 
                            #             k+12, 
                            #             i )
                            cr.fill()
                            #-------------------------------------------
                            
                            cr.set_source_rgb(self.bg_color[0], 
                                            self.bg_color[1],
                                            self.bg_color[2],)
                            
                            
                            
                            #cr.set_source_rgb(1, 0, 0)
                            
                            
                            #cr.set_source_rgb(self.rtype_color)
                            text =  residue.rcode
                            # Posição do texto
                            cr.move_to(x, y)
                            cr.show_text(text)
                        else:
                            if residue.rcode  == '-':
                                cr.set_source_rgb(self.gap_color[0],
                                                  self.gap_color[1],
                                                  self.gap_color[2], )
                            else:
                                cr.set_source_rgb(residue.rtype_color[0],residue.rtype_color[1],residue.rtype_color[2] )
                            
                            text =  residue.rcode
                            x, y = k, i*3 +i*index# Posição do texto
                            cr.move_to(x, y)
                            cr.show_text(text)
                    else:
                        if residue.rcode  == '-':
                            cr.set_source_rgb(self.gap_color[0],
                                              self.gap_color[1],
                                              self.gap_color[2], )
                        else:
                            cr.set_source_rgb(residue.rtype_color[0],residue.rtype_color[1],residue.rtype_color[2] )
                        
                        text =  residue.rcode
                        x, y = k, i*3 +i*index# Posição do texto
                        cr.move_to(x, y)
                        cr.show_text(text)
                    
                    #---------------------------------------------------
                    k = k+self.char_width 
                    if y > ymax:
                        ymax = y
                    if x > xmax:
                        xmax = x
                    #---------------------------------------------------

            
            
            
            
            elif self.main.mode == 1:
                
                k = 0
                for marker in range(10, len(sequence), 10):
                    
                    k =   (12*marker -12)*4 
                    
                    cr.set_source_rgb(1, 1, 1)
                    text =  str(marker)
                    x, y = k, i  # Posição do texto
                    cr.move_to(x, y)
                    cr.show_text(text)
                    
                    
                    x, y = k, i*2
                    cr.move_to(x, y)
                    cr.show_text('|')
                
                k = 0
                for residue in sequence:
                    if residue.is_selected:
                        cr.set_source_rgb(1, 0, 0)
                        text =  residue.rname+' '
                        x, y = k, i*3  # Posição do texto
                        cr.move_to(x, y)
                        cr.show_text(text)
                    
                    else:
                        #print ('')
                        cr.set_source_rgb(0, 0, 1)
                        text =  residue.rname+' '
                        x, y = k, i*3  # Posição do texto
                        cr.move_to(x, y)
                        cr.show_text(text)
                        #print (x, y, text)
                    k = k+self.char_width*4 
            
            #----------------------------------------------------------    
            # It is necessary to change the canvas size 
            # so that the scrollbars work properly.
            self.set_size_request(xmax+50, ymax)
            #----------------------------------------------------------

class GtkSequenceViewer(Gtk.ScrolledWindow):
    """ Class doc """
    
    def __init__ (self, main):
        """ Class initialiser """
        super(GtkSequenceViewer, self).__init__()
        
        self.main = main
        self.label = False
        self.canvas_width  = 800
        self.canvas_height = 200
        self.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        
        
        self.text_drawing_area = TextDrawingArea(main = self )
        self.add(self.text_drawing_area)
        #self.connect("motion-notify-event", self.on_motion)
        #self.connect("button_press_event",  self.button_press)


        self.is_button_1_pressed = False # mouse left button
        self.is_button_2_pressed = False # mouse middle button
        self.is_button_3_pressed = False # mouse right button
        
        
        self.select_in_motion = True
        
        self.mode = 0 # 0 - one letter code / 1 - three letter code

        self.click_x = 0.0
        self.click_y = 0.0
        self.motion_x = 0.0
        self.motion_y = 0.0
        self.sequences = []
        '''
        #-------------------------------------------------------------------------------
        self.sequence  ='---MISYLASIFLLATVSA-VPSGRVEVVFPSVETSRSGVK--TVKFTARVEVVFPSVETSRSGVK--TVKFTA' 
        self.sequence2 ='---MISYLASIFLLTTTTT-VPSGRVEVVFPSVETSRSGVK--TVKFTA-------------------------' 
        self.sequence3 ='---MISYLASIFLLTTTTT-VPSGRVEVVFPSVETSRSGVK--TVKFTA-------TVKFTA------------' 
        
        if self.mode ==1:
            self.set_canvas_width_and_height( width = (len(self.sequence)*12*4)+12  , height = 800)
        elif self.mode ==0:
            self.set_canvas_width_and_height( width = (len(self.sequence)*12)+12  , height = 800)



        self.sequences = []
        
        sequence       = []
        
        counter = 0
        for rnum, res in enumerate(self.sequence):
            
            if res == '-':
                pass
                res3l = one_letter_res_dict[res]
                residue = SeqResidue(rname = res3l , rcode = res, rnum = None)
                sequence.append(residue)
            else:
                res3l = one_letter_res_dict[res]
                
                residue = SeqResidue(rname = res3l , rcode = res, rnum = counter)
                sequence.append(residue)
                counter+=1

        self.sequences.append(sequence)

        
        
        
        sequence       = []
        counter = 0
        for rnum, res in enumerate(self.sequence2):
            
            if res == '-':
                pass
                res3l = one_letter_res_dict[res]
                residue = SeqResidue(rname = res3l , rcode = res, rnum = None)
                sequence.append(residue)
            else:
                res3l = one_letter_res_dict[res]
                
                residue = SeqResidue(rname = res3l , rcode = res, rnum = counter)
                sequence.append(residue)
                counter+=1
        self.sequences.append(sequence)



        sequence       = []
        counter = 0
        for rnum, res in enumerate(self.sequence3):
            
            if res == '-':
                pass
                res3l = one_letter_res_dict[res]
                residue = SeqResidue(rname = res3l , rcode = res, rnum = None)
                sequence.append(residue)
            else:
                res3l = one_letter_res_dict[res]
                
                residue = SeqResidue(rname = res3l , rcode = res, rnum = counter)
                sequence.append(residue)
                counter+=1
        self.sequences.append(sequence)




        self.sequence ='MISYLASIFTTTTHHHHHHHHHHHEVVFPSVETSRSGVKTVKFTA' 
        sequence       = []
        for rnum, res in enumerate(self.sequence):
            res3l = one_letter_res_dict[res]
            residue = Residue(rname = res3l , rcode = res, rnum = rnum)
            sequence.append(residue)
            
        self.sequences.append(sequence)
        #'''
        #print(self.sequences)
        #-------------------------------------------------------------------------------
    def add_sequence_from_vobject (self, vobject):
        """ Function doc """
        #print('aqui', vobject)
        
        if getattr(vobject, 'is_surface', False):
            pass
            #vobject.is_surface = False

        
        else:
            if getattr(vobject, 'e_sequence', False):
                sequence = vobject.e_sequence
            
            else:
                #if vobject.sequence
                sequence = []
                chains = vobject.chains.keys()
                for key in chains:
                    #----------------------------------------------------
                    # - - - - - - - - Finding the color - - - - - - - - -
                    #----------------------------------------------------
                    e_id = vobject.e_id
                    system = self.main.p_session.get_system(index = e_id)
                    color = system.e_color_palette['C']
                    #----------------------------------------------------
    
                    chain = vobject.chains[key]
                    max_index = max(chain.residues.keys())
                    
                    #print('max_index1',max_index)
                    resi_keys = list(chain.residues.keys())
                    resi_keys.sort()
                    
                    
                    index_max = max(resi_keys)
                    max_resi  = chain.residues[index_max]
                    empty_seq = '-'*index_max
                    
                    for index, gap in enumerate(empty_seq):
                        if index+1 in resi_keys:
                            #if the residue is on three_letter_res_dict
                            try:
                                vresidue = chain.residues[index+1]
                                resn     = vresidue.name 
                                res1l    = three_letter_res_dict[resn]
                                residue  = SeqResidue(rname = resn, 
                                                    rcode = res1l, 
                                                    rnum  = index+1,
                                                    color = color,
                                                    chain = chain.name)
                                residue.vres = vresidue
                                sequence.append(residue)
                            
                            #if the residue is not on three_letter_res_dict
                            except:
                                vresidue = chain.residues[index+1]
                                resn     = vresidue.name 
                                rcode    = '%'
                                
                                if vresidue.is_solvent:
                                    pass
                                else:
                                    residue  = SeqResidue(rname = resn, 
                                                        rcode = rcode, 
                                                        rnum  = index+1,
                                                        color = color,
                                                        chain = chain.name)
                                    residue.vres = vresidue
                                    sequence.append(residue)
                        # if it is a 'gap' = '-'
                        else:
                            residue = SeqResidue(rname = 'gap', 
                                                rcode = '-', 
                                                rnum  = index,
                                                chain = chain.name)
                            sequence.append(residue)
                
                delete = True
                pos = -1
                #'''
                while delete:
                    if sequence[-1].rcode == '-':
                        sequence.pop(-1)
                    else:
                        delete = False
                #'''
                vobject.e_sequence = sequence  
            
            self.sequences.append(sequence)
            self.text_drawing_area.queue_draw()
        
    def add_new_3lettercode_sequence (self, inputseq):
        """ Function doc """
        sequence = []
        counter  = 0
        
        for rnum, res in enumerate( inputseq):
            try:
                res1l = three_letter_res_dict[res] 
                residue = SeqResidue(rname = res  , rcode = res1l, rnum = counter)
                sequence.append(residue)
                counter+=1
            
            except:
                pass
        
        self.sequences.append(sequence)
        self.text_drawing_area.queue_draw()
  
    def add_new_sequence (self, inputseq):
        """ Function doc """
        sequence = []
        counter  = 0

        for rnum, res in enumerate( inputseq):
            #print(rnum, res)
            if res == '-':
                pass
                res3l = one_letter_res_dict[res]
                residue = SeqResidue(rname = res3l , rcode = res, rnum = None)
                sequence.append(residue)
            
            else:
                res3l = one_letter_res_dict[res]                
                residue = SeqResidue(rname = res3l , rcode = res, rnum = counter)
                sequence.append(residue)
                counter+=1
        
        self.sequences.append(sequence)
        self.text_drawing_area.queue_draw()
    
    def clear_sequence_list (self):
        """ Function doc """
        self.sequences = []
    
    def refresh_vobject_sequence_list (self, vm_objects_dic):
        """ Function doc """
        self.clear_sequence_list()
        for index, vobject in vm_objects_dic.items():
            if vobject.active:
                self.add_sequence_from_vobject(vobject)
            else:
                pass
        
    def set_font_size (self, size = 20):
        """ Function doc """
        self.text_drawing_area.font_size = size
    
    def set_font_type (self, font_type = None):
        """ Function doc """
        self.text_drawing_area.font_type = font_type
    
    def set_canvas_width_and_height(self, width = 800, height = 800):
        """ Function doc """
        self.text_drawing_area.set_size_request(width, height)

    def get_canvas_xy(self, event):
        """ Function doc """
        char_width  = self.text_drawing_area.char_width
        char_height = self.text_drawing_area.char_height
        gap_height  = self.text_drawing_area.gap_height 
        
        (x, y) = event.x, event.y
        canvas_x, canvas_y = self.text_drawing_area.translate_coordinates(self.text_drawing_area, x, y)
        
        # Obter os objetos de ajuste horizontal e vertical do ScrolledWindow
        hadjustment = self.get_hadjustment()
        vadjustment = self.get_vadjustment()

        # Converter coordenadas da tela para coordenadas do widget
        canvas_x = x + hadjustment.get_value()
        canvas_y = y + vadjustment.get_value()
        
        
        if self.mode == 0:
            canvas_x,canvas_y= (int(canvas_x/char_width), 
                                int(canvas_y/(char_height + gap_height)/1 -2  )
                                )
        elif self.mode == 1:
            canvas_x,canvas_y= (int(canvas_x/(char_width*4)), 
                                int((canvas_y/(char_height + gap_height)  )/3)
                                )
        
        
        return canvas_x, canvas_y
        
    def on_motion(self, widget, event):
        '''(i/self.x_major_ticks)*self.deltaX + self.Xmin'''
        
        char_width  = self.text_drawing_area.char_width
        char_height = self.text_drawing_area.char_height
        gap_height  = self.text_drawing_area.gap_height
        
        self.motion_x = event.x
        self.motion_y = event.y
        
        canvas_x,canvas_y = self.get_canvas_xy(event)

        #if self.mode == 0:
        #    canvas_x,canvas_y= (int(canvas_x/char_width), 
        #                        int(canvas_y/(char_height + gap_height)  )
        #                        )
        #elif self.mode == 1:
        #    canvas_x,canvas_y= (int(canvas_x/(char_width*4)), 
        #                        int(canvas_y/(char_height + gap_height)  )
        #                        )
        try:
            pass
            if self.label:
                
                text = '{} , {} , {}'.format(self.sequences[canvas_y][canvas_x].chain,
                                             self.sequences[canvas_y][canvas_x].rname,
                                             self.sequences[canvas_y][canvas_x].rnum   
                                            )
                self.label.set_text(text)
            #print(canvas_x,canvas_y)
            #print(
            #      self.sequences[canvas_y][canvas_x].chain, 
            #      self.sequences[canvas_y][canvas_x].rname, 
            #      self.sequences[canvas_y][canvas_x].rnum)
        except:
            pass
            
        if self.is_button_1_pressed:
            
            #print(self.sequences[canvas_y][canvas_x].vres.name,
            #      self.sequences[canvas_y][canvas_x].rnum)
            '''
            if self.select_in_motion:
                if self.sequences[canvas_y][canvas_x].rcode == '-':
                    pass
                else:
                    self.sequences[canvas_y][canvas_x].set_active(True)
                    #atom_keys = list(self.sequences[canvas_y][canvas_x].vres.atoms.values())
                    #self.main.vm_session._selection_function_set({atom_keys[0]})
            else:
                self.sequences[canvas_y][canvas_x].set_active(False)
            #'''

            self.main.vm_session.vm_glcore.queue_draw()
            self.text_drawing_area.queue_draw()

    def button_press (self, button, event):
        """ Function doc """
        #print("here!")
        canvas_x, canvas_y= self.get_canvas_xy(event)
        
        if event.button == 1:
            self.is_button_1_pressed = True # mouse left button

            # Checking if the canvas position (y) is smaller than the size of the sequence list
            size = len(self.sequences)
            #print(canvas_x, canvas_y, size)
            if canvas_y <= size-1:
                atom_keys = list(self.sequences[canvas_y][canvas_x].vres.atoms.values())
                self.main.vm_session._selection_function_set({atom_keys[0]})
                self.main.vm_session.vm_glcore.queue_draw()
                
                if self.sequences[canvas_y][canvas_x].is_selected:
                    #print(canvas_x,canvas_y, self.sequences[canvas_y][canvas_x].rname)
                    self.sequences[canvas_y][canvas_x].set_active(False)
                    self.select_in_motion = False
                    
                    #print('if dentro',self.sequences[0][canvas_x].is_selected) 
                    #print('if dentro',self.sequences[1][canvas_x].is_selected)
                else:
                    #print('else dentro',self.sequences[0][canvas_x].is_selected) 
                    #print('else dentro',self.sequences[1][canvas_x].is_selected)
                    if self.sequences[canvas_y][canvas_x].rcode == '-':
                        pass
                    else:
                        self.sequences[canvas_y][canvas_x].set_active(True)
                        self.select_in_motion = True
                    
                        
                    
                #    print('else dentro',self.sequences[0][canvas_x].is_selected) 
                #    print('else dentro',self.sequences[1][canvas_x].is_selected)
                #    
                #print('fora1',self.sequences[0][canvas_x].is_selected) 
                #print('fora1',self.sequences[1][canvas_x].is_selected) 
            
            #self.text_drawing_area.queue_draw()
            self.text_drawing_area.queue_draw()
        if event.button == 2:

            
            self.is_button_2_pressed = True # mouse middle button
            size = len(self.sequences)
            if canvas_y <= size-1:

                frame = self.main.vm_session.get_frame ()
                res   = self.sequences[canvas_y][canvas_x].vres
                res.get_center_of_mass(frame = frame)
                atoms = list(res.atoms.values())
                #print(atoms)
                #atom_keys = list(res.atoms.values())
                self.main.vm_session.vm_glcore.center_on_atom(atoms[0])
                #self.main.vm_session._selection_function_set({atom_keys[0]})
                self.main.vm_session.vm_glcore.queue_draw()
        
        
        if event.button == 3:
            self.is_button_3_pressed = True # mouse right button
            
    def button_release (self, button, event):
        """ Function doc """
        char_width  = self.text_drawing_area.char_width
        char_height = self.text_drawing_area.char_height
        gap_height  = self.text_drawing_area.gap_height 
        
        if event.button == 1:
            #print(self.click_x, self.motion_x, self.click_y,self.motion_y)
            '''
            if self.click_x == self.motion_x and self.click_y == self.motion_y:
                canvas_x,canvas_y= self.get_canvas_xy(event)
                canvas_x,canvas_y= (int(canvas_x/char_width), 
                                    int(canvas_y/(char_height + gap_height)  )
                                    )
                
                if self.sequences[canvas_y][canvas_x].is_selected:
                    self.sequences[canvas_y][canvas_x].is_selected = False
                else:
                    self.sequences[canvas_y][canvas_x].is_selected = True
                self.text_drawing_area.queue_draw()
            '''
            self.is_button_1_pressed = False # mouse left button
        
        if event.button == 2:          
            self.is_button_2_pressed = False # mouse middle button
            
            
            
            
        if event.button == 3:          
            self.is_button_3_pressed = False # mouse right button
    
    def insert_sequence_from_vobject (self, vobject, chain):
        """ Function doc """
        
        
    
            
if __name__ == '__main__':
    #'''        
    main  = Gtk.Window()
    main.set_default_size(400, 200)
    seqview =  GtkSequenceViewer(main = main)
    notebook = Gtk.Notebook()
    tab_label = Gtk.Label("sequence")
    notebook.append_page(seqview, tab_label)



    #seqview.text_drawing_area.add_events(Gdk.EventMask.POINTER_MOTION_MASK)
    #seqview.text_drawing_area.add_events(Gdk.EventMask.BUTTON_PRESS_MASK)
    #seqview.text_drawing_area.add_events(Gdk.EventMask.BUTTON_RELEASE_MASK)
    #seqview.text_drawing_area.add_events(Gdk.EventMask.KEY_PRESS_MASK)
    #seqview.text_drawing_area.add_events(Gdk.EventMask.KEY_RELEASE_MASK)
    #seqview.text_drawing_area.add_events(Gdk.EventMask.SCROLL_MASK)

    #-----------------------------------------------------------------------
    seqview.add_events(Gdk.EventMask.POINTER_MOTION_MASK)
    seqview.add_events(Gdk.EventMask.BUTTON_PRESS_MASK)
    seqview.add_events(Gdk.EventMask.BUTTON_RELEASE_MASK)

    seqview.connect("motion-notify-event" , seqview.on_motion)
    seqview.connect("button_press_event"  , seqview.button_press)
    seqview.connect("button-release-event", seqview.button_release)
    #-----------------------------------------------------------------------
    
    sequence  ='---MISYLASIFLLATVSA-VPSGRVEVVFPSVETSRSGVK--TVKFTARVEVVFPSVETSRSGVK--TVKFTA' 
    sequence2 ='---MISYLASIFLLTTTTT-VPSGRVEVVFPSVETSRSGVK--TVKFTA-------------------------' 
    sequence3 ='---MISYLASIFLLTTTTT-VPSGRVEVVFPSVETSRSGVK--TVKFTA-------TVKFTA------------' 
    #print(sequence)
    seqview.add_new_sequence(sequence)
    seqview.add_new_sequence(sequence2)
    seqview.add_new_sequence(sequence3)
    print(seqview.sequences)
    
    
    main.add(notebook) 
    main.show_all()
    Gtk.main()
    #print ('uhuuu')
    #'''
