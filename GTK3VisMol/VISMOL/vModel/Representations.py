import numpy as np
import math
import ctypes
from OpenGL import GL
import time
import VISMOL.glCore.sphere_data as sphd
import VISMOL.glCore.cylinder_data as cyd
import VISMOL.glCore.matrix_operations as mop
#from   VISMOL.glCore.sphere_representation import _create_frame_sphere_data

#import EDTSurf.edtsurf as  edtsurf
'''
def _create_frame_sphere_data (frame, atoms ,offset, elems, scale, level,qtty ):
    """ Function doc """
    coords  = sphd.sphere_vertices[level]*qtty
    centers = sphd.sphere_vertices[level]*qtty
    for a, atom in enumerate(atoms):
        pos = atom.coords (frame)
        centers[a*offset:(a+1)*offset] = [pos[0],pos[1],pos[2]]*elems
        
        for i in range(elems):
            coords[a*offset+i*3]   *= atom.radius * scale
            coords[a*offset+i*3+1] *= atom.radius * scale
            coords[a*offset+i*3+2] *= atom.radius * scale
            coords[a*offset+i*3]   += pos[0]
            coords[a*offset+i*3+1] += pos[1]
            coords[a*offset+i*3+2] += pos[2]
    
    coords  = np.array(coords, dtype=np.float32)
    centers = np.array(centers, dtype=np.float32)
    
    return coords, centers
    
    
    #self.centers_list.append(self.centers)
    #self.frames.append(self.coords)
'''










class Representation:
    """ Class doc """

    def _make_gl_VAO (self):
        """ Function doc """
        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)
        return vao

    def _make_gl_index_buffer(self, indices):
        """ Function doc """
        ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)
        return ind_vbo

    def _make_gl_coord_buffer(self, coords, program):
        """ Function doc """
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        att_position = GL.glGetAttribLocation(program, 'vert_coord')
        GL.glEnableVertexAttribArray(att_position)
        GL.glVertexAttribPointer(att_position, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        #GL.glDisableVertexAttribArray(att_position)
        return coord_vbo

    def _make_gl_color_buffer(self, colors, program):
        """ Function doc """
        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(program, 'vert_color')
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
        #GL.glDisableVertexAttribArray(att_colors)
        return col_vbo

    def _make_gl_size_buffer (self, dot_sizes, program):
        """ Function doc """
        size_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, size_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, dot_sizes.nbytes, dot_sizes, GL.GL_STATIC_DRAW)
        att_size = GL.glGetAttribLocation(program, 'vert_dot_size')
        GL.glEnableVertexAttribArray(att_size)
        GL.glVertexAttribPointer(att_size, 1, GL.GL_FLOAT, GL.GL_FALSE, dot_sizes.itemsize, ctypes.c_void_p(0))
        #GL.glDisableVertexAttribArray(att_size)
        return size_vbo


    def _set_coordinates_to_buffer (self, coord_vbo = True, sel_coord_vbo = True):
        '''Esta fun????o atribui  as coordenadas que dele ser desenhadas  pela fun????o  "draw_representation"'''
        
        frame = self.glCore._safe_frame_exchange(self.visObj)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, visObj.line_buffers[1])
        
        if coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 
                            self.coord_vbo    )
            
            GL.glBufferData(GL.GL_ARRAY_BUFFER, 
                            frame.nbytes      ,
                            frame             , 
                            GL.GL_STATIC_DRAW)   
        else:
            pass
        
        
        
        if sel_coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 
                            self.sel_coord_vbo )
            
            
            GL.glBufferData(GL.GL_ARRAY_BUFFER, 
                            frame.nbytes      ,
                            frame             , 
                            GL.GL_STATIC_DRAW)    

        else:
            pass
            
    def _make_gl_representation_vao_and_vbos (self, 
                                              indices    = None,
                                              coords     = None,
                                              colors     = None,
                                              dot_sizes  = None,
                                              ):
        """ Function doc """
        print ('building', self.name,' VAO  and VBOs')    
        self.vao        =   self._make_gl_VAO()
        self.ind_vbo    =   self._make_gl_index_buffer( indices                        )
        self.coord_vbo  =   self._make_gl_coord_buffer( coords   , self.shader_program )
        self.col_vbo    =   self._make_gl_color_buffer( colors   , self.shader_program )
        if dot_sizes is not None:
            self.sel_size_vbo   =   self._make_gl_size_buffer ( dot_sizes , self.sel_shader_program )
        else:
            pass
        
    
    def _make_gl_sel_representation_vao_and_vbos (self, 
                                                  indices    = None,
                                                  coords     = None,
                                                  colors     = None,
                                                  dot_sizes  = None,
                                                  ):
        """ Function doc """
        print ('building', self.name,'background selection  VAO  and VBOs')    
        self.sel_vao        =   self._make_gl_VAO()
        self.sel_ind_vbo    =   self._make_gl_index_buffer( indices                             )
        self.sel_coord_vbo  =   self._make_gl_coord_buffer( coords    , self.sel_shader_program )
        self.sel_col_vbo    =   self._make_gl_color_buffer( colors    , self.sel_shader_program )
        if dot_sizes is not None:
            self.sel_size_vbo   =   self._make_gl_size_buffer ( dot_sizes , self.sel_shader_program )
        else:
            pass
        
        

    def _check_VAO_and_VBOs (self, indices = None):
        """ Function doc """
        if self.sel_vao is None:
            self._make_gl_vao_and_vbos ()
        else:
            pass
    
    def _enable_anti_alis_to_lines (self):
        """ Function doc """
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glEnable(GL.GL_LINE_SMOOTH)
        GL.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST)

    
    def define_new_indices_to_VBO ( self, input_indices = []):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        indices = input_indices
        indices = np.array(indices,dtype=np.uint32)
        #ind_vbo = self.ind_vbo
        
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)
        
        #ind_vbo = self.sel_ind_vbo
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.sel_ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)

    def change_vbo_colors  (self,  colors = []):
        """ Function doc """
        colors = np.array(colors,dtype=np.float32)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(self.shader_program, 'vert_color')
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))






class LinesRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'lines', active = True, _type = 'mol', visObj = None, glCore = None):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore

        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None


    def _make_gl_vao_and_vbos (self, indices = None):
        """ Function doc """
        
        #if indices is not None:
        #    print ('_make_gl_vao_and_vbos',indices)
        #    pass
        #else:
        #    print ('_make_gl_vao_and_vbos',indices)
        #    indices = np.array(self.visObj.index_bonds,dtype=np.uint32)

        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        indices = np.array(self.visObj.index_bonds,dtype=np.uint32)
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors
        #colors  = self.visObj.colors_rainbow

        self._make_gl_representation_vao_and_vbos (indices    = indices,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indices
        self._make_gl_sel_representation_vao_and_vbos (indices    = indices    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )



    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()
        GL.glUseProgram(self.shader_program)
        LineWidth = (80/abs(self.glCore.dist_cam_zrp)/2)**0.5  #40/abs(self.glCore.dist_cam_zrp)
        GL.glLineWidth(LineWidth)


        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*2), GL.GL_UNSIGNED_INT, None)

        GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        #if self.sel_vao is None:
        #    print ('_make_gl_vao_and_vbos')    
        #    self._make_gl_vao_and_vbos ()
        #else:
        #    pass
        line_width = self.visObj.vismol_session.vConfig.gl_parameters['line_width_selection'] 
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width)

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)

            #frame = self.glCore._safe_frame_exchange(self.visObj)
            #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_coord_vbo)
            #
            #GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes,
            #                frame, 
            #                GL.GL_STATIC_DRAW)              

            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*2), GL.GL_UNSIGNED_INT, None)  
        GL.glBindVertexArray(0)
        GL.glLineWidth(1)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)



class SticksRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'sticks', active = True, _type = 'mol', visObj = None, glCore = None,  indices = [] ):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore

        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        
        if indices == []:
            self.indices = np.array(self.visObj.index_bonds, dtype=np.uint32)
        else:
            self.indices = np.array(indices, dtype=np.uint32)
            
    def _make_gl_vao_and_vbos (self, indices = []):
        """ Function doc """
        #if indices == []:
        #    self.indices = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #else:
        #    self.indices = np.array(indices, dtype=np.uint32)
        
        #self.indices = np.array([0,1,0,2,1,2], dtype=np.uint32)
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        
        #indices = np.array(self.visObj.index_bonds, dtype=np.uint32)
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors

        self._make_gl_representation_vao_and_vbos (indices    = self.indices,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indices
        self._make_gl_sel_representation_vao_and_vbos (indices    = self.indices    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()

        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))

        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        self.glCore.load_lights(self.shader_program)

        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()

        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(20)
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_DEPTH_TEST)




class RibbonsRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'ribbon', active = True, _type = 'mol', visObj = None, glCore = None):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        
        
        if self.visObj.c_alpha_bonds == []:
            self.active  = False
        else:
            pass
        
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None


    def _make_gl_vao_and_vbos (self, indices = None):
        """ Function doc """
        
        #if indices is not None:
        #    print ('_make_gl_vao_and_vbos',indices)
        #    pass
        #else:
        #    print ('_make_gl_vao_and_vbos',indices)
        #    indices = np.array(self.visObj.index_bonds,dtype=np.uint32)

        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        indices = []
        
        # c_alpha_bonds
        for bond in self.visObj.c_alpha_bonds:
            i = bond.atom_index_i
            j = bond.atom_index_j
            #print (i, j)
            indices.append(i)
            indices.append(j)
            
        if indices == []:
            self.activate = False
        else:
            
            indices = np.array(indices,dtype=np.uint32)

            coords  = self.visObj.frames[0]
            #colors  = self.visObj.colors
            colors  = self.visObj.colors_rainbow

            self._make_gl_representation_vao_and_vbos (indices    = indices,
                                                       coords     = coords ,
                                                       colors     = colors ,
                                                       dot_sizes  = None   ,
                                                       )
            colors_idx = self.visObj.color_indices
            self._make_gl_sel_representation_vao_and_vbos (indices    = indices    ,
                                                           coords     = coords     ,
                                                           colors     = colors_idx ,
                                                           dot_sizes  = None       ,
                                                           )



    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()
        GL.glUseProgram(self.shader_program)
        
        ribbon_width = self.visObj.vismol_session.vConfig.gl_parameters['ribbon_width']
        LineWidth = ((ribbon_width*20)/abs(self.glCore.dist_cam_zrp)/2)  #40/abs(self.glCore.dist_cam_zrp)
        #print(LineWidth)
        GL.glLineWidth(LineWidth)


        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*2), GL.GL_UNSIGNED_INT, None)

        GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        #if self.sel_vao is None:
        #    print ('_make_gl_vao_and_vbos')    
        #    self._make_gl_vao_and_vbos ()
        #else:
        #    pass
        line_width = self.visObj.vismol_session.vConfig.gl_parameters['line_width_selection'] 
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width)

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)

            #frame = self.glCore._safe_frame_exchange(self.visObj)
            #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_coord_vbo)
            #
            #GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes,
            #                frame, 
            #                GL.GL_STATIC_DRAW)              

            GL.glDrawElements(GL.GL_LINES, int(len(self.visObj.index_bonds)*2), GL.GL_UNSIGNED_INT, None)  
        GL.glBindVertexArray(0)
        GL.glLineWidth(1)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)




class NonBondedRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'nonbonded', active = True, _type = 'mol', visObj = None, glCore = None):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore

        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None


    def _make_gl_vao_and_vbos (self, indices = None):
        """ Function doc """
        if indices is not None:
            pass
        else:
            indices = np.array(self.visObj.non_bonded_atoms, dtype=np.uint32)
        
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        #indices = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #indices = np.array(self.visObj.non_bonded_atoms, dtype=np.uint32)
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors

        self._make_gl_representation_vao_and_vbos (indices    = indices,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indices
        self._make_gl_sel_representation_vao_and_vbos (indices    = indices    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()


        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))

        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.non_bonded_atoms)), GL.GL_UNSIGNED_INT, None)

        GL.glBindVertexArray(0)
        GL.glLineWidth(1)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(20)

        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.non_bonded_atoms)), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glLineWidth(1)
        GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)



class DotsRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'dots', active = True, _type = 'mol', visObj = None, glCore = None, indices = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        
        

        
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        
        if indices == []:
            self.indices = np.array(range(len(self.visObj.atoms)), dtype=np.uint32)
        else:
            self.indices = np.array(indices, dtype=np.uint32)

    def _make_gl_vao_and_vbos (self, indices = None):
        """ Function doc """
        #if indices is not None:
        #    pass
        #else:
        
        #dot_qtty  = int(len(self.visObj.frames[0])/3)
        #indices = []
        #for i in range(dot_qtty):
        #    indices.append(i)
        
        indices = self.indices       
        
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        #indices = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #indices = np.array(self.visObj.idex, dtype=np.uint32)
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors


        self._make_gl_representation_vao_and_vbos (indices    = indices,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indices
        self._make_gl_sel_representation_vao_and_vbos (indices    = indices    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        self._enable_anti_alis_to_lines()
        #print ('DotsRepresentation')
        height = self.visObj.vismol_session.glwidget.vm_widget.height

        GL.glUseProgram(self.shader_program)
        #1*self.height dot_size
        #GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))
        GL.glPointSize(0.3*height/abs(self.glCore.dist_cam_zrp)) # dot size not included yet
        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)

        #GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        #GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        GL.glPointSize(200/abs(self.glCore.dist_cam_zrp))
        #GL.glLineWidth(20)
        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)

        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            #print(self.name,'draw_background_sel_representation')
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)



class SpheresRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'spheres', 
                      active = True, 
                       _type = 'mol', 
                      visObj = None, 
                      glCore = None,  
                       atoms = None,
                       #level = 'level_1', 
                       ):
        
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        
        # --------------------------------
        #self.level              = level
        self.level              = self.visObj.vismol_session.vConfig.gl_parameters['sphere_quality']
        self.scale              = self.visObj.vismol_session.vConfig.gl_parameters['sphere_scale']
        
        
        if atoms is None:
            self.atoms          = self.visObj.atoms
        else:
            self.atoms          = atoms
        
        
        self.coords             = None
        self.colors             = None
        self.centers            = None
        self.indices            = None
        self.triangles          = None
        self.frames             = []
        # --------------------------------

        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
     
    
    def _update_sphere_data_to_VBOs (self):
        """ Function doc """
        
        #GL.glDeleteVertexArrays( 1, self.vao)        
        #GL.glDeleteBuffers(1, self.ind_vbo)
        #GL.glDeleteBuffers(1, self.coord_vbo)
        #GL.glDeleteBuffers(1, self.centr_vbo)
        #GL.glDeleteBuffers(1, self.col_vbo)
        #self._make_gl_vao_and_vbos ()
        
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indices.nbytes, self.indices, GL.GL_DYNAMIC_DRAW)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.nbytes, self.coords, GL.GL_STATIC_DRAW)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.colors.itemsize*len(self.colors), self.colors, GL.GL_STATIC_DRAW)

    
    def _create_sphere_data(self ):
        
        """ Function doc """
        init = time.time()
        #cdef Py_ssize_t a, i, qtty, elems, offset, inds_e
        qtty = int(len( self.atoms))
        nucleus = [0.0, 0.0, 0.0]*qtty
        colores = [0.0, 0.0, 0.0]*qtty
        coords  = sphd.sphere_vertices[self.level]*qtty
        centers = sphd.sphere_vertices[self.level]*qtty
        colors  = sphd.sphere_vertices[self.level]*qtty
        indices = np.array(sphd.sphere_triangles[self.level]*qtty, dtype=np.uint32)
        elems  = int(len(sphd.sphere_vertices[self.level])/3)
        offset = int(len(sphd.sphere_vertices[self.level]))
        inds_e = int(len(sphd.sphere_triangles[self.level]))
        
        self.centers_list = []
        self.frames = []
        frame =0
        #for frame in range(len(self.visObj.frames)-1):
        for a, atom in enumerate( self.atoms ):
            pos = atom.coords (frame)
            #print (pos, atom.index, frame)

            colors[a*offset:(a+1)*offset]  = [atom.color[0],atom.color[1],atom.color[2]]*elems
            centers[a*offset:(a+1)*offset] = [pos[0],pos[1],pos[2]]*elems

            for i in range(elems):
                coords[a*offset+i*3]   *= atom.radius * self.scale
                coords[a*offset+i*3+1] *= atom.radius * self.scale
                coords[a*offset+i*3+2] *= atom.radius * self.scale
                coords[a*offset+i*3]   += pos[0]
                coords[a*offset+i*3+1] += pos[1]
                coords[a*offset+i*3+2] += pos[2]
            indices[a*inds_e:(a+1)*inds_e] += a*elems
        end = time.time()
        print('Time used creating nucleus, vertices and colors:', end-init)



        self.coords  = np.array(coords, dtype=np.float32)
        self.frames.append(self.coords)
        self.centers = np.array(centers, dtype=np.float32)
        self.centers_list.append(self.centers)
        self.colors  = np.array(colors, dtype=np.float32)
        self.indices = indices
        
        self.triangles = int(len(self.indices))

        init = time.time()

        if len(self.visObj.frames) > 1:
            '''
            import concurrent.futures
            with concurrent.futures.ProcessPoolExecutor() as executor:
                
                frames  = self.visObj.frames[1:]
                nframes = len(frames)
                
                atoms   = [self.atoms]*nframes
                offsets = [offset]*nframes
                Elems   = [elems]*nframes                
                
                scales  = [self.scale]*nframes
                levels  = [self.level]*nframes
                qttys   = [qtty]*nframes
                
                results = executor.map(_create_frame_sphere_data, frames, atoms, offsets, Elems, scales, levels, qttys)
            
            for result in results:
                self.frames.append(result[0])
                self.centers_list.append(result[1])
            
            #'''

            '''
            for frame in range(1,len(self.visObj.frames)-1):
                
                coords, centers = _create_frame_sphere_data (frame      , 
                                                             self.atoms ,
                                                             offset     , 
                                                             elems      , 
                                                             self.scale , 
                                                             self.level ,
                                                             qtty       )
                self.frames.append(coords)
                self.centers = np.array(centers, dtype=np.float32)
                self.centers_list.append(centers)
            
            
            
            
            
            #'''
            
            #'''
            for frame in range(1,len(self.visObj.frames)-1):
                coords  = sphd.sphere_vertices[self.level]*qtty
                centers = sphd.sphere_vertices[self.level]*qtty
                for a, atom in enumerate( self.atoms ):
                    pos = atom.coords (frame)
                    centers[a*offset:(a+1)*offset] = [pos[0],pos[1],pos[2]]*elems
                    
                    for i in range(elems):
                        coords[a*offset+i*3]   *= atom.radius * self.scale
                        coords[a*offset+i*3+1] *= atom.radius * self.scale
                        coords[a*offset+i*3+2] *= atom.radius * self.scale
                        coords[a*offset+i*3]   += pos[0]
                        coords[a*offset+i*3+1] += pos[1]
                        coords[a*offset+i*3+2] += pos[2]
                self.coords  = np.array(coords, dtype=np.float32)
                self.frames.append(self.coords)
                self.centers = np.array(centers, dtype=np.float32)
                self.centers_list.append(self.centers)
            #'''
        end = time.time()
        print('Time used creating nucleus, vertices and colors:', end-init)
    
        return True
    
    '''
    def _create_sel_sphere_data(self, level):
        """ Function doc """
        init = time.time()
        #cdef Py_ssize_t a, i, qtty, elems, offset, inds_e
        qtty = int(len(self.visObj.atoms))
        nucleus = [0.0, 0.0, 0.0]*qtty
        colores = [0.0, 0.0, 0.0]*qtty
        coords = sphd.sphere_vertices[level]*qtty
        colors = sphd.sphere_vertices[level]*qtty
        indices = np.array(sphd.sphere_triangles[level]*qtty, dtype=np.uint32)
        elems = int(len(sphd.sphere_vertices[level])/3)
        offset = int(len(sphd.sphere_vertices[level]))
        inds_e = int(len(sphd.sphere_triangles[level]))
        for a,atom in enumerate(self.visObj.atoms):
            colors[a*offset:(a+1)*offset] = [atom.color_id[0],atom.color_id[1],atom.color_id[2]]*elems
            for i in range(elems):
                coords[a*offset+i*3] *= atom.radius * self.scale
                coords[a*offset+i*3+1] *= atom.radius * self.scale
                coords[a*offset+i*3+2] *= atom.radius * self.scale
                coords[a*offset+i*3] += atom.pos[0]
                coords[a*offset+i*3+1] += atom.pos[1]
                coords[a*offset+i*3+2] += atom.pos[2]
            indices[a*inds_e:(a+1)*inds_e] += a*elems
        end = time.time()
        print('Time used creating nucleus, vertices and colors for selection:', end-init)
        self.sel_coords = np.array(coords, dtype=np.float32)
        self.sel_colors = np.array(colors, dtype=np.float32)
        self.sel_indices = indices
        return True
        '''
    '''
    def _make_gl_spheres(self, program):
        """ Function doc """
        vertex_array_object = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vertex_array_object)
        
        ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indices.itemsize*int(len(self.indices)), self.indices, GL.GL_DYNAMIC_DRAW)
        
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.itemsize*len(self.coords), self.coords, GL.GL_STATIC_DRAW)
        gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
        GL.glEnableVertexAttribArray(gl_coord)
        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.coords.itemsize, ctypes.c_void_p(0))
        
        centr_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, centr_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
        gl_center = GL.glGetAttribLocation(program, 'vert_centr')
        GL.glEnableVertexAttribArray(gl_center)
        GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.centers.itemsize, ctypes.c_void_p(0))
        
        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.colors.itemsize*len(self.colors), self.colors, GL.GL_STATIC_DRAW)
        gl_colors = GL.glGetAttribLocation(program, 'vert_color')
        GL.glEnableVertexAttribArray(gl_colors)
        GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.colors.itemsize, ctypes.c_void_p(0))
        
        GL.glBindVertexArray(0)
        GL.glDisableVertexAttribArray(gl_coord)
        GL.glDisableVertexAttribArray(gl_center)
        GL.glDisableVertexAttribArray(gl_colors)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
        self.spheres_vao = vertex_array_object
        self.spheres_buffers = (ind_vbo, coord_vbo, col_vbo)
        self.triangles = int(len(self.indices))
        return True
        '''
    '''
    def _make_sel_gl_spheres(self, program):
        """ Function doc """
        vertex_array_object = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vertex_array_object)
        
        ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.sel_indices.itemsize*int(len(self.sel_indices)), self.sel_indices, GL.GL_DYNAMIC_DRAW)
        
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.sel_coords.itemsize*len(self.sel_coords), self.sel_coords, GL.GL_STATIC_DRAW)
        gl_coord = GL.glGetAttribLocation(program, 'vert_coord')
        GL.glEnableVertexAttribArray(gl_coord)
        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, self.sel_coords.nbytes, ctypes.c_void_p(0))
        
        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.sel_colors.itemsize*len(self.sel_colors), self.sel_colors, GL.GL_STATIC_DRAW)
        gl_colors = GL.glGetAttribLocation(program, 'vert_color')
        GL.glEnableVertexAttribArray(gl_colors)
        GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, self.sel_colors.nbytes, ctypes.c_void_p(0))
        
        GL.glBindVertexArray(0)
        GL.glDisableVertexAttribArray(gl_coord)
        GL.glDisableVertexAttribArray(gl_colors)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
        self.sel_spheres_vao = vertex_array_object
        self.sel_spheres_buffers = (ind_vbo, coord_vbo, col_vbo)
        self.sel_triangles = int(len(self.indices))
        return True
    

    '''
    def _make_gl_vao_and_vbos (self):
        """ Function doc """
        #self._create_sphere_data()

        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        #indices = np.array(self.visObj.index_bonds,dtype=np.uint32)
        #coords  = self.visObj.frames[0]
        #colors  = self.visObj.colors

        #self._make_gl_representation_vao_and_vbos (indices    = indices,
        #                                           coords     = coords ,
        #                                           colors     = colors ,
        #                                           dot_sizes  = None   ,
        #                                           )

        self.vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(self.vao)

        self.ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indices.nbytes, self.indices, GL.GL_DYNAMIC_DRAW)
        
        #print (self.coords)
        self.coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.nbytes, self.coords, GL.GL_STATIC_DRAW)
        gl_coord = GL.glGetAttribLocation(self.shader_program , 'vert_coord')
        GL.glEnableVertexAttribArray(gl_coord)
        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.coords.itemsize, ctypes.c_void_p(0))
        
        self.centr_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)
        gl_center = GL.glGetAttribLocation(self.shader_program , 'vert_centr')
        GL.glEnableVertexAttribArray(gl_center)
        GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.centers.itemsize, ctypes.c_void_p(0))

        self.col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.colors.itemsize*len(self.colors), self.colors, GL.GL_STATIC_DRAW)
        gl_colors = GL.glGetAttribLocation(self.shader_program, 'vert_color')
        GL.glEnableVertexAttribArray(gl_colors)
        GL.glVertexAttribPointer(gl_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.colors.itemsize, ctypes.c_void_p(0))
        
        self.triangles = int(len(self.indices))
        
        
        
        colors_idx = self.visObj.color_indices
        self._make_gl_sel_representation_vao_and_vbos (indices    = self.indices    ,
                                                       coords     = self.coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        
        GL.glUseProgram          (self.shader_program )
        self.glCore.load_matrices(self.shader_program , self.visObj.model_mat)
        self.glCore.load_lights  (self.shader_program )
        self.glCore.load_fog     (self.shader_program )
        
        if self.vao is not None:
            GL.glBindVertexArray (self.vao)
            if self.glCore.modified_view:
                pass
            
            else:
                #self.centers = self.glCore._safe_frame_exchange(self.visObj)
                
                frame = self.glCore._get_visObj_frame (self.visObj)
                self.coords = self.frames[frame]
                self.centers =self.centers_list[frame]
                
                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
                GL.glBufferData(GL.GL_ARRAY_BUFFER, self.coords.itemsize*int(len(self.coords)), 
                                                    self.coords, GL.GL_STATIC_DRAW)
                

                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.centr_vbo)
                GL.glBufferData(GL.GL_ARRAY_BUFFER, self.centers.itemsize*len(self.centers), self.centers, GL.GL_STATIC_DRAW)



                GL.glDrawElements(GL.GL_TRIANGLES,  self.triangles , GL.GL_UNSIGNED_INT, None)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass















class GlumpyRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'glumpy', active = True, _type = 'mol', visObj = None, glCore = None, scale=1.0):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type
        self.visObj             = visObj
        self.glCore             = glCore
        self.scale              = scale

        # representation    
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None

        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None


    def _make_gl_vao_and_vbos (self, indices = None):
        """ Function doc """
        
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        coords  = self.visObj.frames[0]
        colors  = self.visObj.colors
        radii   = [self.scale] * len(self.visObj.frames[0])
        radii   = np.array(radii)
        # print ('radii', radii)

        dot_qtty  = int(len(coords)/3)
        indices = np.arange(dot_qtty, dtype=np.uint32)

        self._make_gl_representation_vao_and_vbos (indices    = indices,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None  ,
                                                   )
        colors_idx = self.visObj.color_indices
        self._make_gl_sel_representation_vao_and_vbos (indices    = indices    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None      ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs()
        self._enable_anti_alis_to_lines()
        GL.glUseProgram(self.shader_program)
        #GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        height = self.visObj.vismol_session.glwidget.vm_widget.height
        dist_cam_zrp = self.visObj.vismol_session.glwidget.vm_widget.dist_cam_zrp
        
        GL.glPointSize((50*height/(abs(dist_cam_zrp)))**0.5)
        #GL.glPointSize(55)
        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        if self.glCore.modified_view:
            pass
        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            self._set_coordinates_to_buffer (coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        if self.glCore.modified_view:
            pass
        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea
            '''
            #print(self.name,'draw_background_sel_representation')
            self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
            GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)









class SurfaceRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'surface', active = True, _type = 'mol', visObj = None, glCore = None, indices = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type

        self.visObj             = visObj
        self.glCore             = glCore
        
        

        
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.read_surface_data()
    
    def read_surface_data(self):
        """ Function doc """
        #from random import random 
        #
        #[verts, tris, verts_gpu, tris_gpu] = edtsurf.calc_surface('/home/fernando/programs/EasyHybrid3/Coords/pdbs/1bx4_H.pdb')
        #self.coords2  = verts_gpu
        #self.indices2 = tris_gpu
        #self.colors2  = []
        #
        #
        #size = len( self.coords2 )
        #for i in range(size):
        #    self.colors2.append(float(i/size) + random())
        
        rawdata = open('../EasyHybrid3/Coords/pdbs/1l2y.ply', 'r')
        lines  = rawdata.readlines()
        
        self.coords2 = []
        self.colors2 = []
        self.indices2 = []
        
        
        
        for line in lines:
            line2 = line.split()
            
            if len(line2) == 6:
                #print (line2)
                self.coords2.append(float(line2[0]))
                self.coords2.append(float(line2[1]))
                self.coords2.append(float(line2[2]))
                                                  
                self.colors2.append(float(line2[3])/255)
                self.colors2.append(float(line2[4])/255)
                self.colors2.append(float(line2[5])/255)
            
            if len(line2) == 7:
                
                self.indices2.append(int(line2[1]))
                self.indices2.append(int(line2[2]))
                self.indices2.append(int(line2[3]))
                
                
                
                

    def _make_gl_vao_and_vbos (self, indices = None):
        """ Function doc """
        #if indices is not None:
        #    pass
        #else:
        
        #dot_qtty  = int(len(self.visObj.frames[0])/3)
        #indices = []
        #for i in range(dot_qtty):
        #    indices.append(i)
        

        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        
        #indices = np.array(self.visObj.index_bonds, dtype=np.uint32)
        #indices = np.array(self.visObj.idex, dtype=np.uint32)

        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.array(self.colors2, dtype=np.float32)
        #indices = range(0, len(self.coords2))     
        #indices = np.array(indices, dtype=np.uint32)
        indices = np.array(self.indices2, dtype=np.uint32)
        #print (indices)


        self._make_gl_representation_vao_and_vbos (indices    = indices,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        
        #self.centr_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coords)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
        #gl_center = GL.glGetAttribLocation(self.shader_program , 'vert_centr')
        #GL.glEnableVertexAttribArray(gl_center)
        #GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        
        
        colors_idx = self.visObj.color_indices
        self._make_gl_sel_representation_vao_and_vbos (indices    = indices    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        #self._enable_anti_alis_to_lines()
        
        
        
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        view = self.glCore.glcamera.view_matrix
        
        GL.glUseProgram(self.shader_program )
        
        #print (self.visObj.model_mat,view)
        
        m_normal = np.array(np.matrix(np.dot(view, self.visObj.model_mat)).I.T)
        
        self.glCore.load_matrices(self.shader_program , self.visObj.model_mat)
        self.glCore.load_lights  (self.shader_program )
        self.glCore.load_fog     (self.shader_program )
        GL.glBindVertexArray(self.vao)

        
        
        
        
        
        
        
        
        '''
        #print ('DotsRepresentation')
        height = self.visObj.vismol_session.glwidget.vm_widget.height

        GL.glUseProgram(self.shader_program)
        #1*self.height dot_size
        #GL.glLineWidth(40/abs(self.glCore.dist_cam_zrp))
        GL.glPointSize(0.1*height/abs(self.glCore.dist_cam_zrp)) # dot size not included yet
        self.glCore.load_matrices(self.shader_program, self.visObj.model_mat)
        self.glCore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        '''
        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
            #GL.glDrawElements(GL.GL_POINTS, int(len(self.indices2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_LOOP, int(len(self.coords2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_STRIP, int(len(self.indices2)), GL.GL_UNSIGNED_INT, None)
            GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indices2)), GL.GL_UNSIGNED_INT, None)

        #GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        #GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass
        #self._check_VAO_and_VBOs ()
        #
        #GL.glEnable(GL.GL_DEPTH_TEST)
        #GL.glUseProgram(self.sel_shader_program)
        #GL.glPointSize(200/abs(self.glCore.dist_cam_zrp))
        ##GL.glLineWidth(20)
        #self.glCore.load_matrices(self.sel_shader_program, self.visObj.model_mat)
        #GL.glBindVertexArray(self.sel_vao)
        #
        #if self.glCore.modified_view:
        #    pass
        #
        #else:
        #    '''
        #    This function checks if the number of the called frame will not exceed 
        #    the limit of frames that each object has. Allowing two objects with 
        #    different trajectory sizes to be manipulated at the same time within the 
        #    glArea
        #    '''
        #    #print(self.name,'draw_background_sel_representation')
        #    self._set_coordinates_to_buffer (coord_vbo = False, sel_coord_vbo = True)
        #    GL.glDrawElements(GL.GL_POINTS, int(len(self.visObj.atoms)), GL.GL_UNSIGNED_INT, None)






class WiresRepresentation (Representation):
    """ Class doc """
    
    def __init__ (self, name = 'wires', active = True, _type = 'mol', visObj = None, glCore = None, indices = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = _type
        self.visObj             = visObj
        self.glCore             = glCore
        
        # representation    
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
        
        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None

        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.read_surface_data()
    
    def read_surface_data(self):
        """ Function doc """
        rawdata = open('../EasyHybrid3/Coords/pdbs/1l2y.ply', 'r')
        lines  = rawdata.readlines()
        
        self.coords2 = []
        self.colors2 = []
        self.indices2 = []
        
        for line in lines:
            line2 = line.split()
            if len(line2) == 6:
                self.coords2.append(float(line2[0]))
                self.coords2.append(float(line2[1]))
                self.coords2.append(float(line2[2]))
                self.colors2.append(float(line2[3])/255)
                self.colors2.append(float(line2[4])/255)
                self.colors2.append(float(line2[5])/255)
            if len(line2) == 7:
                self.indices2.append(int(line2[1]))
                self.indices2.append(int(line2[2]))
                self.indices2.append(int(line2[3]))

    def _make_gl_vao_and_vbos (self, indices = None):
        """ Function doc """
        self.shader_program     = self.glCore.shader_programs[self.name]
        self.sel_shader_program = self.glCore.shader_programs[self.name+'_sel']
        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.zeros(len(self.colors2))
        indices = np.array(self.indices2, dtype=np.uint32)
        self._make_gl_representation_vao_and_vbos (indices    = indices,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.visObj.color_indices
        self._make_gl_sel_representation_vao_and_vbos (indices    = indices    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_VAO_and_VBOs ()
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        view = self.glCore.glcamera.view_matrix
        GL.glUseProgram(self.shader_program )
        m_normal = np.array(np.matrix(np.dot(view, self.visObj.model_mat)).I.T)
        self.glCore.load_matrices(self.shader_program , self.visObj.model_mat)
        self.glCore.load_lights  (self.shader_program )
        self.glCore.load_fog     (self.shader_program )
        GL.glBindVertexArray(self.vao)
        if self.glCore.modified_view:
            pass

        else:
            '''
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea'''
            # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
            GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indices2)), GL.GL_UNSIGNED_INT, None)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass