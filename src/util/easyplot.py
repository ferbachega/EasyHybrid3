import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
import cairo
import random
import numpy as np
import sys

try:
    from util.colormaps import COLOR_MAPS  
except:
    from colormaps import COLOR_MAPS
from pprint import pprint

def interpolate_color(color1, color2, fraction):
    return [ ( (x1 + fraction * (x2 - x1))) for (x1, x2) in zip(color1, color2)]

def get_color(value, color_map):#, factor = 2):   
    
    if value == float('inf'):
        return [1,1,1]
    #value = (value+1)/factor   
    
    cutoffs = list(color_map.keys())
    color_list = list(color_map.values())
    counter = 0 
    
    #print(value)
    
    
    
    for cutoff, color in color_map.items():
        
        if value <= cutoff :
            
            fraction = (value -  cutoffs[counter-1]) / ( cutoff - cutoffs[counter-1])
            
            
            color1  = color_list[counter-1]
            color2  = color_list[counter]
            #
            #
            #
            #
            #
            color = interpolate_color(color1, color2, fraction)
            #color = [0,0,0]
            
            #color[0] = color[0] 
            #color[1] = color[1] 
            #color[2] = color[2] 
            #print(value, cutoff, fraction,color )

            return color
        else:
            pass
            #print(value, cutoff)#, fraction,color )
        counter +=1


def bilinear_interpolation(c1, c2, c3, c4, tx, ty):
    return (
        int((1 - tx) * (1 - ty) * c1[0] + tx * (1 - ty) * c2[0] + (1 - tx) * ty * c4[0] + tx * ty * c3[0]),
        int((1 - tx) * (1 - ty) * c1[1] + tx * (1 - ty) * c2[1] + (1 - tx) * ty * c4[1] + tx * ty * c3[1]),
        int((1 - tx) * (1 - ty) * c1[2] + tx * (1 - ty) * c2[2] + (1 - tx) * ty * c4[2] + tx * ty * c3[2]),
        int((1 - tx) * (1 - ty) * c1[3] + tx * (1 - ty) * c2[3] + (1 - tx) * ty * c4[3] + tx * ty * c3[3])
    )
 
 
def expand_array_size (norm_data):
    """ 
    this function receives an i x j matrix and transforms it into an 
    i+2 x j+2 matrix. the original matrix is ​​in the middle of the new 
    matrix, so the added elements are on the edges. they will be 
    important for the interpolation of points when the smooth function 
    is used in ImagePlot.
    
    
       extra column /extra row
         |
       |---|---|---|---|---|
     --| x | x | x | x | x |
       |---|---|---|---|---|
       | x |###|###|###| x |
       |---|---|---|---|---|
       | x |###|###|###| x |
       |---|---|---|---|---|
       | x | x | x | x | x |--
       |---|---|---|---|---|
                         |
                       extra column /extra row
    
    ### = elements in original matrix
    x   = new elements
    
    -It's not a very efficient algorithm, but it works well
    """
    size_x = len(norm_data)
    size_y = len(norm_data[0])
    #print('expanded matrix')
    bigmatrix = []#[None]*size_x+2
    for k in range (0, size_x+2):
        bigmatrix.append([])
        for l in range(0,  size_y+2):
            bigmatrix[k].append(0)
    #pprint (bigmatrix )    
    for k in range (0, size_x+2):    
        for l in range(0,  size_y+2):
           
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|XXX|   |   |   |   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |   |   |   |   |--
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            if   k == 0 and l == 0:
                bigmatrix[k][l] = norm_data[0][0]
            
            elif k == 0 and l == 1:
                bigmatrix[k][l] = norm_data[0][0]
            
            elif k == 1 and l == 0:
                bigmatrix[k][l] = norm_data[0][0]

            
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |   |   |   |   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |   |   |   |XXX|--
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif k == size_x+1   and l == size_y+1 :
                bigmatrix[k][l] = norm_data[size_x-1][size_y-1]
            
            elif k == size_x   and l == size_y+1 :
                bigmatrix[k][l] = norm_data[size_x-1][size_y-1]
            
            elif k == size_x+1    and l == size_y :
                bigmatrix[k][l] = norm_data[size_x-1][size_y-1]

                
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |XXX|XXX|XXX|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |   |   |   |   |--
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif k == 0 and   size_y+1 > l > 0  :
                bigmatrix[k][l] = norm_data[0][l-1]
            
            
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |   |   |   |   |
            #    |---|---|---|---|---|
            #    |XXX|###|###|###|   |
            #    |---|---|---|---|---|
            #    |XXX|###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |   |   |   |   |-- 
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif size_x+1 > k > 0 and l == 0:
                bigmatrix[k][l] = norm_data[k-1][0]
            
            
            #---------------------------------------------------
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |   |   |   |   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|xxx|
            #    |---|---|---|---|---|
            #    |   |###|###|###|xxx|
            #    |---|---|---|---|---|
            #    |   |   |   |   |   |-- 
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif size_x+1 > k > 0 and size_y+1 == l:
                bigmatrix[k][l] = norm_data[k-1][size_y-1]
            #---------------------------------------------------
               
            #---------------------------------------------------
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |   |   |   |   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |xxx|xxx|xxx|   |-- 
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif k == size_x+1 and size_y  > l> 0:
                bigmatrix[k][l] = norm_data[k-2][l-1]
            #---------------------------------------------------


            elif   k == 0        and l == size_y+1:
                bigmatrix[k][l] = norm_data[0][l-2]
            
            elif k == size_x+1 and l == 0:
                bigmatrix[k][l] = norm_data[k-2][0]
            
            else:
                bigmatrix[k][l] = norm_data[k-1][l-1]
    return bigmatrix
 
    
class ColorSquareBilinearInterpolation:
    def __init__(self, surface, x, y, size, c1, c2, c3, c4):
        self.ctx = cairo.Context(surface)
        self.x, self.y, self.size = x, y, size
        self.c1, self.c2, self.c3, self.c4 = c1, c2, c3, c4

    def draw(self):
        for i in range(self.size):
            for j in range(self.size):
                tx = j / (self.size - 1)
                ty = i / (self.size - 1)
                color = bilinear_interpolation(self.c1, self.c2, self.c3, self.c4, tx, ty)

                self.ctx.set_source_rgba(*[comp / 255.0 for comp in color])
                self.ctx.rectangle(self.x + j, self.y + i, 1, 1)
                self.ctx.fill()
                
                cr.rectangle(self.bx+i*self.factor_x -1,
                             self.by+j*self.factor_y -1, 
                             round(self.factor_x)+1, 
                             round(self.factor_y)+1)
                cr.fill()


class Canvas(Gtk.DrawingArea):
    """ Class doc """
    
    def __init__ (self, bg_color = [1,1,1], cmap =  'jet' ):
        """ Class initialiser """
        super().__init__( )
        #Gtk.DrawingArea.__init__(main)
        
        self.color_map = COLOR_MAPS[cmap]
        
        self.bg_color = None
        self.set_bg_color (bg_color[0], bg_color[1], bg_color[2])
        
        self.set_events(Gdk.EventMask.ALL_EVENTS_MASK)
        self.connect("draw", self.on_draw)

        
    def set_bg_color (self, r = 1 , g = 1 , b = 1):
        """ Function doc """
        self.bg_color = [r, g, b]
    
    
    def _make_bg (self, cr):
        """ Function doc """
        cr.set_source_rgb( self.bg_color[0],self.bg_color[1], self.bg_color[2])
        cr.paint()
    
    
    def get_i_and_j_from_click_event (self, event):
        """ Function doc """
        (x, y) = int(event.x), int(event.y)
        
        x_on_plot = int((x-self.bx)/self.factor_x)
        y_on_plot = int((y-self.by)/self.factor_y)
        
        return x_on_plot, y_on_plot, x, y

    
    def check_clicked_points(self, x_on_plot, y_on_plot):
        """ Function doc """
        
        if int(x_on_plot) < 0 or int(x_on_plot) >= self.size_x :
            print('canvas')
            self.points = []
            self.selected_dot = None
            return False
        elif int(y_on_plot) < 0 or int(y_on_plot) >= self.size_y :
            print('canvas')
            self.points = []
            self.selected_dot = None #set no None the redot at IJ matrix
            return False
        
        else:
            return True
            
            self.points.append((x_on_plot, y_on_plot))
        
        print("Mouse clicker at:",  x, y, int(x_on_plot), int(y_on_plot), 
                                    (x-self.bx)/self.factor_x, #/(self.x_final-self.bx), 
                                    (y-self.by)/self.factor_y) #/(self.y_final-self.by) ) #, x - self.x_final , y - self.y_final , self.x_final ,  self.y_final ,   (x-self.bx) /self.factor_x ,    (y-self.by)/self.factor_y )
        print(self.points)

    
    def on_mouse_button_press(self, widget, event):
        #(x, y) = int(event.x), int(event.y)
        ##x, y = device.get_position(widget)
        #print("Mouse moved to:", x, y)
        (x, y) = int(event.x), int(event.y)
        
        x_on_plot = int((x-self.bx)/self.factor_x)
        y_on_plot = int((y-self.by)/self.factor_y)
        

        
        if int(x_on_plot) < 0 or int(x_on_plot) >= self.size_x :
            #print('canvas2222')
            self.points = []
            self.selected_dot = None
        
        elif int(y_on_plot) < 0 or int(y_on_plot) >= self.size_y :
            #print('canvas1111')
            self.points = []
            self.selected_dot = None
        
        else:
            self.points.append((x_on_plot, y_on_plot))
        
        print("Mouse clicker at:",  x, y, int(x_on_plot), int(y_on_plot), 
                                    (x-self.bx)/self.factor_x, #/(self.x_final-self.bx), 
                                    (y-self.by)/self.factor_y) #/(self.y_final-self.by) ) #, x - self.x_final , y - self.y_final , self.x_final ,  self.y_final ,   (x-self.bx) /self.factor_x ,    (y-self.by)/self.factor_y )
        print(self.points)
        self.queue_draw()

    
    def on_motion(self, widget, event):
        (x, y) = int(event.x), int(event.y)
        x_on_plot = int((x-self.bx)/self.factor_x)
        y_on_plot = int((y-self.by)/self.factor_y)
        #print("Mouse moved to:", x, y, int(x_on_plot), int(y_on_plot))

    def get_pixel_rgb(self, x, y):
        # Converte os dados da superfície para um array de bytes
        buf = self.surface.get_data()

        # Cairo usa um formato ARGB32, onde cada pixel é representado por 4 bytes (A, R, G, B)
        array = np.frombuffer(buf, dtype=np.uint8)
        array = array.reshape((self.height, self.width, 4))  # Reshape para uma matriz 2D com 4 valores por pixel

        # Pegar o valor do pixel em x, y
        pixel = array[y, x]

        # O formato é ARGB, então extraímos apenas R, G e B
        a, r, g, b = pixel
        return (r, g, b)


class ImagePlot(Canvas):
    """ Class doc """
    
    def __init__ (self, data = None):
        """ Class initialiser 
        
        Now the class can receive an energy matrix "Z" with 
        infinite value "inf". To do this it was necessary 
        to modify the get_color method.
        
        """
        super().__init__( )

        self.bx = 50 
        self.by = 50 
        self.points = []
        
        self.data = data
        self.dataRC1 = None
        self.dataRC2 = None   
        
        self.cmap = 'jet'
        
        self.norm_data = None
        self.data_min  = None
        self.data_max  = None
        self.connect("motion-notify-event", self.on_motion)
        
        self.RC_label = None
        self.label_mode = 0
        self.is_discrete = True
        
        self.selected_dot = None #it is the i and j coordinates at the energey matrix
        self.sel_dot_rgb  = [1,0,0]
        
        
    def define_datanorm (self):
        """ Function doc 
        When you have a 2D umbrella sampling data, 
        some points have energy estimated with inf (infinity).
        
        In python, "float('inf') is allowed - accepted as a float.
        
        When the value in the matrix is ​​given by inf, IMAGEPLOT 
        should draw a white square (showing that there is no energy 
        calculated at the point)
        
        In this case it is necessary to temporarily replace the inf 
        values ​​with 0 (so that it is possible to extract the maximum 
        and minimum values ​​of the matrix)
        """
        
        
        temp_data = self.data.copy()
        a = len(temp_data)
        b = len(temp_data[0])
        inf = float('inf')
        
        #replacing the inf values ​​of the matrix  
        inf_list = []
        for i in range(a):
            for j in range(b):
                if temp_data[i][j] == inf:
                    temp_data[i][j] = 0
                    inf_list.append([i, j])
        #print(temp_data)
        self.norm_data = temp_data - np.min(temp_data)
        #self.norm_data = self.data - np.min(self.data)

        self.data_min = np.min(self.norm_data)
        self.data_max = np.max(self.norm_data)
        delta = (self.data_max - self.data_min)
        self.norm_data = self.norm_data/delta
        
        #resetting the matrix inf values
        for inf_item in inf_list:
            i = inf_item[0]
            j = inf_item[1]
            self.norm_data[i][j] = float('inf')
        
        return self.data_min, self.data_max, delta , self.norm_data
        
        
    def set_threshold_color (self, _min = 0, _max = None, cmap = 'jet'):
        """ Function doc """
        self.color_map = COLOR_MAPS[self.cmap]
        fraction = _max / self.data_max #number between 0 and 1
        
        new_colormap = {}
        
        for i in range(101):
            value = i/100
            
            color = get_color(value, self.color_map)
            new_colormap[fraction*value] = color
        
        last_color = get_color(1, self.color_map)
        new_colormap[1.1] = last_color
        #pprint(self.color_map)
        #pprint( new_colormap)
        self.color_map = new_colormap
    
        
    def set_cmap (self, cmap = 'jet'):
        """ Function doc """
        self.color_map = COLOR_MAPS[cmap]
    
    
    def set_label_mode (self, mode = 0):
        """ Function doc """
        self.label_mode = mode
        
    
    def draw_color_bar (self, cr, res, gap = 20, label_spacing = 2, font_size = 12, num_labels = 10):
        """ Function doc """
         
       
        cr.set_line_width (1.0)
        #print('\n\n',res,'\n\n')
        for j in range(0, res*10, 1):
            
            #--------------------------------------------------------------------------------------------
            #color bar
            #--------------------------------------------------------------------------------------------
            z = j  / float(res*10)
            color = get_color(z, self.color_map)
            cr.set_source_rgb( color[0], color[1], color[2])
            
            j = res*10 -j-1
            
            x = round(self.size_x * self.factor_x)+self.bx + gap
            y = self.by + j*self.factor_y/10
            
            size_x  = round(20)+1 
            size_y  = round((self.factor_y/10)+1)
            
            
            
            cr.rectangle(
                         x,
                         y, 
                         size_x, 
                         size_y)
            cr.fill()
        
        
        
        factor = self.size_y/num_labels
        if factor ==  0:
            factor = 1
        else:
            pass
        
        
        
        cr.set_source_rgb(0,0,0)
        _min = self.data_min
        _max = self.data_max
        factor2 = (_max - _min)/(num_labels-1)
        #print(_min, _max, factor, self.size_y, self.size_y/factor, self.factor_y)
        counter = 0
        
        '''
           The position of the labels depends on the number of elements 
        in y of the data matrix (self.size_y). It depends on the factor 
        in y (self.factor_y), which in turn takes into account the 
        number of squares in the matrix and the size of the window.
        
        
        '''
        for i in range(0, num_labels):    
            
            j = (num_labels-1)-i
            
            x1 = round(self.size_x * self.factor_x)+self.bx*1.2 +37
            y1 = self.by + (j*self.factor_y*(self.size_y/(num_labels-1)))
            
            cr.move_to (x1,y1)
            
            x2 = round(self.size_x * self.factor_x)+self.bx*1.2 +30  
            cr.line_to (x2,y1) 
            #--------------------------------------------------------------------------------------------
            cr.stroke ()
            
            # . text
            cr.set_font_size(font_size)
            text = '{:<6.3f}'.format(i*factor2)
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            
            x3 = round(self.size_x * self.factor_x)+self.bx*1.2 + 45 
            cr.move_to(x3,y1 + height/2)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1
        
        
        
        '''
        for j in range(0, self.size_y+1, factor ):    

            j = self.size_y -j-1
            
            #z = ((j  / float(self.size_y)) - 0.5)*2
            #print(z)
            cr.set_source_rgb(0,0,0)
            cr.move_to (round(self.size_x * self.factor_x)+self.bx*1.2 +37      , 
                              self.by + j*self.factor_y + self.factor_y  )
            
            
            cr.line_to (round(self.size_x * self.factor_x)+self.bx*1.2 +30  , 
                              self.by + j*self.factor_y + self.factor_y   ) 
            #--------------------------------------------------------------------------------------------
            
            cr.stroke ()
            #--------------------------------------------------------------------------------------------

            #--------------------------------------------------------------------------------------------
            # texto em y
            cr.set_font_size(font_size)
            #print(_min, factor2, counter, _min + factor2*counter)
            text = '{:<6.3f}'.format(_min + factor2*counter )
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            cr.move_to(round(self.size_x * self.factor_x)+self.bx*1.2 + 45 ,  
                             self.by + j*self.factor_y+ self.factor_y + height/2)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1
        '''







        cr.set_line_width (2.0)
        cr.rectangle(round(self.size_x * self.factor_x)+self.bx +gap,
                     self.by   , 
                     round(20+1), 
                     round(self.size_y * self.factor_y))
        cr.stroke ()


    def draw_discrete_square (self, i, j, cr):
        """ Function doc """
        #z is the nergey value
        z = self.norm_data[i][j]
        color = get_color(z, self.color_map)
        
        
        if color:
            cr.set_source_rgb( color[0], color[1], color[2]   )
        else:
            #.Associates the color white when color is null 
            #.(due to matrix points with value inf)
            cr.set_source_rgb( 1, 1, 1   )
        
        cr.rectangle(self.bx+i*self.factor_x -1,
                     self.by+j*self.factor_y -1, 
                     round(self.factor_x)+1, 
                     round(self.factor_y)+1)
        cr.fill()
    
    
    def draw_smooth_square (self, cr = None, square = None, colors = None):
        """ Function doc 
        
        
        """
        x0 , y0 = square[0] # is tuple of two coordinates
        x05, y0 = square[1] # is tuple of two coordinates
        x05,y05 = square[2] # is tuple of two coordinates
        x0 ,y05 = square[3] # is tuple of two coordinates
        
        pattern = cairo.MeshPattern()
        # Começa a definição do patch de padrão de malha
        pattern.begin_patch()
        pattern.move_to(x0 , y0)
        pattern.line_to(x05, y0)
        pattern.line_to(x05,y05)
        pattern.line_to(x0 ,y05)
        
        color1 = colors[0]
        color2 = colors[1]
        color3 = colors[2]
        color4 = colors[3]
        
        ### Define os vértices do patch
        ##pattern.move_to(i,     j)
        ##pattern.line_to(i+100, j)
        ##pattern.line_to(i+100, j+100)
        ##pattern.line_to(i    , j+100)

        #color1 = get_color(z1, self.color_map)
        #color2 = get_color(z2, self.color_map)
        #color3 = get_color(z3, self.color_map)
        #color4 = get_color(z4, self.color_map)
        #
        pattern.set_corner_color_rgb(0, color1[0], color1[1], color1[2]) 
        pattern.set_corner_color_rgb(1, color2[0], color2[1], color2[2]) 
        pattern.set_corner_color_rgb(2, color3[0], color3[1], color3[2]) 
        pattern.set_corner_color_rgb(3, color4[0], color4[1], color4[2]) 
        
        # Finaliza a definição do patch de padrão de malha
        #cr.scale(1, 1)
        pattern.end_patch()
        cr.set_source(pattern)
        cr.paint()          

    
    def draw_inter_square (self, surface, x, y, size, c1, c2, c3, c4, cr):
        """ Function doc """
        #self.ctx = cr
        #self.x, self.y, self.size = x, y, size
        #self.c1, self.c2, self.c3, self.c4 = c1, c2, c3, c4

        for i in range(size):
            for j in range(size):
                tx = j / (size - 1)
                ty = i / (size - 1)
                color = bilinear_interpolation(c1, c2, c3, c4, tx, ty)

                cr.set_source_rgba(*[comp / 255.0 for comp in color])
                cr.rectangle(x + j, y + i, 1, 1)
                cr.fill()


    def draw_image (self, cr):
        """ Function doc """
        if self.is_discrete == True:
            #discrete square
            #'''
            n = 0 
            for i in range(0,self.size_x, ):
                for j in range(0,self.size_y):
                    self.draw_discrete_square (i, j, cr)


        else:
            
            data = expand_array_size(self.norm_data)
        
            for i in range(1, self.size_x+1):
                for j in range(1,  self.size_y+1):
                    
                    x0, y0 = self.bx+(i-1)*self.factor_x -1   , self.by+(j-1)*self.factor_y -1
                    x1, y1 = self.bx+(i )*self.factor_x -1 , self.by+ (j )*self.factor_y -1
                    
                    x05   = (x1 + x0)/2
                    y05   = (y1 + y0)/2
                    '''
                    print(
                    x0 , y0   ,
                    x0 , y05 ,
                    x0 , y1   ,

                    x05 , y0  ,
                    x05 , y05,
                    x05 , y1  ,

                    x1 , y0   ,
                    x1 , y05 ,
                    x1 , y1   ,
                    )
                    '''
                    # calculating energies / avarege 
                    E_x05_y05 = data[i][j]
                    C_x05_y05 = get_color(E_x05_y05, self.color_map)
                    
                    E_x0_y0 = (data[i-1][j ] + data[i-1][j-1] + data[i ][j-1] +  data[i][j])/4
                    C_x0_y0 = get_color(E_x0_y0, self.color_map)
                    
                    E_x1_y0 = (data[i][j-1] + data[i+1][j] + data[i+1][j-1] +  data[i][j])/4
                    C_x1_y0 = get_color(E_x1_y0, self.color_map)
                    
                    E_x0_y1 = (data[i-1][j+1]+ data[i-1][j]+ data[i][j+1] +  data[i][j])/4
                    C_x0_y1 = get_color(E_x0_y1, self.color_map)
                    
                    E_x05_y1 = (data[i][j+1] +  data[i][j])/2
                    C_x05_y1 = get_color(E_x05_y1, self.color_map)
                    
                    E_x05_y0 = (data[i][j-1] +  data[i][j])/2
                    C_x05_y0 = get_color(E_x05_y0, self.color_map)
                    
                    E_x0_y05 = (data[i-1][j] +  data[i][j])/2
                    C_x0_y05 = get_color(E_x0_y05, self.color_map)
                    
                    
                    E_x1_y05 = (data[i+1][j] +  data[i][j])/2
                    C_x1_y05 = get_color(E_x1_y05, self.color_map)
                    
                    E_x1_y1 = (data[i+1][j+1] +data[i][j+1] +data[i+1][j] +  data[i ][j ])/4
                    C_x1_y1 = get_color(E_x1_y1, self.color_map)
                    #C_x1_y1 = [0 ,0 ,0]
                    
                    
                    
                    '''each square is actually 4 squares that will be drawn next.'''
                    #---------------------------------------------------
                    square = [(x0 , y0),   # 
                              (x05, y0),   #  o------o------o
                              (x05,y05),   #  |XXXXXX|      |
                              (x0 ,y05)]   #  |XXXXXX|      |
                                           #  o------o------o
                    colors = [C_x0_y0  ,   #  |      |      |
                              C_x05_y0 ,   #  |      |      |
                              C_x05_y05,   #  o------o------o
                              C_x0_y05 ]   #
                    
                    self.draw_smooth_square(cr = cr, square = square, colors = colors)
                    #---------------------------------------------------

            
                    #---------------------------------------------------
                    square = [(x05 , y0),   # 
                              (x1  , y0),   #  o------o------o
                              (x1  , y05),  #  |      |      |
                              (x05 , y05)]  #  |      |      |
                                            #  o------o------o
                    colors = [C_x05_y0,     #  |XXXXXX|      |
                              C_x1_y0,      #  |XXXXXX|      |
                              C_x1_y05,     #  o------o------o
                              C_x05_y05]    #

                    self.draw_smooth_square(cr = cr, square = square, colors = colors)
                    #---------------------------------------------------
                    
                    
                    #---------------------------------------------------
                    square = [(x0  , y05),  # 
                              (x05 ,y05),   #  o------o------o
                              (x05 , y1 ),  #  |      |XXXXXX|
                              (x0  , y1)]   #  |      |XXXXXX|
                                            #  o------o------o
                    colors = [C_x0_y05,     #  |      |      |
                              C_x05_y05,    #  |      |      |
                              C_x05_y1,     #  o------o------o
                              C_x0_y1]      #

                    self.draw_smooth_square(cr = cr, square = square, colors = colors)
                    #---------------------------------------------------

                    
                    #---------------------------------------------------
                    square = [(x05 , y05),  # 
                              (x1  ,y05),   #  o------o------o
                              (x1  , y1 ),  #  |      |      |
                              (x05  , y1)]  #  |      |      |
                                            #  o------o------o
                    colors = [C_x05_y05,    #  |      |XXXXXX|
                              C_x1_y05,     #  |      |XXXXXX|
                              C_x1_y1,      #  o------o------o
                              C_x05_y1]     #

                    self.draw_smooth_square(cr = cr, square = square, colors = colors)
                    #---------------------------------------------------
                             
    def draw_image_box (self, cr, line_width = 1, color = [0,0,0]):
        cr.set_line_width (line_width)
        cr.set_source_rgb(color[0], color[1], color[2])
        cr.rectangle(self.bx,
                     self.by, 
                     round(self.size_x * self.factor_x), 
                     round(self.size_y * self.factor_y))
        cr.stroke ()
        
    
    def draw_scale (self, cr, font_size = 12, x_major_ticks = 3 , y_major_ticks = 3 ):
        """ Function doc """        
        
        round(self.size_x/x_major_ticks)
        
        for i in range(0,self.size_x, x_major_ticks):
            #--------------------------------------------------------------------------------------------
            # marcacoes em x
            #print(self.bx+ i*self.factor_x + (self.factor_x)/2.0 , self.bx+ i*self.factor_x + (self.factor_x)/2.0)
            cr.move_to (self.bx+ i*self.factor_x + (self.factor_x)/2.0      , self.by )
            cr.line_to (self.bx+ i*self.factor_x + (self.factor_x)/2.0      , self.by -10 ) 
            #--------------------------------------------------------------------------------------------
            cr.stroke ()
            
            
            #--------------------------------------------------------------------------------------------
            cr.set_source_rgb(0, 0, 0)
            cr.set_font_size(font_size)
            
            if self.label_mode == 0:
                text = str(i)
            elif self.label_mode == 1:
                '''
                   The label of reaction coordinate distance list is wrong  
                   self.dataRC2[i][0] - means frames on the vertical (index i)
                '''
                text = '{:3.2f}'.format(self.dataRC2[i][0])
            else:
                text = ''
                pass
            
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            #print(x_bearing, y_bearing, width, height, x_advance, y_advance)
            #cr.move_to(self.bx + i*self.factor_x - width/2.0    ,  self.by -15 )   
            cr.move_to(self.bx + i*self.factor_x - x_advance/2.0 +(self.factor_x)/2   ,  self.by -15 )   
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
        
        for j in range(0,self.size_y, y_major_ticks):
            #--------------------------------------------------------------------------------------------
            # marcacoes em y
            cr.move_to (self.bx      , self.by + j*self.factor_y + (self.factor_y)/2.0 )
            cr.line_to (self.bx -10  , self.by + j*self.factor_y + (self.factor_y)/2.0 ) 
            #--------------------------------------------------------------------------------------------
            cr.stroke ()
            
            
            #--------------------------------------------------------------------------------------------
            # texto em y
            cr.set_source_rgb(0, 0, 0)
            cr.set_font_size(font_size)

            if self.label_mode == 0:
                text = str(j)
            elif self.label_mode == 1:
                '''
                The label of reaction coordinate distance list is wrong  
                self.dataRC1[j][0] - means frames on the vertical (index j)
                '''
                text = '{:3.2f}'.format(self.dataRC1[0][j])
            else:
                text = ''
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)

            cr.move_to(self.bx  -15 -width     , self.by + j*self.factor_y + height/2  + (self.factor_y)/2)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
        
    def draw_dots (self, cr, color = [0,0,0]):
        """ Function doc """
        cr.set_source_rgb( color[0], color[1], color[2])
        #print(self.points)
        for x,y in self.points:
            #else:
            cr.set_source_rgb( color[0], color[1], color[2])
            cr.arc(((x+0.5)*self.factor_x)+ self.bx, ((y+0.5)*self.factor_y)+self.by , 5, 0, 2 * 3.14)
            #cr.stroke ()
            cr.fill()
        # draws the red dot at the energey matrix
        if len(self.points) and self.selected_dot:
            x = self.selected_dot[0]
            y = self.selected_dot[1]
            cr.set_source_rgb( self.sel_dot_rgb[0], self.sel_dot_rgb[1], self.sel_dot_rgb[2])
            cr.arc(((x+0.5)*self.factor_x)+ self.bx, ((y+0.5)*self.factor_y)+self.by , 5, 0, 2 * 3.14)
            cr.fill()
   
    
    def draw_lines (self, cr, line_width = 2, color = [0,0,0]):
        """ Function doc """
        cr.set_source_rgb( color[0], color[1], color[2])
        cr.set_line_width (line_width)
        for x, y in self.points:
            #cr.line_to(x , y)
            cr.line_to(((x+0.5)*self.factor_x)+ self.bx, ((y+0.5)*self.factor_y)+self.by)
        cr.stroke ()


    def on_draw(self, widget, cr, ):
        self.cr =  cr
        self.width = widget.get_allocated_width()   #widget dimentions
        self.height = widget.get_allocated_height() #widget dimentions   
        self._make_bg(cr)

        if self.data:
            pass
        else:
            return None
        self.size_x = len(self.data)
        self.size_y = len(self.data[0])

        
        self.x_plot_edge = self.width  -(self.bx+70)
        self.y_plot_edge = self.height -(self.by)
        
        '''(self.x_plot_edge - self.bx) = number of pixels the plot has at x 
        We divide by the dimension of the data matrix to know how much each 
        rectangle of the plot will be in length'''
        
        self.factor_x = (self.x_plot_edge - self.bx)/self.size_x 
        self.factor_y = (self.y_plot_edge - self.by)/self.size_y
        
        self.x_plot_edge = self.size_x*self.factor_x-self.factor_x
        self.y_plot_edge = self.size_y*self.factor_y-self.factor_y
        
        if self.data_max is None:
            self.define_datanorm()
        
        '''
        self.norm_data = self.data - np.min(self.data)
        
        
        self.data_min = np.min(self.norm_data)
        self.data_max = np.max(self.norm_data)
        #self.data_max = 100
        delta = (self.data_max - self.data_min)
        self.norm_data = self.norm_data/delta
        
        self.set_threshold_color ( _min = 0, _max = 200)
        '''
        
        # Criar uma superfície de imagem temporária
        self.temp_surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)
        self.temp_cr = cairo.Context(self.temp_surface)
        
        #print(np.min(self.norm_data), np.max(self.norm_data) )
        
        self.draw_color_bar (self.temp_cr, res = self.size_y)#(cr, res = self.size_y)

        self.draw_image (self.temp_cr)#(cr)
        
        self.draw_image_box (self.temp_cr, line_width = 1, color = [0,0,0])#(cr, line_width = 1, color = [0,0,0])

        self.draw_scale(self.temp_cr)#(cr)

        self.draw_dots (self.temp_cr)#(cr)

        self.draw_lines (self.temp_cr, line_width = 2, color = [0,0,0])#(cr, line_width = 2, color = [0,0,0])
    
        cr.set_source_surface(self.temp_surface, 0, 0)
        cr.paint()
    
        #rgb = self.get_pixel_rgb_from_surface(self.temp_surface, 0, 0)
        #print(f"RGB at (0, 0): {rgb}")
    
        #width, height = 300, 300
        #surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        #surface.write_to_png("mesh_pattern_example.png")
        
    def on_motion (self, widget, event):
        """ Function doc """
        if self.norm_data is not None:
            pass
        else:
            return False
        x_on_plot, y_on_plot, x, y = self.get_i_and_j_from_click_event (event)
        i = y_on_plot
        j = x_on_plot
        
        j_size = len(self.norm_data)
        i_size = len(self.norm_data[0])

        #rgb = self.get_pixel_rgb_from_surface(self.temp_surface, x, y)
        #print(f"RGB at ({x}, {y}): {rgb}")

        if i < 0 or j < 0:
            pass
        else:
            if i < i_size and j < j_size:
                text = 'i {}    |    j {}    |    rc1 {:4.2f}    |    rc2 {:4.2f}    |    E = {:6.3f}'.format(i, j,  self.dataRC1[j][i], self.dataRC2[j][i],  self.data[j][i])

                if self.RC_label:
                    self.RC_label.set_text(text)
                else:
                    print(text)


    def get_pixel_rgb_from_surface(self, surface, x, y):
        """Lê o valor RGB de um pixel dentro de uma cairo.ImageSurface."""
        
        # Pega os dados da superfície de Cairo
        buf = surface.get_data()

        # Converte os dados em uma array numpy (ARGB formato, 4 bytes por pixel)
        array = np.frombuffer(buf, dtype=np.uint8)
        array = array.reshape((surface.get_height(), surface.get_width(), 4))

        # Pega o valor do pixel na posição especificada (x, y) (formato ARGB)
        pixel = array[y, x]

        # Extraímos os valores de Red, Green e Blue
        r, g, b = pixel[1], pixel[2], pixel[3]

        return (r, g, b)


class XYPlot(Gtk.DrawingArea):
    """ Class doc """
    
    def __init__ (self, bg_color = [1,1,1], mode = None ):
        """ Class initialiser """
        super().__init__( )

        
        self.data = []
        
        self.bg_color = bg_color
        self.decimal = 2

        
        self.x_minor_ticks = 6
        self.x_major_ticks = 5
        
        self.y_minor_ticks = 5
        self.y_major_ticks = 10
        
        self.Xmin_list = []
        self.Ymin_list = []
        
        self.Xmax_list = []
        self.Ymax_list = []





        self.Y_space = 10


        #self.set_bg_color (bg_color[0], bg_color[1], bg_color[2])
        self.set_events(Gdk.EventMask.ALL_EVENTS_MASK)
        
        self.connect("draw", self.on_draw)
        self.connect("button_press_event", self.on_mouse_button_press)
        self.connect("motion-notify-event", self.on_motion)


        self.line_color     = [0,0,0]
        #self.bg_color       = [1,1,1]
        self.sel_color      = [1,0,0]
        
        #self.line_color     = [1,1,1]
        #self.bg_color       = [0,0,0]
        #self.sel_color      = [1,0,0]
        
        self.bglines_color  = [0.5,0.5,0.5]
        self.bx = 100#80
        self.by = 50#50 

        self.y_botton = -1
        self.y_top    =  1

    
    def data_update (self, data):
        """ Function doc """
        X = data['X']
        Y = data['Y']
        
        self.Xmin_list.append(min(X))
        self.Ymin_list.append(min(Y))
        self.Xmax_list.append(max(X))
        self.Ymax_list.append(max(Y))
        
        
        data['Xmin'] =  min(X)
        data['Ymin'] =  min(Y)
        data['Xmax'] =  max(X)
        data['Ymax'] =  max(Y)
        
        return data


    def add (self, X = None, Y = None, 
              symbol = 'dot', sym_color = [1,1,1], sym_fill = True, 
              line = 'solid', line_color = [1,1,1], energy_label = None):
        
        """ Function doc """

        data = {
                'X'          : X     ,
                'Y'          : Y     ,
                'Xmin'       : None,
                'Ymin'       : None,
                'Xmax'       : None,
                'Ymax'       : None,
                'sym'        : symbol,
                'sym_color'  : sym_color,
                'sym_fill'   : sym_fill,
                'line'       : line  ,
                'line_color' : line_color
               }

        data = self.data_update( data)
        

        #self.Xmin_list.append(min(X))
        #self.Ymin_list.append(min(Y))
        #self.Xmax_list.append(max(X))
        #self.Ymax_list.append(max(Y))
        #
        #
        #data = {
        #        'X'          : X     ,
        #        'Y'          : Y     ,
        #        'Xmin'       : min(X),
        #        'Ymin'       : min(Y),
        #        'Xmax'       : max(X),
        #        'Ymax'       : max(Y),
        #        'sym'        : symbol,
        #        'sym_color'  : sym_color,
        #        'sym_fill'   : sym_fill,
        #        'line'       : line  ,
        #        'line_color' : line_color
        #       }
        
        self.energy_label = energy_label
        
        self.data.append(data)
        self.define_xy_limits()
    
    
    def define_xy_limits (self, y_botton = -3, y_top = 2):
        """ Function doc """
        
        self.Ymin = min(self.Ymin_list)
        self.Ymax = max(self.Ymax_list)
        

        
        self.Ymin = self.Ymin + self.Ymin/self.Y_space
        self.Ymax = self.Ymax + self.Ymax/self.Y_space
        self.Ymin = int(self.Ymin)# + self.y_botton #- 3 #round(delta/20)
        self.Ymax = int(self.Ymax) + self.y_top    #+ 2 #round(delta/20)
        delta = self.Ymax - self.Ymin

        
        #self.Ymin = -2.5
        #self.Ymax =  2.5
        
        self.Xmax = max(self.Xmax_list)
        self.Xmin = min(self.Xmin_list)
        #self.Xmax = self.Xmax - self.Xmax/20
        #self.Xmin = self.Xmin - self.Xmin/20
        
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin
        
        if self.deltaX == 0:
            self.deltaX = 1
        
        if self.deltaY == 0:
            self.deltaY = 1
        

    def draw_box (self, cr, line_width = 3.0 , color = [1, 1, 1] ):
        """ Function doc """
        cr.set_line_width (3.0)
        #----------------------------------------------------------------------
        #retangle
        cr.set_source_rgb( color[0], color[1], color[2])
        cr.rectangle(self.bx,
                     self.by, 
                     self.x_box_size,
                     self.y_box_size)
        cr.stroke ()
        #----------------------------------------------------------------------
    
    def draw_XY_axes (self, cr, line_width = 3.0 , color = [1, 1, 1], font_size = 12 ):
        
        
        cr.set_line_width (1.0)
        cr.set_source_rgb( color[0], color[1], color[2])
        #---------------------------------------------------------------------------------------------------
        # x axy
        cr.move_to (self.bx, self.y_box_size+ self.by)
        cr.line_to (self.bx +self.x_box_size , self.y_box_size+ self.by)
        cr.stroke ()
        
        # y axy
        cr.move_to (self.bx, self.y_box_size+ self.by)
        cr.line_to (self.bx,  self.by)
        cr.stroke ()
        #---------------------------------------------------------------------------------------------------

    def draw_XY_scale (self, cr, line_width = 3.0 , color = [1, 1, 1], font_size =12 ):
        """ Function doc """
        
        #---------------------------------------------------------------------------------------------------
        x_major_factor = (self.x_box_size) / float(self.x_major_ticks)
        x_minor_factor = x_major_factor / float(self.x_minor_ticks)
        counter = 1
        for i in range(0, self.x_major_ticks+1):
            # marcacoes em x
            cr.move_to (self.bx+ i*x_major_factor, self.y_box_size+self.by )
            cr.line_to (self.bx+ i*x_major_factor, self.y_box_size+self.by + 15 ) 
            if i == self.x_major_ticks:
                pass
            else:
                for j in range(0, self.x_minor_ticks ):
                    cr.move_to (self.bx+ i*x_major_factor+j*x_minor_factor, self.y_box_size+self.by )
                    cr.line_to (self.bx+ i*x_major_factor+j*x_minor_factor, self.y_box_size+self.by + 10 )
            counter +=1
        #---------------------------------------------------------------------------------------------------
            #--------------------------------------------------------------------------------------------
            # texto em X
            
            cr.set_font_size(font_size)
            #text = '{:^6.2f}'.format((i/self.x_major_ticks)*self.deltaX + self.Xmin)
            text = '{:^6.'+str(self.decimal)+'f}'
            text = text.format((i/self.x_major_ticks)*self.deltaX + self.Xmin)
            
            if len(text) >7:
                num = (i/self.x_major_ticks)*self.deltaX -self.Xmax 
                text = (f"{num:.2e}")
            
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            
            cr.move_to(  
                        self.bx + i*x_major_factor - x_advance/2, 
                        self.by+self.y_box_size + height + 20)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1


            #cr.set_font_size(font_size)
            #
            #text = str(i)
            #x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            ##print(x_bearing, y_bearing, width, height, x_advance, y_advance)
            ##cr.move_to(self.bx + i*self.factor_x - width/2.0    ,  self.by -15 )   
            #cr.move_to(self.bx + i*self.factor_x - x_advance/2.0 +(self.factor_x)/2   ,  self.by -15 )   
            #cr.show_text(text)







        
        #---------------------------------------------------------------------------------------------------
        y_major_factor = (self.y_box_size) / float(self.y_major_ticks)
        y_minor_factor =  y_major_factor   / float(self.y_minor_ticks)
        
        counter = 1
        self.deltaY2 =  self.deltaY
        for i in range(0, self.y_major_ticks+1):
            # marcacoes em x
            cr.move_to (self.bx, self.by+ i*y_major_factor )
            cr.line_to (self.bx - 15, self.by+ i*y_major_factor ) 
            #--------------------------------------------------------------------------------------------
            if i == self.y_major_ticks:
                pass
            else:
                for j in range(0, self.y_minor_ticks ):
                    cr.move_to  (self.bx, 
                                 self.by+ i*y_major_factor+j*y_minor_factor)
                    cr.line_to (self.bx -  10, 
                                self.by+ i*y_major_factor+j*y_minor_factor )
        
        
            #--------------------------------------------------------------------------------------------
            # texto em y
            cr.set_font_size(font_size)
            text = '{:>6.3f}'.format(self.Ymax - (i/self.y_major_ticks)*self.deltaY )
            if len(text) >6:
                num = self.Ymax - (i/self.y_major_ticks)*self.deltaY
                text = (f"{num:.2e}")
            
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            
            cr.move_to(  
                        self.bx -x_advance -20  , 
                        self.by+ i*y_major_factor )
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1
        
        
        
        
        
        
        
        
        
        
        
        
        
            counter +=1
        cr.stroke ()
        #---------------------------------------------------------------------------------------------------
        
        
        
        
        
        ''' backgound lines '''
        cr.set_source_rgb( 0.3, 0.3, 0.3)
        #---------------------------------------------------------------------------------------------------
        x_major_factor = (self.x_box_size) / float(self.x_major_ticks)
        x_minor_factor = x_major_factor / float(self.x_minor_ticks)
        counter = 1
        for i in range(0, self.x_major_ticks+1):
            # marcacoes em x
            cr.move_to (self.bx+ i*x_major_factor, self.y_box_size+self.by )
            cr.line_to (self.bx+ i*x_major_factor,  self.by  ) 
        #---------------------------------------------------------------------------------------------------
        
        
        #---------------------------------------------------------------------------------------------------
        y_major_factor = (self.y_box_size) / float(self.y_major_ticks)
        y_minor_factor =  y_major_factor   / float(self.y_minor_ticks)
        counter = 1
        for i in range(0, self.y_major_ticks+1):
            # marcacoes em x
            cr.move_to (self.bx, self.by+ i*y_major_factor )
            cr.line_to (self.bx +self.x_box_size, self.by+ i*y_major_factor ) 
        cr.stroke ()
        #---------------------------------------------------------------------------------------------------
        
    def on_draw(self, widget, cr):
        self.cr =  cr
        self.width = widget.get_allocated_width()
        self.height = widget.get_allocated_height()

        line_color = self.line_color
        bg_color   = self.bg_color  
        cr.set_source_rgb( bg_color[0], bg_color[1], bg_color[2])
        #line_color = [0,0,0]
        cr.paint()

        self.x = self.width
        self.y = self.height





        self.x_box_size = self.width  - (self.bx*1.5)
        self.y_box_size = self.height - (self.by*1.8)

        self.draw_box    (cr, line_width = 3.0 , color = line_color) #color = [1, 1, 1])
        self.draw_XY_axes(cr, line_width = 3.0 , color = line_color) #color = [1, 1, 1])
            
        if self.data == []:
            return False
        else:

            if self.deltaY > 0:
                for data in self.data:
                    data['Ynorm'] = [(y-self.Ymin)/self.deltaY for y in data['Y']]
                    data['Xnorm'] = [(x-self.Xmin)/self.deltaX for x in data['X']]
            else:
                for data in self.data:
                    data['Ynorm'] = data['Y']
                    data['Xnorm'] = data['X']
                
                
                
            #self.Ynorm_max = max(self.norm_data)
            
            self.draw_XY_scale ( cr, line_width = 3.0 , color = line_color)#color = [1, 1, 1] )
            
            #---------------------------------------------------------------- 
            #                          LINES                                  
            #---------------------------------------------------------------- 
            cr.set_line_width (2.0)                                           
            #---------------------------------------------------------------- 
            
            
  
            cx = 0#self.bx - self.norm_X[0]                                   
            #---------------------------------------------------------------- 
            cr.set_source_rgb( line_color[0], line_color[1], line_color[2])
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                
                #r = random.random()
                #g = random.random()
                #b = random.random()
                #rgb = [r,g,b,]
                
                
                line_color = data['line_color']
                #line_color = rgb
                #print(line_color)
                cr.set_source_rgb( line_color[0], line_color[1], line_color[2]) 
                
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    
                    #color = data['line_color']
                    #color = [0,0,0]
                    #cr.set_source_rgb( color[0], color[1], color[2])
                    
                    cr.line_to (x*self.x_box_size + cx, 
                            #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                            (1*self.y_box_size + self.by) - (y*self.y_box_size))
                cr.stroke ()
            #---------------------------------------------------------------- 
            #                          DOTS                                  
            #---------------------------------------------------------------- 
            #'''
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    color = data['line_color']
                    color = [0,0,0]
                    
                    if data['sym'] is not None:
                        color = data['sym_color']
                        
                        #cr.set_source_rgb( color[0], color[1], color[2])
                        cr.set_source_rgb( self.sel_color[0], self.sel_color[1], self.sel_color[2])
                        
                        cr.arc (x*self.x_box_size + cx, 
                                #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                                (1*self.y_box_size + self.by) - (y*self.y_box_size), 5, 0, 2*3.14)
                    
                    if data['sym_fill']:
                        cr.set_source_rgb( self.sel_color[0], self.sel_color[1], self.sel_color[2])
                        cr.fill()
                    else:
                        cr.set_source_rgb( self.line_color[0], self.line_color[1], self.line_color[2])
                        cr.stroke ()
    
    
    def on_motion(self, widget, event):
        '''(i/self.x_major_ticks)*self.deltaX + self.Xmin'''
        if self.data == []:
            return False 
        
        else:
            (x, y) = int(event.x), int(event.y)
            
            print("Mouse moved to:", 'x = {:10.5f} y = {:10.5f}'.format( 
                
                ((((x-self.bx)) / self.x_box_size)  *  self.deltaX + self.Xmin),
                
                ((self.y_box_size-(y-self.by)) / self.y_box_size)  *  self.deltaY + self.Ymin)
                )
    
    
    def on_mouse_button_press (self, widget, event):
        """ Function doc """
        pass


class XYScatterPlot(Gtk.DrawingArea) :
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
        super().__init__( )
        self.data = []
        
        self.bg_color = [0,0,0]
        self.decimal = 2

        
        self.x_minor_ticks = 6
        self.x_major_ticks = 5
        
        self.y_minor_ticks = 5
        self.y_major_ticks = 10
        
        self.Xmin_list = []
        self.Ymin_list = []
        
        self.Xmax_list = []
        self.Ymax_list = []





        self.Y_space = 10


        #self.set_bg_color (bg_color[0], bg_color[1], bg_color[2])
        self.set_events(Gdk.EventMask.ALL_EVENTS_MASK)
        
        self.connect("draw", self.on_draw)
        #self.connect("button_press_event", self.on_mouse_button_press)
        #self.connect("motion-notify-event", self.on_motion)


        self.line_color     = [0,0,0]
        #self.bg_color       = [1,1,1]
        self.sel_color      = [1,0,0]
        
        #self.line_color     = [1,1,1]
        #self.bg_color       = [0,0,0]
        #self.sel_color      = [1,0,0]
        
        self.bglines_color  = [0.5,0.5,0.5]
        self.bx = 100#80
        self.by = 50#50 
        
        self.y_botton = 0
        self.y_top    = 0     

    def define_xy_limits (self, y_botton = -3, y_top = 2):
        """ Function doc """
        
        self.Ymin = min(self.Ymin_list)
        self.Ymax = max(self.Ymax_list)
        

        
        self.Ymin = self.Ymin + self.Ymin/self.Y_space
        self.Ymax = self.Ymax + self.Ymax/self.Y_space
        self.Ymin = int(self.Ymin) + self.y_botton #- 3 #round(delta/20)
        self.Ymax = int(self.Ymax) + self.y_top    #+ 2 #round(delta/20)
        delta = self.Ymax - self.Ymin

        
        #self.Ymin = -2.5
        #self.Ymax =  2.5
        
        self.Xmax = max(self.Xmax_list)
        self.Xmin = min(self.Xmin_list)
        #self.Xmax = self.Xmax - self.Xmax/20
        #self.Xmin = self.Xmin - self.Xmin/20
        
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin
        
        if self.deltaX == 0:
            self.deltaX = 1
        
        if self.deltaY == 0:
            self.deltaY = 1

    def add (self, X = None, Y = None, Z = None,
              symbol = 'dot', sym_color = [1,1,1], sym_fill = True, 
              line = 'solid', line_color = [1,1,1], energy_label = None):
        
        """ Function doc """

        data = {
                'X'          : X     ,
                'Y'          : Y     ,
                'Z'          : Z     ,
                'Xmin'       : None,
                'Ymin'       : None,
                'Xmax'       : None,
                'Ymax'       : None,
                'sym'        : symbol,
                'sym_color'  : sym_color,
                'sym_fill'   : sym_fill,
                'line'       : line  ,
                'line_color' : line_color
               }

        data = self.data_update( data)
        
        
        tempZ = Z.copy()
        n = tempZ.count(None)
        for i in range(n):
            tempZ.remove(None)
        
        data['Zmax'] = max(tempZ)
        data['Zmin'] = min(tempZ)
        
        self.c_factor =  1/(data['Zmax']-data['Zmin'])
        
        
        self.energy_label = energy_label
        self.data.append(data)
        self.define_xy_limits()
  
    def data_update (self, data):
        """ Function doc """
        X = data['X']
        Y = data['Y']
        
        self.Xmin_list.append(min(X))
        self.Ymin_list.append(min(Y))
        self.Xmax_list.append(max(X))
        self.Ymax_list.append(max(Y))
        
        
        data['Xmin'] =  min(X)
        data['Ymin'] =  min(Y)
        data['Xmax'] =  max(X)
        data['Ymax'] =  max(Y)
        
        return data

    def draw_box (self, cr, line_width = 3.0 , color = [1, 1, 1] ):
        """ Function doc """
        cr.set_line_width (3.0)
        #----------------------------------------------------------------------
        #retangle
        print('draw_box')
        cr.set_source_rgb( color[0], color[1], color[2])
        cr.rectangle(self.bx,
                     self.by, 
                     self.x_box_size,
                     self.y_box_size)
        cr.stroke ()
        #----------------------------------------------------------------------

    def on_draw(self, widget, cr):
        self.cr =  cr
        self.width = widget.get_allocated_width()
        self.height = widget.get_allocated_height()

        line_color = self.line_color
        #bg_color   = self.bg_color  
        bg_color   = [1,1,1]  
        #print(bg_color)
        cr.set_source_rgb( bg_color[0], bg_color[1], bg_color[2])
        #line_color = [0,0,0]
        cr.paint()

        self.x = self.width
        self.y = self.height

        #self.data['sym_fill'] = True



        self.x_box_size = self.width  - (self.bx*1.5)
        self.y_box_size = self.height - (self.by*1.8)

        self.draw_box    (cr, line_width = 3.0 , color = line_color) #color = [1, 1, 1])
        #self.draw_XY_axes(cr, line_width = 3.0 , color = line_color) #color = [1, 1, 1])
            
        if self.data == []:
            return False
        else:
            pass
            #print(
            #'\nself.Ymax     ', self.Ymax      ,
            #'\nself.Ymin     ', self.Ymin      ,
            #'\nself.Xmax     ', self.Xmax      ,
            #'\nself.Xmin     ', self.Xmin      ,
            #'\nself.deltaY   ', self.deltaY    ,
            #'\nself.deltaX   ', self.deltaX    ,
            #)
            
            
    
            if self.deltaY > 0:
                for data in self.data:
                    data['Ynorm'] = [(y-self.Ymin)/self.deltaY for y in data['Y']]
                    data['Xnorm'] = [(x-self.Xmin)/self.deltaX for x in data['X']]
            else:
                for data in self.data:
                    data['Ynorm'] = data['Y']
                    data['Xnorm'] = data['X']
                
                
                
            #self.Ynorm_max = max(self.norm_data)
            
            #self.draw_XY_scale ( cr, line_width = 3.0 , color = line_color)#color = [1, 1, 1] )
            
            #---------------------------------------------------------------- 
            #                          LINES                                  
            #---------------------------------------------------------------- 
            cr.set_line_width (2.0)                                           
            #---------------------------------------------------------------- 
            
            
  
            cx = 0 #self.bx - self.norm_X[0]                                   
            #---------------------------------------------------------------- 
            cr.set_source_rgb( line_color[0], line_color[1], line_color[2])
            '''
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                
                #r = random.random()
                #g = random.random()
                #b = random.random()
                #rgb = [r,g,b,]
                
                
                line_color = data['line_color']
                #line_color = rgb
                #print(line_color)
                cr.set_source_rgb( line_color[0], line_color[1], line_color[2]) 
                
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    
                    #color = data['line_color']
                    #color = [0,0,0]
                    #cr.set_source_rgb( color[0], color[1], color[2])
                    
                    cr.line_to (x*self.x_box_size + cx, 
                            #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                            (1*self.y_box_size + self.by) - (y*self.y_box_size))
                cr.stroke ()
            #'''
            #---------------------------------------------------------------- 
            #                          DOTS                                  
            #---------------------------------------------------------------- 
            #'''
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    color = data['line_color']
                    color = [0,0,0]
                    data['sym_fill'] = True
                    
                    if data['sym'] is not None:
                        color = data['sym_color']
                        
                        #cr.set_source_rgb( color[0], color[1], color[2])
                        #cr.set_source_rgb( self.sel_color[0], self.sel_color[1], self.sel_color[2])
                        
                        cr.arc (x*self.x_box_size + cx, 
                                #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                                (1*self.y_box_size + self.by) - (y*self.y_box_size), 5, 0, 2*3.14)
                    
                        if data['sym_fill']:
                            
                            if data['Z'][i]:
                                #print(data['Z'][i]* 1/42.0)
                                f = data['Z'][i] * 1/42.0
                                #self.sel_color[0]
                                cr.set_source_rgb( f, 0, f)
                                cr.fill()
                            else:
                                cr.set_source_rgb( 1, 1, 1)
                            
                                cr.fill()
                    cr.stroke ()
                    #else:
                    #    cr.set_source_rgb( self.line_color[0], self.line_color[1], self.line_color[2])
                        



