import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, cairo
import numpy as np
import sys


from util.colormaps import COLOR_MAPS  

def interpolate_color(color1, color2, fraction):
    return [ ( (x1 + fraction * (x2 - x1))) for (x1, x2) in zip(color1, color2)]

def get_color(valor, color_map):#, factor = 2):   
    
    
    #valor = (valor+1)/factor   
    
    cutoffs = list(color_map.keys())
    color_list = list(color_map.values())
    counter = 0 
    
    #print(valor)
    
    
    
    for cutoff, color in color_map.items():
        
        if valor <= cutoff :
            
            fraction = (valor -  cutoffs[counter-1]) / ( cutoff - cutoffs[counter-1])
            
            
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
            #print(valor, cutoff, fraction,color )

            return color
        else:
            pass
            #print(valor, cutoff)#, fraction,color )
        counter +=1

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
        #self.connect("button_press_event", self.on_mouse_button_press)
        #self.connect("motion-notify-event", self.on_motion)
    
    def set_cmap (self, cmap = 'jet'):
        """ Function doc """
        
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
            return False
        elif int(y_on_plot) < 0 or int(y_on_plot) >= self.size_y :
            print('canvas')
            self.points = []
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
            print('canvas')
            self.points = []
        elif int(y_on_plot) < 0 or int(y_on_plot) >= self.size_y :
            print('canvas')
            self.points = []
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


class ImagePlot(Canvas):
    """ Class doc """
    
    def __init__ (self, data = None):
        """ Class initialiser """
        super().__init__( )

        self.bx = 50 
        self.by = 50 
        self.points = []
        
        self.data = data
        
    
    def set_cmap (self, cmap = 'jet'):
        """ Function doc """
        
    
    def draw_color_bar (self, cr, res, gap = 20, label_spacing = 2, font_size = 12):
        """ Function doc """
         
       
        cr.set_line_width (1.0)
        
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
        
        
        
        factor = int( self.size_y/10 )
        
        if factor ==  0:
            factor = 1
        else:
            pass
        
        
        _min = np.min(self.data)
        _max = np.max(self.data)
        factor2 = (_max - _min)/10
        print(_min, _max)
        counter = 0
        for j in range(0, self.size_y+1, factor ):    
            
            
            
            j = self.size_y -j-1
            
            #z = ((j  / float(self.size_y)) - 0.5)*2
            #print(z)
            cr.set_source_rgb( 0, 0,0)
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
            text = '{:<6.3f}'.format(_min + factor2*counter )
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            cr.move_to(round(self.size_x * self.factor_x)+self.bx*1.2 + 45 ,  
                             self.by + j*self.factor_y+ self.factor_y + height/2)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1








        cr.set_line_width (2.0)
        cr.rectangle(round(self.size_x * self.factor_x)+self.bx +gap,
                     self.by   , 
                     round(20+1), 
                     round(self.size_y * self.factor_y))
        cr.stroke ()


    def draw_image (self, cr):
        """ Function doc """
        n = 0 
        for i in range(0,self.size_x, ):
            for j in range(0,self.size_y):
                
                z = self.norm_data[i][j]

                color = get_color(z, self.color_map)

                cr.set_source_rgb( color[0], color[1], color[2]   )
                
                cr.rectangle(self.bx+i*self.factor_x -1,
                             self.by+j*self.factor_y -1, 
                             round(self.factor_x)+1, 
                             round(self.factor_y)+1)
                cr.fill()


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
    
            text = str(i)
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
    
            text = str(j)
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)

            cr.move_to(self.bx  -15 -width     , self.by + j*self.factor_y + height/2  + (self.factor_y)/2)
            cr.show_text(str(j))
            #--------------------------------------------------------------------------------------------
        
    def draw_dots (self, cr, color = [0,0,0]):
        """ Function doc """
        cr.set_source_rgb( color[0], color[1], color[2])
        
        for x, y in self.points:
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
        
        
        
        
        self.norm_data = self.data - np.min(self.data)
        
        
        self.data_min = np.min(self.norm_data)
        self.data_max = np.max(self.norm_data)
        delta = (self.data_max - self.data_min)
        #print(self.data_min, self.data_max )
        self.norm_data = self.norm_data/delta

        #print (self.norm_data)

        #print(np.min(self.norm_data), np.max(self.norm_data) )
        
        self.draw_color_bar (cr, res = self.size_y)

        self.draw_image (cr)
        
        self.draw_image_box (cr, line_width = 1, color = [0,0,0])

        self.draw_scale(cr)

        self.draw_dots (cr)

        self.draw_lines (cr, line_width = 2, color = [0,0,0])



class XYPlot(Gtk.DrawingArea):
    """ Class doc """
    
    def __init__ (self, bg_color = [1,1,1] ):
        """ Class initialiser """
        super().__init__( )

        
        self.data = []
        
        self.bg_color = None
        self.decimal = 2

        
        self.x_minor_ticks = 6
        self.x_major_ticks = 5
        
        self.y_minor_ticks = 5
        self.y_major_ticks = 10
        
        self.Xmin_list = []
        self.Ymin_list = []
        
        self.Xmax_list = []
        self.Ymax_list = []





        self.Y_space = 20


        #self.set_bg_color (bg_color[0], bg_color[1], bg_color[2])
        self.set_events(Gdk.EventMask.ALL_EVENTS_MASK)
        
        self.connect("draw", self.on_draw)
        self.connect("button_press_event", self.on_mouse_button_press)
        self.connect("motion-notify-event", self.on_motion)


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
              line = 'solid', line_color = [1,1,1]):
        
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
        
 
        
        self.data.append(data)
        self.define_xy_limits()
    

    def define_xy_limits (self):
        """ Function doc """
        
        self.Ymin = min(self.Ymin_list)
        self.Ymax = max(self.Ymax_list)
        

        
        self.Ymin = self.Ymin + self.Ymin/self.Y_space
        self.Ymax = self.Ymax + self.Ymax/self.Y_space
        self.Ymin = int(self.Ymin) - 3 #round(delta/20)
        self.Ymax = int(self.Ymax) + 2 #round(delta/20)
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
        cr.set_source_rgb( 1, 1, 1)
        cr.rectangle(self.bx,
                     self.by, 
                     self.x_box_size,
                     self.y_box_size)
        cr.stroke ()
        #----------------------------------------------------------------------
    
    def draw_XY_axes (self, cr, line_width = 3.0 , color = [1, 1, 1], font_size = 12 ):
        
        
        cr.set_line_width (1.0)
        cr.set_source_rgb( 1, 1, 1)
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


        cr.set_source_rgb( 0, 0, 0)
        cr.paint()

        self.x = self.width
        self.y = self.height



        self.bx = 100 
        self.by = 50 

        self.x_box_size = self.width  - (self.bx*1.5)
        self.y_box_size = self.height - (self.by*1.8)

        self.draw_box    (cr, line_width = 3.0 , color = [1, 1, 1])
        self.draw_XY_axes(cr, line_width = 3.0 , color = [1, 1, 1])
            
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
            
            self.draw_XY_scale ( cr, line_width = 3.0 , color = [1, 1, 1] )
            
            
            
            #----------------------------------------------------------------------
            cr.set_line_width (2.0)
            #----------------------------------------------------------------------
            cr.set_source_rgb( 1, 1, 1)
            cx = 0#self.bx - self.norm_X[0]
            #----------------------------------------------------------------------
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    color = data['line_color']
                    
                    cr.set_source_rgb( color[0], color[1], color[2])
                    
                    cr.line_to (x*self.x_box_size + cx, 
                            #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                            (1*self.y_box_size + self.by) - (y*self.y_box_size))
                cr.stroke ()
                
            
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    color = data['line_color']
                    
                    if data['sym'] is not None:
                        color = data['sym_color']
                        cr.set_source_rgb( color[0], color[1], color[2])
                        cr.arc (x*self.x_box_size + cx, 
                                #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                                (1*self.y_box_size + self.by) - (y*self.y_box_size), 5, 0, 2*3.14)
                    
                    if data['sym_fill']:
                        cr.fill()
                    else:
                        cr.stroke ()
            #print(
            #'\nself.Ymax     ', self.Ymax      ,
            #'\nself.Ymin     ', self.Ymin      ,
            #'\nself.Xmax     ', self.Xmax      ,
            #'\nself.Xmin     ', self.Xmin      ,
            #'\nself.deltaY   ', self.deltaY    ,
            #'\nself.deltaX   ', self.deltaX    ,
            #)
    
    def on_motion(self, widget, event):
        '''(i/self.x_major_ticks)*self.deltaX + self.Xmin'''
        if self.data == []:
            return False 
        
        else:
            (x, y) = int(event.x), int(event.y)
            
            print("Mouse moved to:", 
                #x, y,  
                #x-self.bx, y-self.by, 
                
                #((((x-self.bx)) / self.x_box_size)  *  self.deltaX + self.Xmin),
                
                ((self.y_box_size-(y-self.by)) / self.y_box_size)  *  self.deltaY + self.Ymin)
            
    
    def on_mouse_button_press (self, widget, event):
        """ Function doc """
        pass
        #self.Ymin = -1.5 -1
        #self.Ymax =  1.5 +1
        ##print(self.Ymin, self.Ymax)
        ##self.define_xy_limits ()
        #self.deltaY  = self.Ymax - self.Ymin
        #self.deltaX  = self.Xmax - self.Xmin
        
        
        
        
