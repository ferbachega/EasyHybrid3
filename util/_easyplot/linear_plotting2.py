import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, cairo
import  random  
import numpy as np
import math





class XYPlot(Gtk.DrawingArea):
    """ Class doc """
    
    def __init__ (self, bg_color = [1,1,1] ):
        """ Class initialiser """
        super().__init__( )

        
        self.data = []
        
        self.bg_color = None
        


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

    def add (self, X = None, Y = None, 
              symbol = 'dot', sym_color = [1,1,1], sym_fill = True, 
              line = 'solid', line_color = [1,1,1]):
        
        """ Function doc """
        
        self.Xmin_list.append(min(X))
        self.Ymin_list.append(min(Y))
        self.Xmax_list.append(max(X))
        self.Ymax_list.append(max(Y))
        
        
        data = {
                'X'          : X     ,
                'Y'          : Y     ,
                'Xmin'       : min(X),
                'Ymin'       : min(Y),
                'Xmax'       : max(X),
                'Ymax'       : max(Y),
                'sym'        : symbol,
                'sym_color'  : sym_color,
                'sym_fill'   : sym_fill,
                'line'       : line  ,
                'line_color' : line_color
               }
        
 
        
        self.data.append(data)
        self.define_xy_limits()
    

    def define_xy_limits (self):
        """ Function doc """
        
        self.Ymin = min(self.Ymin_list)
        self.Ymax = max(self.Ymax_list)
        self.Ymin = self.Ymin + self.Ymin/self.Y_space
        self.Ymax = self.Ymax + self.Ymax/self.Y_space
        
        #self.Ymin = int(self.Ymin)-1
        #self.Ymax = int(self.Ymax)+1
        
        self.Ymin = -2.5
        self.Ymax =  2.5
        
        self.Xmax = max(self.Xmax_list)
        self.Xmin = min(self.Xmin_list)
        #self.Xmax = self.Xmax - self.Xmax/20
        #self.Xmin = self.Xmin - self.Xmin/20
        
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin


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
            text = '{:^6.2f}'.format((i/self.x_major_ticks)*self.deltaX + self.Xmin)
            
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
        
        self.x_box_size = self.width  -(self.bx*1.5)
        self.y_box_size = self.height -(self.by*1.8)
        
        

        
        self.draw_box (cr, line_width = 3.0 , color = [1, 1, 1] )
        self.draw_XY_axes(cr, line_width = 3.0 , color = [1, 1, 1] )
        
        
        print(
        '\nself.Ymax     ', self.Ymax      ,
        '\nself.Ymin     ', self.Ymin      ,
        '\nself.Xmax     ', self.Xmax      ,
        '\nself.Xmin     ', self.Xmin      ,
        '\nself.deltaY   ', self.deltaY    ,
        '\nself.deltaX   ', self.deltaX    ,
        )
        
        

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
        
    def on_motion(self, widget, event):
        '''(i/self.x_major_ticks)*self.deltaX + self.Xmin'''
        
        (x, y) = int(event.x), int(event.y)
        
        print("Mouse moved to:", 
               x, y,  
               x-self.bx, y-self.by, 
               
               ((((x-self.bx)) / self.x_box_size)  *  self.deltaX + self.Xmin),
               
               ((self.y_box_size-(y-self.by)) / self.y_box_size)  *  self.deltaY + self.Ymin)
        
    def on_mouse_button_press (self, widget, event):
        """ Function doc """
        self.Ymin = -1.5 -1
        self.Ymax =  1.5 +1
        #print(self.Ymin, self.Ymax)
        #self.define_xy_limits ()
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin
        
        
        
        
        
        
class MyWindow(Gtk.Window):


    def __init__(self):
        Gtk.Window.__init__(self, title="Cairo and GTK Example")
        x = 700
        y = 500
        self.set_default_size(x, y)
        self.plot = XYPlot()
        
        '''        
        x1 = np.linspace(0.0, 5.0)
        y1 = np.cos(2 * np.pi * x1) * np.exp(-x1)
        x2 = np.linspace(0.0, 2.0)
        y2 = np.cos(2 * np.pi * x2)
        '''
        
        
        
        #'''
        x1    = []
        y1 = []
        for i in range( -100 , 100, 1 ):
            x1.append(i/10)
            y1.append(math.sin(i/10 ))
        #'''
        
        self.plot.add( X = x1, Y = y1, symbol = 'dot', line = 'solid', sym_color = [1,0,1], line_color = [1,0,1])
        
        #'''
        self.raw_X    = []
        self.raw_data = []
        for i in range( -100 , 100, 1 ):
            self.raw_X.append(i/10)
            self.raw_data.append(2*math.cos(i/10 ))

        self.plot.add( X = self.raw_X, Y = self.raw_data, symbol = 'dot', sym_fill = False, sym_color = [0,1,1], line = 'solid', line_color = [0,1,1] )
        #'''
        
        self.plot.x_minor_ticks = 6
        self.plot.x_major_ticks = 10
        self.plot.y_minor_ticks = 6
        self.plot.y_major_ticks = 5
        
        
        
        
        
        
        
        
        #self.drawingarea.connect("draw", self.on_draw)
        #self.drawingarea.connect("button_press_event", self.on_motion2)
        #self.drawingarea.connect("motion-notify-event", self.on_motion)

        self.add(self.plot)
    




class MyWindow2(Gtk.Window):


    def __init__(self):
        Gtk.Window.__init__(self, title="Cairo and GTK Example")

        self.counter = 3
        self.raw_X    = []
        self.raw_data = []
        

        for i in range( -100 , 100, 1 ):
            self.raw_X.append(i/10)
            self.raw_data.append(math.sin(i/10 ))
            #self.raw_data.append((i**3))
            #print(math.sin(i/50))
            #print (i/10, math.sin(i/100))
            #print (i , i**3/10)
        
 
        x = 700
        y = 500
        self.set_default_size(x, y)


        self.points = []
        
        self.x_minor_ticks = 6
        self.x_major_ticks = 5
        
        self.y_minor_ticks = 4
        self.y_major_ticks = 10
        
        
        

        self.X    = self.raw_X   .copy() 
        self.data = self.raw_data.copy() 
        #----------------------------------------------------------------------
        self.Ymax    = max(self.data)
        self.Ymin    = min(self.data)
        
        self.Xmax    = max(self.X)
        self.Xmin    = min(self.X)
        
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin
        #----------------------------------------------------------------------
        
        self.data.append(self.Ymax + self.deltaY/20)
        self.data.append(self.Ymin - self.deltaY/20)
        
        self.X.append(1)
        self.X.append(1)
        
        
         
         
        
        
        
        
        self.drawingarea = Gtk.DrawingArea()
        self.drawingarea.set_events(Gdk.EventMask.ALL_EVENTS_MASK)

        self.drawingarea.connect("draw", self.on_draw)
        self.drawingarea.connect("button_press_event", self.on_motion2)
        self.drawingarea.connect("motion-notify-event", self.on_motion)

        self.add(self.drawingarea)


        
    def on_motion2(self, widget, event):
        (x, y) = int(event.x), int(event.y)
        
        #x_on_plot = (x-self.bx)/self.factor_x
        #y_on_plot = (y-self.by)/self.factor_y
        #self.data.append(x)
        #self.X.append(x)
        
        print (self.raw_X)
        print (self.raw_data)
        self.raw_X   .append(self.counter )
        self.raw_data.append(self.counter**3 )
        print (self.raw_X)
        print (self.raw_data)
        self.counter +=1
        
        
        
        self.X    = self.raw_X   .copy()
        self.data = self.raw_data.copy()
        #----------------------------------------------------------------------
        self.Ymax    = max(self.data)
        self.Ymin    = min(self.data)
        
        self.Xmax    = max(self.X)
        self.Xmin    = min(self.X)
        
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin
        #----------------------------------------------------------------------
        
        self.data.append(self.Ymax + self.deltaY/20)
        self.data.append(self.Ymin - self.deltaY/20)
        
        #self.X.append(1)
        #self.X.append(1)
        
    
    def on_motion(self, widget, event):
        '''(i/self.x_major_ticks)*self.deltaX + self.Xmin'''
        
        (x, y) = int(event.x), int(event.y)
        
        print("Mouse moved to:", 
               x, y,  
               x-self.bx, y-self.by, 
               
               ((((x-self.bx)) / self.x_box_size)  *  self.deltaX + self.Xmin),
               
               ((self.y_box_size-(y-self.by)) / self.y_box_size)  *  self.deltaY + self.Ymin)



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
            text = '{:^6.2f}'.format((i/self.x_major_ticks)*self.deltaX + self.Xmin)
            
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
        
        self.x_box_size = self.width  -(self.bx*1.5)
        self.y_box_size = self.height -(self.by*1.8)
        
        

        
        self.draw_box (cr, line_width = 3.0 , color = [1, 1, 1] )
        
        
        self.draw_XY_axes(cr, line_width = 3.0 , color = [1, 1, 1] )
        
        



        #----------------------------------------------------------------------
        self.Ymax    = max(self.data)
        self.Ymin    = min(self.data)
        
        self.Xmax    = max(self.X)
        self.Xmin    = min(self.X)
        

        #----------------------------------------------------------------------
        print(
        '\nself.Ymax     ', self.Ymax      ,
        '\nself.Ymin     ', self.Ymin      ,
        '\nself.Xmax     ', self.Xmax      ,
        '\nself.Xmin     ', self.Xmin      ,
        '\nself.deltaY   ', self.deltaY    ,
        '\nself.deltaX   ', self.deltaX    ,
        )



        #self.Yfactor = (self.y_box_size )/self.deltaY
        #self.Xfactor = (self.x_box_size )/self.deltaX
        if self.deltaY > 0:
            self.norm_data = [(y-self.Ymin)/self.deltaY for y in self.data]
            self.norm_X    = [(x-self.Xmin)/self.deltaX for x in self.X]
        else:
            self.norm_data = self.data
            self.norm_X    = self.X 
            
            
        self.Ynorm_max = max(self.norm_data)



        print(
        '\nself.Ymax     ', self.Ymax      ,
        '\nself.Ymin     ', self.Ymin      ,
        '\nself.Xmax     ', self.Xmax      ,
        '\nself.Xmin     ', self.Xmin      ,
        '\nself.deltaY   ', self.deltaY    ,
        '\nself.deltaX   ', self.deltaX    ,
        #'\nself.norm_data', self.norm_data ,
        #'\nself.norm_X   ', self.norm_X    ,
        '\nself.Ynorm_max', self.Ynorm_max ,
        )


        self.draw_XY_scale ( cr, line_width = 3.0 , color = [1, 1, 1] )





        #----------------------------------------------------------------------
        cr.set_line_width (2.0)
        #----------------------------------------------------------------------
        cr.set_source_rgb( 1, 1, 1)
        cx = self.bx - self.norm_X[0]
        #----------------------------------------------------------------------
        
        for i in range(len(self.norm_data)-2):           
            x = self.norm_X[i]    
            y = self.norm_data[i]

            cr.line_to (x*self.x_box_size + cx, 
                       #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                       (1*self.y_box_size + self.by) - (y*self.y_box_size))
        cr.stroke ()
        
        cr.set_source_rgb( 0, 1, 1)
        for i in range(len(self.norm_data)-2):           
            x = self.norm_X[i]    
            y = self.norm_data[i]

            cr.arc (x*self.x_box_size + cx, 
                       #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                       (1*self.y_box_size + self.by) - (y*self.y_box_size), 5, 0, 2*3.14)
            cr.fill()

        #cr.stroke ()

win = MyWindow()
win.connect("delete-event", Gtk.main_quit)
win.show_all()
Gtk.main()

