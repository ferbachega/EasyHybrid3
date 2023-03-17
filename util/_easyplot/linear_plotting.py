import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, cairo
import  random  
import numpy as np
import math


    
    
def interpolate_color(color1, color2, fraction):
    return [ ( (x1 + fraction * (x2 - x1))) for (x1, x2) in zip(color1, color2)]


def get_color(valor, factor = 2):   
    
    
    valor = (valor+1)/factor   
    
    cutoffs = list(color_map.keys())
    color_list = list(color_map.values())
    counter = 0 
    
    for cutoff, color in color_map.items():
        if valor <= cutoff :
            fraction = (valor -  cutoffs[counter-1]) / ( cutoff - cutoffs[counter-1])
            color1  = color_list[counter-1]
            color2  = color_list[counter]
            color = interpolate_color(color1, color2, fraction)
            return color
        else:
            pass
        counter +=1




class MyWindow(Gtk.Window):


    def __init__(self):
        Gtk.Window.__init__(self, title="Cairo and GTK Example")


        self.X =  []
        self.data = []
        
        '''
        for i in range(-50, 100):
            self.X.append(i/10)
            self.data.append(math.sin(i/10))
            #print(math.sin(i/50))
            print (i/10, math.sin(i/10))
        '''
        for i in range(-550, 550, 50 ):
            self.X.append(i/10)
            #self.data.append(math.sin(i/10))
            self.data.append(i**3/10)
            #print(math.sin(i/50))
            print (i/10, math.sin(i/10))
        
 
        self.x = 700
        self.y = 500
        self.set_default_size(self.x, self.y)


        self.points = []
        
        self.x_minor_ticks = 5
        self.x_major_ticks = 10
        
        self.y_minor_ticks = 5
        self.y_major_ticks = 10
        
        self.drawingarea = Gtk.DrawingArea()
        self.drawingarea.set_events(Gdk.EventMask.ALL_EVENTS_MASK)

        self.drawingarea.connect("draw", self.on_draw)
        self.drawingarea.connect("button_press_event", self.on_motion2)
        self.add(self.drawingarea)


        
    def on_motion2(self, widget, event):
        (x, y) = int(event.x), int(event.y)
        
        x_on_plot = (x-self.bx)/self.factor_x
        y_on_plot = (y-self.by)/self.factor_y
        

        
        if int(x_on_plot) < 0 or int(x_on_plot) >= self.size_x :
            print('canvas')
            self.points = []
        elif int(y_on_plot) < 0 or int(y_on_plot) >= self.size_y :
            print('canvas')
            self.points = []
        else:
            self.points.append((x, y))
        
        print("Mouse clicker at:",  x, y, int(x_on_plot), int(y_on_plot), 
                                    (x-self.bx)/self.factor_x, #/(self.x_final-self.bx), 
                                    (y-self.by)/self.factor_y) #/(self.y_final-self.by) ) #, x - self.x_final , y - self.y_final , self.x_final ,  self.y_final ,   (x-self.bx) /self.factor_x ,    (y-self.by)/self.factor_y )
        
        #self.cr.move_to (x, y)
        #self.cr.rectangle(x, y, 3, 3)
        #self.on_draw()
        self.drawingarea.queue_draw()
    
    def on_motion(self, widget, event):
        (x, y) = int(event.x), int(event.y)
        #x, y = device.get_position(widget)
        print("Mouse moved to:", x, y)



    def draw_box (self, cr):
        """ Function doc """
        borda = 20 
        

        cr.set_line_width (1.0)
        
        for j in range(0,self.size_y, 1):
            
            #--------------------------------------------------------------------------------------------
            #color bar
            #--------------------------------------------------------------------------------------------
            z = ((j  / float(self.size_y)) - 0.5)*2
            #print (z)
            color = get_color(z)
            cr.set_source_rgb( color[0], color[1], color[2])
            
            j = self.size_y -j-1
            
            cr.rectangle(
                         round(self.size_x * self.factor_x)+self.bx +borda     ,
                         self.by + j*self.factor_y  , 
                         
                         round(20)+1, 
                         
                         round(self.factor_y)+1)
            cr.fill()
        
        for j in range(0,self.size_y, 6):    
            
            j = self.size_y -j-1
            
            z = ((j  / float(self.size_y)) - 0.5)*2
            #print(z)
            cr.set_source_rgb( 0, 0,0)
            cr.move_to (round(self.size_x * self.factor_x)+self.bx*1.2 +40      , 
                              self.by + j*self.factor_y + (self.factor_y)/2 )
            
            
            cr.line_to (round(self.size_x * self.factor_x)+self.bx*1.2 +25  , self.by + j*self.factor_y + (self.factor_y)/2 ) 
            #--------------------------------------------------------------------------------------------
            
            cr.stroke ()
            #--------------------------------------------------------------------------------------------





        cr.set_line_width (2.0)
        cr.rectangle(round(self.size_x * self.factor_x)+self.bx +borda,
                     self.by   , 
                     round(20+1), 
                     round(self.size_y * self.factor_y))
        cr.stroke ()


    def on_draw(self, widget, cr):
        self.cr =  cr
        #print('aqui', cr)
        self.width = widget.get_allocated_width()
        self.height = widget.get_allocated_height()
        
        
        #cr.set_source_rgb( 1, 1, 1)
        cr.set_source_rgb( 0, 0, 0)
        cr.paint()

        self.x = self.width
        self.y = self.height


        
        self.bx = 60 
        self.by = 50 
        


        
        
        
        
        
        #bordas inferiores - garante que todo o plot cabe na janela 
        self.x_final = self.width  -(self.bx*2)
        self.y_final = self.height -(self.by*2)
        
        
        
        self.size_y = len(self.data)
        self.size_x = len(self.X)
        
        
        self.Ymax    = max(self.data)
        self.Ymin    = min(self.data)
        
        self.Xmax    = max(self.X)
        self.Xmin    = min(self.X)
        
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin
        
        
        
        self.Yfactor = (self.y_final - self.by)/self.deltaY
        self.Xfactor = (self.x_final - self.bx)/self.deltaX
        
        
        self.new_data = [x*self.Yfactor for x in self.data]
        self.new_X    = [x*self.Xfactor for x in self.X]
        
        print(
        '\n self.size_x :', self.size_x ,
        '\n self.size_y :', self.size_y ,
        
        '\n self.Ymax   :', self.Ymax   ,
        '\n self.Ymin   :', self.Ymin   ,

        '\n self.deltaX :', self.deltaX ,
        '\n self.deltaY :', self.deltaY ,

        '\n self.Xfactor:', self.Xfactor,
        '\n self.Yfactor:', self.Yfactor,

        '\n self.Xmax_new:', max(self.X) ,
        '\n self.Xmin_new:', min(self.X) ,
        
        '\n self.Ymax_new:', max(self.new_data) ,
        '\n self.Ymin_new:', min(self.new_data) ,
        )
        
        
        self.new_Ymax    = max(self.new_data)
        self.new_Ymin    = min(self.new_data)
        self.delta = self.new_Ymax - self.new_Ymin
        

     
        
        
        
        
        #'''
        
        #----------------------------------------------------------------------
        cr.set_source_rgb( 1, 1, 1)
        cx = self.bx-self.new_X[0]

        cy = self.deltaY/2 +self.by
        #----------------------------------------------------------------------
        for i in range(len(self.new_data)):           
            x = self.new_X[i]    
            y = self.new_data[i]

            cr.line_to (x + cx + self.bx/2, 
                       #(self.new_Ymax + self.y_final + self.by ) - (y+self.y_final))
                       (self.new_Ymax + self.y_final + self.by + self.by/2) - (y+self.y_final))
        cr.stroke ()
        #----------------------------------------------------------------------
        
        
        
        
        #----------------------------------------------------------------------
        for i in range(len(self.new_data)):
            x = self.new_X[i]    
            y = self.new_data[i]
            
            
            cr.arc(x + cx + self.bx/2, 
                   (self.new_Ymax+self.y_final+self.by+self.by/2) - (y+self.y_final), 5, 0, 2*3.14)
            cr.fill()
            
            
            
            cr.arc(cx + self.bx/2 , 
                   (self.new_Ymax+self.y_final+self.by+self.by/2) - (y+self.y_final), 5, 0, 2*3.14)
            cr.fill()
        
            
            
            cr.arc(x+cx + self.bx/2 , 
                  (self.y_final+self.by), 5, 0, 2*3.14)
            cr.fill()
        #----------------------------------------------------------------------
        
        
        
        
        #'''
        #cr.stroke ()
        #----------------------------------------------------------------------
            
        cr.set_source_rgb( 1, 1, 1)
        '''
        #----------------------------------------------------------------------
        #retangle
        cr.set_source_rgb( 1, 1, 1)
        cr.set_line_width (1.0)
        cr.rectangle(self.bx,
                     self.by, 
                     self.x_final,
                     self.y_final)

        cr.stroke ()
        #----------------------------------------------------------------------
        '''      
        #self.x_minor_ticks = 5
        #self.x_major_ticks = 10
        
        
        cr.move_to (self.bx, self.y_final+self.by )
        cr.line_to (self.bx+ self.x_final, self.y_final+self.by)
        cr.stroke ()
        
        x_major_factor = (self.bx+ self.x_final) / float(self.x_major_ticks)
        x_minor_factor = x_major_factor / float(self.x_minor_ticks)
        
        
        for i in range(0, self.x_major_ticks):
            # marcacoes em x
            #print(self.bx+ i*self.factor_x + (self.factor_x)/2.0 , self.bx+ i*self.factor_x + (self.factor_x)/2.0)
            cr.move_to (self.bx+ i*x_major_factor, self.y_final+self.by )
            cr.line_to (self.bx+ i*x_major_factor, self.y_final+self.by + 15 ) 
            #--------------------------------------------------------------------------------------------
            for j in range(0, self.x_minor_ticks):
                cr.move_to (self.bx+ i*x_major_factor+j*x_minor_factor, self.y_final+self.by )
                cr.line_to (self.bx+ i*x_major_factor+j*x_minor_factor, self.y_final+self.by + 10 )
            
        cr.stroke ()
            
        '''
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
        '''
        

win = MyWindow()
win.connect("delete-event", Gtk.main_quit)
win.show_all()
Gtk.main()

