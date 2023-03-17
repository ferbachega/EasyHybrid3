import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, cairo
import  random  
import numpy as np
import math
'''
color_map = {
        0	                  : [0                 	,0	                ,0                ],
        0.142857142857143     : [0.02269028242686	,0.049634617798597	,0.469850365435695],
        0.285714285714286     : [0.014712149007199	,0.304056764250118	,0.208741430834425],
        0.428571428571429     : [0.189397609083558	,0.454376493558188	,0.021637600293833],
        0.571428571428571     : [0.891586421514143	,0.389697937913429	,0.042785205791797],
        0.714285714285714     : [0.978214236714738	,0.544449995205519	,0.791623070953407],
        0.857142857142857     : [0.884026030284464	,0.807586693400741	,0.990677803425028],
        1.0                   : [0.990677803425028	,0.990677803425028	,0.990677803425028],
}
#'''
'''
color_map = {
        0	                    : [0        ,0       ,100/255 ],
        0.175	                : [0        ,113/255 ,197/255 ],       
        0.25	                : [0        ,217/255 ,236/255 ],       
        0.5                     : [0.00000  ,255/255 ,0.0     ],
        0.75                    : [255/255  ,219/255 ,0.0     ],
        0.87                    : [255/255  ,0       ,0.0     ],
        1.0                     : [100/255  ,0.0	 ,0.0     ],                  
}

#'''
#jet

#'''
color_map = {
        0	                    : [0        ,0       , 100/255 ],
        0.2  	                : [0        ,100/255  ,191/255 ],       
        0.4 	                : [0        ,222/255 ,171/255 ],       
        0.6                     : [173/255  ,220/255 ,0.0     ],
        0.8                     : [191/255  , 100/255 ,0.0     ],
        1.0                     : [100/255  ,0.0	 ,0.0     ],                  
}

#'''

'''
color_map = {
 0      :[  0.5199999415509024        ,  0.23784035734339465       ,    0.014399781237605947     ],
 0.125  :[      0.8199999347340485    ,      0.33784033341279346   ,        0.016399661791829116 ],
 0.25   :[     0.8999999297090052     ,     0.43965031326255594    ,       0.00899952568157918   ],
 0.375  :[      0.929999720711983     ,     0.5508316161641296     ,      0.04649971339115546    ],
 0.5    :[    0.9699995785589836      ,    0.6534571258568723      ,     0.10669978434506405     ],
 0.625  :[      0.999999440760131     ,     0.7300005801808864     ,      0.18999979351538043    ],
 0.75   :[     0.9999996007480701     ,     0.8110004602574017     ,      0.3699998211951647     ],
 0.875  :[      0.9999998170281835    ,      0.8433336609371632    ,       0.5299999027260401    ],
 1      :[  0.9999995530516227        ,  0.9631665995878685        ,   0.8699996387296692        ],
}
#'''

def interpolate_color(color1, color2, fraction):
    return [ ( (x1 + fraction * (x2 - x1))) for (x1, x2) in zip(color1, color2)]


def get_color(valor, factor = 2):   
    
    
    valor = (valor+1)/factor   
    
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




class MyWindow(Gtk.Window):


    def __init__(self):
        Gtk.Window.__init__(self, title="Cairo and GTK Example")
        

        self.x = 700
        self.y =  500
        self.set_default_size(self.x, self.y)

        self.bx = 50 
        self.by = 50 
        
        self.size_x = 14 
        self.size_y = 16 
        
        #bordas inferior - garante que todo o plot cabe na janela 
        self.x_final = self.x -self.bx
        self.y_final = self.y -self.by
        
        self.factor_x = int((self.x_final - self.bx)/self.size_x)
        self.factor_y = int((self.y_final - self.by)/self.size_y)
        
        print (
        
        '\n bx       ',  self.bx      ,  
        '\n by       ',  self.by      ,  
        '\n size_x   ',  self.size_x  ,  
        '\n size_y   ',  self.size_y  ,  
        '\n x_final  ',  self.x_final ,  
        '\n y_final  ',  self.y_final ,  
        '\n factor_x ',  self.factor_x,  
        '\n factor_y ',  self.factor_y,   
        )
        
        self.points = []
        
        
        self.drawingarea = Gtk.DrawingArea()
        self.drawingarea.set_events(Gdk.EventMask.ALL_EVENTS_MASK)

        self.drawingarea.connect("draw", self.on_draw)
        #self.drawingarea.connect("motion-notify-event", self.on_motion)
        self.drawingarea.connect("button_press_event", self.on_motion2)
        #hbox.pack_start(self.drawingarea, True, True, 0)
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
        print('aqui', cr)
        width = widget.get_allocated_width()
        height = widget.get_allocated_height()
        
        
        

        
        
        
        
        
        
        
        
        
        
        cr.set_source_rgb( 0.5,0.5  , 0.5 )
        #cr.set_source_rgb( 0, 0, 0)
        cr.paint()
        #cr.set_source_rgb( 1, 1, 1)
        #cr.rectangle(0,
        #             0, 
        #             width, 
        #             height)
        #cr.fill()

        
        #print(width, height)
        self.x = width
        self.y = height


        self.bx = 50 
        self.by = 50 
        
        self.size_x = 100
        self.size_y = 100
        
        #bordas inferior - garante que todo o plot cabe na janela 
        self.x_final = self.x -(self.bx+40)
        self.y_final = self.y -(self.by)
        
        self.factor_x = ((self.x_final - self.bx)/self.size_x)
        self.factor_y = ((self.y_final - self.by)/self.size_y)

        self.x_final = self.size_x*self.factor_x-self.factor_x
        self.y_final = self.size_y*self.factor_y-self.factor_y



        
        
        self.draw_box (cr)
        
        
        '''
        for j in range(0,self.size_y, 1):
            #--------------------------------------------------------------------------------------------
            z = (j  / float(self.size_y)) - 0.5
            
            color = get_color(z)
            cr.set_source_rgb( color[0], color[1], color[2])
            #print(z, color)
            # marcacoes em y
            cr.rectangle(
                         round(self.size_x * self.factor_x)+self.bx*1.2      ,
                         self.by + j*self.factor_y  , 
                         
                         round(20)+1, 
                         
                         round(self.factor_y)+2)
            cr.fill()
        
        for j in range(0,self.size_y, 5):    
            
            cr.set_source_rgb( 0, 0,0)
            cr.move_to (round(self.size_x * self.factor_x)+self.bx*1.2 +40      , 
                              self.by + j*self.factor_y + (self.factor_y)/2 )
            
            
            cr.line_to (round(self.size_x * self.factor_x)+self.bx*1.2 +25  , self.by + j*self.factor_y + (self.factor_y)/2 ) 
            #--------------------------------------------------------------------------------------------
            
            cr.stroke ()
            #--------------------------------------------------------------------------------------------
        '''










        n = 0 
        for i in range(0,self.size_x, ):
            for j in range(0,self.size_y):
                
                z = math.sin(i*3.1415/(self.size_x/3)) * math.sin(j*3.1415/(self.size_y/2))
                color = get_color(z)
                #print(z, color)
                
                
                #cr.set_source_rgb( i/float(self.size_x), j/float(self.size_x) , ((i*j)**0.5) /float(self.size_x) )
                
                #cr.set_source_rgb( z, z*-1 , z  )
                #cr.set_source_rgb( z*-1, z  , z-1   )
                #cr.set_source_rgb( 1, 1  , 1   )
                cr.set_source_rgb( color[0], color[1], color[2]   )
                
                cr.rectangle(self.bx+i*self.factor_x -1,
                             self.by+j*self.factor_y -1, 
                             round(self.factor_x)+1, 
                             round(self.factor_y)+1)
                cr.fill()
                

        cr.set_line_width (3.0)
        cr.rectangle(self.bx,
                     self.by, 
                     round(self.size_x * self.factor_x), 
                     round(self.size_y * self.factor_y))
        
        #drawing axes x ,y 
        cr.set_source_rgb( 0, 0, 0)
        #cr.set_line_width (2.0)
        #
        #cr.move_to (self.bx      , self.by    )        
        #cr.line_to (self.x_final , self.by ) 
        #
        #cr.move_to (self.bx , self.by    )
        #cr.line_to (self.bx , self.y_final ) 
        
        
        font_size = 18
        
        
        
        for i in range(0,self.size_x, 5 ):
            #--------------------------------------------------------------------------------------------
            # marcacoes em x
            cr.move_to (self.bx+ i*self.factor_x + (self.factor_x)/2      , self.by )
            cr.line_to (self.bx+ i*self.factor_x + (self.factor_x)/2      , self.by -10 ) 
            #--------------------------------------------------------------------------------------------
            
            
            
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
        
        for j in range(0,self.size_y, 10):
            #--------------------------------------------------------------------------------------------
            # marcacoes em y
            cr.move_to (self.bx      , self.by + j*self.factor_y + (self.factor_y)/2 )
            cr.line_to (self.bx -10  , self.by + j*self.factor_y + (self.factor_y)/2 ) 
            #--------------------------------------------------------------------------------------------

            
            
            #--------------------------------------------------------------------------------------------
            # texto em y
            cr.set_source_rgb(0, 0, 0)
            cr.set_font_size(font_size)
    
            text = str(j)
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)

            cr.move_to(self.bx  -15 -width     , self.by + j*self.factor_y + height/2  + (self.factor_y)/2)
            cr.show_text(str(j))
            #--------------------------------------------------------------------------------------------
        


        
        
        

            
        



        cr.stroke ()
        
        for x, y in self.points:
            cr.arc(x , y , 5, 0, 2 * 3.14)
            cr.fill()

win = MyWindow()
win.connect("delete-event", Gtk.main_quit)
win.show_all()
Gtk.main()

#import math
#import cairo
#
#WIDTH, HEIGHT = 256, 256
#
#surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
#ctx = cairo.Context(surface)
#
#ctx.scale(WIDTH, HEIGHT)  # Normalizing the canvas
#
#pat = cairo.LinearGradient(0.0, 0.0, 0.0, 1.0)
#pat.add_color_stop_rgba(1, 0.7, 0, 0, 0.5)  # First stop, 50% opacity
#pat.add_color_stop_rgba(0, 0.9, 0.7, 0.2, 1)  # Last stop, 100% opacity
#
#ctx.rectangle(0, 0, 1, 1)  # Rectangle(x0, y0, x1, y1)
#ctx.set_source(pat)
#ctx.fill()
#
#ctx.translate(0.1, 0.1)  # Changing the current transformation matrix
#
#ctx.move_to(0, 0)
## Arc(cx, cy, radius, start_angle, stop_angle)
#ctx.arc(0.2, 0.1, 0.1, -math.pi / 2, 0)
#ctx.line_to(0.5, 0.1)  # Line to (x,y)
## Curve(x1, y1, x2, y2, x3, y3)
#ctx.curve_to(0.5, 0.2, 0.5, 0.4, 0.2, 0.8)
#ctx.close_path()
#
#ctx.set_source_rgb(0.3, 0.2, 0.5)  # Solid color
#ctx.set_line_width(0.02)
#ctx.stroke()
#
#surface.write_to_png("example.png")  # Output to PNG
