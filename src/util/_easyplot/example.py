import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf
import cairo
import numpy as np

class PlotArea(Gtk.DrawingArea):
    def __init__(self):
        Gtk.DrawingArea.__init__(self)
        self.connect("draw", self.on_draw)
        self.connect("configure-event", self.on_configure)
        
        self.surface = None
        
    def on_configure(self, widget, event):
        width = self.get_allocated_width()
        height = self.get_allocated_height()
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        cr = cairo.Context(self.surface)
        
        # clear the background
        cr.set_source_rgb(1, 1, 1)
        cr.paint()

        # set up the plot area
        cr.set_line_width(2)
        cr.rectangle(0, 0, width, height)
        cr.clip()
        cr.translate(40, height - 40)
        cr.scale(width/100, -height/100)

        # plot x and y axes
        cr.set_source_rgb(0, 0, 0)
        cr.move_to(0, 0)
        cr.line_to(2*np.pi*100, 0)
        cr.move_to(0, -1)
        cr.line_to(0, 1)
        cr.stroke()

        # plot a simple line graph
        x = np.linspace(0, 2*np.pi, 100)
        y = np.sin(x)
        for i in range(100):
            cr.line_to(x[i]*10, y[i]*10)
        cr.stroke()

        self.queue_draw()
        
    def on_draw(self, widget, cr):
        if self.surface is not None:
            cr.set_source_surface(self.surface, 0, 0)
            cr.paint()
        
class MainWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="Line Plot Example")
        self.set_default_size(800, 600)

        plot_area = PlotArea()
        self.add(plot_area)
        
win = MainWindow()
win.connect("delete-event", Gtk.main_quit)
win.show_all()
Gtk.main()
