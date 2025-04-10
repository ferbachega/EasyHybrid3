import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
import cairo
import math

class RamachandranWindow(Gtk.Window):
    def __init__(self):
        super().__init__(title="Diagrama de Ramachandran")
        self.set_default_size(500, 500)
        self.set_border_width(10)

        self.drawing_area = Gtk.DrawingArea()
        self.drawing_area.connect("draw", self.on_draw)
        self.add(self.drawing_area)
        
        self.connect("destroy", Gtk.main_quit)
        self.show_all()
        
    def to_canvas_coords(self, phi, psi, origin_x, origin_y, size):
        """Converte coordenadas phi, psi (-180 a 180) para canvas (x, y)"""
        x = origin_x + ((phi + 180) / 360) * size
        y = origin_y + size - ((psi + 180) / 360) * size
        return x, y

    def draw_region(self, cr, region_points, origin_x, origin_y, size, color):
        cr.set_source_rgba(*color)
        cr.set_line_width(1.0)
        first = True
        for phi, psi in region_points:
            x, y = self.to_canvas_coords(phi, psi, origin_x, origin_y, size)
            if first:
                cr.move_to(x, y)
                first = False
            else:
                cr.line_to(x, y)
        cr.close_path()
        cr.fill()

    def on_draw(self, widget, cr):
        width = widget.get_allocated_width()
        height = widget.get_allocated_height()
        cr.set_source_rgb(1, 1, 1)  # fundo branco
        cr.paint()

        # Área de desenho
        margin = 40
        size = min(width, height) - 2 * margin
        origin_x = margin
        origin_y = margin
        step = size / 12  # 30° por passo

        # Grade
        cr.set_line_width(1)
        cr.set_source_rgb(0.8, 0.8, 0.8)
        for i in range(13):
            x = origin_x + i * step
            y = origin_y + i * step
            cr.move_to(x, origin_y)
            cr.line_to(x, origin_y + size)
            cr.move_to(origin_x, y)
            cr.line_to(origin_x + size, y)
        cr.stroke()

        # Eixos φ = 0 e ψ = 0
        cr.set_source_rgb(0, 0, 0)
        cr.set_line_width(2)
        cr.move_to(origin_x + size / 2, origin_y)
        cr.line_to(origin_x + size / 2, origin_y + size)
        cr.move_to(origin_x, origin_y + size / 2)
        cr.line_to(origin_x + size, origin_y + size / 2)
        cr.stroke()

        # Rótulos
        cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        cr.set_font_size(10)
        for i in range(13):
            label = str(-180 + i * 30)
            x = origin_x + i * step
            y = origin_y + size + 15
            cr.move_to(x - 10, y)
            cr.show_text(label)
            cr.move_to(origin_x - 30, origin_y + size - i * step + 3)
            cr.show_text(label)

        # Regiões permitidas (dados aproximados)
        alpha_region = [
            (-80, -60), (-70, -40), (-60, -30), (-50, -40), (-40, -60),
            (-50, -70), (-60, -80), (-70, -70), (-80, -60),# ( 100,  111), ( 150,  151)
        ]
        beta_region = [
            (-160, 120), (-150, 130), (-140, 140), (-130, 135),
            (-120, 120), (-130, 110), (-140, 110), (-150, 115), (-160, 120)
        ]
        l_alpha_region = [
            (50, 30), (60, 40), (70, 30), (70, 20), (60, 10), (50, 20), (50, 30)
        ]

        self.draw_region(cr, alpha_region, origin_x, origin_y, size, (0.5, 1.0, 0.5, 0.5))  # verde
        self.draw_region(cr, beta_region, origin_x, origin_y, size, (0.5, 0.5, 1.0, 0.5))   # azul
        self.draw_region(cr, l_alpha_region, origin_x, origin_y, size, (1.0, 0.5, 0.5, 0.5))  # vermelho claro

        # Pontos de exemplo
        #pontos = [(-60, -45), (60, 120), (-135, 135), (0, 0)]
        pontos = [ ]
        cr.set_source_rgb(1, 0, 0)
        for phi, psi in pontos:
            x, y = self.to_canvas_coords(phi, psi, origin_x, origin_y, size)
            cr.arc(x, y, 4, 0, 2 * math.pi)
            cr.fill()

#win = RamachandranWindow()
#win.connect("destroy", Gtk.main_quit)
#win.show_all()
#Gtk.main()
