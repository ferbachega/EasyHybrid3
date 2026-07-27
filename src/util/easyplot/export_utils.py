"""
Exportacao de figuras em alta resolucao para os plots de
util/easyplot/ (image_plot.ImagePlot, xy_plots.XYPlot).

[EN] Ambos os widgets sao Gtk.DrawingArea com um metodo
on_draw(self, widget, cr, export_scale=1) que desenha usando Cairo puro
-- nada de matplotlib. export_plot_to_png() reusa esse MESMO codigo de
desenho (nao tira "screenshot" do widget na tela) para renderizar
diretamente numa superficie Cairo maior, com o parametro export_scale
props (ver o docstring de cada on_draw() para o motivo de cada classe
precisar de um tratamento levemente diferente internamente: ImagePlot
desenha primeiro numa superficie raster interna de tamanho fixo,
XYPlot desenha direto no cr recebido). O resultado e' um redesenho
NITIDO na resolucao pedida, nao um upscale borrado de uma imagem
pequena.
"""
import os
import cairo


def export_plot_to_png(plot_widget, filepath, scale=4):
    """ Renderiza plot_widget (um ImagePlot ou XYPlot ja' visivel na
    tela, com layout/dados definidos) para um arquivo PNG em
    'scale' vezes a resolucao normal de tela.

    plot_widget : instancia de ImagePlot ou XYPlot (qualquer
                  Gtk.DrawingArea cujo on_draw aceite export_scale=).
    filepath    : caminho de saida (.png). O diretorio precisa existir.
    scale       : fator de escala em relacao ao tamanho ATUAL do widget
                  na tela (ex.: se o widget esta' com 500x400 pixels na
                  janela, scale=4 gera um PNG de 2000x1600 -- da' pra
                  pensar nisso como "4x a resolucao da tela", parecido
                  com tirar um print em um monitor 4x mais denso).
                  4 costuma dar um resultado adequado para a maioria
                  dos usos (slides, relatorios); para figuras de
                  publicacao/impressao, prefira 6-8.

    Levanta ValueError se o widget ainda nao tem tamanho alocado (ex.:
    a janela nunca foi mostrada/redimensionada) -- nesse caso nao ha'
    largura/altura de referencia para escalar.
    """
    width = plot_widget.get_allocated_width()
    height = plot_widget.get_allocated_height()
    if width <= 0 or height <= 0:
        raise ValueError(
            "export_plot_to_png: o widget do grafico ainda nao tem "
            "tamanho alocado (width={}, height={}) -- mostre a janela "
            "antes de exportar.".format(width, height)
        )

    out_dir = os.path.dirname(os.path.abspath(filepath))
    if out_dir and not os.path.isdir(out_dir):
        raise ValueError("export_plot_to_png: diretorio nao existe: {}".format(out_dir))

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32,
                                  int(width * scale), int(height * scale))
    cr = cairo.Context(surface)

    # widget real (nao um dublê): get_allocated_width/height() continuam
    # retornando o tamanho ORIGINAL da tela -- e' isso que faz as
    # proporcoes (margens, fontes, grade de dados) ficarem identicas as
    # da tela, so' que desenhadas numa superficie fisicamente maior.
    plot_widget.on_draw(plot_widget, cr, export_scale=scale)

    surface.write_to_png(filepath)
    return filepath
