#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#
#  The "Plots" tabs of the log window (2026-10-28): what util/log_analysis.py
#  finds in a job log, drawn with util/analysis_plot.AnalysisPlot -- a summary
#  on top, the plots below (two columns), Save Image... under each plot and
#  Refresh for a job that is still running.

import os

import gi

gi.require_version("Gtk", "3.0")
from gi.repository import GLib, Gtk  # noqa: E402

_KIND_TITLES = {"optimization": "Geometry optimization", "md": "Molecular dynamics", "scan1d": "Relaxed scan (1D)",
                "scan2d": "Relaxed scan (2D)", "neb": "Nudged elastic band", "umbrella": "Umbrella sampling",
                "normal_modes": "Normal modes", "refinement": "Energy refinement", "energy": "Energy"}


def _plot_widget(panel, folder):
    from util.analysis_plot import AnalysisPlot
    plot = AnalysisPlot()
    plot.set_size_request(480, 300)
    if panel["kind"] == "heatmap":
        plot.set_heatmap(panel["matrix"], panel.get("xlabel", ""), panel.get("ylabel", ""), panel["title"],
                         row_labels=panel.get("row_labels"), col_labels=panel.get("col_labels"),
                         value_label=panel.get("value_label", "value"), overlays=panel.get("overlays"))
    else:
        plot.set_lines(panel["series"], panel.get("xlabel", ""), panel.get("ylabel", ""), panel["title"],
                       hline=panel.get("hline"), xlabels=panel.get("xlabels"))
    box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=2)
    box.pack_start(plot, True, True, 0)
    save = Gtk.Button(label="Save Image...")
    save.set_halign(Gtk.Align.END)
    save.set_relief(Gtk.ReliefStyle.NONE)

    def on_save(_button):
        chooser = Gtk.FileChooserDialog(title="Save Image", action=Gtk.FileChooserAction.SAVE,
                                        transient_for=box.get_toplevel())
        chooser.add_buttons("Cancel", Gtk.ResponseType.CANCEL, "Save", Gtk.ResponseType.OK)
        chooser.set_do_overwrite_confirmation(True)
        if folder and os.path.isdir(folder):
            chooser.set_current_folder(folder)
        safe = "".join(c if c.isalnum() else "_" for c in panel["title"])[:40].strip("_") or "plot"
        chooser.set_current_name(safe + ".png")
        if chooser.run() == Gtk.ResponseType.OK:
            plot.export(chooser.get_filename())
        chooser.destroy()

    save.connect("clicked", on_save)
    box.pack_start(save, False, False, 0)
    return box


def _summary_widget(title, rows):
    frame = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
    head = Gtk.Label(xalign=0)
    head.set_markup("<b>{}</b>".format(GLib.markup_escape_text(title)))
    frame.pack_start(head, False, False, 0)
    grid = Gtk.Grid(column_spacing=16, row_spacing=2)
    for i, (label, value) in enumerate(rows):
        a = Gtk.Label(label=label, xalign=1)
        a.get_style_context().add_class("dim-label")
        b = Gtk.Label(label=value, xalign=0, selectable=True)
        b.set_line_wrap(True)
        grid.attach(a, 0, i, 1, 1)
        grid.attach(b, 1, i, 1, 1)
    frame.pack_start(grid, False, False, 0)
    return frame


def build_plots_page(title, summary, panels, folder=None, on_refresh=None, empty_text=None):
    """ A scrolled page: summary, then the panels two per row. """
    outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8, border_width=10)
    top = Gtk.Box(spacing=8)
    top.pack_start(_summary_widget(title, summary), True, True, 0)
    if on_refresh is not None:
        refresh = Gtk.Button(label="Refresh")
        refresh.set_tooltip_text("Read the log again (a job that is still running).")
        refresh.set_valign(Gtk.Align.START)
        refresh.connect("clicked", lambda *_: on_refresh())
        top.pack_end(refresh, False, False, 0)
    outer.pack_start(top, False, False, 0)
    if not panels:
        label = Gtk.Label(label=empty_text or "Nothing to plot in this log.", xalign=0)
        label.get_style_context().add_class("dim-label")
        outer.pack_start(label, False, False, 0)
    grid = Gtk.Grid(row_spacing=10, column_spacing=10, column_homogeneous=True)
    for k, panel in enumerate(panels):
        grid.attach(_plot_widget(panel, folder), k % 2, k // 2, 1, 1)
    outer.pack_start(grid, True, True, 0)
    scrolled = Gtk.ScrolledWindow()
    scrolled.set_hexpand(True)
    scrolled.set_vexpand(True)
    scrolled.add(outer)
    return scrolled


def analysis_tabs(logfile, simulation_type=None):
    """ [(tab label, make_page)] for TabbedLogWindow.add_lazy_tab: the pDynamo
        plots, and the plots of the QC program when there is one. The
        analysis runs once, when the first of these tabs is opened. """
    from util.log_analysis import analyze
    folder = os.path.dirname(logfile)
    cache = {}

    def result():
        if "r" not in cache:
            cache["r"] = analyze(logfile, simulation_type)
        return cache["r"]

    def make_pdynamo(refresh):
        r = result()
        title = _KIND_TITLES.get(r["kind"], "Log") + " — " + os.path.basename(logfile)
        return build_plots_page(title, r["summary"] or [("log", "no table recognized")], r["panels"], folder, refresh,
                                "No table that can be plotted was found in this log.")

    def make_qc(refresh):
        qc = result().get("qc")
        return build_plots_page("{} files next to the log".format(qc["program"]), qc["summary"], qc["panels"],
                                folder, refresh)

    def reset():
        cache.clear()

    tabs = [("Plots", make_pdynamo, reset)]
    # . a quick look at the folder decides whether there is a QC-program tab
    from util.log_analysis import _orca_files, _xtb_files
    if _orca_files(folder):
        tabs.append(("ORCA plots", make_qc, reset))
    elif _xtb_files(folder):
        tabs.append(("xTB plots", make_qc, reset))
    return tabs
