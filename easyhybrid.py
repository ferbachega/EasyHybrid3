#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2025 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      Provides functions for selecting atoms and residues in pDynamo systems
#      to facilitate QM/MM partitioning and molecular simulations.
#
#import gc
#gc.disable()          # teste de diagnostico: desliga o coletor de lixo automatico
#gc.set_threshold(0)   # garante que nao ha coleta automatica por contagem

from _version import EASYHYBRID_VERSION

import os, sys, time, re, json
import logging
import traceback

# --- Fix: engasgo de rotacao em GPUs integradas Intel (driver Mesa/GLX) ---
# Em GPUs integradas Intel com driver Mesa, um frame que atrasa levemente
# (jitter normal do loop do GTK/Python) perde a janela de vblank e o
# driver trava o frame seguinte ate o PROXIMO vblank, dobrando a latencia
# daquele frame. Isso aparece como engasgo aleatorio durante rotacao/pan/
# zoom do mouse, mesmo com o custo de render() baixo e estavel (medido e
# confirmado em testes). Desligar o vsync (vblank_mode=0) resolve, mas so
# fazemos isso quando detectamos uma GPU Intel: em GPUs NVIDIA o problema
# nao ocorre, entao nao vale a pena aceitar o tearing la sem necessidade.
# CRITICO: essa variavel de ambiente so tem efeito se for definida ANTES
# do contexto GL ser criado pelo Mesa - por isso essa deteccao roda aqui,
# no topo do arquivo, antes de qualquer import do GTK/OpenGL.
#
# O usuario pode sobrepor esse auto-detect via Preferences > Startup >
# "V-Sync (vblank_mode)" (ver setup_interface.py/.glade), que grava a
# preferencia ('auto'/'on'/'off') em gl_parameters['vblank_mode'] e persiste
# no .config.json (gui/config.py). Lemos esse arquivo aqui, como JSON cru
# (sem importar gui.config), porque isso precisa rodar antes do GTK ser
# importado -- so tem efeito na PROXIMA vez que o EasyHybrid for aberto.
def _load_vblank_mode_preference():
    """ Le a preferencia 'vblank_mode' salva em .config.json ('auto', 'on'
        ou 'off'). Retorna 'auto' (comportamento padrao/legado) se o
        arquivo nao existir, nao tiver a chave, ou nao puder ser lido.
    """
    try:
        config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".config.json")
        with open(config_path, "r", encoding="utf-8") as config_file:
            data = json.load(config_file)
        return data.get("vblank_mode", "auto")
    except (OSError, ValueError):
        return "auto"

def _maybe_disable_vsync_for_intel_igpu():
    if "vblank_mode" in os.environ:
        return  # usuario ou launcher ja definiu explicitamente; respeita

    preference = _load_vblank_mode_preference()
    if preference == "off":
        os.environ["vblank_mode"] = "0"
        return
    if preference == "on":
        os.environ["vblank_mode"] = "1"
        return
    # preference == "auto" (default): mantem o auto-detect de GPU Intel abaixo.

    if sys.platform != "linux":
        return
    try:
        import glob
        card_pattern = re.compile(r"^card\d+$")
        for vendor_path in glob.glob("/sys/class/drm/card*/device/vendor"):
            card_name = vendor_path.split("/")[-3]
            if not card_pattern.match(card_name):
                continue
            with open(vendor_path) as f:
                vendor_id = f.read().strip()
            if vendor_id == "0x8086":  # Intel
                os.environ["vblank_mode"] = "0"
                return
    except OSError:
        pass  # falha na deteccao nao deve impedir o programa de iniciar

_maybe_disable_vsync_for_intel_igpu()

'''
O módulo re é nativo do Python (faz parte da biblioteca padrão).  
Você não precisa instalar nada.Ele implementa expressões regulares 
(regular expressions) para busca e manipulação de texto.
'''


'''
os.environ["PDYNAMO3_HOME"]         = "/home/fernando/programs/pDynamo3"
os.environ["PDYNAMO3_PARAMETERS"]   = "/home/fernando/programs/pDynamo3/parameters"
os.environ["PDYNAMO3_SCRATCH"]      = "/home/fernando/programs/pDynamo3/scratch"
os.environ["PDYNAMO3_ORCACOMMAND"]  = "/home/fernando/programs/orca_6_1_0_linux_x86-64_shared_openmpi418/orca"
os.environ["PDYNAMO3_XTBCOMMAND"]   = "/home/fernando/programs/xtb-6.6.1/bin/xtb"
os.environ["PDYNAMO3_MOPACCOMMAND"] = "/home/fernando/programs/MOPAC2016A/bin/mopac"
os.environ["PDYNAMO3_DFTBCOMMAND"]  = "/home/fernando/programs/dftbplus-24.1.x86_64-linux/bin/dftb+"
sys.path.append("/home/fernando/programs/pDynamo3")
#sys.path.append("/home/fernando/programs/amber24_src/lib/python3.12/site-packages")
'''

#filepath = '/home/fernando/programs/pDynamo3/installation/shellScripts/environment_bash.com'

def parse_bash_env_file(filepath):
    """
    Parse a bash environment file and return a dictionary
    with environment variables and values.
    """
    env_vars = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            # Ignore empty lines and comments
            if not line or line.startswith("#"):
                continue

            # Remove inline comments
            line = line.split("#")[0].strip()

            # Remove "; export VAR"
            line = re.sub(r";\s*export\s+\w+", "", line)

            if "=" in line:
                var, value = line.split("=", 1)
                var = var.strip()
                value = value.strip()

                env_vars[var] = value

    # Expand variables
    for var, value in env_vars.items():
        # Expand variables like $HOME, $PDYNAMO3_HOME etc
        value = os.path.expandvars(value)
        # Expand ~
        value = os.path.expanduser(value)
        env_vars[var] = value
        # Add to environment
        os.environ[var] = value

    # Add pDynamo to python path
    if "PDYNAMO3_HOME" in env_vars:
        sys.path.append(env_vars["PDYNAMO3_HOME"])

    return env_vars


try:
    from paths import PDYNAMO_HOME
    shell_scripts = os.path.join(PDYNAMO_HOME, 'installation/shellScripts/environment_bash.com')
    parse_bash_env_file(shell_scripts)
    logging.info('Importing environment variables from: %s', shell_scripts)
except Exception as e:
    logging.warning('paths.py / PDYNAMO_HOME not found -- pDynamo3 environment '
                     'variables were not imported (%s).', e)





import gi
# [EN] macOS/Quartz fix: PyGObject's automatic Gdk.init_check() (run as a
# side effect of the first "from gi.repository import Gtk/Gdk..." below)
# used to run before the Quartz backend was fully loaded, causing issues
# specific to that backend. gi.disable_legacy_autoinit() skips that
# automatic init; Gtk.init([]) right after the import does it explicitly,
# once the backend is actually ready. Gated to macOS only -- Linux/
# Windows keep using PyGObject's normal automatic init, unchanged.
if sys.platform == "darwin":
    gi.disable_legacy_autoinit()
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, GdkPixbuf, GLib
import cairo
import ctypes
import ctypes.util
if sys.platform == "darwin":
    Gtk.init([])


#               Installation is not necessary anymore.
#.This retrieves the absolute path of the script file that is currently being executed.
easy_main_file = os.path.abspath(__file__)

#.This extracts the directory (folder) from the absolute path obtained in the previous step.
EASYHYBRID_HOME = os.path.dirname(easy_main_file)

#.Adding GRAPHIC ENGINE LIB
sys.path.append(os.path.join(EASYHYBRID_HOME,"src/graphics_engine/src"))
sys.path.append(os.path.join(EASYHYBRID_HOME,"src/"))


from gui.main import MainWindow
from gui.eSession import EasyHybridSession
from gui.config   import VismolConfig
import time
import threading









# --- Bundled brand font for the splash screen -------------------------
# cairo's "select_font_face(name, ...)" is a *toy* API: it just hands
# "name" to the platform's font backend (fontconfig on Linux/macOS) and
# asks for whatever matches. "Sans" therefore renders as a DIFFERENT
# real font depending on what happens to be installed on the machine
# (typically DejaVu Sans on Linux, Helvetica/Arial elsewhere) -- not
# something we control or can guarantee looks the same everywhere.
#
# To get a font we actually ship (and that renders identically on any
# machine), we register the .ttf files bundled with vismol's own
# renderer (src/graphics_engine/.../libgl/fonts/, already used for the
# 3D viewport's atom labels) with fontconfig's *in-process* font list
# via FcConfigAppFontAddFile. This does NOT install anything on the
# system -- it only makes the font available to this process, for as
# long as it runs. Once registered, "select_font_face("Amiko", ...)"
# resolves to our own bundled file instead of whatever "Sans" happens
# to mean locally.
_BRAND_FONT_FAMILY   = "Amiko"
_BRAND_FONT_REGISTERED = False

def _register_bundled_fonts():
    """ Registers the bundled Amiko .ttf files with fontconfig for this
        process only. Safe to call more than once. On any failure (no
        fontconfig, unexpected platform, missing files, ...) it just
        logs a warning and leaves things as they were -- callers should
        always be prepared to fall back to a generic family name
        ("Sans") if this doesn't succeed.
    """
    global _BRAND_FONT_REGISTERED
    if _BRAND_FONT_REGISTERED:
        return True
    try:
        libname = ctypes.util.find_library("fontconfig") or "libfontconfig.so.1"
        fontconfig = ctypes.CDLL(libname)
        fontconfig.FcConfigAppFontAddFile.restype  = ctypes.c_int
        fontconfig.FcConfigAppFontAddFile.argtypes = [ctypes.c_void_p, ctypes.c_char_p]

        fonts_dir = os.path.join(EASYHYBRID_HOME,
                                  "src", "graphics_engine", "src",
                                  "vismol", "libgl", "fonts")
        ok = True
        for filename in ("Amiko-Regular.ttf", "Amiko-Bold.ttf"):
            path = os.path.join(fonts_dir, filename)
            ok = bool(fontconfig.FcConfigAppFontAddFile(None, path.encode("utf-8"))) and ok
        if not ok:
            raise RuntimeError("FcConfigAppFontAddFile reported failure for one or more files")
        _BRAND_FONT_REGISTERED = True
    except Exception as e:
        logging.warning('Could not register the bundled Amiko font (%s) -- '
                         'the splash screen will fall back to the system '
                         '"Sans" font instead.', e)
    return _BRAND_FONT_REGISTERED


# Splash Screen
class SplashScreen(Gtk.Window):
    """ Splash screen shown while modules load (see load_modules()/main()
        below). The background artwork lives in splash.png; the title,
        version and tagline text are drawn on top at runtime with Cairo
        instead of being baked into the image's own pixels (as the
        previous splash.png did) -- so the version string always matches
        EASYHYBRID_VERSION (_version.py) with no need to regenerate the
        image on every release. Based on tests/splash.py, integrated
        here with: the version pulled from EASYHYBRID_VERSION instead of
        a second hardcoded copy, the image decoded once and cached
        instead of on every "draw" event, the window sized from the
        image's own real dimensions instead of a hardcoded guess, kept
        centered (set_position), and the same defensive fallback the
        previous version had if splash.png is missing/unreadable (log a
        warning, show an empty window, instead of crashing on the very
        first frame of startup).
    """
    _BOX_HEIGHT = 92          # height, in pixels, of the semi-transparent text band
    _TEXT_X     = 24          # left margin for all three lines of text

    def __init__(self):
        super().__init__(title="EasyHybrid")
        self.set_decorated(False)  # Sem bordas
        self.set_position(Gtk.WindowPosition.CENTER)

        self._font_family = (_BRAND_FONT_FAMILY
                              if _register_bundled_fonts() else "Sans")

        self._image_surface = None
        width, height = 600, 492  # fallback size if splash.png can't be read at all
        try:
            self._image_surface = cairo.ImageSurface.create_from_png(
                os.path.join(EASYHYBRID_HOME, "splash.png"))
            width  = self._image_surface.get_width()
            height = self._image_surface.get_height()
        except Exception as e:
            logging.warning('splash.png could not be loaded (%s) -- '
                             'showing an empty splash window instead.', e)

        self.set_default_size(width, height)

        area = Gtk.DrawingArea()
        area.connect("draw", self._on_draw)
        self.add(area)

    def _on_draw(self, widget, cr):
        width  = widget.get_allocated_width()
        height = widget.get_allocated_height()

        if self._image_surface is not None:
            cr.set_source_surface(self._image_surface, 0, 0)
            cr.paint()

        # . Semi-transparent band along the bottom edge, holding the
        # text -- readable regardless of what's directly behind it in
        # the artwork itself.
        box_top = height - self._BOX_HEIGHT
        cr.set_source_rgba(0.0, 0.0, 0.0, 0.65)
        cr.rectangle(0, box_top, width, self._BOX_HEIGHT)
        cr.fill()

        cr.set_source_rgb(1.0, 1.0, 1.0)

        cr.select_font_face(self._font_family, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        cr.set_font_size(24)
        cr.move_to(self._TEXT_X, box_top + 30)
        cr.show_text("EasyHybrid")

        cr.select_font_face(self._font_family, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        cr.set_font_size(13)
        cr.move_to(self._TEXT_X, box_top + 52)
        cr.show_text("Version {}".format(EASYHYBRID_VERSION))

        cr.set_font_size(14)
        cr.move_to(self._TEXT_X, box_top + 72)
        cr.show_text("A Graphical Environment for QC/MM Simulations")

        return False

def load_modules(callback_final):
    def _load():
        logging.info("Starting module loading...")
        time.sleep(1.5)   
        GLib.idle_add(callback_final)  #Call the finalize function in the main loop
    threading.Thread(target=_load).start()


def main():
    splash = SplashScreen()
    splash.show_all()
    
    def on_finalizado():
        logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
                            datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)

        vconfig = VismolConfig(home = EASYHYBRID_HOME)
        # [REVERTIDO] Uma chamada explicita a vconfig.load_easyhybrid_config(...)
        # foi adicionada aqui numa correcao anterior, mas era redundante: o
        # proprio __init__ do VismolConfig ja chama _check_config_file()
        # (que carrega o .config.json salvo) seguido de _check_startup_path()
        # (que valida/corrige o startup_path). Chamar load_easyhybrid_config
        # de novo DEPOIS disso fazia merge com o conteudo BRUTO do arquivo
        # salvo, desfazendo a correcao que _check_startup_path() ja tinha
        # aplicado -- causa direta de um AttributeError/TypeError com
        # startup_path=None mais adiante, em EasyHybridPreferencesWindow.
        
        
        vm_session = EasyHybridSession(vm_config = vconfig)
        vm_session.vm_widget.insert_glmenu()
        main_window = MainWindow(vm_session = vm_session,
                                 home       =  EASYHYBRID_HOME,
                                 version    = EASYHYBRID_VERSION)
        vm_session.main_session = main_window                  
        #print(vm_session.vm_config.gl_parameters)
        #main_window.window.connect('destroy', Gtk.main_quit)
        
        
        # do now show these itens:
        
        #main_window.builder.get_object('toolbutton_export_img').hide()
        #main_window.builder.get_object('button_task_list')           .hide()
        
        main_window.builder.get_object('test_item')       .hide() # IR spectrum
        
        #This is the editor
        #This is the editor
        
        
        main_window.builder.get_object('_show_cell')      .hide() # IR spectrum
        
        
        
        
        #main_window.builder.get_object('toolbutton_terminal')   .hide()
        #main_window.builder.get_object('menuitem_reimaging')   .hide()
        #main_window.builder.get_object('menuitem_RMSD_tool')   .hide()
        
        
        main_window.builder.get_object('menuitem_rama')   .hide()
        main_window.builder.get_object('menuitem_transition_state_search').hide() # Baker method
        main_window.builder.get_object('menuitem_reaction_path')   .hide() # IRC / Reaction Path...
        main_window.builder.get_object('menuitem_cpr')   .hide() # Conjugate Peak Refinement...
        main_window.builder.get_object('menuitem_extras') .hide() # IR spectrum
        main_window.builder.get_object('menuitem_rdf_analysis')   .hide() # RDF Analysis (g(r))
        
        
        
        #main_window.builder.get_object('menuitem_advanced_rc_scans').hide()
        
        splash.destroy()
        # . BUG FIX: this used to unconditionally read sys.argv[-1] (the
        #   script's own path when no file was actually passed) and call
        #   vm_session.load_molecule() -- vismol's GENERIC loader, which
        #   only recognises .aux/.gro/.mol2/.pdb/.psf/.top/.prmtop/.xyz
        #   by extension and does nothing for a ".easy" EasyHybrid
        #   session file (parse_file() falls through every branch,
        #   leaving its `vismol_object` local unset -> UnboundLocalError
        #   -- confirmed against a real .easy that reproduced exactly
        #   "app opens, main treeview stays empty, no error shown").
        #   The bare `except: pass` was load-bearing for the "no file
        #   given" case (argv[-1] is always something), which is why it
        #   swallowed real load failures too. Fixed on both counts: only
        #   attempt a load when a file was actually passed, and use
        #   vm_session.load() -- EasyHybridSession's own override that
        #   actually dispatches ".easy" correctly (see gui/eSession.py) --
        #   instead of vm_session.load_molecule(); a genuine failure is
        #   now printed instead of silently disappearing.
        if len(sys.argv) > 1:
            filein = sys.argv[-1]
            try:
                vm_session.load(filein)
            except Exception:
                traceback.print_exc()
        #Gtk.main()
        return 0

    load_modules(on_finalizado)

    Gtk.main()


if __name__ == "__main__":
    main()

