#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  easyhybrid_terminal.py
#
#  Copyright 2022-2025 Fernando Bachega <ferbachega@gmail.com>
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
# ============================================================================
#  GUIA RAPIDO DO TERMINAL  (sintaxe DSL:  comando arg=valor arg2=valor2)
# ----------------------------------------------------------------------------
#  Descoberta
#     help                         lista todos os comandos
#     list                         lista os objetos carregados (com indices)
#
#  Representacoes  (rep = lines|sticks|spheres|dash|...)
#     show rep=sticks              aplica na SELECAO ATIVA
#     show rep=sticks obj=1        mira so o objeto 1 (nao muda a selecao)
#     show rep=spheres obj=0 chain=A
#     show rep=sticks  obj=0 resn=HIS name=CA
#     hide rep=lines   obj=0 resi=10-25
#
#  Selecao ativa (persistente)  -- depois, comandos sem obj= agem sobre ela
#     select obj=0                 objeto inteiro
#     select obj=0 chain=A,B       cadeias A e B
#     select obj=0 resi=10-30      faixa de residuos
#     select obj=0 chain=A resn=HIS name=CA
#     deselect                     limpa a selecao
#
#  Filtros (combinaveis, logica AND), aceitos por show/hide/select:
#     chain = A   | A,B                 (nome da cadeia)
#     resi  = 45  | 10-20 | 1-5,8,12    (indice de residuo; faixas e listas)
#     resn  = HIS                       (nome do residuo)
#     name  = CA                        (nome do atomo)
#
#  Trajetoria (ex.: 2000 frames)
#     frame                        mostra o frame atual
#     frame n=1000                 pula para o frame 1000
#     next  /  prev                avanca / retrocede um frame
#
#  Camera
#     center                       centraliza na selecao ativa
#     center obj=0                 centro de massa do objeto 0
#     center obj=0 chain=A resn=HIS
#     zoom dir=in                  aproxima (in) ou afasta (out)
#     zoom dir=out steps=10        varios passos de uma vez
#
#  Cena / arquivo
#     axes show=true | false
#     load file=/caminho/sistema.pdb
#     load file="/com espaco/sis.xyz"   (use aspas para caminhos com espaco)
#
#  Teclas: setas Up/Down navegam o historico; Tab completa o comando;
#          duplo-Tab lista as opcoes que casam com o prefixo.
#
#  API em scripts (mesma logica, sem a string DSL):
#     cmd.show(rep="sticks", obj=1, chain="A")
#     cmd.select(obj=0, resi="10-30")
#     for i in range(0, 2000, 100): cmd.frame(n=i)
# ============================================================================
import gi
import sys
import io
import time
import shlex
import inspect
import traceback
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, GLib

import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


# ============================================================================
#  Command — API unificada de comandos do EasyHybrid
# ----------------------------------------------------------------------------
#  Esta classe e a FONTE UNICA de comandos. Substitui as tres peças antes
#  fragmentadas (a Command stub do terminal, o CommandLine de eSession e o
#  eval/exec cru). O modelo de execucao e o DSL solicitado:
#
#        comando arg1=valor1 arg2=valor2
#
#  Dupla via de acesso (por isso "ambos" para API):
#    - No terminal (texto):     cmd.run("show rep=sticks")
#    - Em scripts Python:       cmd.show(rep="sticks")
#  As duas rotas caem no MESMO metodo cmd_*, entao nao ha duplicacao de logica.
#
#  Cada comando publico e um metodo nomeado cmd_<nome>. O prefixo cmd_ serve a
#  tres propositos: (1) o parser sabe o que e despachavel, (2) o autocomplete
#  enumera os comandos por introspecao (resolvendo o antigo command_list vazio),
#  (3) metodos auxiliares internos nao viram comandos por acidente.
# ============================================================================
class Command:
    def __init__(self, console, vm_session=None):
        self.console = console
        self.vm_session = vm_session

    # ----------------------------------------------------------------- parser
    def _parse(self, cmd_text):
        """
        Converte 'comando arg=val arg2=val2' em (nome, kwargs).

        Robustez em relacao ao parser antigo:
          - usa shlex para respeitar aspas ("nome com espaco");
          - NAO usa str.replace(func, '') (que corrompia args contendo o nome
            do comando);
          - tolera entrada vazia, espacos extras e '=' coladas ou separadas;
          - coage tipos: true/false -> bool, inteiros, floats, resto string.
        """
        cmd_text = (cmd_text or "").strip()
        if not cmd_text:
            return None, {}
        try:
            tokens = shlex.split(cmd_text)
        except ValueError:
            # aspas desbalanceadas etc.: cai para split simples
            tokens = cmd_text.split()
        name = tokens[0]
        kwargs = {}
        for tok in tokens[1:]:
            if "=" in tok:
                k, v = tok.split("=", 1)
                kwargs[k.strip()] = self._coerce(v.strip())
            else:
                # argumento posicional solto -> acumula em _args
                kwargs.setdefault("_args", []).append(self._coerce(tok))
        return name, kwargs

    @staticmethod
    def _coerce(value):
        """ Converte string do DSL para o tipo Python provavel. """
        low = value.lower()
        if low in ("true", "yes", "on"):
            return True
        if low in ("false", "no", "off"):
            return False
        if low in ("none", "null"):
            return None
        for cast in (int, float):
            try:
                return cast(value)
            except (ValueError, TypeError):
                pass
        return value

    # ------------------------------------------------------------- dispatcher
    def run(self, cmd_text):
        """
        Ponto de entrada do TERMINAL (texto DSL). Faz parse, localiza o metodo
        cmd_<nome> e o executa com os kwargs. Retorna uma string de log (ou
        None) que o terminal imprime. Erros viram log legivel, nunca derrubam
        a janela.
        """
        name, kwargs = self._parse(cmd_text)
        if name is None:
            return None
        method = getattr(self, "cmd_" + name, None)
        if method is None or not callable(method):
            return "Comando '{}' nao reconhecido. Use 'help'.".format(name)
        try:
            return method(**kwargs)
        except TypeError as te:
            # tipicamente argumento errado: mostra a assinatura para ajudar
            sig = inspect.signature(method)
            return "Uso: {} {}\n  ({})".format(name, self._sig_hint(sig), te)
        except Exception as exc:
            return "Erro em '{}': {}".format(name, exc)

    @staticmethod
    def _sig_hint(sig):
        parts = []
        for pname, p in sig.parameters.items():
            if pname == "_args":
                continue
            if p.default is inspect.Parameter.empty:
                parts.append(pname)
            else:
                parts.append("{}={!r}".format(pname, p.default))
        return " ".join(parts)

    # ---------------------------------------------------- introspeccao p/ TAB
    def command_names(self):
        """ Lista os nomes de comando disponiveis (sem o prefixo cmd_).
            Usado pelo autocomplete e pelo 'help'. """
        return sorted(
            attr[len("cmd_"):]
            for attr in dir(self)
            if attr.startswith("cmd_") and callable(getattr(self, attr))
        )

    # ---------------------------------------------------- helpers de alvo
    def _resolve_object(self, obj):
        """ Recebe um indice (int) e devolve o vm_object correspondente, ou
            None se nao existir. Aceita tambem None (sem alvo). """
        if obj is None or self.vm_session is None:
            return None
        try:
            return self.vm_session.vm_objects_dic.get(int(obj))
        except (ValueError, TypeError):
            return None

    def _selection_for_atoms(self, vm_object, atoms):
        """ Monta uma selecao temporaria (VMSele) contendo `vm_object` e o
            conjunto `atoms` dado, SEM tocar na selecao ativa do usuario.

            Usada por show/hide quando recebem filtros inline (obj=, chain=,
            ...), para que a acao seja pontual e nao mude o estado de selecao.
        """
        active = self.vm_session.selections[self.vm_session.current_selection]
        temp = active.__class__(self.vm_session)          # mesmo tipo VMSele
        temp.selected_objects = {vm_object}
        temp.selected_atoms = set(atoms)
        temp.selected_atom_ids = set(
            getattr(a, "atom_id", getattr(a, "index", None)) for a in atoms)
        return temp

    @staticmethod
    def _match_field(value, spec):
        """ Verdadeiro se `value` casa com `spec`, onde spec pode ser:
              - valor unico:        "A"        -> value == "A"
              - lista por virgula:  "A,B,C"    -> value em {A,B,C}
              - faixa numerica:     "10-20"    -> 10 <= int(value) <= 20
              - lista de faixas:    "1-5,8,12" -> combina os anteriores
            Comparacao textual e case-sensitive (nomes de cadeia/residuo
            costumam ser); faixas exigem que value seja inteiro.

            Exemplos:
              _match_field("A",  "A,B")   -> True
              _match_field("15", "10-20") -> True
              _match_field("CA", "CA")    -> True
        """
        value = str(value)
        for part in str(spec).split(","):
            part = part.strip()
            if "-" in part and part.replace("-", "").isdigit():
                lo, hi = part.split("-", 1)
                try:
                    if int(lo) <= int(value) <= int(hi):
                        return True
                except ValueError:
                    pass
            elif part == value:
                return True
        return False

    def _filter_atoms(self, vm_object, chain=None, resi=None, resn=None,
                      name=None):
        """ Devolve a lista de atomos de `vm_object` que satisfazem TODOS os
            filtros fornecidos (logica AND). Filtro None = ignorado.

            Os atributos comparados existem no modelo de atomo do VisMol:
              chain -> atom.chain.name     resn -> atom.residue.name
              resi  -> atom.residue.index  name -> atom.name

            NOTA: os metodos selecting_by_* da sessao NAO filtram por nome --
            eles expandem a partir de um atomo clicado. Por isso filtramos
            aqui diretamente.

            Exemplos de specs aceitas (via _match_field):
              chain="A"        chain="A,B"
              resi="45"        resi="10-20"      resi="1-5,8,12"
              resn="HIS"       name="CA"
        """
        def keep(atom):
            if chain is not None:
                cn = getattr(getattr(atom, "chain", None), "name", None)
                if cn is None or not self._match_field(cn, chain):
                    return False
            if resn is not None:
                rn = getattr(getattr(atom, "residue", None), "name", None)
                if rn is None or not self._match_field(rn, resn):
                    return False
            if resi is not None:
                ri = getattr(getattr(atom, "residue", None), "index", None)
                if ri is None or not self._match_field(ri, resi):
                    return False
            if name is not None:
                an = getattr(atom, "name", None)
                if an is None or not self._match_field(an, name):
                    return False
            return True
        return [a for a in vm_object.atoms.values() if keep(a)]

    # ========================================================================
    #  COMANDOS REAIS  (cada um e tambem chamavel direto via API Python)
    #  Ancorados na API publica de vm_session, nao em stubs.
    # ========================================================================
    def cmd_help(self, **_):
        """ Lista os comandos disponiveis e a sintaxe de cada um. """
        lines = ["Comandos disponiveis:"]
        for name in self.command_names():
            method = getattr(self, "cmd_" + name)
            doc = (method.__doc__ or "").strip().split("\n")[0]
            lines.append("  {:<12} {}".format(name, doc))
        return "\n".join(lines)

    def cmd_list(self, **_):
        """ Lista os objetos moleculares carregados na sessao. """
        if self.vm_session is None:
            return "Sessao indisponivel."
        objs = self.vm_session.vm_objects_dic
        if not objs:
            return "Nenhum objeto carregado."
        return "\n".join("  [{}] {}".format(i, getattr(o, "name", "?"))
                         for i, o in objs.items())

    def cmd_show(self, rep="lines", obj=None, chain=None, resi=None,
                 resn=None, name=None, **_):
        """ Mostra uma representacao.

            Sem filtros, age sobre a SELECAO ATIVA. Com obj= (e opcionalmente
            chain/resi/resn/name) mira so esses atomos, SEM alterar a selecao
            ativa -- a acao e pontual.

            rep  : lines | sticks | spheres | dash | ...
            obj  : indice do objeto (veja 'list')
            chain: nome da cadeia      ex: A   ou A,B
            resi : indice de residuo   ex: 45  ou 10-20  ou 1-5,8
            resn : nome do residuo     ex: HIS
            name : nome do atomo       ex: CA

            Exemplos:
              show rep=sticks
              show rep=sticks obj=1
              show rep=spheres obj=0 chain=A
              show rep=sticks  obj=0 resn=HIS name=CA
              show rep=lines   obj=0 resi=10-25
        """
        return self._show_or_hide(True, rep, obj, chain, resi, resn, name)

    def cmd_hide(self, rep="lines", obj=None, chain=None, resi=None,
                 resn=None, name=None, **_):
        """ Oculta uma representacao. Mesmos filtros de 'show'.

            Exemplos:
              hide rep=lines
              hide rep=lines   obj=0
              hide rep=spheres obj=0 chain=B
              hide rep=sticks  obj=0 resn=HOH
        """
        return self._show_or_hide(False, rep, obj, chain, resi, resn, name)

    def _show_or_hide(self, show, rep, obj, chain, resi, resn, name):
        """ Implementacao compartilhada de show/hide. Decide o alvo:
              - nenhum filtro          -> selecao ativa (selection=None);
              - obj= [+ filtros finos] -> selecao temporaria so com o alvo.
        """
        if self.vm_session is None:
            return "Sessao indisponivel."
        verbo = "Exibindo" if show else "Ocultando"

        # Caso 1: sem alvo explicito -> usa a selecao ativa do usuario.
        if obj is None:
            self.vm_session.show_or_hide(rep_type=rep, selection=None, show=show)
            return "{} '{}' na selecao ativa.".format(verbo, rep)

        # Caso 2: alvo explicito por objeto (e talvez filtros finos).
        target = self._resolve_object(obj)
        if target is None:
            return "Objeto {} nao encontrado. Use 'list'.".format(obj)

        has_fine = any(f is not None for f in (chain, resi, resn, name))
        if has_fine:
            atoms = self._filter_atoms(target, chain, resi, resn, name)
            if not atoms:
                return "Nenhum atomo casa com os filtros em obj={}.".format(obj)
            selection = self._selection_for_atoms(target, atoms)
            escopo = "obj={} ({} atomos filtrados)".format(obj, len(atoms))
        else:
            selection = self._selection_for_atoms(target, target.atoms.values())
            escopo = "obj={} (inteiro)".format(obj)

        self.vm_session.show_or_hide(rep_type=rep, selection=selection, show=show)
        return "{} '{}' em {}.".format(verbo, rep, escopo)

    def cmd_select(self, obj=None, chain=None, resi=None, resn=None,
                   name=None, **_):
        """ Define a SELECAO ATIVA (persistente). Comandos seguintes sem
            obj= passam a agir sobre ela ate um novo 'select' ou 'deselect'.

            Filtros combinaveis (AND), iguais aos de show/hide:
              obj  : indice do objeto (obrigatorio)
              chain: A   ou A,B
              resi : 45  ou 10-20  ou 1-5,8
              resn : HIS
              name : CA

            Exemplos:
              select obj=0
              select obj=0 chain=A
              select obj=0 chain=A resn=HIS
              select obj=0 resi=10-30
              select obj=0 chain=A resn=HIS name=CA
        """
        if self.vm_session is None:
            return "Sessao indisponivel."
        target = self._resolve_object(obj)
        if target is None:
            return "Uso: select obj=N [chain= resi= resn= name=]  ('list')"

        atoms = self._filter_atoms(target, chain, resi, resn, name)
        if not atoms:
            return "Nenhum atomo casa com os filtros em obj={}.".format(obj)

        # Grava na selecao ativa (persistente).
        active = self.vm_session.selections[self.vm_session.current_selection]
        active.selected_objects = {target}
        active.selected_atoms = set(atoms)
        active.selected_atom_ids = set(
            getattr(a, "atom_id", getattr(a, "index", None)) for a in atoms)
        # Sincroniza o flag .selected de cada atomo do objeto.
        for a in target.atoms.values():
            a.selected = False
        for a in atoms:
            a.selected = True

        filtros = ", ".join(
            "{}={}".format(k, v) for k, v in
            (("chain", chain), ("resi", resi), ("resn", resn), ("name", name))
            if v is not None) or "objeto inteiro"
        return "Selecionado obj={} ({}): {} atomos.".format(obj, filtros, len(atoms))

    def cmd_deselect(self, **_):
        """ Limpa a selecao ativa. """
        if self.vm_session is None:
            return "Sessao indisponivel."
        active = self.vm_session.selections[self.vm_session.current_selection]
        active.selected_objects = set()
        active.selected_atoms = set()
        active.selected_atom_ids = set()
        return "Selecao limpa."

    def cmd_frame(self, n=None, **_):
        """ Sem n: mostra o frame atual. Com n=int: pula para esse frame. """
        if self.vm_session is None:
            return "Sessao indisponivel."
        if n is None:
            return "Frame atual: {}".format(self.vm_session.get_frame())
        self.vm_session.set_frame(frame=int(n))
        return "Frame -> {}".format(int(n))

    def cmd_next(self, **_):
        """ Avanca um frame da trajetoria. """
        if self.vm_session is None:
            return "Sessao indisponivel."
        self.vm_session.forward_frame()
        return "Frame -> {}".format(self.vm_session.get_frame())

    def cmd_prev(self, **_):
        """ Retrocede um frame da trajetoria. """
        if self.vm_session is None:
            return "Sessao indisponivel."
        self.vm_session.reverse_frame()
        return "Frame -> {}".format(self.vm_session.get_frame())

    def cmd_axes(self, show=True, **_):
        """ Mostra/oculta os eixos (show=true|false). """
        if self.vm_session is None:
            return "Sessao indisponivel."
        if show:
            self.vm_session.show_axes()
            return "Eixos visiveis."
        self.vm_session.hide_axes()
        return "Eixos ocultos."

    def cmd_center(self, obj=None, chain=None, resi=None, resn=None,
                   name=None, **_):
        """ Centraliza a camera num alvo (anima a translacao ate ele).

            Sem obj=: centraliza na SELECAO ATIVA (centroide dos atomos).
            Com obj= e sem filtros finos: centro de massa do objeto.
            Com obj= e filtros: centroide dos atomos filtrados.

            obj=N  chain=A  resi=45|10-20  resn=HIS  name=CA

            Exemplos:
              center                       (na selecao ativa)
              center obj=0                 (centro de massa do objeto 0)
              center obj=0 chain=A         (centroide da cadeia A)
              center obj=0 resn=HIS name=CA
        """
        if self.vm_session is None:
            return "Sessao indisponivel."
        glcore = self.vm_session.vm_glcore

        # Caso A: sem objeto -> centraliza na selecao ativa.
        if obj is None:
            active = self.vm_session.selections[self.vm_session.current_selection]
            atoms = list(active.selected_atoms)
            if not atoms:
                return "Nada selecionado. Use 'select' ou passe obj=N."
            vm_object = next(iter(active.selected_objects), None)
            target = self._centroid(atoms, vm_object)
            glcore.center_on_coordinates(vm_object, target)
            return "Centralizado na selecao ativa ({} atomos).".format(len(atoms))

        # Caso B: objeto explicito.
        vm_object = self._resolve_object(obj)
        if vm_object is None:
            return "Objeto {} nao encontrado. Use 'list'.".format(obj)

        has_fine = any(f is not None for f in (chain, resi, resn, name))
        if has_fine:
            atoms = self._filter_atoms(vm_object, chain, resi, resn, name)
            if not atoms:
                return "Nenhum atomo casa com os filtros em obj={}.".format(obj)
            target = self._centroid(atoms, vm_object)
            glcore.center_on_coordinates(vm_object, target)
            return "Centralizado em obj={} ({} atomos).".format(obj, len(atoms))
        else:
            # objeto inteiro -> centro de massa (mesmo padrao do autocenter)
            glcore.center_on_coordinates(vm_object, vm_object.mass_center)
            return "Centralizado no centro de massa de obj={}.".format(obj)

    def _centroid(self, atoms, vm_object):
        """ Centroide (media das coordenadas) de uma lista de atomos no frame
            atual. Usa atom.coords(frame) quando disponivel; recai em
            atom.coords se for atributo. Retorna np.array float32 [x,y,z]. """
        import numpy as np
        frame = 0
        try:
            frame = self.vm_session.get_frame()
        except Exception:
            pass
        pts = []
        for a in atoms:
            c = getattr(a, "coords", None)
            if callable(c):
                pts.append(np.asarray(c(frame), dtype=np.float32))
            elif c is not None:
                pts.append(np.asarray(c, dtype=np.float32))
        if not pts:
            # fallback: centro de massa do objeto
            return vm_object.mass_center
        return np.mean(np.vstack(pts), axis=0).astype(np.float32)

    def cmd_zoom(self, dir="in", steps=5, **_):
        """ Aproxima ou afasta a camera (equivale a girar o scroll do mouse).

            dir  : in | out          (in aproxima, out afasta)
            steps: quantos 'cliques' de scroll aplicar (default 5)

            Exemplos:
              zoom dir=in
              zoom dir=out steps=10
        """
        if self.vm_session is None:
            return "Sessao indisponivel."
        glcore = self.vm_session.vm_glcore
        direction = 1 if str(dir).lower() in ("in", "+", "1", "up") else -1
        try:
            n = max(1, int(steps))
        except (ValueError, TypeError):
            n = 5
        for _i in range(n):
            glcore.mouse_scroll(direction)
        return "Zoom {} x{}.".format("in" if direction == 1 else "out", n)

    def cmd_load(self, file=None, **_):
        """ Carrega uma molecula de arquivo (file=/caminho/arquivo). """
        if self.vm_session is None:
            return "Sessao indisponivel."
        if file is None:
            return "Uso: load file=/caminho/do/arquivo"
        self.vm_session.load_molecule(file)
        return "Carregado: {}".format(file)

    # --- atalhos de API direta (mesma logica, nomes "pythonicos") -----------
    #  Permitem cmd.show(rep="sticks") em scripts externos sem o prefixo cmd_.
    # --- atalhos de API direta (mesma logica, nomes "pythonicos") -----------
    #  Permitem uso em scripts externos sem o prefixo cmd_ nem a string DSL:
    #      cmd.show(rep="sticks", obj=1, chain="A")
    #      cmd.select(obj=0, resi="10-30")
    def show(self, rep="lines", obj=None, chain=None, resi=None, resn=None, name=None):
        return self.cmd_show(rep=rep, obj=obj, chain=chain, resi=resi, resn=resn, name=name)
    def hide(self, rep="lines", obj=None, chain=None, resi=None, resn=None, name=None):
        return self.cmd_hide(rep=rep, obj=obj, chain=chain, resi=resi, resn=resn, name=name)
    def select(self, obj=None, chain=None, resi=None, resn=None, name=None):
        return self.cmd_select(obj=obj, chain=chain, resi=resi, resn=resn, name=name)
    def center(self, obj=None, chain=None, resi=None, resn=None, name=None):
        return self.cmd_center(obj=obj, chain=chain, resi=resi, resn=resn, name=name)
    def zoom(self, dir="in", steps=5):      return self.cmd_zoom(dir=dir, steps=steps)
    def frame(self, n=None):                return self.cmd_frame(n=n)
    def list(self):                         return self.cmd_list()


class TerminalWindow():
    """ Janela do terminal do EasyHybrid. """

    def open_window(self):
        """ """
        if self.visible == False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.main.home, 'src/gui/windows/setup/easyhybrid_terminal.glade'))
            self.builder.connect_signals(self)

            self.window = self.builder.get_object('window')
            self.window.set_default_size(650, 350)
            self.window.set_title('EasyHybrid Terminal')
            self.window.connect('destroy-event', self.close_window)
            self.window.set_keep_above(True)

            self.tag_table = self.textbuffer.get_tag_table()
            self.textbuffer.create_tag("green", foreground="green")
            self.textbuffer.create_tag("red", foreground="red")
            self.textbuffer.create_tag("blue", foreground="blue")

            self.locals = {}
            # Contexto de execucao: a Command unificada
            self.cmd = Command(self, self.vm_session)
            self.locals['cmd'] = self.cmd
            # command_list agora e preenchido de verdade (resolvia o vazio)
            self.command_list = self.cmd.command_names()

            self.entry_terminal = self.builder.get_object('entry_terminal')
            self.textview = self.builder.get_object('entry_text_buffer')
            self.textview.set_buffer(self.textbuffer)
            self.textview.set_wrap_mode(Gtk.WrapMode.WORD)
            self.textview.get_style_context().add_provider(self.css_provider, Gtk.STYLE_PROVIDER_PRIORITY_USER)
            self.entry_terminal.get_style_context().add_provider(self.css_provider, Gtk.STYLE_PROVIDER_PRIORITY_USER)

            self.window.show_all()
            self.visible = True
        else:
            self.window.present()

    def close_window(self, button, data=None):
        """ Function doc """
        self.window.destroy()
        self.visible = False
        self.main.cmd_terminal_button.set_active(False)

    def __init__(self, main=None):
        """ Class initialiser """
        self.main = main
        self.visible = False
        self.p_session = main.p_session
        self.vm_session = main.vm_session

        self.cmd_history = []
        self.cmd_history_counter = 0
        self.textbuffer = self.main.terminal_text_buffer
        self.command_list = []

        text = '''
        --------------------------------------------------------
               Welcome to the EasyHybrid Terminal.
        --------------------------------------------------------

                   Created by J.F.R Bachega


        '''
        end_iter = self.textbuffer.get_end_iter()
        self.textbuffer.insert(end_iter, text)

        self.last_tab_time = 0
        self.tab_timer_id = None

        self.css_provider = Gtk.CssProvider()
        self.css_provider.load_from_data(b"""
        textview {
            font-family: Monospace;
            font-size: 12pt;
        }

        entry {
            font-family: Monospace;
            font-size: 12pt;
        }
        """)

    def run_cmd(self, cmd):
        """ Executa um comando DSL via Command.run e imprime o log. """
        self.write_output(">" + cmd, "normal")
        log = self.cmd.run(cmd)
        if log is not None:
            self.write_output(log, "normal")

    def write_output(self, text: str, color: str = "normal"):
        end_iter = self.textbuffer.get_end_iter()
        tag = self.tag_table.lookup(color)
        if not tag:
            tag = self.textbuffer.create_tag(color, foreground=color)
        self.textbuffer.insert_with_tags(end_iter, text + "\n", tag)
        self.textview.scroll_to_iter(end_iter, 0.0, False, 0, 0)

    def on_entry_terminal(self, widget):
        """
        Enter no campo de entrada. Agora roteia pelo DSL (Command.run) em vez
        de eval/exec cru -- mais seguro e consistente com o autocomplete.
        """
        command = self.entry_terminal.get_text()
        self.entry_terminal.set_text("")
        if not command.strip():
            return
        self.cmd_history.append(command)
        self.cmd_history_counter = 0
        self.write_output(">>> " + command, "normal")
        try:
            log = self.cmd.run(command)
            if log is not None:
                self.write_output(str(log), "normal")
        except Exception:
            self.write_output(traceback.format_exc().strip(), "red")

    # ------------------------------------------------------------ navegacao
    def on_entry_terminal_backspace(self, widget):
        pass

    def on_entry_terminal_move_cursor(self, widget, data, data2, data3):
        pass

    def on_entry_terminal_change(self, widget):
        pass

    def update_window(self, system_names=True, coordinates=False, selections=True):
        pass

    def on_key_press_event(self, widget, event):
        """ Setas: navega no historico. Tab: autocomplete. Duplo-Tab: lista. """
        k_name = Gdk.keyval_name(event.keyval)
        size = -len(self.cmd_history)

        if k_name in ['Down', 'Up']:
            if k_name == 'Up':
                self.cmd_history_counter += -1
            if k_name == 'Down':
                self.cmd_history_counter += 1

            if self.cmd_history_counter >= 0:
                self.entry_terminal.set_text('')
                self.cmd_history_counter = 0
            elif self.cmd_history_counter < 0 and self.cmd_history_counter > size:
                self.entry_terminal.set_text(self.cmd_history[self.cmd_history_counter])
            else:
                self.cmd_history_counter = size
                self.entry_terminal.set_text(self.cmd_history[self.cmd_history_counter])
            # leva o cursor para o fim do texto recuperado
            self.entry_terminal.set_position(-1)
            return True

        if event.keyval == Gdk.KEY_Tab:
            now = time.time()
            if now - self.last_tab_time < 0.2:
                # duplo-Tab: cancela o timer e LISTA todas as opcoes
                if self.tab_timer_id is not None:
                    GLib.source_remove(self.tab_timer_id)
                    self.tab_timer_id = None
                self.last_tab_time = 0
                self._list_completions()
            else:
                # Tab simples (agendado): tenta completar o prefixo atual
                self.last_tab_time = now
                self.tab_timer_id = GLib.timeout_add(200, self._simple_Tab)
            # retornar True impede o GTK de mover o foco para fora do entry
            return True
        return False

    # ----------------------------------------------------------- autocomplete
    def _current_word(self):
        """ A primeira palavra digitada (o nome do comando em construcao). """
        return self.entry_terminal.get_text().strip().split(" ")[0]

    def _matches(self, prefix):
        return [c for c in self.command_list if c.startswith(prefix)]

    def _simple_Tab(self):
        """
        Tab simples: completa o nome do comando.
          - 1 match  -> completa e adiciona espaco;
          - varios   -> completa ate o maior prefixo comum;
          - nenhum   -> nada.
        """
        self.tab_timer_id = None
        full = self.entry_terminal.get_text()
        # so completa enquanto o usuario ainda digita o NOME (sem espaco)
        if " " in full.strip():
            return False
        prefix = full.strip()
        matches = self._matches(prefix)
        if not matches:
            return False
        if len(matches) == 1:
            self.entry_terminal.set_text(matches[0] + " ")
        else:
            common = os.path.commonprefix(matches)
            if len(common) > len(prefix):
                self.entry_terminal.set_text(common)
            else:
                self._list_completions()
        self.entry_terminal.set_position(-1)
        return False  # remove o timer

    def _list_completions(self):
        """ Mostra no buffer as opcoes que casam com o prefixo atual. """
        prefix = self._current_word()
        matches = self._matches(prefix) if prefix else list(self.command_list)
        if matches:
            self.write_output("  ".join(matches), "blue")

    def _simple_Tab_legacy(self):
        # mantido apenas como referencia; nao usado
        self.tab_timer_id = None
        return False
