#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: umbrella sampling runs read back (Umbrella Sampling Browser, 2026-10-23)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      The folder written by Simulation > Umbrella Sampling
#      (pdynamo/p_methods/umbrella_sampling.py):
#
#        <name>_umb_sam/
#          output.log                 parameters (a python dict after the system summary)
#          optmization.ptGeo/         frame{i}.pkl / frame{i}_{j}.pkl  initial structure of each window
#          equilibration/             window{i}[_{j}].ptGeo           equilibration frames
#          data_collection/           window{i}[_{j}].ptRes           restraint values, every step:
#                                         header.pkl (labels, restraintEnergyModels: target, force)
#                                         block*.pkl (RealArray1D, RC1 RC2 RC1 RC2 ... interleaved)
#          data_collection_traj/      window{i}[_{j}].ptGeo           data collection frames
#          wham/                      *_pmf.log (x [y] F), WHAM output
#
#      GTK-free: the window (gui/windows/analysis/umbrella_browser_window.py)
#      draws what this module computes. Window keys are tuples: ( i, ) for 1D
#      and ( i, j ) for 2D. RC values in the units of the restraint (A or degrees);
#      energies in kJ/mol (the restraint E = k ( x - x0 )^2 of pDynamo, k in
#      kJ/mol/unit^2).
#

import ast
import glob
import heapq
import math
import os
import pickle
import re

import numpy as np

KB = 0.0083144626          # kJ/mol/K
SOURCES = ( ( "dc",  "Data collection", "data_collection_traj" ),
            ( "eq",  "Equilibration",   "equilibration" ),
            ( "opt", "Initial structure (optimization)", "optmization.ptGeo" ) )


def _key_of ( name ):
    m = re.match ( r"(?:window|frame)(\d+)(?:_(\d+))?", name )
    if not m: return None
    return ( int ( m.group ( 1 ) ), ) if m.group ( 2 ) is None else ( int ( m.group ( 1 ) ), int ( m.group ( 2 ) ) )


def _numeric_frames ( folder ):
    files = [ f for f in os.listdir ( folder ) if f.startswith ( "frame" ) and f.endswith ( ".pkl" ) ]
    return [ os.path.join ( folder, f ) for f in sorted ( files, key = lambda f: int ( re.sub ( r"\D", "", f ) or 0 ) ) ]


def read_parameters ( path ):
    """ The parameters dict printed in output.log (objects replaced by None); {} if absent. """
    try:
        text = open ( path, errors = "replace" ).read ( )
    except OSError:
        return { }
    start = text.find ( "\n{'" )
    if start < 0: return { }
    depth, end = 0, None
    for k in range ( start + 1, len ( text ) ):
        c = text[k]
        if c == "{": depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0: end = k + 1; break
    if end is None: return { }
    body = re.sub ( r"<[^<>]*? object at 0x[0-9a-f]+>", "None", text[start + 1:end] )
    try:
        return ast.literal_eval ( body )
    except Exception:
        return { }


class USRun:
    """ One umbrella sampling folder. """

    def __init__ ( self, folder ):
        self.folder = os.path.abspath ( folder )
        if not os.path.isdir ( os.path.join ( self.folder, "data_collection" ) ):
            raise ValueError ( "{} has no data_collection folder: not an EasyHybrid umbrella sampling run.".format ( folder ) )
        self.params = read_parameters ( os.path.join ( self.folder, "output.log" ) )
        md = self.params.get ( "MD_parm" ) or { }
        self.temperature = float ( md.get ( "temperature" ) or 300.0 )
        self.freq_res = int ( md.get ( "trajectory_frequency_dc_ptRes" ) or 1 )
        self.freq_geo = int ( md.get ( "trajectory_frequency_dc_ptGeo" ) or md.get ( "trajectory_frequency" ) or 1 )
        self.windows = { }
        for path in glob.glob ( os.path.join ( self.folder, "data_collection", "window*.ptRes" ) ):
            key = _key_of ( os.path.basename ( path ) )
            if key is not None: self.windows.setdefault ( key, { } )["res"] = path
        for code, _label, sub in SOURCES:
            base = os.path.join ( self.folder, sub )
            if not os.path.isdir ( base ): continue
            for name in os.listdir ( base ):
                key = _key_of ( name )
                if key is None: continue
                path = os.path.join ( base, name )
                if code == "opt" and name.endswith ( ".pkl" ): self.windows.setdefault ( key, { } )["opt"] = path
                elif code != "opt" and name.endswith ( ".ptGeo" ): self.windows.setdefault ( key, { } )[code] = path
        if not self.windows:
            raise ValueError ( "No umbrella windows found in {}.".format ( folder ) )
        dims = { len ( k ) for k in self.windows }
        if len ( dims ) != 1: raise ValueError ( "Mixed 1D and 2D window names in {}.".format ( folder ) )
        self.dims = dims.pop ( )
        self.shape = tuple ( max ( k[d] for k in self.windows ) + 1 for d in range ( self.dims ) )
        self.keys = sorted ( self.windows )
        self.labels = [ "RC{}".format ( d + 1 ) for d in range ( self.dims ) ]
        self._rc, self._header = { }, { }

    # ------------------------------------------------------------------ restraint data
    def header ( self, key ):
        if key not in self._header:
            path = self.windows.get ( key, { } ).get ( "res" )
            h = None
            if path and os.path.exists ( os.path.join ( path, "header.pkl" ) ):
                with open ( os.path.join ( path, "header.pkl" ), "rb" ) as handle: h = pickle.load ( handle )
            self._header[key] = h
        return self._header[key]

    def target ( self, key ):
        """ ( x0 per dimension ), from the restraint models of the window; None if unknown. """
        h = self.header ( key )
        if not h: return None
        return tuple ( float ( m.lowEquilibriumValue ) for m in h["restraintEnergyModels"] )

    def force ( self, key ):
        h = self.header ( key )
        if not h: return None
        return tuple ( float ( m.lowForceConstant ) for m in h["restraintEnergyModels"] )

    def rc ( self, key ):
        """ ( n, dims ) restraint values, every ptRes step (cached). """
        if key not in self._rc:
            path = self.windows.get ( key, { } ).get ( "res" )
            values = [ ]
            if path:
                blocks = sorted ( glob.glob ( os.path.join ( path, "block*.pkl" ) ), key = lambda f: int ( re.sub ( r"\D", "", os.path.basename ( f ) ) or 0 ) )
                for b in blocks:
                    with open ( b, "rb" ) as handle: values.extend ( list ( pickle.load ( handle ) ) )
            nd = len ( ( self.header ( key ) or { } ).get ( "labels" ) or self.labels )
            a = np.array ( values, dtype = float )
            self._rc[key] = a[:len ( a ) // nd * nd].reshape ( -1, nd ) if len ( a ) else np.zeros ( ( 0, self.dims ) )
        return self._rc[key]

    def stats ( self ):
        """ { key: { n, mean, std, target, deviation (|mean - x0| / std of the restraint), drift } }. """
        out = { }
        for key in self.keys:
            x = self.rc ( key ); t = self.target ( key ); k = self.force ( key )
            if not len ( x ):
                out[key] = { "n": 0 }; continue
            mean = x.mean ( axis = 0 ); std = x.std ( axis = 0 )
            half = len ( x ) // 2
            drift = ( x[half:].mean ( axis = 0 ) - x[:half].mean ( axis = 0 ) ) if half else np.zeros ( self.dims )
            # . the width of the restraint alone: sigma = sqrt ( kT / 2k )
            sigma = np.array ( [ math.sqrt ( KB * self.temperature / ( 2.0 * f ) ) for f in k ] ) if k else std
            out[key] = { "n": int ( len ( x ) ), "mean": mean, "std": std, "target": np.array ( t ) if t else None, "drift": drift,
                         "deviation": float ( np.linalg.norm ( ( mean - np.array ( t ) ) / sigma ) ) if t else float ( "nan" ) }
        return out

    # ------------------------------------------------------------------ overlap
    def overlap_1d ( self, bins = 120 ):
        """ ( pairs [ ( key_a, key_b, overlap ) ] of neighbouring windows, histogram data ). Overlap
            coefficient sum min ( p, q ) of the normalized histograms on a common grid (0 to 1). """
        lo = min ( self.rc ( k )[:, 0].min ( ) for k in self.keys if len ( self.rc ( k ) ) )
        hi = max ( self.rc ( k )[:, 0].max ( ) for k in self.keys if len ( self.rc ( k ) ) )
        edges = np.linspace ( lo, hi, bins + 1 )
        hist = { k: np.histogram ( self.rc ( k )[:, 0], bins = edges )[0] / max ( 1, len ( self.rc ( k ) ) ) for k in self.keys }
        order = sorted ( self.keys, key = lambda k: ( self.target ( k ) or ( k[0], ) )[0] )
        pairs = [ ( a, b, float ( np.minimum ( hist[a], hist[b] ).sum ( ) ) ) for a, b in zip ( order, order[1:] ) ]
        return pairs, { "edges": edges, "hist": hist, "order": order }

    def overlap_2d ( self, bins = 60 ):
        """ Neighbour pairs along i and along j with their 2D overlap coefficient. """
        allx = np.vstack ( [ self.rc ( k ) for k in self.keys if len ( self.rc ( k ) ) ] )
        ex = np.linspace ( allx[:, 0].min ( ), allx[:, 0].max ( ), bins + 1 )
        ey = np.linspace ( allx[:, 1].min ( ), allx[:, 1].max ( ), bins + 1 )
        hist = { }
        for k in self.keys:
            x = self.rc ( k )
            hist[k] = np.histogram2d ( x[:, 0], x[:, 1], bins = [ ex, ey ] )[0] / max ( 1, len ( x ) )
        pairs = [ ]
        for ( i, j ) in self.keys:
            for nb in ( ( i + 1, j ), ( i, j + 1 ) ):
                if nb in hist:
                    pairs.append ( ( ( i, j ), nb, float ( np.minimum ( hist[( i, j )], hist[nb] ).sum ( ) ) ) )
        return pairs

    def coverage_2d ( self, bins = 80 ):
        allx = np.vstack ( [ self.rc ( k ) for k in self.keys if len ( self.rc ( k ) ) ] )
        h, ex, ey = np.histogram2d ( allx[:, 0], allx[:, 1], bins = bins )
        return h, ex, ey

    # ------------------------------------------------------------------ WHAM (own, for comparisons)
    def wham ( self, part = None, bins = None, tolerance = 1e-6, max_iterations = 20000 ):
        """ PMF (kJ/mol, minimum 0) from the RC samples with the standard WHAM equations.
            part: None (all), ( first, last ) fractions of each window's samples (e.g. ( 0, 0.5 )).
            1D: ( centers, F ); 2D: ( ( cx, cy ), F[nx, ny] ). Empty bins: nan. """
        beta = 1.0 / ( KB * self.temperature )
        data, targets, forces = [ ], [ ], [ ]
        for k in self.keys:
            x = self.rc ( k ); t = self.target ( k ); f = self.force ( k )
            if not len ( x ) or t is None: continue
            if part is not None:
                x = x[int ( part[0] * len ( x ) ):int ( part[1] * len ( x ) )]
            data.append ( x ); targets.append ( np.array ( t ) ); forces.append ( np.array ( f ) )
        allx = np.vstack ( data )
        if self.dims == 1:
            nb = bins or 100
            edges = [ np.linspace ( allx[:, 0].min ( ), allx[:, 0].max ( ), nb + 1 ) ]
        else:
            nb = bins or 60
            edges = [ np.linspace ( allx[:, d].min ( ), allx[:, d].max ( ), nb + 1 ) for d in range ( 2 ) ]
        centers = [ 0.5 * ( e[1:] + e[:-1] ) for e in edges ]
        grid = np.stack ( np.meshgrid ( *centers, indexing = "ij" ), axis = -1 ).reshape ( -1, self.dims )
        counts = [ ]
        for x in data:
            h = np.histogramdd ( x, bins = edges )[0].ravel ( )
            counts.append ( h )
        counts = np.array ( counts ); N = counts.sum ( axis = 1 )
        bias = np.array ( [ ( ( grid - t ) ** 2 * f ).sum ( axis = 1 ) for t, f in zip ( targets, forces ) ] )   # ( windows, bins )
        expb = np.exp ( -beta * ( bias - bias.min ( axis = 1, keepdims = True ) ) )
        shift = bias.min ( axis = 1 )
        numerator = counts.sum ( axis = 0 )
        f = np.zeros ( len ( data ) )
        for _ in range ( max_iterations ):
            denom = ( N[:, None] * np.exp ( beta * ( f - shift ) )[:, None] * expb ).sum ( axis = 0 )
            P = np.where ( denom > 0, numerator / np.where ( denom > 0, denom, 1 ), 0.0 )
            P /= P.sum ( )
            new = -np.log ( ( expb * P[None, :] ).sum ( axis = 1 ) ) / beta + shift
            new -= new[0]
            if np.max ( np.abs ( new - f ) ) < tolerance: f = new; break
            f = new
        with np.errstate ( divide = "ignore" ):
            F = -np.log ( P ) / beta
        F[~np.isfinite ( F )] = np.nan
        F -= np.nanmin ( F )
        if self.dims == 1: return centers[0], F
        return ( centers[0], centers[1] ), F.reshape ( len ( centers[0] ), len ( centers[1] ) )

    # ------------------------------------------------------------------ PMF files of pDynamo's WHAM
    def pmf_files ( self ):
        found = glob.glob ( os.path.join ( self.folder, "*pmf*.log" ) ) + glob.glob ( os.path.join ( self.folder, "*", "*pmf*.log" ) )
        # . *_pmf_errors.log (bootstrap errors, extra columns) is not a PMF file
        found = [ f for f in found if not f.endswith ( "_pmf_errors.log" ) ]
        return sorted ( set ( found ), key = lambda f: ( os.path.basename ( os.path.dirname ( f ) ) != "wham", f ) )

    # ------------------------------------------------------------------ coordinates
    def frame_files ( self, key, source = "dc" ):
        w = self.windows.get ( key, { } )
        if source == "opt":
            return [ w["opt"] ] if w.get ( "opt" ) else [ ]
        path = w.get ( source )
        return _numeric_frames ( path ) if path and os.path.isdir ( path ) else [ ]

    @staticmethod
    def read_frame ( path ):
        with open ( path, "rb" ) as handle: c3 = pickle.load ( handle )
        return np.array ( list ( c3 ), dtype = np.float32 ).reshape ( -1, 3 )

    def frames ( self, key, source = "dc", stride = 1, only = None ):
        files = self.frame_files ( key, source )
        idx = list ( range ( 0, len ( files ), max ( 1, int ( stride ) ) ) ) if only is None else [ i for i in only if 0 <= i < len ( files ) ]
        return [ self.read_frame ( files[i] ) for i in idx ], idx

    def representative ( self, key, source = "dc" ):
        """ Index of the frame whose restraint value is closest to the window's target
            (each dimension scaled by the window's own spread); ( index, distance ). """
        files = self.frame_files ( key, source )
        if not files: return None, None
        if source != "dc": return 0, None
        x = self.rc ( key ); t = self.target ( key )
        if not len ( x ) or t is None: return 0, None
        step = max ( 1, self.freq_geo // max ( 1, self.freq_res ) )
        scale = np.where ( x.std ( axis = 0 ) > 0, x.std ( axis = 0 ), 1.0 )
        best, dist = 0, None
        # . frame 0 is the starting structure, placed on the target: skipped when there are others
        for m in range ( 1 if len ( files ) > 1 else 0, len ( files ) ):
            r = m * step
            if r >= len ( x ): break
            d = float ( np.linalg.norm ( ( x[r] - np.array ( t ) ) / scale ) )
            if dist is None or d < dist: best, dist = m, d
        return best, dist

    def frame_rc ( self, key, frame_index ):
        x = self.rc ( key )
        r = frame_index * max ( 1, self.freq_geo // max ( 1, self.freq_res ) )
        return x[r] if r < len ( x ) else None


def read_pmf ( path ):
    """ pDynamo WHAM PMF file: 1D ( x, F ) or 2D ( ( xs, ys ), F[nx, ny] ) with inf -> nan. """
    rows = [ ]
    for line in open ( path ):
        parts = line.split ( )
        if not parts: continue
        try:
            rows.append ( [ float ( v ) for v in parts ] )
        except ValueError:
            continue
    a = np.array ( rows, dtype = float )
    a[~np.isfinite ( a )] = np.nan
    if a.shape[1] == 2:
        return 1, a[:, 0], a[:, 1]
    xs = np.unique ( a[:, 0] ); ys = np.unique ( a[:, 1] )
    F = np.full ( ( len ( xs ), len ( ys ) ), np.nan )
    ix = { v: n for n, v in enumerate ( xs ) }; iy = { v: n for n, v in enumerate ( ys ) }
    for x, y, f in a[:, :3]: F[ix[x], iy[y]] = f
    return 2, ( xs, ys ), F


def lowest_barrier_path ( F, start, end ):
    """ Path on a 2D grid (8 neighbours) that minimizes the highest F met, then the sum of F
        (the minimum energy path through the lowest saddle). nan cells are walls.
        Returns [ ( i, j ) ]. """
    nx, ny = F.shape
    G = np.where ( np.isfinite ( F ), F, np.inf )
    best = { start: ( G[start], G[start] ) }
    heap = [ ( G[start], G[start], start ) ]
    prev = { }
    while heap:
        peak, total, cell = heapq.heappop ( heap )
        if cell == end: break
        if best.get ( cell, ( np.inf, np.inf ) ) < ( peak, total ): continue
        i, j = cell
        for di in ( -1, 0, 1 ):
            for dj in ( -1, 0, 1 ):
                if not di and not dj: continue
                a, b = i + di, j + dj
                if not ( 0 <= a < nx and 0 <= b < ny ) or not np.isfinite ( G[a, b] ): continue
                cand = ( max ( peak, G[a, b] ), total + G[a, b] * math.hypot ( di, dj ) )
                if cand < best.get ( ( a, b ), ( np.inf, np.inf ) ):
                    best[( a, b )] = cand; prev[( a, b )] = cell
                    heapq.heappush ( heap, ( cand[0], cand[1], ( a, b ) ) )
    if end not in prev and end != start: return [ ]
    path, cell = [ end ], end
    while cell != start:
        cell = prev[cell]; path.append ( cell )
    return path[::-1]


def default_endpoints ( F ):
    """ The global minimum and the lowest cell at least 30 % of the grid diagonal away from it. """
    G = np.where ( np.isfinite ( F ), F, np.inf )
    a = np.unravel_index ( int ( np.argmin ( G ) ), G.shape )
    ii, jj = np.indices ( G.shape )
    far = np.hypot ( ii - a[0], jj - a[1] ) >= 0.3 * math.hypot ( *G.shape )
    H = np.where ( far, G, np.inf )
    b = np.unravel_index ( int ( np.argmin ( H ) ), G.shape )
    return tuple ( int ( v ) for v in a ), tuple ( int ( v ) for v in b )


def export_frames ( path, system, frames ):
    """ Frames ( list of N x 3 arrays ) written with pDynamo (.dcd, .xyz, .pdb by extension;
        .pdb / .xyz: one file per frame when the format holds one structure). """
    from pBabel import ExportTrajectory, ExportSystem
    from pScientific.Geometry3 import Coordinates3
    c3 = system.coordinates3
    saved = np.array ( list ( c3 ), dtype = float ).reshape ( -1, 3 )
    def put ( X ):
        for i in range ( len ( X ) ):
            for k in range ( 3 ): system.coordinates3[i, k] = float ( X[i][k] )
    try:
        ext = os.path.splitext ( path )[1].lower ( )
        if ext == ".dcd":
            trajectory = ExportTrajectory ( path, system, log = None )
            trajectory.WriteHeader ( )
            for X in frames:
                put ( X ); trajectory.WriteOwnerData ( )
            trajectory.WriteFooter ( ); trajectory.Close ( )
            return [ path ]
        if ext == ".xyz":
            from pScientific import PeriodicTable
            symbols = [ PeriodicTable.Symbol ( a.atomicNumber ) for a in system.atoms ]
            with open ( path, "w" ) as handle:
                for m, X in enumerate ( frames ):
                    handle.write ( "{}\nframe {}\n".format ( len ( symbols ), m ) )
                    for s, p in zip ( symbols, X ): handle.write ( "{:2s} {:12.5f} {:12.5f} {:12.5f}\n".format ( s, *p ) )
            return [ path ]
        out = [ ]
        root, ext = os.path.splitext ( path )
        for m, X in enumerate ( frames ):
            put ( X )
            name = path if len ( frames ) == 1 else "{}_{:03d}{}".format ( root, m, ext or ".pdb" )
            ExportSystem ( name, system, log = None ); out.append ( name )
        return out
    finally:
        put ( saved )
