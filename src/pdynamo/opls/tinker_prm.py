#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Tinker OPLS parameter file reader
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Reads Tinker OPLS parameter files (oplsaal.prm = OPLS-AA/L,
#      oplsaa08.prm = OPLS-AA 2008 / BOSS 4.8) and converts every term to
#      pDynamo's conventions, keyed by atom CLASS SYMBOL (CT, HC, CA, ...)
#      -- the names pDynamo's OPLS parameter sets use.
#
#      Unit/convention conversions (checked against pDynamo's own
#      opls/protein set, see parameters.compare_with_pdynamo()):
#        bond     Tinker  k (b - b0)^2, "bond c1 c2 k b0"     -> ( b0, k )      same form
#        angle    Tinker  k (t - t0)^2, "angle c1 c2 c3 k t0" -> ( t0, k )      same form
#        torsion  Tinker  torsionunit 0.5:  0.5 V (1 + cos ( n phi - d ))
#                 pDynamo k ( cos ( n phi - d ) + 1 )                         -> k = V / 2
#        imptors  Tinker  imptorunit 0.5, "imptors a b CENTRAL d V 180 2"
#                 pDynamo central + sorted peripherals, K                     -> K = V / 2
#        vdw      sigma (radiustype SIGMA, radiussize DIAMETER), epsilon      same
#        ureybrad "ureybrad c1 c2 c3 k r0"                                    -> ( r0, k )
#      Class 0 is Tinker's wild card -> "X".
#

import re
from collections import OrderedDict

WILD = "X"


class TinkerPRM:
    """ The content of one Tinker .prm file. """
    def __init__ ( self, path ):
        self.path      = path
        self.forcefield = None
        self.vdwindex  = "CLASS"
        self.torsionunit = 1.0
        self.imptorunit  = 1.0
        self.atoms     = OrderedDict ( )   # type -> { class, symbol, description, z, mass, valence }
        self.vdw       = { }               # type or class -> ( sigma, epsilon )
        self.charges   = { }               # type -> charge
        self.bonds     = [ ]               # ( c1, c2, b0, k )
        self.angles    = [ ]               # ( c1, c2, c3, t0, k )
        self.torsions  = [ ]               # ( c1, c2, c3, c4, [ ( V, phase, n ) ] )
        self.imptors   = [ ]               # ( c1, c2, c3, c4, V, phase, n )
        self.ureybrad  = [ ]               # ( c1, c2, c3, r0, k )
        self.read ( )

    def read ( self ):
        quoted = re.compile ( r'"[^"]*"' )
        with open ( self.path, errors = "replace" ) as handle:
            for line in handle:
                tokens = line.split ( )
                if not tokens or tokens[0].startswith ( "#" ): continue
                key = tokens[0].lower ( )
                try:
                    if   key == "forcefield":  self.forcefield  = " ".join ( tokens[1:] )
                    elif key == "vdwindex":    self.vdwindex    = tokens[1].upper ( )
                    elif key == "torsionunit": self.torsionunit = float ( tokens[1] )
                    elif key == "imptorunit":  self.imptorunit  = float ( tokens[1] )
                    elif key == "atom":
                        match = quoted.search ( line )
                        description = match.group ( 0 ).strip ( '"' ).strip ( ) if match else ""
                        head = line[:match.start ( )].split ( ) if match else tokens[:4]
                        tail = line[match.end ( ):].split ( ) if match else tokens[5:]
                        atom_type, atom_class, symbol = int ( head[1] ), int ( head[2] ), head[3]
                        self.atoms[atom_type] = { "class": atom_class, "symbol": symbol, "description": description,
                                                  "z": int ( tail[0] ), "mass": float ( tail[1] ),
                                                  "valence": int ( tail[2] ) if len ( tail ) > 2 else None }
                    elif key == "vdw":      self.vdw[int ( tokens[1] )] = ( float ( tokens[2] ), float ( tokens[3] ) )
                    elif key == "charge":   self.charges[int ( tokens[1] )] = float ( tokens[2] )
                    elif key == "bond":     self.bonds.append ( ( int ( tokens[1] ), int ( tokens[2] ), float ( tokens[4] ), float ( tokens[3] ) ) )
                    elif key == "angle":    self.angles.append ( ( int ( tokens[1] ), int ( tokens[2] ), int ( tokens[3] ), float ( tokens[5] ), float ( tokens[4] ) ) )
                    elif key == "ureybrad": self.ureybrad.append ( ( int ( tokens[1] ), int ( tokens[2] ), int ( tokens[3] ), float ( tokens[5] ), float ( tokens[4] ) ) )
                    elif key == "torsion":
                        classes = tuple ( int ( t ) for t in tokens[1:5] )
                        values  = tokens[5:]
                        terms   = [ ]
                        for k in range ( 0, len ( values ) - 2, 3 ):
                            terms.append ( ( float ( values[k] ), float ( values[k + 1] ), int ( values[k + 2] ) ) )
                        self.torsions.append ( classes + ( terms, ) )
                    elif key == "imptors":
                        self.imptors.append ( tuple ( int ( t ) for t in tokens[1:5] ) +
                                              ( float ( tokens[5] ), float ( tokens[6] ), int ( tokens[7] ) ) )
                except ( IndexError, ValueError ):
                    continue          # malformed line: ignored

    # . Class symbols ---------------------------------------------------------
    def class_symbols ( self ):
        """ { class number: symbol } (0 -> wild card). A class used with
            more than one symbol keeps the first one met. """
        symbols = { 0: WILD }
        for data in self.atoms.values ( ):
            symbols.setdefault ( data["class"], data["symbol"] )
        return symbols

    def type_lj ( self, atom_type ):
        """ ( sigma, epsilon ) of a TYPE, whatever the vdw index. """
        data = self.atoms.get ( atom_type )
        if data is None: return None
        key = atom_type if self.vdwindex == "TYPE" else data["class"]
        return self.vdw.get ( key )

    # . Terms in pDynamo conventions, keyed by class symbols -------------------
    def converted_terms ( self ):
        """ { "bond": [ ( (s1, s2), (b0, k) ) ], "angle": [ ( (s1, s2, s3), (t0, k) ) ],
              "dihedral": [ ( (s1..s4), [ (k, n, phase) ] ) ],
              "outofplane": [ ( (central, p1, p2, p3), K ) ],
              "ureybradley": [ ( (s1, s2, s3), (r0, k) ) ] }
            in file order (several class numbers may give the same symbols). """
        sym = self.class_symbols ( )
        get = lambda c: sym.get ( c, "?{}".format ( c ) )
        out = { "bond": [ ], "angle": [ ], "dihedral": [ ], "outofplane": [ ], "ureybradley": [ ] }
        for ( c1, c2, b0, k ) in self.bonds:
            out["bond"].append ( ( ( get ( c1 ), get ( c2 ) ), ( b0, k ) ) )
        for ( c1, c2, c3, t0, k ) in self.angles:
            out["angle"].append ( ( ( get ( c1 ), get ( c2 ), get ( c3 ) ), ( t0, k ) ) )
        for ( c1, c2, c3, c4, terms ) in self.torsions:
            converted = [ ( self.torsionunit * V, n, phase ) for ( V, phase, n ) in terms if V != 0.0 ]
            out["dihedral"].append ( ( ( get ( c1 ), get ( c2 ), get ( c3 ), get ( c4 ) ), converted ) )
        for ( c1, c2, c3, c4, V, phase, n ) in self.imptors:
            # . Tinker: the third atom is the central one.
            out["outofplane"].append ( ( ( get ( c3 ), get ( c1 ), get ( c2 ), get ( c4 ) ), self.imptorunit * V ) )
        for ( c1, c2, c3, r0, k ) in self.ureybrad:
            out["ureybradley"].append ( ( ( get ( c1 ), get ( c2 ), get ( c3 ) ), ( r0, k ) ) )
        return out
