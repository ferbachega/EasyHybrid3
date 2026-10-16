#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: DYFF with user-edited parameters
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      DYFF generates every bonded parameter with "factories" (pMolecule/
#      MMModel/DYFFUtilities.py) that return ( value, key ), the key being
#      the atom types plus the bond tags (e.g. "C:Res:C:Res:C:Res:C:Res:S").
#      Its parameter files are empty, and pDynamo's reader of explicit
#      cosine rows has a bug, so edits cannot go through files. Instead,
#      MMModelDYFFEdited wraps the factories: when the generated key has a
#      user edit, the edited value is returned instead of the rule-based one.
#
#      edits = { "bond": { key: ( b0, k ) }, "angle" | "dihedral" | "outofplane":
#                { key: [ ( force constant, period ), ... ] } }   (kJ/mol units)
#
#      The edits live on the Builder object (vismol_object.dyff_parameter_
#      edits, kept by Undo) and are used by the DYFF Parameters window,
#      define_MMModel('DYFF') and the Builder's Optimize.
#

from pMolecule.MMModel import MMModelDYFF

_FACTORY_TERM = { "harmonicBondParameters"     : "bond",
                  "cosineAngleParameters"      : "angle",
                  "cosineDihedralParameters"   : "dihedral",
                  "cosineOutOfPlaneParameters" : "outofplane" }
KEY_SEPARATOR = ":"


class MMModelDYFFEdited ( MMModelDYFF ):
    """ DYFF with parameters overridden by key (see the module docstring). """

    _attributable = dict ( MMModelDYFF._attributable )
    _attributable.update ( { "parameterEdits" : None } )

    def ParameterFactories ( self, atomTypes ):
        factories = super ( MMModelDYFFEdited, self ).ParameterFactories ( atomTypes )
        edits = self.parameterEdits or { }
        wrapped = { }
        for name, factory in factories.items ( ):
            term  = _FACTORY_TERM.get ( name )
            table = edits.get ( term ) or { }
            if not table:
                wrapped[name] = factory
                continue
            def edited ( types, indices, connectivity, key, factory = factory, table = table, term = term ):
                value, new_key = factory ( types, indices, connectivity, key )
                override = table.get ( KEY_SEPARATOR.join ( new_key ) )
                if override is None:
                    return value, new_key
                if term == "bond":
                    return ( float ( override[0] ), float ( override[1] ) ), new_key
                return [ ( float ( fc ), int ( n ) ) for ( fc, n ) in override ], new_key
            wrapped[name] = edited
        return wrapped


def dyff_model ( parameter_set = "dyff-1.0", edits = None ):
    """ MMModelDYFF, or MMModelDYFFEdited when there are edits. """
    edits = { term: dict ( table ) for term, table in ( edits or { } ).items ( ) if table }
    if not edits:
        return MMModelDYFF.WithParameterSet ( parameter_set )
    model = MMModelDYFFEdited.WithParameterSet ( parameter_set )
    model.parameterEdits = edits
    return model
