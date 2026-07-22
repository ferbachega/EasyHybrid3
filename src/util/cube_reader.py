#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Gaussian Cube file reader (orbitals/density/potential from external
#          QC programs, e.g. ORCA via orca_plot)
#
#  Description:
#      Reads volumetric scalar grid data in the Gaussian Cube format (the
#      same format written by ORCA's orca_plot utility, Gaussian, and most
#      other quantum chemistry packages). Returns plain numpy arrays --
#      no dependency on pDynamo3 in this module, so it can be tested and
#      reused independently of the rest of EasyHybrid.
#
#      Format reference: the file has a 2-line title/comment, then:
#        line 3: NATOMS  X0 Y0 Z0                 (grid origin, Bohr)
#        line 4: NX  dX1 dY1 dZ1                  (voxel vector along axis 1)
#        line 5: NY  dX2 dY2 dZ2                  (voxel vector along axis 2)
#        line 6: NZ  dX3 dY3 dZ3                  (voxel vector along axis 3)
#        NATOMS lines: atomic_number  charge  x y z   (Bohr)
#        then NX*NY*NZ scalar values, whitespace separated (any number of
#        values per line -- most writers use 6, but readers should not
#        assume a fixed count per line), value ordering: Z varies fastest,
#        then Y, then X slowest (row-major / C order for shape (NX,NY,NZ)).
#
#      NATOMS on line 3 can be negative in some variants (signals extra
#      data such as orbital index/count on the same line -- "multi-cube"
#      files with several values per grid point, e.g. several MOs in one
#      file). This reader does NOT support that variant -- only the
#      standard single-value-per-point cube, which is what orca_plot
#      writes for one orbital/density/potential at a time.
#
from util.debug import dprint
import numpy as np


BOHR_PER_ANGSTROM = 1.889725989


class CubeFileError ( Exception ):
    """ Erro ao ler ou interpretar um arquivo .cube. """
    pass


class CubeGrid:
    """ Representa os dados lidos de um arquivo .cube.

    Atributos:
        title, comment : str            -- as duas primeiras linhas do arquivo
        natoms          : int
        atoms           : list de (atomic_number:int, charge:float, x:float, y:float, z:float), em Bohr
        origin          : np.ndarray shape (3,)          -- origem do grid, em Bohr
        voxel_vectors   : np.ndarray shape (3,3)          -- vetor de cada eixo do grid (linha i = eixo i), em Bohr
        dims            : tuple (nx, ny, nz)
        values          : np.ndarray shape (nx, ny, nz)   -- dado escalar (densidade/orbital/potencial)
        is_orthogonal   : bool  -- True se voxel_vectors for diagonal (o caso comum;
                                    RegularGrid do pDynamo so aceita esse caso)
    """
    __slots__ = ( "title", "comment", "natoms", "atoms", "origin",
                  "voxel_vectors", "dims", "values", "is_orthogonal" )

    def __init__ ( self, title, comment, natoms, atoms, origin, voxel_vectors, dims, values ):
        self.title         = title
        self.comment       = comment
        self.natoms        = natoms
        self.atoms         = atoms
        self.origin        = origin
        self.voxel_vectors = voxel_vectors
        self.dims          = dims
        self.values        = values
        # so a diagonal e nao-zero (dentro de uma tolerancia numerica pequena)
        off_diag = voxel_vectors - np.diag ( np.diag ( voxel_vectors ) )
        self.is_orthogonal = bool ( np.allclose ( off_diag, 0.0, atol = 1e-8 ) )

    @property
    def spacing ( self ):
        """ Espacamento (dx, dy, dz) em Bohr -- só faz sentido se is_orthogonal. """
        return ( self.voxel_vectors[0,0], self.voxel_vectors[1,1], self.voxel_vectors[2,2] )

    def value_range ( self ):
        return ( float ( self.values.min ( ) ), float ( self.values.max ( ) ) )


def read_cube_file ( path ):
    """ Le um arquivo .cube e devolve um CubeGrid.

    Levanta CubeFileError se o arquivo nao seguir o formato esperado, ou
    se o numero de valores lidos nao bater com nx*ny*nz declarado no
    cabecalho (arquivo truncado/corrompido). """
    with open ( path, "r" ) as f:
        lines = f.readlines ( )

    if len ( lines ) < 6:
        raise CubeFileError ( "The .cube file is too short (fewer than 6 header lines): {}".format ( path ) )

    title   = lines[0].rstrip ( "\n" )
    comment = lines[1].rstrip ( "\n" )

    try:
        line3 = lines[2].split ( )
        natoms_raw = int ( line3[0] )
        origin     = np.array ( [ float ( line3[1] ), float ( line3[2] ), float ( line3[3] ) ], dtype = np.float64 )
    except ( IndexError, ValueError ) as error:
        raise CubeFileError ( "Could not parse line 3 (NATOMS/origin) of {}: {}".format ( path, error ) )

    if natoms_raw < 0:
        raise CubeFileError (
            "A negative NATOMS (line 3) indicates a cube with multiple values per "
            "grid point (e.g. several orbitals in the same file) -- a format not "
            "supported by this reader. Generate one cube per orbital/density/"
            "potential at a time instead (which is what orca_plot does by default)."
        )
    natoms = natoms_raw

    dims          = []
    voxel_vectors = np.zeros ( (3,3), dtype = np.float64 )
    for axis in range ( 3 ):
        try:
            parts = lines[3+axis].split ( )
            n     = int ( parts[0] )
            vec   = [ float ( parts[1] ), float ( parts[2] ), float ( parts[3] ) ]
        except ( IndexError, ValueError ) as error:
            raise CubeFileError ( "Could not parse line {} (axis {} dimension) of {}: {}".format ( 4+axis, axis, path, error ) )
        if n <= 0:
            raise CubeFileError (
                "Negative or zero voxel count on axis {} -- this reader only "
                "supports the standard convention (units always in Bohr, voxel "
                "count always positive). File: {}".format ( axis, path )
            )
        dims.append ( n )
        voxel_vectors[axis, :] = vec
    nx, ny, nz = dims

    atom_start = 6
    atoms = []
    for i in range ( natoms ):
        try:
            parts = lines[atom_start + i].split ( )
            atomic_number = int   ( parts[0] )
            charge        = float ( parts[1] )
            x, y, z       = float ( parts[2] ), float ( parts[3] ), float ( parts[4] )
        except ( IndexError, ValueError ) as error:
            raise CubeFileError ( "Could not parse atom line {} of {}: {}".format ( i, path, error ) )
        atoms.append ( ( atomic_number, charge, x, y, z ) )

    data_start = atom_start + natoms
    flat_values = []
    for line in lines[data_start:]:
        line = line.strip ( )
        if line == "":
            continue
        flat_values.extend ( line.split ( ) )

    expected = nx * ny * nz
    if len ( flat_values ) != expected:
        raise CubeFileError (
            "The number of values read ({}) does not match the expected nx*ny*nz ({}) "
            "in {} -- is the file truncated or corrupted?".format ( len ( flat_values ), expected, path )
        )

    values = np.array ( flat_values, dtype = np.float64 ).reshape ( ( nx, ny, nz ) )  # Z mais rapido (C order)

    return CubeGrid ( title, comment, natoms, atoms, origin, voxel_vectors, ( nx, ny, nz ), values )


if __name__ == "__main__":
    import sys
    if len ( sys.argv ) < 2:
        dprint ( "usage: python3 cube_reader.py file.cube" )
        sys.exit ( 1 )
    grid = read_cube_file ( sys.argv[1] )
    dprint ( "title   :", grid.title )
    dprint ( "comment :", grid.comment )
    dprint ( "natoms  :", grid.natoms )
    dprint ( "atoms   :", grid.atoms )
    dprint ( "origin (Bohr):", grid.origin )
    dprint ( "dims (nx,ny,nz):", grid.dims )
    dprint ( "is_orthogonal:", grid.is_orthogonal )
    dprint ( "spacing (Bohr):", grid.spacing )
    dprint ( "value range:", grid.value_range ( ) )
