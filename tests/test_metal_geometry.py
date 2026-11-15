"""
Transition-metal geometry perception for DYFF (pdynamo.metal_geometry, dyff_edits; 2026-10-13).
Synthetic complexes with ideal polyhedra: the coordinate geometry must reach the DYFF metal type.
"""
import contextlib, io, math, os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from pdynamo import metal_geometry as MG

SQ = [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)]
TET = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
TBP = [(0, 0, 1), (0, 0, -1), (1, 0, 0), (-0.5, math.sqrt(3) / 2, 0), (-0.5, -math.sqrt(3) / 2, 0)]
SPY = [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0), (0, 0, 1)]
OCT = SQ + [(0, 0, 1), (0, 0, -1)]


def complex_xyz(directions, d=2.0):
    v = np.array(directions, float); v = d * v / np.linalg.norm(v, axis=1)[:, None]
    return np.vstack([[0.0, 0.0, 0.0], v])


@pytest.mark.parametrize("dirs,expected", [(SQ, "SPl"), (TET, "Tet"), (TBP, "TBp"), (SPY, "SPy"), (OCT, "Oct"), (SQ[:2] + [(-1, 0, 0)], "TSh")])
def test_classify_ideal_polyhedra(dirs, expected):
    xyz = complex_xyz(dirs)
    bonds = [(0, i) for i in range(1, len(xyz))]
    out = MG.perceive([78] + [17] * (len(xyz) - 1), xyz, bonds)
    assert out[0]["class"] == expected


def test_hapto_ligand_is_one_site():
    # . eta5 ring (5 C bonded to each other and to the metal) + 3 monodentate ligands -> 4 sites (piano stool)
    ring = [(1.2 * math.cos(2 * math.pi * k / 5), 1.2 * math.sin(2 * math.pi * k / 5), 1.7) for k in range(5)]
    legs = [(1.9 * math.cos(2 * math.pi * k / 3), 1.9 * math.sin(2 * math.pi * k / 3), -0.9) for k in range(3)]
    xyz = np.vstack([[0, 0, 0], ring, legs])
    bonds = [(0, i) for i in range(1, 9)] + [(1 + k, 1 + (k + 1) % 5) for k in range(5)]
    info = MG.perceive([26] + [6] * 8, xyz, bonds)[0]
    assert info["hapto"] == [5] and len(info["sites"]) == 4 and info["tag"] == "Tet"


def test_tbp_ideal_angles():
    xyz = complex_xyz(TBP)
    out = MG.perceive([26] + [6] * 5, xyz, [(0, i) for i in range(1, 6)])[0]
    vals = sorted(set(out["angles"].values()))
    assert vals == [90.0, 120.0, 180.0]


def _dyff_metal_type(smiles_free_metal_complex):
    """ DYFF type of the metal of a small RDKit-built complex (coordinates embedded by hand). """
    pytest.importorskip("pBabel")                 # . pBabel before pMolecule (pDynamo circular import)
    from pMolecule import Atom, Bond, BondType, Connectivity, ConvertInputConnectivity, System
    from pScientific.Geometry3 import Coordinates3
    from pdynamo.dyff_edits import dyff_model
    z, xyz, bonds, fc = smiles_free_metal_complex
    c = Connectivity()
    for k, (zz, q) in enumerate(zip(z, fc)): c.AddNode(Atom.WithOptions(atomicNumber=zz, label="A%d" % (k + 1), formalCharge=q))
    for i, j in bonds: c.AddEdge(Bond.WithNodes(c.nodes[i], c.nodes[j], type=BondType.Single))
    with contextlib.redirect_stdout(io.StringIO()):
        ConvertInputConnectivity(c, {})
        s = System.FromConnectivity(connectivity=c)
        c3 = Coordinates3.WithExtent(len(z))
        for i, p in enumerate(xyz):
            for k in range(3): c3[i, k] = float(p[k])
        s.coordinates3 = c3
        s.DefineMMModel(dyff_model("dyff-1.0", None), log=None)
    st = s.mmState
    return st.atomTypes[st.atomTypeIndices[0]]


def test_square_planar_pt_chloride_gets_spl_type():
    # . [PtCl4]2- in the covalent model: Cl neutral (single bond), Pt formal -2 (ox +2)
    xyz = complex_xyz(SQ, 2.3)
    label = _dyff_metal_type(([78] + [17] * 4, xyz, [(0, i) for i in range(1, 5)], [-2, 0, 0, 0, 0]))
    assert label.startswith("Pt:SPl")


def test_tetrahedral_zn_chloride_gets_tet_type():
    xyz = complex_xyz(TET, 2.25)
    label = _dyff_metal_type(([30] + [17] * 4, xyz, [(0, i) for i in range(1, 5)], [-2, 0, 0, 0, 0]))
    assert label.startswith("Zn:Tet")


def test_pyramidal_zn_gets_spy_type():
    xyz = complex_xyz(SPY, 2.25)
    label = _dyff_metal_type(([30] + [17] * 5, xyz, [(0, i) for i in range(1, 6)], [-3, 0, 0, 0, 0, 0]))
    assert label.startswith("Zn:SPy")


def test_switch_restores_pdynamo_geometry(monkeypatch):
    # . pDynamo alone (switch on) sees 4 neighbours as Tet; with the perceiver the square-planar Pd is SPl
    xyz = complex_xyz(SQ, 2.3)
    spec = ([46] + [17] * 4, xyz, [(0, i) for i in range(1, 5)], [-2, 0, 0, 0, 0])
    monkeypatch.setenv("EASYHYBRID_DYFF_NO_METAL_GEOMETRY", "1")
    try:
        off = _dyff_metal_type(spec)
    except Exception:
        off = None                                   # . untyped is also pDynamo's behaviour
    assert off is None or off.split(":")[1] == "Tet"
    monkeypatch.delenv("EASYHYBRID_DYFF_NO_METAL_GEOMETRY")
    assert _dyff_metal_type(spec).startswith("Pd:SPl")


def test_high_spin_flag_lengthens_metal_bonds():
    # . [FeCl6]4- octahedron, covalent model (Cl neutral, Fe -4 -> ox +2): DYFF optimum with and without the
    #   high-spin flag (dyff_edits.HIGH_SPIN_BOND_OFFSETS[26])
    pytest.importorskip("pBabel")
    import pBabel  # noqa: F401
    from pMolecule import Atom, Bond, BondType, Connectivity, ConvertInputConnectivity, System
    from pMolecule.NBModel import NBModelFull
    from pScientific.Geometry3 import Coordinates3
    from pSimulation import LBFGSMinimize_SystemGeometry
    from pdynamo.dyff_edits import dyff_model, HIGH_SPIN_BOND_OFFSETS
    xyz = complex_xyz(OCT, 2.4)

    def fe_cl(high_spin):
        c = Connectivity()
        for k, (z, q) in enumerate(zip([26] + [17] * 6, [-4] + [0] * 6)):
            c.AddNode(Atom.WithOptions(atomicNumber=z, label="A%d" % (k + 1), formalCharge=q))
        for i in range(1, 7): c.AddEdge(Bond.WithNodes(c.nodes[0], c.nodes[i], type=BondType.Single))
        with contextlib.redirect_stdout(io.StringIO()):
            ConvertInputConnectivity(c, {}); s = System.FromConnectivity(connectivity=c)
            c3 = Coordinates3.WithExtent(7)
            for i, p in enumerate(xyz):
                for k in range(3): c3[i, k] = float(p[k])
            s.coordinates3 = c3
            s.DefineMMModel(dyff_model("dyff-1.0", None, high_spin_atoms={0} if high_spin else None), log=None)
            for i in range(7): s.mmState.charges[i] = 0.0
            s.DefineNBModel(NBModelFull.WithDefaults())
            LBFGSMinimize_SystemGeometry(s, log=None, maximumIterations=500, rmsGradientTolerance=0.01)
        return float(np.mean([np.linalg.norm([s.coordinates3[i, k] - s.coordinates3[0, k] for k in range(3)]) for i in range(1, 7)]))
    low, high = fe_cl(False), fe_cl(True)
    assert abs((high - low) - HIGH_SPIN_BOND_OFFSETS[26]) < 0.05


def test_unbonded_coordinated_water_counts_as_a_site():
    # . Zn with 3 bonded donors + 1 water kept WITHOUT a metal-O bond (protein preparation): tetrahedral, not trigonal
    xyz = complex_xyz(TET, 2.05)
    o = xyz[4]; u = o / np.linalg.norm(o); t = np.cross(u, [1.0, 0.0, 0.0]); t /= np.linalg.norm(t)
    h1 = o + 0.96 * (0.6 * u + 0.8 * t); h2 = o + 0.96 * (0.6 * u - 0.8 * t)
    xyz = np.vstack([xyz, h1, h2])
    z = [30, 7, 7, 7, 8, 1, 1]
    bonds = [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6)]
    info = MG.perceive(z, xyz, bonds)[0]
    assert info["class"] == "Tet" and info["unbonded_waters"] == [4]
    assert MG.perceive(z, xyz, bonds, include_waters=False)[0]["class"] in ("Tri", "Pyr")
