"""
util/pdb_tools.dedupe_pdb_atom_names (run on every plain-PDB import):

  1) atoms sharing a name in one residue (ligand written as "C", "C", "O")
     are renamed element + count, so pDynamo does not drop them;
  2) the same name under different alternate locations (ALTLOC A/B) is NOT a
     duplicate: such residues keep their PDB names (3KS3 HIS 4, HOH 353 used to
     become N1, C1, ... / O1, O2, which broke the Prepare DYFF templates);
  3) insertion codes separate residues (52 and 52A are not one residue).
"""
import os
import sys

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_TESTS_DIR), "src"))
from util.pdb_tools import dedupe_pdb_atom_names  # noqa: E402


def _line(serial, name, resname, resseq, x, element, altloc=" ", icode=" "):
    return "HETATM{:5d} {:<4s}{:1s}{:>3s} A{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}  1.00  0.00          {:>2s}\n".format(
        serial, (" " + name) if len(name) < 4 else name, altloc, resname, resseq, icode, x, 0.0, 0.0, element)


def _names(path):
    return [l[12:16].strip() for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def test_real_duplicates_are_renamed(tmp_path):
    p = tmp_path / "lig.pdb"
    p.write_text(_line(1, "C", "UNL", 1, 0.0, "C") + _line(2, "C", "UNL", 1, 1.5, "C") + _line(3, "O", "UNL", 1, 3.0, "O"))
    assert _names(dedupe_pdb_atom_names(str(p))) == ["C1", "C2", "O1"]


def test_alternate_locations_are_kept(tmp_path):
    text = (_line(1, "O", "HOH", 353, 0.0, "O", altloc="A") + _line(2, "O", "HOH", 353, 0.5, "O", altloc="B") +
            _line(3, "O", "HOH", 52, 9.0, "O") + _line(4, "O", "HOH", 52, 12.0, "O", icode="A"))
    p = tmp_path / "alt.pdb"
    p.write_text(text)
    dedupe_pdb_atom_names(str(p))
    assert p.read_text() == text


def test_3ks3_is_untouched(tmp_path):
    src = os.path.join(_TESTS_DIR, "data", "metalloproteins", "3KS3.pdb")
    out = tmp_path / "3KS3.pdb"
    dedupe_pdb_atom_names(src, str(out))
    assert out.read_text() == open(src).read()
