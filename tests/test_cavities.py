"""
Cavity detection (src/util/cavities.py) on geometries with a known answer:
a closed shell of atoms around a spherical void, the same shell with a
narrow or a wide opening, a box region, and a pocket followed along frames.
"""

import numpy as np
import pytest

pytest.importorskip("pScientific")
import pBabel  # noqa: F401  (pDynamo import order)
from util import cavities as CV


def _sphere(radius, spacing):
    n = int(4 * np.pi * radius * radius / (spacing * spacing * 0.87))
    k = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * k / n)
    theta = np.pi * (1 + 5 ** 0.5) * k
    return radius * np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


R_SHELL, R_ATOM = 9.0, 1.7


@pytest.fixture(scope="module")
def shell():
    # . two layers so that no probe passes between the atoms
    coords = np.vstack([_sphere(R_SHELL, 1.6), _sphere(R_SHELL + 1.8, 1.6)])
    return coords, np.full(len(coords), R_ATOM)


def test_closed_shell_is_one_buried_cavity(shell):
    coords, radii = shell
    result = CV.detect_cavities(coords, radii)
    assert len(result["cavities"]) == 1
    cavity = result["cavities"][0]
    ball = 4.0 / 3.0 * np.pi * (R_SHELL - R_ATOM) ** 3
    # . the atom bumps take a little of the ball
    assert 0.90 * ball < cavity["volume"] < ball
    assert cavity["buried"]
    assert np.linalg.norm(cavity["centre"]) < 0.5
    # . every lining atom is on the inner layer
    assert np.all(np.linalg.norm(coords[cavity["lining_atoms"]], axis=1) < R_SHELL + 0.5)
    mesh = CV.cavity_mesh(result, 1)
    assert mesh is not None
    sphere_area = 4 * np.pi * (R_SHELL - R_ATOM) ** 2
    assert 0.95 * sphere_area < cavity["area"] < 1.15 * sphere_area


def test_narrow_opening_gives_open_pocket(shell):
    coords, radii = shell
    keep = coords[:, 2] < 8.3        # . free opening smaller than Probe Out
    result = CV.detect_cavities(coords[keep], radii[keep])
    assert len(result["cavities"]) == 1
    cavity = result["cavities"][0]
    assert not cavity["buried"]
    assert 10.0 < cavity["depth_max"] < 2 * (R_SHELL - R_ATOM)


def test_wide_opening_is_bulk(shell):
    coords, radii = shell
    keep = coords[:, 2] < 6.5        # . Probe Out gets in: no cavity
    result = CV.detect_cavities(coords[keep], radii[keep])
    assert result["cavities"] == []


def test_box_keeps_only_points_inside(shell):
    coords, radii = shell
    whole = CV.detect_cavities(coords, radii, lining=False)["cavities"][0]["volume"]
    box = (np.zeros(3), np.full(3, 3.0))
    result = CV.detect_cavities(coords, radii, box=box, lining=False)
    assert len(result["cavities"]) == 1
    points = CV.grid_points(result, result["cavities"][0]["points"])
    assert np.all(np.abs(points) <= 3.0 + 1e-6)
    assert result["cavities"][0]["volume"] < whole


def test_label_components():
    mask = np.zeros((6, 6, 6), dtype=bool)
    mask[0:2, 0:2, 0:2] = True
    mask[4:6, 4:6, 4:6] = True
    mask[3, 0, 0] = True             # . touches the first block only diagonally: separate
    labels, n = CV.label_components(mask)
    assert n == 3
    assert len(np.unique(labels[mask])) == 3


def test_pocket_followed_along_frames(shell):
    coords, radii = shell
    reference = CV.detect_cavities(coords, radii)["cavities"][0]["lining_atoms"]
    shifted = coords + np.array([5.0, -3.0, 2.0])          # . the molecule drifts
    rows = CV.pocket_volume_trajectory([coords, shifted], radii, reference, processes=1)
    assert rows[0]["cavities"] == [1] and rows[1]["cavities"] == [1]
    assert abs(rows[0]["volume"] - rows[1]["volume"]) < 0.05 * rows[0]["volume"]


def test_points_mesh_and_files(shell, tmp_path):
    coords, radii = shell
    result = CV.detect_cavities(coords, radii)
    vertices, colors, indexes, normals = CV.points_mesh(result, 1)
    n_points = len(result["cavities"][0]["points"])
    assert len(vertices) == n_points * 32 * 3 * 3 and len(vertices) == len(normals) == len(colors)
    CV.write_points_pdb(result, str(tmp_path / "cav.pdb"))
    assert sum(1 for line in open(tmp_path / "cav.pdb") if line.startswith("HETATM")) == n_points
    CV.write_table_csv(result, str(tmp_path / "cav.csv"))
    assert open(tmp_path / "cav.csv").read().splitlines()[1].startswith("1,")
