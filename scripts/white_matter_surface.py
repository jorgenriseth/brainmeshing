from pathlib import Path

import click
import numpy as np
import pyvista as pv
import tempfile
from typing import Optional

from gonzo.simple_mri import load_mri
from brainmeshing.surfaces import (
    fs_surf_to_stl,
    hemisphere_surface_refinement,
    grow_white_connective_tissue,
    surface_union,
    svmtk2pyvista,
    pyvista2svmtk,
)


@click.command(name="ventricle-surface")
@click.option("--fs_dir", type=Path, required=True)
@click.option("--seg", type=Path, required=True)
@click.option("--ventricles", "--ventricle_surf", type=Path, required=True)
@click.option("--output", type=Path, required=True)
@click.option("--edge_length", type=float, default=0.5)
@click.option("--remesh_iter", type=int, default=1)
@click.option("--gapsize", type=float, default=0.2)
@click.option("--tmpdir", type=Path)
def main(
    fs_dir: Path,
    seg: Path,
    ventricles: Path,
    output: Path,
    edge_length: float,
    remesh_iter: int,
    gapsize: float,
    tmpdir: Optional[Path] = None,
):
    tempdir = tempfile.TemporaryDirectory()
    tmppath = Path(tempdir.name) if tmpdir is None else tmpdir
    print(tmppath)
    (output.parent).mkdir(exist_ok=True)
    ventricle_surf = pv.read(ventricles)
    fs_surf_to_stl(fs_dir / "surf", tmppath)
    seg_mri = load_mri(seg, dtype=np.int16)
    hemisphere_surface_refinement(
        tmppath / "lh_white.stl",
        tmppath / "rh_white.stl",
        tmppath,
        max_edge_length=edge_length,
        remesh_iter=remesh_iter,
        fix_boundaries=False,
        gapsize=gapsize,
        overlapping=True,
    )
    rh_white_refined = pv.read(tmppath / "rh_white_refined.stl")
    lh_white_refined = pv.read(tmppath / "lh_white_refined.stl")
    connective = grow_white_connective_tissue(seg_mri)
    white_svm = surface_union(lh_white_refined, connective, rh_white_refined)
    white_svm.difference(pyvista2svmtk(ventricle_surf))
    white = svmtk2pyvista(white_svm)
    pv.save_meshio(output, white)


if __name__ == "__main__":
    main()
