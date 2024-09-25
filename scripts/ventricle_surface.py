from pathlib import Path

import click
import numpy as np
import pyvista as pv

from gonzo.simple_mri import load_mri
from brainmeshing.ventricles import extract_ventricle_surface


@click.command(name="ventricle-surface")
@click.option("-i", "--input", type=Path, required=True)
@click.option("-o", "--output", type=Path, required=True)
@click.option("--min_radius", type=int, default=2)
@click.option("--initial_smoothing", type=float, default=0)
@click.option("--surface_smoothing", type=float, default=0)
@click.option("--taubin_iter", type=int, default=0)
@click.option("--voxelized", type=bool, is_flag=True)
def main(input: Path, output: Path, **kwargs):
    Path(output).parent.mkdir(exist_ok=True)
    seg_mri = load_mri(input, dtype=np.int16)
    surf = extract_ventricle_surface(seg_mri, **kwargs)
    pv.save_meshio(output, surf)


if __name__ == "__main__":
    main()
