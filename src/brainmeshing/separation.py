from pathlib import Path

import pyvista as pv
import SVMTK as svmtk
import numpy as np

from brainmeshing.surfaces import pyvista2svmtk, svmtk2pyvista


def separate_white_and_gray_surfaces(
    pial_in: Path, white_in: Path, pial_out: Path, white_out: Path
):
    pial = svmtk.Surface(pial_in)
    white = svmtk.Surface(white_in)

    pial_grid = svmtk2pyvista(pial)
    white_grid = svmtk2pyvista(white)
    pial_grid.compute_normals(
        inplace=True, auto_orient_normals=True, consistent_normals=True
    )
    white_grid.compute_normals(
        inplace=True, auto_orient_normals=True, consistent_normals=True
    )
    pial_grid.compute_implicit_distance(white_grid, inplace=True)
    white_grid.compute_implicit_distance(pial_grid, inplace=True)

    pial_degenerate = pial_grid.extract_points(
        np.where(pial_grid.point_data["implicit_distance"] < 0)
    )
    white_degenerate = white_grid.extract_points(
        np.where(white_grid.point_data["implicit_distance"] > 0)
    )

    _, closest_points = (pial_degenerate + white_degenerate).find_closest_cell(
        pial_grid.points, return_closest_point=True
    )
    d_exact = np.linalg.norm(pial_grid.points - closest_points, axis=1)

    pial_degenerate_neighbourhood = pial_grid.extract_points(np.where(d_exact < 0.5))
    all_regions = pial_degenerate_neighbourhood.connectivity("all")
    region_ids = np.unique(all_regions["RegionId"])
    regions = [
        all_regions.extract_points(
            np.where(all_regions.point_data["RegionId"] == label)
        )
        for label in region_ids
    ]

    new_pial_svm = pyvista2svmtk(pial_grid)
    for region in regions:
        patch = pyvista2svmtk(region).convex_hull()
        patch.isotropic_remeshing(1.0, 1, True)
        new_pial_svm.union(patch)
        new_pial_svm.fill_holes()
    try:
        print("Pia num self intersections", new_pial_svm.num_self_intersections())
        print(
            "Repair self intersections, (success, num_left): ",
            new_pial_svm.repair_self_intersections(),
        )
    except svmtk.PreconditionError as e:
        print(e)
    new_pial = svmtk2pyvista(new_pial_svm)
    pv.save_meshio(pial_out, new_pial)

    # White surface fixing
    all_idcs = white_grid.point_data["vtkOriginalPointIds"][
        white_degenerate.point_data["vtkOriginalPointIds"]
    ]
    neighbor_idcs = np.unique(
        np.concatenate(
            [
                sum(white_grid.point_neighbors_levels(idx, 1), start=[])
                for idx in all_idcs
            ]
        )
    )
    white_degenerate_neighbourhood = white_grid.extract_points(neighbor_idcs)
    all_regions = white_degenerate_neighbourhood.connectivity("all")
    region_ids = np.unique(all_regions.point_data["RegionId"])
    regions = [
        all_regions.extract_points(
            np.where(all_regions.point_data["RegionId"] == label)
        )
        for label in region_ids
    ]
    holey_white, ridx = white_grid.remove_points(
        white_degenerate_neighbourhood.point_data["vtkOriginalPointIds"]
    )
    new_white_svm = pyvista2svmtk(holey_white)
    new_white_svm.fill_holes()
    try:
        print("White num self intersections", new_white_svm.num_self_intersections())
        print(
            "Repair self intersections, (success, num_left): ",
            new_white_svm.repair_self_intersections(),
        )
    except svmtk.PreconditionError as e:
        print(e)
    new_white = svmtk2pyvista(new_white_svm)
    pv.save_meshio(white_out, new_white)
