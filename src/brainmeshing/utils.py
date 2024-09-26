import pyvista as pv


def print_edge_info(mesh: pv.UnstructuredGrid):
    edges = mesh.extract_all_edges()
    edge_lengths = edges.compute_cell_sizes(length=True)["Length"]  # type: ignore
    print(f"Min edge length: {edge_lengths.min()}")
    print(f"Max edge length: {edge_lengths.max()}")
    print(f"Mean edge length: {edge_lengths.mean()}")
