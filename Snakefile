configfile: "snakeconfig.yaml"
wildcard_constraints:
  resolution = r"\d+"

rule extract_ventricles:
  input:
    f"{config['FS_DIR']}/mri/aseg.mgz"
  output:
    "surfaces/ventricles.stl"
  params:
    presmooth = config["presmooth"]
  shell:
    "brainmeshing ventricle-surf"
    " -i {input}"
    " -o {output}"
    " --min_radius 2"
    " --initial_smoothing {params.presmooth}"
    " --surface_smoothing 1"
    " --taubin_iter 20"

rule separate_white_and_gray:
  input:
    lh_pial=f"{config['FS_DIR']}/surf/lh.pial",
    lh_white=f"{config['FS_DIR']}/surf/lh.white",
    rh_pial=f"{config['FS_DIR']}/surf/rh.pial",
    rh_white=f"{config['FS_DIR']}/surf/rh.white",
  output:
    lh_pial=f"surfaces/lh_pial.stl",
    lh_white=f"surfaces/lh_white.stl",
    rh_pial=f"surfaces/rh_pial.stl",
    rh_white=f"surfaces/rh_white.stl",
  params:
    fs_dir = config["FS_DIR"]
  shell:
    "brainmeshing separate-surfaces"
    " --fs_dir {params.fs_dir}"
    " --outputdir $(dirname {output[0]})"

rule preprocess_white_surface:
  input:
    seg=f"{config['FS_DIR']}/mri/aseg.mgz",
    lh_white=f"surfaces/lh_white.stl",
    rh_white=f"surfaces/rh_white.stl",
    ventricles="surfaces/ventricles.stl"
  output:
    "surfaces/white.stl"
  shell:
    "brainmeshing wm-surfaces"
    " --inputdir $(dirname {input.lh_white})"
    " --seg {input.seg}"
    " --ventricles {input.ventricles}"
    " --output {output}"
    " --tmpdir surfaces/"


rule preprocess_gray_surfaces:
  input:
    seg=f"{config['FS_DIR']}/mri/aseg.mgz",
    ventricles="surfaces/ventricles.stl",
    lh_pial=f"surfaces/lh_pial.stl",
    rh_pial=f"surfaces/rh_pial.stl",
  output:
    lh="surfaces/lh_pial_novent.stl",
    rh="surfaces/rh_pial_novent.stl",
    subcort="surfaces/subcortical_gm.stl"
  shell:
    "brainmeshing gm-surfaces"
    " --inputdir $(dirname {input.lh_pial})"
    " --seg {input.seg}"
    " --ventricles {input.ventricles}"
    " --outputdir $(dirname {output[0]})"

rule generate_mesh:
  input:
    expand(
      "surfaces/{surf}.stl",
      surf=["lh_pial_novent", "rh_pial_novent", "subcortical_gm", "white"]
    )
  output:
    hdf="meshes/mesh{resolution}.hdf",
    xdmf="meshes/mesh{resolution}_xdmfs/subdomains.xdmf"
  shell:
    "brainmeshing meshgen"
    " --surfacedir $(dirname {input[0]})"
    " --resolution {wildcards.resolution}"
    " --output {output.hdf}"


rule mesh_segmentation:
  input:
    seg=f"{config['FS_DIR']}/mri/aparc+aseg.mgz",
    mesh="meshes/mesh{resolution}.hdf",
  output:
    "meshes/mesh{resolution}_aparc.hdf"
  params: config["FS_DIR"]
  shell:
    "brainmeshing subdomains"
    " {input} {output}"
