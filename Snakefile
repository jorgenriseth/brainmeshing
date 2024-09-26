configfile: "snakeconfig.yaml"
rule extract_ventricles:
  input:
    f"{config['FS_DIR']}/mri/aseg.mgz"
  output:
    "surfaces/ventricles.stl"
  params:
    presmooth = config["presmooth"]
  shell:
    "python scripts/ventricle_surface.py"
    " -i {input}"
    " -o {output}"
    " --min_radius 2"
    " --initial_smoothing {params.presmooth}"
    " --surface_smoothing 1"
    " --taubin_iter 20"


rule preprocess_white_surface:
  input:
    seg=f"{config['FS_DIR']}/mri/aseg.mgz",
    ventricles="surfaces/ventricles.stl"
  output:
    "surfaces/white.stl"
  params:
    fs_dir = config["FS_DIR"]
  shell:
    "python scripts/white_matter_surface.py"
    " --fs_dir {params.fs_dir}"
    " --seg {input.seg}"
    " --ventricles {input.ventricles}"
    " --output {output}"
