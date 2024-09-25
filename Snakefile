configfile: "snakeconfig.yaml"
rule extract_ventricles:
  input:
    f"{config['FS_DIR']}/mri/aseg.mgz"
  output:
    "surfaces/ventricles.stl"
  shell:
    "python scripts/ventricle_surface.py"
    " -i {input}"
    " -o {output}"
    " --min_radius 2"
    " --initial_smoothing 1"
    " --surface_smoothing 1"
    " --taubin_iter 20"
