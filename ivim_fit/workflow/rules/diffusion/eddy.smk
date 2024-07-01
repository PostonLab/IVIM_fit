def get_dwi_for_eddy(filename):
    # Replace target_string in the filename with replacement_string
    new_filename = filename.replace("{dir}", "AP")
    new_filename = new_filename.replace("{run}", "1")
    return new_filename


rule run_eddy_correct:
    input:
        nii=get_dwi_for_eddy(rules.apply_topup_jac.output.nii),
    output:
        dwi=bids(
            root=work,
            suffix="dwi.nii.gz",
            datatype="dwi",
            desc="eddy",
            **subj_wildcards
        ),
    threads: 16  #needs to be set to avoid multiple gpus from executing
    resources:
        runtime=360,  #6 hours (this is a conservative estimate, may be shorter)
        mem_mb=32000,
    container:
        config["singularity"]["fsl"]
    group:
        "subj"
    shell:
        "eddy_correct"
        " {input.nii} {output.dwi} 0"


rule cp_eddy_outputs:
    input:
        mask=get_b0_mask(),
    output:
        bids(
            root=root,
            suffix="mask.nii.gz",
            desc="brain",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
    group:
        "subj"
    run:
        for in_file, out_file in zip(input, output):
            shell("cp -v {in_file} {out_file}")
