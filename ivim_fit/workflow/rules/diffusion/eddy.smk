

rule run_eddy_correct:
    input:
        nii=rules.apply_topup_jac.output.nii,
    output:
        dwi=bids(
            root=work,
            suffix="dwi.nii.gz",
            datatype="dwi",
            desc="eddy",
            **input_wildcards["dwi"]
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
