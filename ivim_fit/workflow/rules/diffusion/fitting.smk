def get_bval_for_fitting(wildcards):

    bvals = bids(
        root=work,
        suffix="dwi.bval",
        desc="denoise",
        datatype="dwi",
        **input_wildcards["dwi"]
    )
    # Replace target_string in the filename with replacement_string
    new_filename = bvals.replace("{dir}", "AP")
    new_filename = new_filename.replace("{run}", "1")
    return new_filename


rule fit_ivim:
    input:
        dwi=rules.run_eddy_correct.output.dwi,
        mask=rules.synthstrip_b0_fix_header.output.mask,
        bval=get_bval_for_fitting,
    params:
        out_dir=bids(root=work, datatype="dwi", **subj_wildcards),
    output:
        fitted=bids(root=work, suffix="ivimfit.csv", **subj_wildcards),
    threads: 8
    resources:
        mem_mb=32000,  # this is going to be dependent on image size
        runtime=720,
    group:
        "subj"
    script:
        "../../scripts/fitting/fit.py"


rule average_dwi:
    input:
        dwi_t1_space=rules.resample_dwi_to_t1w.output.dwi,
        mask=rules.reg_dwi_to_t1.output.warped_avgb0,
        t1_seg=rules.segment_t1.output.t1_seg,
        bval=get_bval_for_fitting,
        label_lookup=os.path.join(
            workflow.basedir, config["synthseg_label_lookup"]
        ),
    params:
        out_dir=bids(root=work, datatype="dwi", **subj_wildcards),
    output:
        avg_dwi=bids(root=work, suffix="average_dwi.csv", **subj_wildcards),
    resources:
        mem_mb=32000,  # this is going to be dependent on image size
    group:
        "subj"
    script:
        "../../scripts/fitting/avg_dwi.py"


rule average_n_fit_ivim:
    input:
        avg_dwi_t1_space=rules.average_dwi.output.avg_dwi,
        bval=get_bval_for_fitting,
        label_lookup=os.path.join(
            workflow.basedir, config["synthseg_label_lookup"]
        ),
    params:
        out_dir=bids(root=work, datatype="dwi", **subj_wildcards),
    output:
        avg_n_fitted=bids(
            root=work, suffix="average_n_fitted_ivim.csv", **subj_wildcards
        ),
    resources:
        mem_mb=32000,  # this is going to be dependent on image size
    group:
        "subj"
    script:
        "../../scripts/fitting/avg_n_fit.py"


rule concat_directions:
    input:
        fitted_dwi=rules.fit_ivim.output.fitted,
    params:
        out_dir=bids(root=work, datatype="dwi", **subj_wildcards),
    output:
        concatted=bids(
            root=work, suffix="concat-ivimfit.csv", **subj_wildcards
        ),
    resources:
        mem_mb=32000,  # this is going to be dependent on image size
    group:
        "subj"
    script:
        "../../scripts/fitting/concat_fitted.py"


rule resample_fitted_dwi_to_t1w:
    input:
        ref=rules.reg_dwi_to_t1.output.warped_avgb0,
        concatted=rules.concat_directions.output.concatted,
        xfm_itk=rules.convert_xfm_ras2itk.output.xfm_itk,
        script=os.path.join(workflow.basedir, f"scripts/fitting/transform.sh"),
    params:
        interpolation="Linear",
        out_dir=bids(root=work, datatype="dwi", **subj_wildcards),
    output:
        resampled=bids(
            root=work, suffix="resampled-ivimfit.csv", **subj_wildcards
        ),
    container:
        config["singularity"]["ants"]
    resources:
        mem_mb=32000,  # this is going to be dependent on image size
    group:
        "subj"
    shell:
        "chmod a+x {input.script} && "
        "{input.script} {input.concatted} {input.ref} {input.xfm_itk} "
        "{params.interpolation} {params.out_dir} {output.resampled}"


rule mean_directions:
    input:
        resampled=rules.resample_fitted_dwi_to_t1w.output.resampled,
    params:
        out_dir=bids(root=work, datatype="dwi", **subj_wildcards),
    output:
        mean=bids(root=work, suffix="mean-ivimfit.csv", **subj_wildcards),
    resources:
        mem_mb=32000,  # this is going to be dependent on image size
    group:
        "subj"
    script:
        "../../scripts/fitting/mean_fitted_loop.py"
