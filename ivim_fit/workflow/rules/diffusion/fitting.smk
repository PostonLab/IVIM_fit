def get_bval_for_fitting(wildcards):

    bvals = bids(
        root=work,
        suffix="dwi.bval",
        desc="denoise",
        datatype="dwi",
        **input_wildcards["dwi"]
    )
    print(bvals)
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
    resources:
        mem_mb=32000,  # this is going to be dependent on image size
    group:
        "subj"
    script:
        "../../scripts/fitting/fit.py"


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
