
wildcard_constraints:
    shell="[0-9]+",


rule import_dwi:
    input:
        dwi=re.sub(".nii.gz", ".{ext}", input_path["dwi"]),
    output:
        dwi=bids(
            root=work,
            suffix="dwi.{ext,nii.gz|bval|bvec|json}",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
    group:
        "subj"
    shell:
        "cp {input.dwi} {output.dwi}"


rule dwidenoise:
    input:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                datatype="dwi",
                **input_wildcards["dwi"],
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    output:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                desc="denoise",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="denoise.log", **input_wildcards["dwi"]),
    group:
        "subj"
    shell:
        "dwidenoise {input[0]} {output[0]} 2> {log} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"


def get_concat_or_cp_cmd(wildcards, input, output):
    """Concatenate (if multiple inputs) or copy"""
    if len(input) > 1:
        cmd = f"mrcat {input} {output}"
    elif len(input) == 1:
        cmd = f"cp {input} {output}"
    else:
        # no inputs
        cmd = None
    return cmd


rule concat_dwi:
    input:
        dwi_niis=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwi.nii.gz",
                    desc="denoise",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        dwi_concat=bids(
            root=work,
            suffix="dwi.nii.gz",
            desc="denoise",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="concat_dwi.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} 2> {log}"


rule concat_runs_bvec:
    input:
        bv_files=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwi.bvec",
                    desc="{{desc}}",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        out_fname=bids(
            root=work,
            suffix="dwi.bvec",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/concat_bv.py"


rule concat_runs_bval:
    input:
        bv_files=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwi.bval",
                    desc="{{desc}}",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        out_fname=bids(
            root=work,
            suffix="dwi.bval",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/concat_bv.py"


# Combine multiple json from multiple scans (currently only copying first)
rule concat_runs_json:
    input:
        jsons=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwi.json",
                    datatype="dwi",
                    desc="denoise",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        json=bids(
            root=work,
            suffix="dwi.json",
            desc="denoise",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input.jsons[0]} {output.json}"


rule get_shells_from_bvals:
    input:
        bval="{dwi_prefix}.bval",
    output:
        json="{dwi_prefix}.shells.json",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_shells_from_bvals.py"


# Write 4D dwi_file with average shells
rule get_shell_avgs:
    input:
        dwi="{dwi_prefix}.nii.gz",
        shells="{dwi_prefix}.shells.json",
    output:
        avgshells="{dwi_prefix}.avgshells.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_shell_avgs.py"


# Extract individual shell (e.g. B0)
rule get_shell_avg:
    input:
        dwi="{dwi_prefix}_dwi.nii.gz",
        shells="{dwi_prefix}_dwi.shells.json",
    params:
        bval="{shell}",
    output:
        avgshell="{dwi_prefix}_b{shell}.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_shell_avg.py"


# Extract vols from particular shell (e.g. B0)
rule get_shell_vols:
    input:
        dwi="{dwi_prefix}_dwi.nii.gz",
        shells="{dwi_prefix}_dwi.shells.json",
    params:
        bval="{shell}",
    output:
        shell_vols="{dwi_prefix}_b{shell}s.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_shell_vols.py"


# now have nii with just the b0's, want to create the topup phase-encoding text files for each one:
rule get_phase_encode_txt:
    input:
        bzero_nii=bids(
            root=work,
            suffix="b0.nii.gz",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
        json=bids(
            root=work,
            suffix="dwi.json",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
    output:
        phenc_txt=bids(
            root=work,
            suffix="phenc.txt",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_phase_encode_txt.py"


rule concat_phase_encode_txt:
    input:
        phenc_txts=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="phenc.txt",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        phenc_concat=bids(
            root=work, suffix="phenc.txt", datatype="dwi", **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cat {input} > {output}"


rule concat_bzeros:
    input:
        bzero_niis=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="b0.nii.gz",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        bzero_concat=bids(
            root=work,
            suffix="concatb0.nii.gz",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="concat_bzeros.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} 2> {log}"


def get_b0_mask():
    # Method options
    methods = {
        "b0_BET": "bet_from-b0",
        "b0_SyN": f"b0SyN_from-{config['template']}",
        "b0_synthstrip": "synthstrip_from-topupb0",
    }

    # Get BIDS name of file
    return bids(
        root=work,
        suffix="mask.nii.gz",
        desc="brain",
        method=methods.get(config["masking_method"]),
        datatype="dwi",
        **subj_wildcards,
    )
