
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
