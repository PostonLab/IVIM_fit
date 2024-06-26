


rule moco_scan_bzeros_4d:
    """ run motion-correction (rigid reg to init volume) on the bzeros from 
    a single acquisition (ie before concatenating across acq/dir entities)"""
    input:
        nii_4d=bids(
            root=work,
            suffix="b0s.nii.gz",
            datatype="dwi",
            desc="denoise",
            **input_wildcards["dwi"]
        ),
    output:
        affine_dir=directory(
            bids(
                root=work,
                suffix="transforms",
                desc="moco",
                datatype="dwi",
                **input_wildcards["dwi"]
            )
        ),
        nii_avg3d=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="moco",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
    threads: 8  #doesn't need to be more than the number of bzeros 
    resources:
        mem_mb=32000,
    shadow:
        "minimal"
    container:
        config["singularity"]["prepdwi"]  #-- this rule needs niftyreg, c3d and mrtrix
    group:
        "subj"
    shell:
        """
        dedent () {{
            xargs -L1 echo
        }}
        c4d {input.nii_4d} -slice w 0:-1 -oo dwi_%03d.nii.gz
        if [ ! -e 'dwi_001.nii.gz' ]; then 
            mkdir -p {output.affine_dir}
            cp dwi_000.nii.gz {output.nii_avg3d}
        else 
            parallel --eta --jobs {threads} \\
                reg_aladin -flo dwi_{{1}}.nii.gz -ref dwi_000.nii.gz \\
                    -res warped_{{1}}.nii.gz -aff affine_xfm_ras_{{1}}.txt \\
                    --rigOnly \\
                ::: $(
                    ls dwi_???.nii.gz | tail -n +2 | grep -Po '(?<=dwi_)[0-9]+'
                )
            mkdir -p {output.affine_dir}
            cp affine_xfm_ras_*.txt {output.affine_dir}
            echo -e '1 0 0 0
                     0 1 0 0
                     0 0 1 0
                     0 0 0 1' |
                dedent > {output.affine_dir}/affine_xfm_ras_000.txt
            mrcat dwi_000.nii.gz warped_*.nii.gz merged_4d.nii.gz
            mrmath merged_4d.nii.gz mean {output.nii_avg3d} -axis 3
        fi
        """


rule moco_bzeros_3d:
    """ motion-correct the avgb0 scans from each dwi acquisition"""
    input:
        b0s=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="b0.nii.gz",
                    datatype="dwi",
                    desc="moco",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        affine_dir=directory(
            bids(
                root=work,
                suffix="transforms",
                desc="mocoavgb0",
                datatype="dwi",
                **subj_wildcards
            )
        ),
        nii_4d=bids(
            root=work,
            suffix="b0s.nii.gz",
            desc="mocoavgb0",
            datatype="dwi",
            **subj_wildcards
        ),
        nii_avg3d=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="mocoavgb0",
            datatype="dwi",
            **subj_wildcards
        ),
    threads: 8  #doesn't need to be more than the number of bzeros 
    resources:
        mem_mb=32000,
    shadow:
        "minimal"
    container:
        config["singularity"]["prepdwi"]  #-- this rule needs niftyreg, c3d and mrtrix
    group:
        "subj"
    params:
        flo_indices=lambda wcards, input: (
            " ".join([f"{i}" for i in range(1, len(input.b0s))])
        ),
        flo_imgs=lambda wcards, input: " ".join(input.b0s[1:]),
    shell:
        """
        dedent () {{
            xargs -L1 echo
        }}
        parallel --eta --jobs {threads} --link \\
            reg_aladin -flo {{2}}  -ref {input.b0s[0]} -res warped_{{1}}.nii \\
                -aff affine_xfm_ras_{{1}}.txt --rigOnly \\
            ::: {params.flo_indices} \\
            ::: {params.flo_imgs}

        mkdir -p {output.affine_dir}
        cp affine_xfm_ras_*.txt {output.affine_dir}
        echo -e '1 0 0 0
                 0 1 0 0
                 0 0 1 0
                 0 0 0 1' |
            dedent > {output.affine_dir}/affine_xfm_ras_000.txt
        mrcat {input.b0s[0]} warped_*.nii {output.nii_4d}
        mrmath {output.nii_4d} mean {output.nii_avg3d} -axis 3
        """
