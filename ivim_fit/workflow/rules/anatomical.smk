# Reading T1w images
rule import_t1:
    input:
        nii=lambda wildcards: expand(
            input_path["T1w"],
            zip,
            **filter_list(input_zip_lists["T1w"], wildcards),
        )[0],
    output:
        nii=bids(
            root=work,
            datatype="anat",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    shell:
        "cp {input.nii} {output.nii}"


rule synthstrip_t1:
    input:
        t1=rules.import_t1.output.nii,
    output:
        mask=temp(
            bids(
                root=work,
                datatype="anat",
                desc="nofixhdrbrain",
                suffix="mask.nii.gz",
                **subj_wildcards,
            )
        ),
    group:
        "subj"
    container:
        config["singularity"]["synthstrip"]
    threads: 8
    shadow:
        "minimal"
    shell:
        "python3 /freesurfer/mri_synthstrip -i {input.t1} -m {output.mask}"


rule fixheader_synthstrip:
    input:
        t1=rules.import_t1.output.nii,
        mask=rules.synthstrip_t1.output.mask,
    output:
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            **subj_wildcards,
        ),
    group:
        "subj"
    container:
        config["singularity"]["itksnap"]
    shell:
        "c3d {input.t1} {input.mask} -copy-transform -o {output.mask}"


rule n4_t1_withmask:
    input:
        t1=rules.import_t1.output.nii,
        mask=rules.fixheader_synthstrip.output.mask,
    output:
        t1=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T1w.nii.gz",
            **subj_wildcards,
        ),
    threads: 8
    container:
        config["singularity"]["ants"]
    group:
        "subj"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -d 3 -i {input.t1} -x {input.mask} -o {output}"


rule segment_t1:
    input:
        t1w=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix="T1wSynthSR.nii.gz",
            **subj_wildcards
        ),
    output:
        t1_seg=bids(
            root=work,
            datatype="anat",
            space="orig",
            label="wholebrain",
            suffix="probseg.nii.gz",
            **subj_wildcards,
        ),
    params:
        out_dir=bids(root=work, **subj_wildcards),
    resources:
        mem_mb=32000,
        runtime=60,
        disk_mb=8000,
    threads: 8
    container:
        config["singularity"]["freesurfer"]
    group:
        "subj"
    shell:
        "mri_synthseg --i {input.t1w} --o {output.t1_seg} --parc --threads {threads} --cpu"
