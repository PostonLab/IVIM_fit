rule run_synthSR:
    input:
        nii="{prefix}.nii.gz",
    output:
        out_nii="{prefix}SynthSR.nii.gz",
    threads: 8
    group:
        "subj"
    container:
        config["singularity"]["synthsr"]
    shadow:
        "minimal"
    shell:
        "python /SynthSR/scripts/predict_command_line.py "
        "--cpu --threads {threads} "
        "{input} {output}"
