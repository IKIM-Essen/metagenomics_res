## read QC
rule fastqc:
    input:
        get_trimmed_fastqs,
    output:
        html=temp("results/{project}/qc/fastqc/{sample}_trimmed.html"),
        zip=temp("results/{project}/qc/fastqc/{sample}_trimmed_fastqc.zip"),
    log:
        "logs/{project}/fastqc/{sample}.log",
    wrapper:
        "v1.23.5/bio/fastqc"


rule multiqc:
    input:
        expand(
            [
                "results/{{project}}/qc/fastqc/{sample}_trimmed_fastqc.zip",
                "results/{{project}}/report_prerequisites/qc/{sample}.fastp.json",
            ],
            sample=get_samples(),
        ),
    output:
        report(
            "results/{project}/output/report/all_samples/multiqc.html",
            htmlindex="multiqc.html",
            category="1. Quality of reads",
            labels={" ": "MultiQC report"},
        ),
        "results/{project}/output/report/all_samples/multiqc_data.zip",
    params:
        extra=(
            "--zip-data-dir "
            "--config config/multiqc_config.yaml "
            "--title 'Results for data from {project} project'"
        ),
        use_input_files_only=True,
    log:
        "logs/{project}/multiqc.log",
    wrapper:
        "v3.3.1/bio/multiqc"


rule qc_summary:
    input:
        # move these to report_prerequistes
        jsons=expand(
            "results/{{project}}/report_prerequisites/qc/{sample}.fastp.json",
            sample=get_samples(),
        ),
        human_logs=expand(
            "results/{{project}}/report_prerequisites/qc/filter_human_{sample}.log",
            sample=get_samples(),
        ),
        host_logs=get_host_map_statistics,
    output:
        csv="results/{project}/output/report/all_samples/read_QC.csv",
        vis_csv=temp("results/{project}/output/report/all_samples/read_QC_visual.csv"),
    params:
        other_host=config["host-filtering"]["do-host-filtering"],
        hostname=config["host-filtering"]["host-name"],
    log:
        "logs/{project}/report/qc_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/read_quality.py"


rule qc_summary_report:
    input:
        rules.qc_summary.output.vis_csv,
    output:
        report(
            directory("results/{project}/output/report/all_samples/read_QC/"),
            htmlindex="index.html",
            category="1. Quality of reads",
            labels={" ": "Summary table"},
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="read_QC",
        header="Quality summary of reads based on fastp report and mapping to host genome(s)",
        pattern=config["tabular-config"],
    log:
        "logs/{project}/report/read_qc_rbt_csv.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report {input} --pin-until {params.pin_until} {output} && "
        "(sed -i '{params.pattern} {params.header}</a>' "
        "{output}/indexes/index1.html && "
        "sed -i 's/report.xlsx/{params.name}_report.xlsx/g' {output}/indexes/index1.html) && "
        "mv {output}/report.xlsx {output}/{params.name}_report.xlsx && "
        "cp {params.styles}* {output}/css/ > {log} 2>&1"
