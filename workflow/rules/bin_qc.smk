## bin QC
rule checkm2_DB_download:
    output:
        dbfile=get_checkm2_db(),  #"{}/{}".format(config["data-handling"]["resources"], config["checkm2"]),
    params:
        direct=lambda wildcards, output: Path(output.dbfile).parent.parent,
    log:
        "logs/checkm2_DB_download.log",
    conda:
        "../envs/checkm2.yaml"
    shell:
        "checkm2 database --download --path {params.direct} > {log} 2>&1"


rule checkm2_run:
    input:
        bins="results/{project}/output/fastas/{sample}/bins/",
        dbfile=get_checkm2_db(),
    output:
        outdir=temp(directory("results/{project}/qc/checkm2/{sample}/")),
        stats="results/{project}/output/report/{sample}/checkm2_quality_report.tsv",
    params:
        outname="quality_report.tsv",
    log:
        "logs/{project}/checkm2/{sample}.log",
    threads: 24
    conda:
        "../envs/checkm2.yaml"
    shell:
        "(checkm2 predict -x fa.gz --threads {threads} --force "
        "--input {input.bins}/ --output-directory {output.outdir}/ && "
        "cp {output.outdir}/{params.outname} {output.stats}) > {log} 2>&1"


rule bin_summary_sample:
    input:
        tool="results/{project}/output/report/{sample}/{sample}_DASTool_summary.tsv",
        checkm="results/{project}/output/report/{sample}/checkm2_quality_report.tsv",
        gtdb="results/{project}/output/classification/bins/{sample}/{sample}.summary.tsv",
    output:
        csv_bins="results/{project}/output/report/{sample}/{sample}_bin_quality.csv",
        csv_tax="results/{project}/output/report/{sample}/{sample}_bin_taxonomy.csv",
    params:
        max_cont=config["MAG-criteria"]["max-contamination"],
        min_comp=config["MAG-criteria"]["min-completeness"],
    log:
        "logs/{project}/bin_quality/{sample}.log",
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/bin_summary_sample.py"


use rule qc_summary_report as bin_sample_report with:
    input:
        "results/{project}/output/report/{sample}/{sample}_bin_quality.csv",
    output:
        report(
            directory("results/{project}/output/report/{sample}/bin/"),
            htmlindex="index.html",
            category="4. Binning",
            subcategory="4.2 Quality of bins",
            labels={"sample": "{sample}"},
        ),
    params:
        pin_until="bin",
        styles="resources/report/tables/",
        name="{sample}_bin_quality",
        header="Quality table of all bins for sample {sample}",
        pattern=config["tabular-config"],
    log:
        "logs/{project}/report/{sample}/bin_rbt_csv.log",


use rule qc_summary_report as taxonomy_report with:
    input:
        "results/{project}/output/report/{sample}/{sample}_bin_taxonomy.csv",
    output:
        report(
            directory("results/{project}/output/report/{sample}/taxonomy/"),
            htmlindex="index.html",
            category="4. Binning",
            subcategory="4.3 Taxonomic classification",
            labels={"sample": "{sample}"},
        ),
    params:
        pin_until="bin",
        styles="resources/report/tables/",
        name="{sample}_bin_taxonomy",
        header="Taxonomy summary for sample {sample}",
        pattern=config["tabular-config"],
    log:
        "logs/{project}/report/{sample}/taxonomy_rbt_csv.log",


rule bin_summary_all:
    input:
        csv_bins=expand(
            "results/{{project}}/output/report/{sample}/{sample}_bin_quality.csv",
            sample=get_samples(),
        ),
    output:
        "results/{project}/output/report/all_samples/binning_quality.csv",
    log:
        "logs/{project}/bin_quality/all.log",
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/bin_summary_all.py"


use rule qc_summary_report as bin_all_report with:
    input:
        "results/{project}/output/report/all_samples/binning_quality.csv",
    output:
        report(
            directory("results/{project}/output/report/all_samples/binning_QC/"),
            htmlindex="index.html",
            category="4. Binning",
            subcategory="4.1 Binning overview",
            labels={" ": "Summary table"},
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="binning_quality",
        header="Bin quality summary for all samples",
        pattern=config["tabular-config"],
    log:
        "logs/{project}/report/all_bin_rbt_csv.log",
