report_input = list()
if config["host-filtering"]["do-host-filtering"]:
    report_input.append("results/{project}/output/report/host_filtering/")


rule snakemake_report:
    input:
        # 1. Read quality
        report_input,
        "results/{project}/output/report/all_samples/multiqc.html",
        "results/{project}/output/report/all_samples/read_QC/",
        "results/{project}/output/report/all_samples/host_contamination.html",
        # 2. read taxonomy
        expand(
            "results/{{project}}/output/report/{sample}/{sample}_read_taxonomy_krona.html",
            sample=get_samples(),
        ),
        # 3. Assembly results
        "results/{project}/output/report/all_samples/assembly_QCy/",
        # 4. Binning results
        "results/{project}/output/report/all_samples/binning_QC/",
        expand(
            "results/{{project}}/output/report/{sample}/bin/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/output/report/{sample}/taxonomy/",
            sample=get_samples(),
        ),
    output:
        "results/{{project}}/output/report/{rundate}_snakemake_report_{{project}}.zip".format(rundate=get_rundate()),
    params:
        style="resources/report/custom-stylesheet.css",
    #    for_testing=get_if_testing("--snakefile ../workflow/Snakefile"),
    log:
        "logs/{project}/snakemake-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake --nolock --report {output} --report-stylesheet {params.style} "
        "> {log} 2>&1"
        #"{params.for_testing} "
