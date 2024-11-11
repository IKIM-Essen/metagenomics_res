report_input = list()
if config["host-filtering"]["do-host-filtering"]:
    report_input.append("results/{project}/output/report/host_filtering/")


rule snakemake_report:
    input:
        # 1. Quality control
        report_input,
        "results/{project}/output/report/all/multiqc.html",
        # 2. Species diversity
        "results/{project}/output/report/all/quality_summary/",
        # 3. Assembly results
        "results/{project}/output/report/all/assembly/",
        # 4. Binning results
        "results/{project}/output/report/all/binning/",
        expand(
            "results/{{project}}/output/report/{sample}/bin/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/output/report/{sample}/checkm2/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/output/report/{sample}/dastool/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/output/report/{sample}/taxonomy/",
            sample=get_samples(),
        ),
        # 5. Taxonomic classification
        expand(
            "results/{{project}}/output/report/{sample}/mags/",
            sample=get_samples(),
        ),
        expand(
            "results/{{project}}/output/report/{sample}/{sample}_kaiju.out.html",
            sample=get_samples(),
        ),
    output:
        "results/{project}/output/report/report_{project}.zip",
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
