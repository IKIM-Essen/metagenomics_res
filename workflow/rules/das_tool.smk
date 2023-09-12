rule postprocess_vamb:
    input:
        "results/{project}/vamb/{sample}/vamb_res/clusters.tsv",
    output:
        "results/{project}/contig2bins/{sample}/vamb_contig2bin.tsv",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_vamb.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(awk '{{print $2 \"\t\" $1}}' {input} > {output} && "
        "sed -i 's/S1C//g' {output}) > {log} 2>&1"


rule postprocess_metabat:
    input:
        "results/{project}/metabat2/{sample}/",
    output:
        "results/{project}/contig2bins/{sample}/metabat2_contig2bin.tsv",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_metabat2.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/binner_postprocess.py"


use rule postprocess_metabat as postprocess_rosella with:
    input:
        "results/{project}/rosella/{sample}/",
    output:
        "results/{project}/contig2bins/{sample}/rosella_contig2bin.tsv",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_rosella.log",


rule postprocess_metacoag:
    input:
        "results/{project}/metacoag/{sample}/contig_to_bin.tsv",
    output:
        "results/{project}/contig2bins/{sample}/metacoag_contig2bin.tsv",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_metacoag.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(awk '{{print $1 \"\t\" $NF}}' {input} | "
        "sed 's/len=[0-9]*,//g' > {output}) > {log} 2>&1"


## move metabinner output to contig2bin folder
rule postprocess_metabinner:
    input:
        "results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
    output:
        "results/{project}/contig2bins/{sample}/metabinner_contig2bin.tsv",
    log:
        "logs/{project}/contig2bins/{sample}/postprocess_metabinner.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cp {input} {output} > {log} 2>&1"


## tests if binner created bins and writes file for DAS Tool
rule binner_control:
    input:
        get_all_contig2bin_files,
    output:
        "results/{project}/das_tool/{sample}/binner_control.csv",
    params:
        get_binners(),
    log:
        "logs/{project}/das_tool/{sample}/binner_control.log",
    conda:
        "../envs/unix.yaml"
    script:
        "../scripts/binner_control.py"


rule dastool_run:
    input:
        binc="results/{project}/das_tool/{sample}/binner_control.csv",
        contigs="results/{project}/megahit/{sample}/final.contigs.fa",
    output:
        outfile="results/{project}/das_tool/{sample}/{sample}_DASTool_summary.tsv",
    params:
        path_bin_list=get_paths_binner,
        outdir=lambda wildcards, output: Path(output.outfile).parent,
    threads: get_DAS_Tool_threads()
    log:
        "logs/{project}/das_tool/{sample}/das_tool_run.log",
    conda:
        "../envs/das_tool.yaml"
    shell:
        "DAS_Tool -t {threads} "
        "--debug --write_bins "
        "-i {params.path_bin_list[0]} "
        "-l {params.path_bin_list[1]} "
        "-c {input.contigs} "
        "-o {params.outdir}/{wildcards.sample} "
        "> {log} 2>&1 "

