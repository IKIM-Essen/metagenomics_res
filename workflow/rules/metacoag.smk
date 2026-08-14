rule coverm_metacoag:
    input:
        bact_reads=get_filtered_gz_fastqs,
        contigs=get_assembly,
    output:
        temp("results/{project}/binning_prep/{sample}/abundance.tsv"),
    log:
        "logs/{project}/coverm/{sample}.log",
    conda:
        "../envs/coverm.yaml"
    threads: 30
    shell:
        "coverm contig -1 {input.bact_reads[0]} -2 {input.bact_reads[1]} "
        "-r {input.contigs} -o {output} -t {threads} > {log} 2>&1"


rule edit_abundance_file:
    input:
        "results/{project}/binning_prep/{sample}/abundance.tsv",
    output:
        temp("results/{project}/binning_prep/{sample}/abundance_metacoag.tsv"),
    log:
        "logs/{project}/coverm/{sample}.log",
    conda:
        "../envs/unix.yaml"
    threads: 1
    shell:
        "cp {input} {output} && sed -i '1d' {output} > {log} 2>&1"


rule fastg_assembly_tree:
    input:
        contigs=get_assembly,
    output:
        fastg=(
            "results/{project}/output/evidence/assembly/{sample}_assembly_tree.fastg"
            if retain_assembly_evidence()
            else temp("results/{project}/binning_prep/{sample}/assembly_tree.fastg")
        ),
    log:
        "logs/{project}/fastg_assembly_tree/{sample}.log",
    conda:
        "../envs/megahit.yaml"
    threads: 2
    script:
        "../scripts/fastg_assembly_tree.py"


rule fastg2gfa:
    input:
        rules.fastg_assembly_tree.output.fastg,
    output:
        gfa=(
            "results/{project}/output/evidence/assembly/{sample}_assembly_tree.gfa"
            if retain_assembly_evidence()
            else temp("results/{project}/binning_prep/{sample}/assembly_tree.gfa")
        ),
    log:
        "logs/{project}/fastg2gfa/{sample}.log",
    conda:
        "../envs/metacoag.yaml"
    threads: 2
    params:
        fastg2gfa_program="workflow/scripts/fastg2gfa",
    shell:
        "{params.fastg2gfa_program} {input} > {output} 2> {log}"


rule metacoag_run:
    input:
        contigs=get_assembly,
        gfa=rules.fastg2gfa.output.gfa,
        abd="results/{project}/binning_prep/{sample}/abundance_metacoag.tsv",
        #assembly folder needs to be there
        folder=rules.megahit.output.outdir,
    output:
        out_tsv=temp("results/{project}/binning/metacoag/{sample}/contig_to_bin.tsv"),
    log:
        "logs/{project}/metacoag/{sample}.log",
    conda:
        "../envs/metacoag.yaml"
    threads: 4
    params:
        outdir=lambda wildcards, output: Path(output.out_tsv).parent,
    shell:
        "metacoag --assembler megahit --graph {input.gfa} "
        "--contigs {input.contigs} --abundance {input.abd} "
        "--output {params.outdir} --min_length 300 "
        "--bin_mg_threshold 0.2 --min_bin_size 100000 "
        "--nthreads {threads} > {log} 2>&1"
