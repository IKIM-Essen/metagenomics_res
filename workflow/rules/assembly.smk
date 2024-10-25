from pathlib import Path


rule megahit:
    input:
        fastqs=get_filtered_gz_fastqs,
    output:
        contigs=("results/{project}/megahit/{sample}/final.contigs.fa"),
        outdir=directory("results/{project}/megahit/{sample}/"),
    params:
        threshold=get_contig_length_threshold(),
    threads: 64
    log:
        "results/{project}/report_prerequisites/assembly/{sample}_megahit.log",
    conda:
        "../envs/megahit.yaml"
    shell:
        #"--cleaning-rounds 10 --merge-level 20,0.90 "
        "megahit -1 {input.fastqs[0]} -2 {input.fastqs[1]} "
        "--min-contig-len {params.threshold} -t {threads} "
        "--out-dir {output.outdir} -f > {log} 2>&1"


rule remove_megahit_intermediates:
    input:
        contigs=get_assembly,
    output:
        touch("logs/{project}/assembly/{sample}_intermediate_removal.done"),
    params:
        outdir=lambda wildcards, input: Path(input.contigs).parent,
    log:
        "logs/{project}/assembly/{sample}_intermediate_removal.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(find {params.outdir}/ -mindepth 1 -type d -exec rm -rf {{}} +) > {log} 2>&1"


rule map_to_assembly:
    input:
        contigs=get_assembly,
        fastqs=get_filtered_gz_fastqs,
    output:
        bam=temp("results/{project}/report_prerequisites/assembly/{sample}_reads_mapped.bam"),
    threads: 64
    log:
        "logs/{project}/assembly/{sample}_mapping_reads.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "(minimap2 -a -xsr -t {threads} {input.contigs} {input.fastqs} | "
        "samtools view -bh | "
        "samtools sort --threads {threads} -o {output.bam}) > {log} 2>&1"


rule index_assembly_alignment:
    input:
        rules.map_to_assembly.output.bam,
    output:
        bai=temp("results/{project}/report_prerequisites/assembly/{sample}_reads_mapped.bam.bai"),
    threads: 20
    log:
        "logs/{project}/assembly/{sample}_mapping_reads_index.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools index {input} > {log} 2>&1"


rule reads_mapped_assembly:
    input:
        bam=rules.map_to_assembly.output.bam,
        bai=rules.index_assembly_alignment.output.bai,
    output:
        bai="results/{project}/report_prerequisites/assembly/{sample}_reads_mapped.txt",
    threads: 20
    log:
        "logs/{project}/assembly/{sample}_mapping_reads.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools view -c -F 4 --threads {threads} "
        "-o {output.bai} {input.bam} > {log} 2>&1"


rule gzip_assembly:
    input:
        contigs=get_assembly,
    output:
        "results/{project}/output/fastas/{sample}/{sample}_final.contigs.fa.gz",
    threads: 64
    log:
        "logs/{project}/assembly/{sample}_gzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "pigz -c {input.contigs} > {output} 2> {log}"


rule assembly_summary:
    input:
        qc_csv=rules.qc_summary.output.csv,
        asbl=expand(
            "results/{{project}}/report_prerequisites/assembly/{sample}_megahit.log",
            sample=get_samples(),
        ),
        mapped=expand(
            "results/{{project}}/report_prerequisites/assembly/{sample}_reads_mapped.txt",
            sample=get_samples(),
        ),
    output:
        csv="results/{project}/output/report/all/assembly_summary.csv",
        vis_csv=temp("results/{project}/output/report/all/assembly_summary_visual.csv"),
    log:
        "logs/{project}/report/assembly_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assembly_summary.py"


use rule qc_summary_report as assembly_report with:
    input:
        rules.assembly_summary.output.vis_csv,
    output:
        report(
            directory("results/{project}/output/report/all/assembly/"),
            htmlindex="index.html",
            category="3. Assembly results",
            labels={
                "sample": "all samples",
            },
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="assembly_summary",
        header="Assembly summary",
        pattern=config["tablular-config"],
    log:
        "logs/{project}/report/assembly_rbt_csv.log",
