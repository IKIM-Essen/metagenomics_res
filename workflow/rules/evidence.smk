"""Evidence exports used by downstream strain/ARG-host inference."""

if culture_free_enabled():

    rule export_culture_free_sample_design:
        output:
            tsv="results/{project}/output/evidence/sample_design.tsv",
            json="results/{project}/output/evidence/sample_design.audit.json",
        log:
            "logs/{project}/evidence/sample_design.log",
        conda:
            "../envs/python.yaml"
        params:
            records=culture_free_sample_records(),
            require_matched_negative=(
                control_aware_enabled()
                and bool(
                    culture_free_config()
                    .get("control-aware", {})
                    .get("require-matched-negative", True)
                )
            ),
            enforce_complete=control_aware_enabled(),
        script:
            "../scripts/export_sample_design.py"

    rule audit_filtered_read_pairs:
        input:
            fastqs=get_filtered_gz_fastqs,
        output:
            json="results/{project}/output/evidence/reads/{sample}.paired_fastq.audit.json",
        log:
            "logs/{project}/evidence/{sample}.paired_fastq.log",
        conda:
            "../envs/python.yaml"
        params:
            sketch_size=culture_free_config()
            .get("read-audit", {})
            .get("sketch-size", 128),
        shell:
            "python workflow/scripts/audit_paired_fastq.py "
            "--r1 {input.fastqs[0]} --r2 {input.fastqs[1]} --output {output.json} "
            "--sketch-size {params.sketch_size} > {log} 2>&1"


if retain_assembly_evidence():

    rule assembly_contig_coverage_evidence:
        input:
            bam=rules.map_to_assembly.output.bam,
            bai=rules.index_assembly_alignment.output.bai,
        output:
            tsv="results/{project}/output/evidence/assembly/{sample}.contig_coverage.tsv",
        log:
            "logs/{project}/evidence/{sample}.contig_coverage.log",
        conda:
            "../envs/minimap2.yaml"
        threads: 4
        params:
            mapq=culture_free_config()["assembly-evidence"]["minimum-mapping-quality"],
            baseq=culture_free_config()["assembly-evidence"]["minimum-base-quality"],
        shell:
            "samtools coverage --min-MQ {params.mapq} --min-BQ {params.baseq} "
            "-o {output.tsv} {input.bam} > {log} 2>&1"

    rule assembly_alignment_statistics_evidence:
        input:
            bam=rules.map_to_assembly.output.bam,
            bai=rules.index_assembly_alignment.output.bai,
        output:
            flagstat="results/{project}/output/evidence/assembly/{sample}.flagstat.tsv",
            stats="results/{project}/output/evidence/assembly/{sample}.samtools_stats.tsv",
            idxstats="results/{project}/output/evidence/assembly/{sample}.idxstats.tsv",
        log:
            "logs/{project}/evidence/{sample}.alignment_stats.log",
        conda:
            "../envs/minimap2.yaml"
        threads: 8
        shell:
            "((samtools flagstat --threads {threads} {input.bam} > {output.flagstat}) && "
            "(samtools stats --threads {threads} {input.bam} > {output.stats}) && "
            "(samtools idxstats {input.bam} > {output.idxstats})) > {log} 2>&1"

    rule assembly_allele_count_evidence:
        input:
            contigs=get_assembly,
            bam=rules.map_to_assembly.output.bam,
            bai=rules.index_assembly_alignment.output.bai,
        output:
            tsv="results/{project}/output/evidence/assembly/{sample}.allele_counts.tsv.gz",
        log:
            "logs/{project}/evidence/{sample}.allele_counts.log",
        conda:
            "../envs/minimap2.yaml"
        threads: 4
        params:
            mapq=culture_free_config()["assembly-evidence"]["minimum-mapping-quality"],
            baseq=culture_free_config()["assembly-evidence"]["minimum-base-quality"],
            max_depth=culture_free_config()["assembly-evidence"]["maximum-pileup-depth"],
        shell:
            "set -o pipefail; samtools mpileup -q {params.mapq} -Q {params.baseq} "
            "-d {params.max_depth} "
            "-f {input.contigs} {input.bam} 2> {log} | "
            "python workflow/scripts/pileup_to_allele_counts.py | gzip -c > {output.tsv}"

    rule assembly_paired_link_evidence:
        input:
            bam=rules.map_to_assembly.output.bam,
            bai=rules.index_assembly_alignment.output.bai,
        output:
            pairs="results/{project}/output/evidence/assembly/{sample}.paired_links.tsv.gz",
            summary="results/{project}/output/evidence/assembly/{sample}.paired_link_summary.tsv",
        log:
            "logs/{project}/evidence/{sample}.paired_links.log",
        conda:
            "../envs/minimap2.yaml"
        threads: 4
        params:
            mapq=culture_free_config()["assembly-evidence"]["minimum-mapping-quality"],
        shell:
            "set -o pipefail; samtools view --threads {threads} -f 65 -F 3852 "
            "-q {params.mapq} {input.bam} "
            "2> {log} | python workflow/scripts/export_paired_links.py "
            "--pairs {output.pairs} --summary {output.summary}"


if retain_assembly_evidence() and direct_unicard_enabled():

    rule lock_culture_free_sample_evidence:
        input:
            fastqs=get_filtered_gz_fastqs,
            read_audit=rules.audit_filtered_read_pairs.output.json,
            sample_design=rules.export_culture_free_sample_design.output.tsv,
            sample_design_audit=rules.export_culture_free_sample_design.output.json,
            fastp="results/{project}/output/report/prerequisites/qc/{sample}.fastp.json",
            contigs=get_assembly,
            bam=rules.map_to_assembly.output.bam,
            bai=rules.index_assembly_alignment.output.bai,
            fastg=rules.fastg_assembly_tree.output.fastg,
            gfa=rules.fastg2gfa.output.gfa,
            coverage=rules.assembly_contig_coverage_evidence.output.tsv,
            flagstat=rules.assembly_alignment_statistics_evidence.output.flagstat,
            stats=rules.assembly_alignment_statistics_evidence.output.stats,
            idxstats=rules.assembly_alignment_statistics_evidence.output.idxstats,
            allele_counts=rules.assembly_allele_count_evidence.output.tsv,
            paired_links=rules.assembly_paired_link_evidence.output.pairs,
            paired_link_summary=rules.assembly_paired_link_evidence.output.summary,
            direct_unicard=rules.direct_uniCARD_reads.output,
            assembly_unicard=rules.rich_uniCARD_assembly_proteins.output.tsv,
            amr_lock=rules.lock_culture_free_AMR_evidence.output.json,
            deeparg=get_deeparg_sample_outputs,
            deeparg_lock=get_deeparg_lock,
        output:
            json="results/{project}/output/evidence/manifests/{sample}.json",
        log:
            "logs/{project}/evidence/{sample}.lock.log",
        conda:
            "../envs/python.yaml"
        params:
            sample_metadata=lambda wildcards: culture_free_sample_metadata(
                wildcards.sample
            ),
            evidence_filters={
                "minimum_mapping_quality": culture_free_config()["assembly-evidence"][
                    "minimum-mapping-quality"
                ],
                "minimum_base_quality": culture_free_config()["assembly-evidence"][
                    "minimum-base-quality"
                ],
                "maximum_pileup_depth": culture_free_config()["assembly-evidence"][
                    "maximum-pileup-depth"
                ],
            },
            tool_versions={
                "diamond": "2.1.12",
                "minimap2": "2.30",
                "samtools": "1.23",
            },
        script:
            "../scripts/lock_sample_evidence.py"
