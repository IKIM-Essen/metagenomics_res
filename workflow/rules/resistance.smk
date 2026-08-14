# Resistance analysis
if not config["card"]["use-shared"]:

    rule download_CARD_data:
        output:
            json=get_card_db_file(),
        log:
            "logs/CARD_data_download.log",
        conda:
            "../envs/unix.yaml"
        threads: 30
        params:
            download=config["card"]["url"],
            folder=lambda wildcards, output: Path(output.json).parent,
            filename=lambda wildcards, output: Path(output.json).name,
        shell:
            "(cd {params.folder} && "
            "wget {params.download} && "
            "tar -xvf data ./{params.filename}) > {log} 2>&1"


rule CARD_load_DB_for_reads:
    input:
        get_card_db_file(),
    output:
        temp(touch("results/CARD_load_DB_for_reads.done")),
    log:
        "logs/CARD_load_DB_for_reads.log",
    conda:
        "../envs/card.yaml"
    shell:
        "rgi clean --local && " "rgi load --card_json {input} --local > {log} 2>&1"


rule CARD_annotation:
    input:
        json=get_card_db_file(),
        load=rules.CARD_load_DB_for_reads.output,
    output:
        done=temp(touch("results/CARD_annotation.done")),
        ann=get_card_annotation_file(),
    log:
        "logs/CARD_annotation.log",
    conda:
        "../envs/card.yaml"
    params:
        folder=lambda wildcards, input: Path(input.json).parent,
        file=lambda wildcards, input: Path(input.json).name,
    shell:
        "(cd {params.folder}/ && "
        "rgi card_annotation -i {params.file}) && "
        "rgi load -i {input.json} --card_annotation {output.ann} --local > {log} 2>&1"


rule CARD_read_run:
    input:
        fa=get_filtered_gz_fastqs,
        db=rules.CARD_annotation.output,
    output:
        txt="results/{project}/output/resistance/CARD/reads/{sample}/{sample}.gene_mapping_data.txt",
    log:
        "logs/{project}/ARGs/reads/{sample}.log",
    conda:
        "../envs/card.yaml"
    threads: 60
    params:
        folder=lambda wildcards, output: Path(output.txt).parent,
    shell:
        "rgi bwt -1 {input.fa[0]} -2 {input.fa[1]} "
        "-o {params.folder}/{wildcards.sample} --local "
        "--clean -n {threads} > {log} 2>&1"


if direct_unicard_enabled():

    rule direct_uniCARD_reads:
        input:
            fastqs=get_filtered_gz_fastqs,
            dmdb=get_unicard_dmnd(),
        output:
            r1="results/{project}/output/resistance/uniCARD/direct/{sample}/{sample}.R1.tsv.gz",
            r2="results/{project}/output/resistance/uniCARD/direct/{sample}/{sample}.R2.tsv.gz",
        log:
            "logs/{project}/uniCARD/direct/{sample}.log",
        conda:
            "../envs/diamond.yaml"
        threads: 60
        resources:
            heavy=1,
        params:
            db_wo_ext=lambda wildcards, input: os.path.splitext(input.dmdb)[0],
            sensitivity=config["culture-free-evidence"]["direct-unicard"]["sensitivity"],
            evalue=config["culture-free-evidence"]["direct-unicard"]["evalue"],
            max_targets=config["culture-free-evidence"]["direct-unicard"][
                "max-target-seqs"
            ],
            outfmt=(
                "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend "
                "sstart send evalue bitscore qlen slen qcovhsp scovhsp qframe"
            ),
        shell:
            "set -o pipefail; ((diamond blastx -q {input.fastqs[0]} -d {params.db_wo_ext} "
            "--{params.sensitivity} --evalue {params.evalue} "
            "--max-target-seqs {params.max_targets} --outfmt {params.outfmt} "
            "--threads {threads} | gzip -c > {output.r1}) && "
            "(diamond blastx -q {input.fastqs[1]} -d {params.db_wo_ext} "
            "--{params.sensitivity} --evalue {params.evalue} "
            "--max-target-seqs {params.max_targets} --outfmt {params.outfmt} "
            "--threads {threads} | gzip -c > {output.r2})) > {log} 2>&1"

    rule rich_uniCARD_assembly_proteins:
        input:
            faa=get_proteins,
            dmdb=get_unicard_dmnd(),
        output:
            tsv="results/{project}/output/resistance/uniCARD/assembly_evidence/{sample}.tsv.gz",
        log:
            "logs/{project}/uniCARD/assembly_evidence/{sample}.log",
        conda:
            "../envs/diamond.yaml"
        threads: 60
        resources:
            heavy=1,
        params:
            db_wo_ext=lambda wildcards, input: os.path.splitext(input.dmdb)[0],
            sensitivity=config["culture-free-evidence"]["direct-unicard"]["sensitivity"],
            evalue=config["culture-free-evidence"]["direct-unicard"]["evalue"],
            max_targets=config["culture-free-evidence"]["direct-unicard"][
                "max-target-seqs"
            ],
            outfmt=(
                "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend "
                "sstart send evalue bitscore qlen slen qcovhsp scovhsp"
            ),
        shell:
            "set -o pipefail; diamond blastp -q {input.faa[0]} -d {params.db_wo_ext} "
            "--{params.sensitivity} --evalue {params.evalue} "
            "--max-target-seqs {params.max_targets} --outfmt {params.outfmt} "
            "--threads {threads} | gzip -c > {output.tsv} 2> {log}"

    rule export_CARD_ARO_categories:
        input:
            card_json=get_card_db_file(),
        output:
            tsv="results/{project}/output/resistance/uniCARD/card_aro_categories.tsv",
        log:
            "logs/{project}/uniCARD/card_aro_categories.log",
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/export_card_aro_categories.py"

    rule lock_culture_free_AMR_evidence:
        input:
            card_json=get_card_db_file(),
            unicard_fasta=get_uniCARD_db(),
            unicard_dmnd=get_unicard_dmnd(),
            hierarchy=get_CARD_hierarchy(),
            categories=rules.export_CARD_ARO_categories.output.tsv,
        output:
            json="results/{project}/output/resistance/uniCARD/evidence-lock.json",
        log:
            "logs/{project}/uniCARD/evidence_lock.log",
        conda:
            "../envs/python.yaml"
        params:
            card_version=config["card"]["version"],
            builder_commit=config["culture-free-evidence"]["direct-unicard"][
                "builder-commit"
            ],
            uniref_release=config["culture-free-evidence"]["direct-unicard"][
                "uniref-release"
            ],
            diamond_version="2.1.12",
            sensitivity=config["culture-free-evidence"]["direct-unicard"]["sensitivity"],
            evalue=config["culture-free-evidence"]["direct-unicard"]["evalue"],
            max_target_seqs=config["culture-free-evidence"]["direct-unicard"][
                "max-target-seqs"
            ],
        script:
            "../scripts/lock_amr_evidence.py"


if deeparg_enabled():

    rule stage_reads_for_DeepARG_SS:
        input:
            fastqs=get_filtered_gz_fastqs,
        output:
            fasta=temp("results/{project}/deeparg_staging/{sample}.paired.fasta"),
        log:
            "logs/{project}/deeparg/{sample}.stage.log",
        conda:
            "../envs/python.yaml"
        shell:
            "python workflow/scripts/paired_fastq_to_fasta.py "
            "--r1 {input.fastqs[0]} --r2 {input.fastqs[1]} --output {output.fasta} "
            "> {log} 2>&1"

    rule direct_DeepARG_SS:
        input:
            fasta=rules.stage_reads_for_DeepARG_SS.output.fasta,
        output:
            arg="results/{project}/output/resistance/deeparg/direct/{sample}/{sample}.mapping.ARG",
            potential="results/{project}/output/resistance/deeparg/direct/{sample}/{sample}.mapping.potential.ARG",
            alignments="results/{project}/output/resistance/deeparg/direct/{sample}/{sample}.align.daa.tsv",
        log:
            "logs/{project}/deeparg/{sample}.log",
        container:
            config["culture-free-evidence"]["deeparg"]["container"]
        threads: 16
        resources:
            heavy=1,
        params:
            prefix=lambda wildcards, output: str(
                Path(output.arg).parent / wildcards.sample
            ),
            data_path=config["culture-free-evidence"]["deeparg"]["data-path"],
            device=config["culture-free-evidence"]["deeparg"]["device"],
            min_probability=config["culture-free-evidence"]["deeparg"][
                "min-probability"
            ],
            identity=config["culture-free-evidence"]["deeparg"]["alignment-identity"],
            evalue=config["culture-free-evidence"]["deeparg"]["alignment-evalue"],
            overlap=config["culture-free-evidence"]["deeparg"]["alignment-overlap"],
            max_alignments=config["culture-free-evidence"]["deeparg"][
                "max-alignments-per-entry"
            ],
        shell:
            "deeparg predict --model SS --type nucl -i {input.fasta} -o {params.prefix} "
            "--data-path {params.data_path} --hf-device {params.device} "
            "--min-prob {params.min_probability} "
            "--arg-alignment-identity {params.identity} "
            "--arg-alignment-evalue {params.evalue} "
            "--arg-alignment-overlap {params.overlap} "
            "--arg-num-alignments-per-entry {params.max_alignments} > {log} 2>&1"

    rule lock_DeepARG_evidence:
        output:
            json="results/{project}/output/resistance/deeparg/evidence-lock.json",
        conda:
            "../envs/python.yaml"
        params:
            runtime_version=config["culture-free-evidence"]["deeparg"][
                "runtime-version"
            ],
            model_release=config["culture-free-evidence"]["deeparg"]["model-release"],
            database_path=config["culture-free-evidence"]["deeparg"]["data-path"],
            database_digest=config["culture-free-evidence"]["deeparg"][
                "database-digest"
            ],
            container=config["culture-free-evidence"]["deeparg"]["container"],
            parameters={
                "min_probability": config["culture-free-evidence"]["deeparg"][
                    "min-probability"
                ],
                "alignment_identity": config["culture-free-evidence"]["deeparg"][
                    "alignment-identity"
                ],
                "alignment_evalue": config["culture-free-evidence"]["deeparg"][
                    "alignment-evalue"
                ],
                "alignment_overlap": config["culture-free-evidence"]["deeparg"][
                    "alignment-overlap"
                ],
                "max_alignments_per_entry": config["culture-free-evidence"]["deeparg"][
                    "max-alignments-per-entry"
                ],
            },
        script:
            "../scripts/write_deeparg_provenance.py"


rule uniCARD_makeDB:
    input:
        db=get_uniCARD_db(),
    output:
        dmdb=get_unicard_dmnd(),
    log:
        "logs/unicard_makeDB.log",
    conda:
        "../envs/diamond.yaml"
    threads: 60
    params:
        filename=lambda wildcards, input: Path(input.db).name,
        filename_wo_ext=Path(get_uniCARD_db_wo_ext()).name,
        folder=lambda wildcards, input: Path(input.db).parent,
    shell:
        "(cd {params.folder} && diamond makedb --in {params.filename} -d {params.filename_wo_ext}) > {log} 2>&1"


rule uniCARD_run:
    input:
        faa=get_proteins,
        dmdb=get_unicard_dmnd(),
    output:
        tsv=temp(
            "results/{project}/output/resistance/uniCARD/assembly/{sample}_out.tsv"
        ),
    log:
        "logs/{project}/uniCARD/{sample}.log",
    # ensure priority over potential contig classification
    priority: 1
    conda:
        "../envs/diamond.yaml"
    threads: 60
    # to ensure only one of these commands are run at the same time
    resources:
        heavy=1,
    params:
        db_wo_ext=lambda w, input: os.path.splitext(input.dmdb)[0],
        max_targets=15,
    shell:
        "diamond blastp -q {input.faa} -d {params.db_wo_ext} "
        "-p {threads} --max-target-seqs {params.max_targets} "
        "-o {output.tsv} > {log} 2>&1"


rule gz_uniCARD_tsv:
    input:
        unicard_in=rules.uniCARD_run.output.tsv,
    output:
        tsv=temp(
            "results/{project}/output/resistance/uniCARD/assembly/{sample}_out.tsv.gz"
        ),
    log:
        "logs/{project}/uniCARD/{sample}_gz.log",
    conda:
        "../envs/unix.yaml"
    threads: 16
    shell:
        "(gzip -k {input.unicard_in}) > {log} 2>&1"


rule filter_uniCARD:
    input:
        tsv=rules.gz_uniCARD_tsv.output.tsv,
        json=get_CARD_hierarchy(),
    output:
        csv="results/{project}/output/resistance/uniCARD/assembly/{sample}.csv",
    log:
        "logs/{project}/uniCARD/{sample}_filter.log",
    conda:
        "../envs/python.yaml"
    threads: 16
    params:
        # filter window for filtering UniCARD results, default: 1%
        filter_window=0.01,
    shell:
        "python workflow/scripts/uniCARD_filtering.py --infile {input.tsv} "
        "--card_hierarchy {input.json} --outfile {output.csv} "
        "--filter_window {params.filter_window} > {log} 2>&1"


rule resistance_abundance_per_sample:
    input:
        csv=rules.filter_uniCARD.output.csv,
        json=get_CARD_hierarchy(),
    output:
        abd_class="results/{project}/output/resistance/uniCARD/per_sample_resistance_abundance/{sample}_drug_class_abundance.csv",
        abd_ARGs="results/{project}/output/resistance/uniCARD/per_sample_resistance_abundance/{sample}_ARGs_abundance.csv",
    log:
        "logs/{project}/uniCARD/{sample}_abundance.log",
    conda:
        "../envs/python.yaml"
    threads: 5
    script:
        "../scripts/uniCARD_abundance_sample.py"


rule resistance_abundance_all:
    input:
        contig_files=expand(
            "results/{{project}}/output/classification/assembly/{sample}/{sample}_classified_contigs.csv",
            sample=get_samples(),
        ),
        abd_class_files=expand(
            "results/{{project}}/output/resistance/uniCARD/per_sample_resistance_abundance/{sample}_drug_class_abundance.csv",
            sample=get_samples(),
        ),
        abd_ARGs_files=expand(
            "results/{{project}}/output/resistance/uniCARD/per_sample_resistance_abundance/{sample}_ARGs_abundance.csv",
            sample=get_samples(),
        ),
    output:
        abd_class_possible="results/{project}/output/resistance/uniCARD/{project}_drug_class_abundance_all_possible.csv",
        abd_class_present="results/{project}/output/resistance/uniCARD/{project}_drug_class_abundance_all_present.csv",
        abd_ARGs_possible="results/{project}/output/resistance/uniCARD/{project}_ARGs_abundance_all_possible.csv",
        abd_ARGs_present="results/{project}/output/resistance/uniCARD/{project}_ARGs_abundance_all_present.csv",
    log:
        "logs/{project}/uniCARD/all_abundances.log",
    conda:
        "../envs/python.yaml"
    threads: 5
    script:
        "../scripts/uniCARD_abundance.py"
