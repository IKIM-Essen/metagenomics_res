import os


configfile: "config/config.yaml"


def get_root():
    return os.getcwd()


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fastqs(wildcards):
    return (
        pep.sample_table.loc[wildcards.sample]["fq1"],
        pep.sample_table.loc[wildcards.sample]["fq2"],
    )


def get_project():
    return config["project-name"]


def get_adapters(wildcards):
    return config["adapter_seqs"]


def get_threshold():
    return config["binning"]["min_contig_length"]


def get_kmersize():
    return config["binning"]["kmer_length"]


def expand_samples_for_project(paths, **kwargs):
    def inner(wildcards):
        return expand(
            paths,
            sample=get_samples(),
            **kwargs,
        )

    return inner


def get_trimmed_fastqs(wildcards):
    return [
        "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
        "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
    ]


def get_bacterial_reads(wildcards):
    return [
        "results/{project}/filtered/bacteria/{sample}_bacteria_1.fastq.gz",
        "results/{project}/filtered/bacteria/{sample}_bacteria_2.fastq.gz",
    ]


def get_kraken_db():
    return config["kraken"]["kraken-db"]


def get_kraken_ref():
    return config["kraken"]["ref-analysis"]


def get_taxID_dict():
    return config["kraken"]["taxIDs-ref"]


def get_taxID(wildcards):
    return config["kraken"]["taxIDs-ref"][wildcards.kraken_ref]


def get_rosella_install():
    folder = get_rosella_folder()
    script = f"{folder}/install.sh"
    return script


def get_rosella_git():
    return config["rosella"]["gitURL"]


def get_rosella_folder():
    return config["rosella"]["rosella_dir"]


def get_root():
    return os.getcwd()


def get_paths_binner(wildcards):
    file = f"results/{wildcards.project}/das_tool/{wildcards.sample}/binner_control.csv"
    lines = open(file).readlines()
    paths = str(lines[0].rstrip('\n'))
    binner= str(lines[1].rstrip('\n'))
    return paths, binner

def get_all_contig2bin(wildcards):
    return [
        "results/{project}/metabinner/{sample}/metabinner_res/metabinner_result.tsv",
        "results/{project}/binning_rev/{sample}/vamb_contig2bin.tsv",
        "results/{project}/binning_rev/{sample}/metabat_contig2bin.tsv",
        "results/{project}/metacoag/{sample}/contig_to_bin.tsv",
        "results/{project}/binning_rev/{sample}/rosella_contig2bin.tsv",
    ]