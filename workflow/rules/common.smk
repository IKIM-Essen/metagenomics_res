import os


configfile: "config/config.yaml"


def get_data_path():
    return config["data-handling"]["data"]


def get_resource_path():
    return config["data-handling"]["resources"]


def get_project():
    return config["project-name"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fastqs(wildcards):
    file_r1 = pep.sample_table.loc[wildcards.sample]["fq1"]
    folder = str(Path(file_r1).parent)
    filename_r1 = Path(file_r1).name
    filename_r2 = Path(pep.sample_table.loc[wildcards.sample]["fq2"]).name
    return [folder, filename_r1, filename_r2]


def get_local_fastqs(wildcards):
    path = get_data_path()
    return (
        "{data}{{project}}/{{sample}}_R1.fastq.gz".format(data=path),
        "{data}{{project}}/{{sample}}_R2.fastq.gz".format(data=path),
    )


def get_adapters(wildcards):
    return config["adapter-seqs"]


def get_trimmed_fastqs(wildcards):
    return [
        "results/{project}/trimmed/fastp/{sample}.1.fastq.gz",
        "results/{project}/trimmed/fastp/{sample}.2.fastq.gz",
    ]


def get_prefiltered_fastqs(wildcards):
    if config["host-filtering"]["do-host-filtering"]:
        return [
            "results/{project}/host_filtering/non_host/{sample}_R1.fastq.gz",
            "results/{project}/host_filtering/non_host/{sample}_R2.fastq.gz",
        ]
    else:
        return get_trimmed_fastqs(wildcards)


def get_host_map_statistics(wildcards):
    if config["host-filtering"]["do-host-filtering"]:
        logs = expand(
            "results/{{project}}/report_prerequisites/qc/filter_host_{sample}.log",
            sample=get_samples(),
        )
        return logs
    else:
        return []


def get_human_ref():
    if config["human-filtering"]["use-local"]:
        path = config["human-filtering"]["local-path"]
    else:
        path = config["human-filtering"]["download-path"]
    filename = path.split("/")[-1]
    local_ref = "".join([get_resource_path(),"ref_genome/", filename])
    return local_ref


def get_human_local_folder():
    path = config["human-filtering"]["local-path"]
    folder = Path(path).parent
    return folder


def get_checkm2_db():
    file = "{}{}".format(get_resource_path(), config["checkm2"])
    return file


def get_filtered_fastqs(wildcards):
    return [
        "results/{project}/filtered/fastqs/{sample}_R1.fastq",
        "results/{project}/filtered/fastqs/{sample}_R2.fastq",
    ]


def get_filtered_gz_fastqs(wildcards):
    return [
        "results/{project}/filtered/fastqs/{sample}_R1.fastq.gz",
        "results/{project}/filtered/fastqs/{sample}_R2.fastq.gz",
    ]


def get_assembly(wildcards):
    return "results/{project}/megahit/{sample}/final.contigs.fa"

    
def get_gz_assembly(wildcards):
    return "results/{project}/output/fastas/{sample}/{sample}.fa.gz"


def get_kaiju_files():
    file = "".join([get_resource_path(), config["kaiju"]["fmi-file"]])
    path = str(Path(file).parent)
    fmi = Path(file).name
    names = ["nodes.dmp", fmi, "names.dmp"]
    files = ["/".join([path,name]) for name in names]
    return files


## binning parameters
def get_contig_length_threshold():
    return config["min-contig-length"]


def get_binners():
    return config["das-tool"]["binner-list"]


def get_all_contig2bin_files(wildcards):
    binners = get_binners()
    file_list = [
        "".join(["results/{project}/output/contig2bins/{sample}/",binner,"_contig2bin.tsv"])
        for binner in binners
    ]
    return file_list


## reads in binner control file and returns list with paths to contig2bin files
## and a list with name of the binners that produced results
def get_paths_binner(wildcards):
    file = "results/{}/das_tool/binner_control_{}.csv".format(wildcards.project,wildcards.sample)
    lines = open(file).readlines()
    paths = str(lines[0].rstrip("\n"))
    binner = str(lines[1].rstrip("\n"))
    return paths, binner


def bins_for_sample(wildcards):
    if len(get_paths_binner[0]) > 0:
        return True
    else:
        return False


def get_mag_fa(wildcards):
    folder = "results/{project}/output/fastas/{sample}/mags/"
    files = [
        os.path.join(folder, binID)
        for binID in os.listdir(folder)
        if binID.endswith("fa.gz")
    ]
    return files


def get_binIDs_for_sample(wildcards):
    folder = "results/{project}/output/fastas/{sample}/mags/"
    binIDs = [binID for binID in os.listdir(folder) if binID.endswith("fa.gz")]
    return binIDs


def get_mag_ARGs(wildcards):
    bin_fastas = (get_mag_fa(wildcards),)
    bin_fastas = list(bin_fastas)[0]
    binIDs = [os.path.basename(binID) for binID in bin_fastas]
    folder = "results/{project}/output/ARGs/mags/{sample}/"
    arg_files = [
        os.path.join(folder, binID.replace(".fa.gz", ".txt")) for binID in binIDs
    ]
    return arg_files


def get_gtdb_folder():
    path = "".join([get_resource_path(), config["gtdb"]["db-folder"]])
    return path


def get_genomad_DB_file():
    path = "".join([get_resource_path(),"genomad_db/names.dmp"])
    return path


def get_card_db_file():
    path = "".join([get_resource_path(),"CARD_db/",config["card"]["dbfile"]])
    return path


def get_card_annotation_file():
    version = config["card"]["version"]
    path = "".join([get_resource_path(),"CARD_db/card_database_",version,".fasta"])
    return path
