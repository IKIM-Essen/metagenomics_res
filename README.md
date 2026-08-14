# ResMAG - name pending

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


ResMAG is a state-of-the-art and user-friendly Snakemake workflow designed for the analysis of metagenomic data. It integrates multiple bioinformatics tools and algorithms to facilitate key steps in metagenome analysis, including quality control, assembling, bin refinement, metagenome-assembled genome (MAG) reconstruction and taxonomic classification. ResMAG has a special focus on highly diverse samples, such as wastewater, and the identification of antibiotic resistance genes.<br />

---

## Key Features

**Binning Techniques**: Employ a collection of two state-of-the-art binning tools to partition metagenomic contigs into individual bins, allowing for comprehensive and accurate analysis.<br />

**MAG Reconstruction**: Utilize cutting-edge algorithms to reconstruct high-quality metagenome-assembled genomes (MAGs), especially from highly diverse microbial communities, like wastewater.<br />

**Taxonomic Classification**: Apply advanced taxonomic classification methods to assign taxonomic labels to reads, contigs and MAGs and identify the microbial community composition within the metagenomic samples.<br />

**Antibiotic Resistance Gene Identification**: Perform in-depth analysis to detect and characterize antibiotic resistance genes within the metagenomic data, providing valuable insights into antimicrobial resistance profiles.<br />

**Culture-free evidence profile**: Optionally retain assembly/linkage evidence, run paired-mate direct
DIAMOND--UniCARD searches before assembly, preserve full competitive hit fields, export typed ARO
category accessions from pinned CARD data, and checksum-lock every AMR database artifact. These outputs
are research evidence and are not resistance, susceptibility, or clinical calls.<br />

**Performance Refinement**: Continuously optimize the pipeline by incorporating the latest advancements in metagenomics research, ensuring the highest accuracy and efficiency in metagenomic data analysis.<br />

---

## Overview

---

## Preparations

Create a snakemake environment using [mamba](https://mamba.readthedocs.io/en/latest/) via:

 ```mamba create -c conda-forge -c bioconda -n snakemake snakemake snakemake-storage-plugin-fs```

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Download GTDB
The GTDB needs to be downloaded and decompressed, it requires about 140 Gb.

1. Change to the directory where the GTDB should be stored
2. Download the latest or your desired version of GTDB
   Please make sure this version is compatible with GTDB-tk version `2.6.1`
   ```
   wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
   ```
3. Decompress the downloaded archive
   ```
   tar xzf gtdbtk_data.tar.gz
   ```
4. After successful step 3: the archive can be removed
5. Please specify the path to your decompressed GTDB in the config file (see [Configuring workflow](#configuring-the-workflow))

### Create UniCARD database
Please follow the steps outlined in the [UniCARD repository](https://github.com/IKIM-Essen/uniCARD) to create a UniCARD database for resistance gene annotation.
And specify the path to your database in the config file (see [Configuring workflow](#configuring-the-workflow))

---

## Usage

Please obtain a copy of this workflow by [cloning](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) this repository.
```
git clone https://github.com/IKIM-Essen/metagenomics_res.git
```

### Configuring the workflow
1. Edit the `config/config.yaml` file:
   - Specify a project name (`project-name`)
   - Specify filtering options for human reads (`human-filtering`)
   - Specify host filtering options, if you have a non-human host (`host-filtering`)
   - Specify options for different databases:
     - GTDB database needs to be downloaded before (see [Download GTDB](#download-gtdb))
     - UniCARD database needs to be create before (see [Create UniCARD database](#create-unicard-database))
     - other databases (kaiju, CheckM2, CARD, genomad) can be given as a local path or downloaded when running the pipeline
2. Provide sample information in the `config/pep/samples.csv` file while keeping the header and the format as:

```
sample_name,fq1,fq2
sample1,path/to/your/fastq/sample1_R1.fastq.gz,path/to/your/fastq/sample1_R2.fastq.gz
```

For culture-free AMR development, enable `culture-free-evidence` in `config/config.yaml`, replace the
UniCARD builder/UniRef placeholders with the exact database provenance, and decide whether to retain the
large assembly evidence. RGI-BWT is optional and disabled in the example profile. Direct
DIAMOND--UniCARD evidence is emitted separately for R1 and R2 with 18 fields so mate identity and
competing CARD wild-type/UniRef hits are not collapsed. Assembly proteins are searched again with a rich
17-field contract for comparison; this output is separate from the legacy filtered summary. The CARD
category export retains canonical `ARO:` identifiers; the legacy UniCARD hierarchy labels are
display-only.

The profile also exports the evidence that downstream strain and ARG-host inference otherwise loses:

- stable FASTG and GFA assembly graphs;
- coordinate-sorted read-to-contig BAM/BAI;
- mapping/base-quality-filtered contig coverage and within-contig allele counts;
- insert-size/alignment statistics plus one first-mate record per pair and aggregated contig links;
- synchronized paired-FASTQ counts, full pair/name digests, and a compact sequence identity sketch;
- a validated sample/control design and a checksum-addressed manifest for every sample.

Set `control-aware.enabled: true` to require the PEP columns `sample_role`, `case_id`,
`material_stage`, and `matched_control`. Raw specimens and plate sweeps must then name a sample whose
role is `negative_control`. This validates the experimental relationship; contamination correction and
clinical thresholds remain downstream, specimen-specific validation tasks.

DeepARG-SS is an optional secondary arm. Before enabling it, provide a tested modern DeepARG container,
runtime/model release, database path, and database digest. ResMAG retains `.mapping.ARG`,
`.mapping.potential.ARG`, and the documented `.align.daa.tsv` alignment table without converting broad
DeepARG classes into exact alleles or clinical calls. See the
[official DeepARG project](https://github.com/gaarangoa/deeparg) for its current runtime and model
distribution.


### Run the workflow
Activate the conda environment:
```conda activate snakemake```

Test your configuration by performing a dry-run via
```snakemake --use-conda -n```

Executing the workflow:
```snakemake --use-conda --cores $N -k```

using `$N` cores. It is recommended to use all available cores.

---

## Output

When the culture-free evidence profile is enabled, the principal additional outputs are:

```text
results/<project>/output/
├── evidence/
│   ├── sample_design.tsv
│   ├── sample_design.audit.json
│   ├── reads/<sample>.paired_fastq.audit.json
│   ├── assembly/<sample>_assembly_tree.{fastg,gfa}
│   ├── assembly/<sample>_{reads_mapped.bam,reads_mapped.bam.bai}
│   ├── assembly/<sample>.{contig_coverage,flagstat,samtools_stats,idxstats}.tsv
│   ├── assembly/<sample>.allele_counts.tsv.gz
│   ├── assembly/<sample>.paired_links.tsv.gz
│   ├── assembly/<sample>.paired_link_summary.tsv
│   └── manifests/<sample>.json
└── resistance/
    ├── uniCARD/direct/<sample>/<sample>.{R1,R2}.tsv.gz
    ├── uniCARD/assembly_evidence/<sample>.tsv.gz
    └── deeparg/direct/<sample>/...                 # optional
```

These are research evidence products. A hit is not a phenotype; an unlinked ARG must not be assigned to
every organism; and no-hit evidence must not be interpreted as susceptibility.

---

## Tools

A list of tools used in the pipeline:

| Tool      | Link                                              |
| --------- | ------------------------------------------------- |
| CARD      | https://doi.org/10.1093/nar&gkac920               |
| CheckM2   | https://doi.org/10.1038/s41592-023-01940-w        |
| CoverM    | https://github.com/wwood/CoverM                   |
| DAS Tool  | https://doi.org/10.1038/s41564-018-0171-1         |
| DIAMOND   | https://doi.org/10.1038/s41592-021-01101-x        |
| fastp     | https://doi.org/10.1093/bioinformatics/bty560     |
| FastQC    | www.bioinformatics.babraham.ac.uk/projects/fastqc |
| geNomad   | https://doi.org/10.1038/s41587-023-01953-y        |
| GTDB-Tk   | https://doi.org/10.1093/bioinformatics/btz848     |
| Kaiju     | https://doi.org/10.1038/ncomms11257               |
| MEGAHIT   | https://doi.org/10.1093/bioinformatics/btv033     |
| MetaBAT   | http://dx.doi.org/10.7717/peerj.1165              |
| MetaCoAG  | https://doi.org/10.1101/2021.09.10.459728         |
| minimap2  | https://doi.org/10.1093/bioinformatics/bty191     |
| MultiQC   | www.doi.org/10.1093/bioinformatics/btw354         |
| pprodigal | https://github.com/sjaenick/pprodigal             |
| Rust-Bio  | https://doi.org/10.1093/bioinformatics/btv573     |
| samtools  | https://doi.org/10.1093/gigascience/giab008       |
| Snakemake | www.doi.org/10.12688/f1000research.29032.1        |
| UniCARD   | https://github.com/IKIM-Essen/uniCARD             |

---

## Contributions

Pull requests and feature suggestions are very welcome! Feel free to fork and submit improvements.
For any other questions, or feedback, please contact the project maintainer Josefa Welling (@josefawelling) at josefa.welling@uk-essen.de.

We appreciate your input and support in using and improving ResMAG.

---

## Acknowledgements

We would like to express our gratitude towards Katharina Block, Adrian Doerr, Miriam Balzer, Alexander Thomas, Johannes Köster, Ann-Kathrin Doerr and the IKIM who have contributed to the development and testing of ResMAG. Their valuable insights and feedback have been helpful throughout the creation of the workflow.

---

## Citation

A paper is on its way. If you use ResMAG in your work, don't forget to give credits to the authors by citing the URL of this repository.

---

## License

ResMAG is licensed under the [BSD-2 Clause](https://www.open-xchange.com/hubfs/2_Clause_BSD_License.pdf?hsLang=en).
