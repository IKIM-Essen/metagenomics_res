# ResMAG - name pending

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


ResMAG is a state-of-the-art and user-friendly Snakemake workflow designed for the analysis of metagenomic data. It integrates multiple bioinformatics tools and algorithms to facilitate key steps in metagenome analysis, including bin refinement, metagenome-assembled genome (MAG) reconstruction, taxonomic classification of MAGs, and identification of antibiotic resistance genes.<br />

### Key Features

**Binning Techniques**: Employ a collection of five state-of-the-art binning tools to partition metagenomic contigs into individual bins, allowing for comprehensive and accurate analysis.<br />

**MAG Reconstruction**: Utilize cutting-edge algorithms to reconstruct high-quality metagenome-assembled genomes (MAGs) from sequencing data.<br />

**Taxonomic Classification**: Apply advanced taxonomic classification methods to assign taxonomic labels to MAGs and identify the microbial community composition within the metagenomic samples.<br />

**Antibiotic Resistance Gene Identification**: Perform in-depth analysis to detect and characterize antibiotic resistance genes within the metagenomic data, providing valuable insights into antimicrobial resistance profiles.<br />

**Performance Refinement**: Continuously optimize the pipeline by incorporating the latest advancements in metagenomics research, ensuring the highest accuracy and efficiency in metagenomic data analysis.<br />

### Overview
```mermaid
%%{init: {
   'theme':'base',
   'themeVariables': {
      'secondaryColor': '#fff',
      'tertiaryColor': '#fff',
      'tertiaryBorderColor' : '#fff'}
   }}%%

flowchart TB;

   subgraph " "
      direction TB

      %% Nodes
      A[/short reads/]
      B["<b>QC</b> <br> <i>fastp<i>"]
      C["<b>Host read filtering</b> <br> <i>Kraken 2<i>"]
      D["<b>Assembly</b> <br> <i>MegaHIT</i>"]
      E["<b>Binning</b>"]
      F["<b>Bin refinement</b> <br> <i>DAS Tool<i>"]
      G[/MAGs/]
      H["<b>Resistance analysis</b> <br> <i>HyDRA<i>"]
      I["<b>taxonomic classification</b> <br> <i>Kaiju and GTDBTk<i>"]
      J[/MultiQC report/]
      K[/Assembly summary/]

      %% input & output node design
      classDef in_output fill:#fff,stroke:#cde498,stroke-width:4px
      class A,G,J,K in_output
      %% rule node design
      classDef rule fill:#cde498,stroke:#000
      class B,C,D,E,F,H,I rule

      %% Node links
      A --> B
      B --> C
      B --- J
      C --> D
      D --> E
      D ---- K
      E --"<i>MetaBAT 2<i>"--> F
      E --"<i>MetaBinner<i>"--> F
      E --"<i>MetaCoAG<i>"--> F
      E --"<i>Rosella<i>"--> F
      E --"<i>Vamb<i>"--> F
      F --> G
      G --- H
      G --- I

   end

```

## Usage
### Preparations
To prepare the workflow
1. Clone it to your desired working folder via git or your preferred IDE
2. Edit the `config/config.yaml` file:
   - Specify a project name (`project-name`)
   - Specify filtering options for human reads (`human-filtering`)
   - Specify host filtering options, if you have a non-human host (`host-filtering`)
   - Specify options for GTDB database (see [Downloading GTDB](#Downloading-GTDB))
3. Provide a sample information in the `config/pep/samples.csv` file with keeping the header and format as `.csv`:

```
sample_name,fq1,fq2
sample1,path/to/your/fastq/sample1_R1.fastq.gz,path/to/your/fastq/sample1_R2.fastq.gz
```

#### Downloading GTDB
The GTDB files need to be downloaded and unarchived, it requires about 110 Gb.
1. Create a new folder `resources/gtdb/` and change to this directory
2. Download the latest version of GTDB
3. Unarchive the downloaded file
4. After successful step 3: the archive can be removed

```
mkdir resources/gtdb/
cd resources/gtdb/
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
rm gtdbtk_data.tar.gz
```

### Install Snakemake
Create a snakemake environment using [mamba](https://mamba.readthedocs.io/en/latest/) via:

 ```mamba create -c conda-forge -c bioconda -n snakemake snakemake```

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Run the workflow
Activate the conda environment:
```conda activate snakemake```

Test your configuration by performing a dry-run via
```snakemake --use-conda -n```

Executing the workflow:
```snakemake --use-conda --cores $N --rerun-incomplete```

using `$N` cores. It is recommended to use all available cores.

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

## Output

## Tools

A list of tools used in the pipeline:

| Tool      | Link                                              |
| --------- | ------------------------------------------------- |
| CoverM    | https://github.com/wwood/CoverM                   |
| DAS Tool  | https://doi.org/10.1038/s41564-018-0171-1         |
| fastp     | https://doi.org/10.1093/bioinformatics/bty560     |
| FastQC    | www.bioinformatics.babraham.ac.uk/projects/fastqc |
| MEGAHIT   | https://doi.org/10.1093/bioinformatics/btv033     |
| MetaBAT   | http://dx.doi.org/10.7717/peerj.1165              |
| MetaCoAG  | https://doi.org/10.1101/2021.09.10.459728         |
| minimap2  | https://doi.org/10.1093/bioinformatics/bty191     |
| MultiQC   | www.doi.org/10.1093/bioinformatics/btw354         |
| samtools  | https://doi.org/10.1093/gigascience/giab008       |
| Snakemake | www.doi.org/10.12688/f1000research.29032.1        |

## Contact Information

For any questions, or feedback, please contact the project maintainer at josefa.welling@uk-essen.de. We appreciate your input and support in using and improving ResMAG.

## Acknowledgements

We would like to express our gratitude towards Katharina Block, Adrian Doerr, Miriam Balzer, Alexander Thomas, Johannes Köster, Ann-Kathrin Doerr and the IKIM who have contributed to the development and testing of ResMAG. Their valuable insights and feedback have been helpful throughout the creation of the workflow.

## Citation

A paper is on its way. If you use ResMAG in your work before the paper, then please consider citing this GitHub.

## License

ResMAG is released under the [BSD-2 Clause](https://www.open-xchange.com/hubfs/2_Clause_BSD_License.pdf?hsLang=en). Please review the license file for more details.