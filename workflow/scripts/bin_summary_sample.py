import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

## input files
in_dastool = (
    snakemake.input.tool
)  # "ResMAG/results/autobrewer/das_tool/ABS_24_07/ABS_24_07_DASTool_summary.tsv"
in_checkm = (
    snakemake.input.checkm
)  # "ResMAG/results/autobrewer/qc/checkm2/ABS_24_07/quality_report.tsv"
in_gtdb = (
    snakemake.input.gtdb
)  # "ResMAG/results/autobrewer/classification/ABS_24_07/ABS_24_07.bac120.summary.tsv"

## output files
### csv
csv_path_mags = (
    snakemake.output.csv_mags
)  # "ResMAG/results/autobrewer/report/ABS_24_07/mags_summary.csv" #
csv_path_bins = (
    snakemake.output.csv_bins
)  # "ResMAG/results/autobrewer/report/ABS_24_07/bin_summary.csv"
csv_path_tax = (
    snakemake.output.csv_tax
)  # "ResMAG/results/autobrewer/report/ABS_24_07/bin_taxonomy.csv"
csv_path_checkm = (
    snakemake.output.csv_checkm
)  # "ResMAG/results/autobrewer/report/ABS_24_07/checkm_summary.csv"
csv_path_dastool = (
    snakemake.output.csv_dastool
)  # "ResMAG/results/autobrewer/report/ABS_24_07/DASTool_summary.csv"

## params
max_cont = snakemake.params.max_cont
min_comp = snakemake.params.min_comp


def save_csv_table(csv_path, summary_df):
    summary_df.to_csv(csv_path)

    summary_df.columns.name = "bin"
    summary_df.index.name = None


tool_df = pd.read_table(in_dastool)
tool_df.drop(["bin_set"], axis=1, inplace=True)
# tool_df.rename({"bin":"bin"},axis=1,inplace=True)
tool_df.set_index("bin", inplace=True)
col_order = [
    "bin_score",
    "contigs",
    "size",
    "N50",
    "SCG_set",
    "unique_SCGs",
    "SCG_completeness",
    "redundant_SCGs",
    "SCG_redundancy",
]
tool_df = tool_df[col_order]

del col_order[:2]
tool_red_df = tool_df.drop(col_order, axis=1)

save_csv_table(csv_path_dastool, tool_red_df)


checkm_df = pd.read_table(in_checkm)
rm_list = ["Translation_Table_Used", "Additional_Notes"]
checkm_df.drop(rm_list, axis=1, inplace=True)
checkm_df.columns = checkm_df.columns.str.lower()
checkm_df.rename(
    {"name": "bin", "gc_content": "GC_content", "contig_n50": "contig_N50"},
    axis=1,
    inplace=True,
)
checkm_df["bin"] = checkm_df["bin"].str.replace(".fa", "", 1, case=True)
checkm_df.set_index("bin", inplace=True)
col_order = [
    "completeness",
    "contamination",
    "genome_size",
    "GC_content",
    "contig_N50",
    "total_coding_sequences",
    "coding_density",
    "average_gene_length",
    "completeness_model_used",
]
checkm_df = checkm_df[col_order]

del col_order[:5]
checkm_red_df = checkm_df.drop(col_order, axis=1)
save_csv_table(csv_path_checkm, checkm_red_df)


gtdb_df = pd.read_table(in_gtdb)
gtdb_df.rename({"user_genome": "bin"}, axis=1, inplace=True)
gtdb_df.set_index("bin", inplace=True)
cols = [
    "classification",
    "closest_genome_reference",
    "closest_genome_reference_radius",
    "closest_genome_ani",
    "closest_genome_af",
    "classification_method",
]
gtdb_red_df = gtdb_df[cols]

save_csv_table(csv_path_tax, gtdb_red_df)


bins_df = pd.concat(
    [
        tool_red_df,
        checkm_red_df,
        gtdb_red_df[
            ["classification", "closest_genome_reference", "closest_genome_ani"]
        ],
    ],
    axis=1,
)
col_order = [
    "completeness",
    "contamination",
    "bin_score",
    "contigs",
    "genome_size",
    "contig_N50",
    "GC_content",
    "classification",
    "closest_genome_reference",
]
bins_df = bins_df[col_order]
bins_df.index.name = "bin"

save_csv_table(csv_path_bins, bins_df)


mags_df = bins_df.loc[
    (bins_df["completeness"] >= min_comp) & (bins_df["contamination"] <= max_cont), :
]
mags_df.index.name = "MAG"
save_csv_table(csv_path_mags, mags_df)


"""
if we want to sort by mean of deviation from min completeness & max contamination

mean_dict={}
for binid in bins_df.index:
    cont=bins_df.at[binid,"contamination"]
    cont_per=round(((max_contamination - cont)/max_contamination), 4)
    comp=bins_df.at[binid,"completeness"]
    comp_per=round((comp/100),4)
    mean_dict[binid]=round(((cont_per + comp_per)/2),4)

to_sort=pd.DataFrame.from_dict(mean_dict,orient="index",columns=['mean'])

# index list sorted by mean used to reindex original df
new_ind=to_sort.sort_values(["mean"],ascending = [False]).index.to_list()
sorted_all=bins_df.reindex(new_ind)

"""
