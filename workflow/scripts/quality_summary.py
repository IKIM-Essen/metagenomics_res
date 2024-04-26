import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")

jsons = snakemake.input.jsons
human_logs = snakemake.input.human_logs

if snakemake.params.other_host:
    hostname=snakemake.params.hostname
    host_logs = snakemake.input.host_logs

outfile = snakemake.output.csv #"/local/work/josefa/ResMAG/results/lanuv/output/report/all/seq_summary.csv"
outfile_vis = snakemake.output.vis_csv #"/local/work/josefa/ResMAG/results/lanuv/output/report/all/seq_summary_vis.csv"

results_dict = {}

for json_file in jsons:
    sample = json_file[json_file.rfind("/") + 1 : json_file.rfind(".fastp.json")]
    sample_results_dict = {}
        
    with open(json_file, "r") as read_file:
        jdata = json.load(read_file)
        # fastp json shows reads, not read pairs
        reads_before=int(jdata["summary"]["before_filtering"]["total_reads"])
        sample_results_dict["#reads"] = reads_before

        total_pre_host_filt = int(jdata["summary"]["after_filtering"]["total_reads"])
        sample_results_dict["#reads_afterQC"] = total_pre_host_filt

        sample_results_dict["%reads_filtered_by_quality"] = round(((1 - (total_pre_host_filt/reads_before)) * 100),4)

        sample_results_dict["#bp_afterQC"] = int(jdata["summary"]["after_filtering"]["total_bases"])

        sample_results_dict["%Q30_reads_afterQC"] = float(jdata["summary"]["after_filtering"]["q30_rate"]) * 100

    if snakemake.params.other_host:
        for host_log in host_logs:
            if host_log.find(sample) >= 0:
                with open(host_log, "r") as file:
                    for line in file.readlines():
                        if line.find("processed") >= 0:
                            non_host_reads = int(line.split()[2]) * 2
                            no_reads = total_pre_host_filt - non_host_reads
                            prct_reads = round(((int(no_reads) / int(total_pre_host_filt)) * 100),4)
                            break

                sample_results_dict[f"#{hostname}_reads"] = no_reads
                sample_results_dict[f"%{hostname}"] = prct_reads
                break


    for human_log in human_logs:
        if human_log.find(sample) >= 0:
            with open(human_log, "r") as file:
                for line in file.readlines():
                    if line.find("processed") >= 0:
                        non_human_reads = int(line.split()[2]) * 2
                        if snakemake.params.other_host:
                            no_reads = non_host_reads - non_human_reads
                        else:
                            no_reads = total_pre_host_filt - non_human_reads
                        prct_reads = round(((int(no_reads) / int(total_pre_host_filt)) * 100),4)
                        break

            sample_results_dict[f"#human_reads"] = no_reads
            sample_results_dict[f"%human"] = prct_reads 
            break
    
    sample_results_dict["#reads_after_filtering"] = non_human_reads
    sample_results_dict["%reads_filtered_total"] = round(((1 - (non_human_reads/reads_before)) * 100),4)

    results_dict[sample] = sample_results_dict


results_df = pd.DataFrame.from_dict(results_dict, orient="index")
results_df.index.name = "sample"

results_df=results_df.sort_index()
results_df.to_csv(outfile)

## edit decimal separator for visual of table

def format_int_with_commas(x):
    #Formats an integer with commas as thousand separators
    return f"{x:,}"

def format_float_with_percent(x):
    #Formats a float with percent sign
    return f"{x:.2f}%"

cols=[col for col in results_df.columns.to_list() if col.find("#")>=0]
results_df[cols]=results_df[cols].map(format_int_with_commas)

cols=[col for col in results_df.columns.to_list() if col.find("%")>=0]
results_df[cols]=results_df[cols].map(format_float_with_percent)

results_df.to_csv(outfile_vis)
