import pandas as pd
import json
import sys

sys.stderr = open(snakemake.log[0], "w")

# one of reads, assembly, mags
in_case = snakemake.params.case

# minimal coverage and identity of reference sequence
min_coverage = snakemake.params.coverage
min_identity = snakemake.params.identity

def filter_RGI_main(df, path_to_json):

    with open(path_to_json, 'r') as f:
        data = json.load(f)

    cov_dict={}
    for orfid in df.index.to_list():
        reference_id=df.at[orfid,"ID"]

        #extract hit start and end from json
        hit_start= data[orfid][reference_id]["hit_start"]
        hit_end= data[orfid][reference_id]["hit_end"]

        #calculate hit length on aminoacids
        hit_len = (hit_end-hit_start)/3

        #extract reference (aminoacids) length from json
        target_str = data[orfid][reference_id]["sequence_from_broadstreet"]
        target_len = len(target_str)

        #calculate coverage
        coverage=round(((hit_len/target_len)*100),2)

        cov_dict[orfid]=coverage

    #remove reference ID from output df 
    df.pop('ID')
    df.insert(2, '%identity', df.pop('Best_Identities'))
    df.insert(3, '%coverage', df.index.map(cov_dict))

    # filter for minimum coverage and identity
    df_filt = df[df['%identity'] >= min_identity]
    df_filt = df_filt[df_filt['%coverage'] >= min_coverage]

    # sort ARGs by %coverage and %identity
    df_filt=df_filt.sort_values(["%coverage","%identity"], ascending=[False,False])

    return(df_filt)


if in_case == "reads":
    infile=snakemake.input.txt
    outfile=snakemake.output.csv

    cols=["ARO Term","Average Percent Coverage","All Mapped Reads","Drug Class","Resistance Mechanism","AMR Gene Family"]

    reads_df=pd.read_table(infile, index_col='ARO Accession')
    reads_df_red=reads_df[cols]
    reads_df_red.rename({"Average Percent Coverage":'%coverage'},axis=1,inplace=True)

    reads_df_filt = reads_df_red[reads_df_red['%coverage'] > min_coverage]
    reads_df_filt=reads_df_filt.sort_values(["%coverage","All Mapped Reads"],ascending=False)
    
    reads_df_filt.to_csv(outfile)


elif in_case == "assembly":
    infile=snakemake.input.txt
    json_f = snakemake.input.json
    outfile=snakemake.output.csv

    cols=["Contig","Best_Hit_ARO","Best_Identities","Drug Class","Antibiotic","Resistance Mechanism","AMR Gene Family","ARO","ID"]

    asbl_df=pd.read_table(infile, index_col="ORF_ID")
    asbl_df_red=asbl_df[cols]
    asbl_df_red.loc[:,'Contig'] = asbl_df_red['Contig'].map(lambda x: x.rsplit('_', 1)[0])

    asbl_df_filt = filter_RGI_main(asbl_df_red,json_f)
    asbl_df_filt.to_csv(outfile)


'''
if in_case == "mags":
    infile=f"/local/work/josefa/ResMAG/results/clinic_june_240903/output/ARGs/mags/{sample}/{sample}.txt"
    json_f= f"/local/work/josefa/ResMAG/results/clinic_june_240903/output/ARGs/mags/{sample}/{sample}.json"
    outfile=f"/local/work/josefa/ResMAG/results/clinic_june_240903/output/ARGs/mags/{sample}/{sample}_mags_ARGs.csv"



    
    for sample in samples:
    sample_path=f"{inpath}{sample}/"
    outfile=f"{sample_path}/ARGs_mags_{sample}.csv"
    list_outfile=f"{sample_path}/ARGs_mags_list_summary_{sample}.csv"

    mags=[f for f in os.listdir(sample_path) if f.endswith(".txt")]
    first_MAG=True
    
    for mag in mags:
        infile=f"{sample_path}{mag}"
        json_f= infile.replace(".txt",".json")

        df=pd.read_table(infile)
        if df.empty:
            continue
        df.set_index("ORF_ID",inplace=True)
        df_red=df[cols]

        with open(json_f, 'r') as f:
            data = json.load(f)
        eval_dict={}
        for orfid in df.index.to_list():
            ident=df.at[orfid,"ID"]
            evalue=data[orfid][ident]["evalue"]
            eval_dict[orfid]=evalue
        eval_df=pd.DataFrame.from_dict(eval_dict,orient='index',columns=["evalue"])
        df_red=pd.concat([eval_df,df_red],axis=1)
        if first_MAG:
            df_all=df_red
            first_MAG=False
        else:
            df_all=pd.concat([df_all,df_red])
    if len(mags)<1:
        df_all=df_red[0:0]
    df_all=df_all.sort_values(["Cut_Off","Best_Hit_Bitscore"], ascending=[True,False])
    df_all.to_csv(outfile)
'''

