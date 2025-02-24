import pandas as pd

# ARG thresholds
min_coverage = 50
min_identity = 50

inpath = '/local/work/josefa/metagenomics_res/'
samples = ['WWLD2401_Mock_red']#, 'WWLD2401_red']

def get_normalized_score(value, maximum, threshold):
    norm = (value - threshold)/(maximum - threshold)
    return(norm)

def get_score(cov, ident, plas, af, ani):
    weight_h = 2/5
    weight_l = 1/5
    score=round(((weight_l*cov + weight_l*ident + weight_h*plas + weight_l*af + weight_l*ani)*100),2)
    return(score)

for sample in samples:

    ## input files
    args_in = f'{inpath}results/test/output/ARGs/assembly/{sample}/{sample}_ARGs_contigs.csv'
    c2b_in = f'{inpath}results/test/output/contig2bins/{sample}/DASTool_contig2bin.tsv'
    clf_in = f'{inpath}results/test/output/classification/bins/{sample}/{sample}.summary.tsv'
    tax_tab = f'{inpath}results/test/output/report/{sample}/{sample}_bin_taxonomy.csv'
    plas_in = f'{inpath}results/test/output/plasmids/{sample}/{sample}_plasmid_summary.tsv'

    arg_df=pd.read_csv(args_in,index_col='ORF_ID')

    cols=['Best_Hit_ARO', 'Drug Class', 'Contig']
    arg_df_red = arg_df[cols]

    arg_df_red['norm_coverage'] = arg_df['%coverage'].map(lambda x: get_normalized_score(x,100,min_coverage))
    arg_df_red['norm_identity'] = arg_df['%identity'].map(lambda x: get_normalized_score(x,100,min_identity))

    # Add plasmid information
    plas_df = pd.read_table(plas_in, usecols=['seq_name','plasmid_score'])
    plas_df.rename({'seq_name':'Contig'}, inplace=True, axis=1)
    arg_df_red.insert(2,'on_plasmid', arg_df_red['Contig'].isin(plas_df['Contig']))
    # Merge arg_df_red with c2b_df based on Contig, keeping all arg_df_red rows (left join)
    arg_df_red = arg_df_red.merge(plas_df, on='Contig', how='left')

    arg_df_red['norm_plas'] = arg_df_red['plasmid_score'].map(lambda x: get_normalized_score(x,1,0.7))

    # Add if ARG is in bin/MAG
    c2b_df = pd.read_table(c2b_in,names=['Contig', 'bin'])

    # Merge arg_df_red with c2b_df based on Contig, keeping all arg_df_red rows (left join)
    arg_df_red = arg_df_red.merge(c2b_df, on='Contig', how='left')

    # Add classification of bin
    tax_df = pd.read_csv(tax_tab,usecols=['bin','species'])

    # Merge arg_df_red with clf_df based on bin, keeping all arg_df_red rows (left join)
    arg_df_red = arg_df_red.merge(tax_df, on='bin', how='left')

    # Add 
    clf_df = pd.read_table(clf_in,usecols=['user_genome','closest_genome_ani','closest_genome_af'])
    clf_df.rename({'user_genome':'bin'}, inplace=True, axis=1)

     # Merge arg_df_red with clf_df based on bin, keeping all arg_df_red rows (left join)
    arg_df_red = arg_df_red.merge(clf_df, on='bin', how='left')
    arg_df_red['norm_ani'] = arg_df_red['closest_genome_ani'].map(lambda x: get_normalized_score(x,100,95))
    arg_df_red['norm_af'] = arg_df_red['closest_genome_af'].map(lambda x: get_normalized_score(x,1,0.5))

    for index, row in arg_df_red.iterrows():

        arg_df_red.loc[index,'confidence'] = get_score(row['norm_coverage'], row['norm_identity'], row['norm_plas'], row['norm_af'], row['norm_ani'])
    print(arg_df_red)
    #doest work because of missing values