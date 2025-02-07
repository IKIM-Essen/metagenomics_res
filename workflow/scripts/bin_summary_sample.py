import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], 'w')

## input files
in_dastool= snakemake.input.tool
in_checkm = snakemake.input.checkm
in_gtdb=snakemake.input.gtdb

## output files
csv_path_bins = snakemake.output.csv_bins
csv_path_tax = snakemake.output.csv_tax 


## parameters for MAG decision
max_cont= snakemake.params.max_cont
min_comp= snakemake.params.min_comp

## helper to split taxonomic assignment
tax_level_encoder = {'d':'domain', 'p':'phylum', 'c':'class', 'o':'order', 'f':'family', 'g':'genus', 's':'species'}

def find_lowest_assigned_level(row):
    # Start from species and go up
    for level in reversed(tax_level_encoder.values()):  
        # If the value is not empty, return the level
        if row[level] != '':  
            return level
    # If all are empty, return empty string
    return ''  


### include information from DASTool report ###

# read in these columns from DASTool report
tool_cols=['bin','contigs','N50']
tool_df = pd.read_table(in_dastool, usecols=tool_cols, index_col='bin')


### include information from checkm2 report ###

# read in these columns from checkm2 report and rename to new column name
checkm_cols = {'Name':'bin', 'Completeness':'completeness', 'Contamination':'contamination', 'Genome_Size':'genome_size', 'GC_Content':'GC_content'}
checkm_df = pd.read_table(in_checkm, usecols=checkm_cols.keys())
checkm_df.rename(checkm_cols,axis=1,inplace=True)

# remove .fa from bin identifier
checkm_df['bin']=checkm_df['bin'].apply(lambda x: x.replace('.fa','') if x.find('.fa') >= 0 else x)
checkm_df.set_index('bin',inplace=True)


### include information from GTDB report ### 

# read in these columns from GTDB report and rename identifier column
gtdb_cols=['user_genome','classification','closest_genome_reference']
gtdb_df=pd.read_table(in_gtdb, usecols=gtdb_cols)
gtdb_df.rename({'user_genome':'bin'},axis=1,inplace=True)

taxonomy_split = gtdb_df['classification'].str.split(';', expand=True)

# go through each level of the taxonomy
for col in taxonomy_split:
    # check on each level if it's not assigned, than set as empty string
    for index, value in taxonomy_split[col].items():
        if value == None or type(value) == float or value.split('__')[-1] == '':
            taxonomy_split.loc[index,col] = ''

taxonomy_split.columns = tax_level_encoder.values()
taxonomy_split.set_index(gtdb_df['bin'], inplace=True)

gtdb_df.set_index('bin',inplace=True)


### write bin report combining different information ###

bins_df = pd.concat([tool_df,checkm_df,gtdb_df], axis=1)
col_order = ['completeness','contamination','genome_size','GC_content','contigs','N50','closest_genome_reference','classification']
bins_df=bins_df[col_order]
bins_df.index.name ='bin'

# insert new column if bin is a MAG
bins_df.insert(0,'is_MAG','bin')
bins_df['is_MAG']= bins_df.apply(lambda row: 'MAG' if row['completeness'] >= min_comp and row['contamination'] <= max_cont else 'bin',axis=1)

# sorted for MAGs first, than highest completeness and after that lowest contamination
bins_df.sort_values(['is_MAG','completeness','contamination'],inplace=True, ascending=[True,False,True])

# save bin report to csv file
bins_df.to_csv(csv_path_bins)

# Combine taxonomy information with MAG information 
tax_df = pd.concat([bins_df['is_MAG'],taxonomy_split], axis=1)


# Apply function to find lowest assigned level
tax_df.insert(1,'lowest_tax_level', tax_df.apply(find_lowest_assigned_level, axis=1))

tax_df.to_csv(csv_path_tax)