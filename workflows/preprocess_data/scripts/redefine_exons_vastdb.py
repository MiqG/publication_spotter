# Script purpose
# --------------
# We would like to be able to map mutations to the genome annotation created by
# vastDB, based on GCRh38.
# The first and last exons are not in VastDB.
# We need to treat the different types of event differently:
# - HsaEX: they are the most common exons in human RNA-Seq with exon skipping. 
# - HsaALTA and HsaALTD: as they have been mapped from the exon-exon junction perspective
#   we are missing coordinates for these when they are at the first or last exon of a gene.
# - HsaIR: we should not consider splice sites.
#
# We ommit events length 0.
#
# Outline
# -------
# For exons, we need to consider the splice sites' margins.
# - for EX, substract/add 2 to each side
# - for ALTD, add only to the end
# - for ALTA, substract only at start
# - for INT, do nothing

from sso_targets import util
import pandas as pd
import numpy as np
import argparse
import multiprocessing as mp

# variables
MARGIN = 2

##### FUNCTIONS ######    
def split_coords(coords):
    """
    turn 'chr2:12345-4567' into
    """
    coords = coords.str.split(':|-',expand=True)
    coords.columns = ['chr','start','end']
    coords[coords == ''] = np.nan
    coords['start'] = coords['start'].astype(float)
    coords['end'] = coords['end'].astype(float)
    return coords


def add_margin(df, margin=MARGIN):
    # sort based on 'start'
    df = df.sort_values('start')
    df = df.reset_index(drop=True)
    if len(df)>1:
        # add margin
        out = pd.DataFrame({
            'EVENT': df['EVENT'],
            'start_wmargin': df['start'].values-margin,
            'end_wmargin': df['end'].values+margin
        })
        # remove margin from first and last exon(s)
        #out.at[0,'start_wmargin']  = out.at[0,'start_wmargin'] + margin
        #out.at[len(out)-1,'end_wmargin']  =  out.at[len(out)-1,'end_wmargin'] - margin
    else:
        out = pd.DataFrame({
            'EVENT': df['EVENT'],                
            'start_wmargin': df['start'].values,
            'end_wmargin': df['end'].values
        })
    return out


def process_exon_info(exon_info, margin=MARGIN):
    # split coordinates
    exon_info[['chr','start','end']] = split_coords(coords = exon_info['CO_A'])
    
    # add margins to EX
    idx = exon_info['EVENT'].str.contains('EX')
    exon_info.loc[idx,'start_wmargin'] = exon_info['start'] - margin
    exon_info.loc[idx,'end_wmargin'] = exon_info['end'] + margin
    
    # add only end margin to ALTD
    idx = exon_info['EVENT'].str.contains('ALTD')
    exon_info.loc[idx,'start_wmargin'] = exon_info['start']
    exon_info.loc[idx,'end_wmargin'] = exon_info['end'] + margin

    # add start margin to ALTA
    idx = exon_info['EVENT'].str.contains('ALTA')
    exon_info.loc[idx,'start_wmargin'] = exon_info['start'] - margin
    exon_info.loc[idx,'end_wmargin'] = exon_info['end']
        
    # maintain coordinates in INT
    idx = exon_info['EVENT'].str.contains('INT')
    exon_info.loc[idx,'start_wmargin'] = exon_info['start']
    exon_info.loc[idx,'end_wmargin'] = exon_info['end'] 
    
    # annotate event types
    exon_info['EVENT_type'] = np.nan
    exon_info.loc[exon_info['EVENT'].str.contains('INT'),'EVENT_type'] = 'IR'
    exon_info.loc[exon_info['EVENT'].str.contains('EX'),'EVENT_type'] = 'AS'
    exon_info.loc[exon_info['EVENT'].str.contains('ALTD'),'EVENT_type'] = 'ALTD'
    exon_info.loc[exon_info['EVENT'].str.contains('ALTA'),'EVENT_type'] = 'ALTA'

    return exon_info


def process_exon_info_bygene(exon_info, n_jobs):
    """
    [Deprecated]
    For every gene, add corresponding exon margins
    """
    # split coordinates
    exon_info[['chr','start','end']] = split_coords(coords = exon_info['CO_A'])
    
    # add margin
    exon_info_by_gene = [x for _, x in exon_info.groupby('GENE')]
    pool = mp.Pool(n_jobs)
    result = pool.map(add_margin, exon_info_by_gene)
    pool.close()
    result = pd.concat(result)

    # join result by EVENT name
    processed = pd.merge(exon_info, result, how='left', on='EVENT')
    
    # for those that did not have a gene, start_wmargin and end_wmargin equals 
    # regular start and end
    idxs = list(processed.index[processed['GENE'].isna()])
    for idx in idxs:
        processed.at[idx,'start_wmargin'] = processed.at[idx,'start']
        processed.at[idx,'end_wmargin'] = processed.at[idx,'end']
    
    return processed


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exon_info_file',type=str)
    parser.add_argument('--output_file',type=str)
    parser.add_argument('--n_jobs', type=int, default=-1)
    args = parser.parse_args()
    return args


def main():
    # read arguments
    args = parse_args()
    exon_info_file = args.exon_info_file
    output_file = args.output_file
    n_jobs = args.n_jobs
    
    # load data
    exon_info = pd.read_table(exon_info_file)
    # for the moment we deal with 'EX', 'IR', and middle ALTA/D
    idx = (~exon_info['CO_A'].isna())
    exon_info = exon_info.loc[idx].reset_index().copy()
    
    # process
    result = process_exon_info(exon_info)
    
    # save
    util.create_required_dirs(output_file)
    util.write_table(result, output_file, index=False)


##### SCRIPT ######
"""
Development
-----------
exon_info_file = '/mnt/0D9213830D921383/Documents/phd/projects/sso_targets/data/references/EVENT_INFO-hg38_noseqs.tsv'

"""
if __name__=='__main__':
    main()
    print('Done!')



