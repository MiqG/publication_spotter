#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# utilities for general scripting and file processing.

import os
import uuid
import shutil
import pickle
import subprocess
from datetime import date, datetime
import numpy as np
import pandas as pd
from io import StringIO

import multiprocessing as mp

import matplotlib.pyplot as plt
import seaborn as sns

##### FILE PROCESSING #####
def create_required_dirs(file_name):
    dir_path = os.path.dirname(file_name)
    if not os.path.exists(dir_path): 
        print('creating directory:',dir_path)
        os.makedirs(dir_path)


def execute_ipynotebook(proj_root, params, report_template_file, out_report_file, timeout=None):
    # create temporary folder
    tmp_dir = os.path.join(proj_root, str(uuid.uuid4()))
    os.mkdir(tmp_dir)
    # save params as pickle at temporary folder
    pickle.dump( params, open(os.path.join(tmp_dir,'params.p'), 'wb') )
    # copy notebook template at temporary folder
    report_file = os.path.join(tmp_dir, os.path.basename(report_template_file))
    shutil.copyfile(report_template_file, report_file)
    # run notebook
    cmd = """jupyter nbconvert --to notebook --ExecutePreprocessor.timeout=%(timeout)s --inplace --execute %(report_file)s""" % {'report_file':report_file, 'timeout':timeout}
    print(cmd)
    shell(cmd)
    # copy executed notebook to reports folder
    shutil.copyfile(report_file, out_report_file)
    # remove temporary folder
    shutil.rmtree(tmp_dir)


def unix_join_files(files, result_file):
    """
    use UNIX core tools to inner join multiple files on their first column.
    """
    tmp = os.path.join(os.path.dirname(result_file),'tmp')
    # init
    cmd = "cp %(file)s %(result_file)s" % {'file':files[0], 'result_file':result_file}
    util.shell(cmd)
    for file in files[1:]:
        # we need to skip the first 6 columns of the second file
        # all combined files have the same length and are sorted
        cmd = """join --header --nocheck-order %(result_file)s <(<%(file)s cut -f2,7-) | tr " " "\t" | awk '{ t=$1; $1=$2; $2=t; print;}' OFS=$'\t' > %(tmp_file)s""" % {'file':file, 'result_file':result_file, 'tmp':tmp}               
        print(cmd)
        util.shell(cmd)
        # overwrite new file
        cmd = "mv %(tmp)s %(result_file)s" % {'result_file':result_file, 'tmp':tmp}
        util.shell(cmd)
        print('Finished joining', file)


def get_file_index(file, usecols=0):
    """returns the column of the file as a list, ommiting the first row"""
    return list(pd.read_table(file,usecols=[usecols],index_col=usecols).index)


def get_file_header(file, nrows=1):
    """returns the first line of a file"""
    return list(pd.read_table(file, nrows=nrows).columns)


##### SCRIPTING #####
def now():
    return datetime.now().strftime('%H:%M:%S')


def today():
    """Return formatted today's date"""
    return date.today().strftime('%Y-%m-%d')


def wget(url, path='.'):
    #cmd = 'wget --user-agent="Chrome" --no-clobber -P '+path+' '+url
    cmd = 'wget --user-agent="Chrome" -P '+path+' '+url
    subprocess.call(cmd, shell=True)
    return True


def shell(cmd):
    process = subprocess.Popen(cmd, 
                               shell=True, 
                               executable='/bin/bash', 
                               stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)

def deconvolve_nested_dicts(keys_of_interest, nested_dicts):
    """"""
    deconvolved = {}
    for key_of_interest in keys_of_interest:
        deconvolved_keys = [key_of_interest]
        for nested_dict in nested_dicts:
            if key_of_interest in nested_dict.keys(): deconvolved_keys = deconvolved_keys + nested_dict[key_of_interest]
        
        deconvolved[key_of_interest] = deconvolved_keys
    
    return deconvolved


def csv2pandas(string):
    return pd.read_csv(StringIO(string))


def DataFrame_to_list_of_Series(df, axis=0):
    if axis == 1:
        df_aslist = list(df.iterrows())
    elif axis == 0:
        df_aslist = list(df.iteritems())
    return [pd.Series(arr, name=name) for name, arr in df_aslist]


def parallel_apply(df, func, axis=0, processes=None):
    """parallelize
    Parameters
    ----------
    func : pd.Series
    """
    # split dataframe over axis
    df_aslist = DataFrame_to_list_of_Series(df, axis)
    
    pool = mp.Pool(processes=processes)
    df = pd.concat(pool.map(func, df_aslist), axis=1)
    pool.close()

    # prepare output
    if axis == 1: df = df.T
    
    return df

def write_table(df, file_name, sep='\t', **kws):
    if file_name.endswith('.gz') and 'compression' not in kws:
        df.to_csv(file_name, sep=sep, compression='gzip', **kws)
    else:
        df.to_csv(file_name, sep=sep, **kws)

def normalise(X):
    return (X - X.mean(1).values.reshape(-1,1))/X.std(1).values.reshape(-1,1)
