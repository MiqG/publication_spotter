#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
# Last Update: 2021-01-13
#
# Purpose
# -------
# Read and preprocess files

import pandas as pd

def auc_split_point(arr):
    r"""
    return the arr point of a vector of values with the maximum value difference 
    from a line from y_min to y_max
    """
    arr = arr[np.argsort(arr)]
    # create line from min to max
    x = [0, len(arr)]
    y = [np.min(arr), np.max(arr)]
    coefficients = np.polyfit(x, y, 1)
    polynomial = np.poly1d(coefficients)

    # get line values from arr
    y_axis = polynomial(np.arange(len(arr)))

    # find where the difference is maximal between the lines
    diff = y_axis - arr
    idx = np.where(diff == np.max(diff))
    x_max = arr[idx][0]
    return x_max


def read_files(data_file, metadata_file):
    """
    Read data matrix (genomic feature x samples) and sample metadata.
    """
    data = pd.read_table(data_file, index_col=0)
    metadata = pd.read_table(metadata_file)
    return data, metadata


def prepare_data(data, metadata, sample_col):
    """
    Keep only samples metadata that are also in data. The output will have the same order.
    """
    common_samples = set(data.columns).intersection(set(metadata[sample_col]))
    
    data = data[common_samples]
    
    metadata = metadata.set_index(sample_col).loc[common_samples]
    metadata = metadata.rename_axis(sample_col).reset_index()
    
    return data, metadata


def subset_data(data, metadata, sample_col, subset_col, subset_values):
    """
    Subset data and metadata based on the vaules in a column from metadata.
    """
    subset_values = subset_values.split(',')
    idx = metadata[subset_col].isin(subset_values)
    metadata = metadata.loc[idx]
    data = data[metadata[sample_col]]
    
    return data, metadata


def filter_data_by_std(data, thresh_std, axis=1):
    """
    Consider only those genomic features (rows) in data that are highly variant.
    
    Parameters
    ----------
    thresh_std: float; keep genomic features with std higher than 'thresh_std'.
                if 'thresh_std'<0, define thresh_std with the elbow method.
    """
    data_std = data.std(axis=axis)
    if thresh_std < 0: thresh_std = auc_split_point(data_std.values)
    rows_oi = data.index[data_std >= thresh_std]
    return data.loc[rows_oi], thresh_std



