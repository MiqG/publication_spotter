    #
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Sample clustering based on their splicing signatures.
# 
# Outline
# -------
# 1. Load splicing matrix PANCAN + CCLE
# 2. Compute Z-scores by cancer
# 3. Embed events into lower dimensions with PCA and UMAP (with Seurat's defaults)
# 4. Cluster

import pandas as pd
import gc
import argparse
import numpy as np
from sklearn.decomposition import PCA
from umap import UMAP
from umap.umap_ import nearest_neighbors
import leidenalg
from scanpy.neighbors import _compute_connectivities_umap
from scanpy._utils import get_igraph_from_adjacency

"""
Development
-----------
import os
ROOT = '~/projects/publication_splicing_dependency'
PREP_DIR = os.path.join(ROOT,'data','prep')
matrix_file = os.path.join(PREP_DIR,'exon_psi_imputed','CCLE-EX.tsv.gz')
features_as_rows = True
log_transform = True
standardize = True
"""

SAVE_PARAMS = {'sep':'\t', 'index':False, 'compression':'gzip'}


##### FUCTIONS #####
def load_data(filename):
    # read
    df = pd.read_table(filename, index_col=0)
    gc.collect()
    return df


def preprocess_matrix(X, log_transform, standardize):
    if log_transform:
        X = np.log2(X + 1)
    if standardize: 
        X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
        
    return X
    

def reduce_dimensions(X, pca_kws={'n_components':50}, umap_kws={'n_neighbors':30, 'metric':'cosine', 'min_dist':0.3}):
    """
    Perform PCA and UMAP dimension reduction.
    """
    # PCA
    pca = PCA(**pca_kws)
    pcs = pca.fit_transform(X)
    
    # UMAP with PCs
    umap = UMAP(**umap_kws)
    umaps = umap.fit_transform(pcs)
    
    # prepare outputs
    #fitted = {'pca':pca,'umap':umap}
    pcs = pd.DataFrame(pcs, 
                       columns = ['PC%s'%n for n in range(pcs.shape[1])], 
                       index = X.index)
    explained_variance = pd.DataFrame([pca.explained_variance_ratio_],
                                      index = pcs.index,
                                      columns = ['explained_variance_ratio-PC%s'%n for n in range(pcs.shape[1])])
    umaps = pd.DataFrame(umaps,
                         columns = ['UMAP%s'%n for n in range(umaps.shape[1])],
                         index = X.index)
    
    return pcs, explained_variance, umaps


def leiden_clustering(X, nn_kws = {'n_neighbors':30, 'metric':'cosine', 'metric_kwds':{}, 'angular':False, 'random_state':np.random}, partition_type = leidenalg.RBConfigurationVertexPartition, leiden_kws = {'n_iterations':-1, 'seed':0}):
    # compute nearest neighbors with UMAP
    knn_indices, knn_dists, forest = nearest_neighbors(X, **nn_kws)
    
    # compute connectivities
    distances, connectivities = _compute_connectivities_umap(knn_indices, knn_dists, X.shape[0], nn_kws['n_neighbors'])
    
    # use connectivites as adjacency matrix to get igraph
    G = get_igraph_from_adjacency(connectivities, directed=True)
    
    # run leiden on graph
    leiden_kws['weights'] = np.array(G.es['weight']).astype(np.float64)
    
    partition = leidenalg.find_partition(G, partition_type, **leiden_kws)
    labels = np.array(partition.membership)
    
    return labels


def cluster_samples(X):
    # reduce dimensions
    print('Reducing dimensions...')
    pcs, explained_variance, umaps = reduce_dimensions(X)
    
    # cluster using PCs
    print('Doing Leiden Clustering...')
    labels = leiden_clustering(pcs.values)
    
    # prepare output
    event_clustering = pd.concat([pcs, explained_variance, umaps],axis=1)
    event_clustering['leiden_labels'] = labels
    event_clustering = event_clustering.reset_index()
    
    return event_clustering
    
    
def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--matrix_file',type=str)
    parser.add_argument('--output_file',type=str)
    parser.add_argument('--features_as_rows',type=bool,default=True)
    parser.add_argument('--standardize',type=bool,default=True)
    parser.add_argument('--log_transform',type=bool,default=False)
    
    args = parser.parse_args()
    return args

    
def main():
    args = parse_args()
    matrix_file = args.matrix_file
    output_file = args.output_file
    features_as_rows = args.features_as_rows
    standardize = args.standardize
    log_transform = args.log_transform
    
    # load
    matrix = load_data(matrix_file)
    
    # preprocess
    if features_as_rows: matrix = matrix.T.copy()
    matrix = preprocess_matrix(matrix, log_transform, standardize)
    
    # cluster
    result = cluster_samples(matrix.dropna(axis=1))
    
    # save
    result.to_csv(output_file, **SAVE_PARAMS)
    

##### SCRIPTS #####
if __name__ == '__main__':
    main()
    print('Done!')