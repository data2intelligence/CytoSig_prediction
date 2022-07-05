#!/usr/bin/env python

import os, pathlib, tables, pandas, numpy, re
import scipy.sparse as sp_sparse

from scipy import io
from glob import glob

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

src_path = os.path.join(base_path, 'src')
data_path = os.path.join(base_path, 'data')

sc_path = os.path.join(data_path, 'single_cell')


def filter_matrix(matrix):
    # jump empty barcodes
    matrix = matrix.loc[~matrix.index.isnull(), matrix.sum() > 0]
    
    # jump all zero genes
    matrix = matrix.loc[(matrix == 0).mean(axis=1) < 1]
    
    # jump ambiguous genes, if any
    cnt_map = matrix.index.value_counts()
    if cnt_map.max() > 1:
        matrix = matrix.loc[cnt_map.loc[matrix.index] == 1]
    
    return matrix



def get_matrix_from_h5(filename):
    fin = tables.open_file(filename, 'r')

    mat_group = fin.get_node(fin.root, 'matrix')
    barcodes = fin.get_node(mat_group, 'barcodes').read()
    data = getattr(mat_group, 'data').read()
    indices = getattr(mat_group, 'indices').read()
    indptr = getattr(mat_group, 'indptr').read()
    shape = getattr(mat_group, 'shape').read()
    matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        
    feature_ref = {}
    feature_group = fin.get_node(mat_group, 'features')
    feature_ids = getattr(feature_group, 'id').read()
    feature_names = getattr(feature_group, 'name').read()
    feature_types = getattr(feature_group, 'feature_type').read()
    feature_ref['id'] = feature_ids
    feature_ref['name'] = feature_names
    feature_ref['feature_type'] = feature_types
        
    tag_keys = getattr(feature_group, '_all_tag_keys').read()
        
    for key in tag_keys:
        feature_ref[key] = getattr(feature_group, key.decode("utf-8")).read()
        
    matrix = pandas.DataFrame.sparse.from_spmatrix(matrix)
    matrix.index = map(lambda v: v.decode("utf-8"), feature_ref['name'])
    matrix.columns = map(lambda v: v.decode("utf-8"), barcodes)
    
    fin.close()
    
    return filter_matrix(matrix)




def load_COVID19_GSE145926():
    """
    Load a COVID19 single-cell RNASeq data from patient lavage samples
    Final column names as cell_type.patient.barcode
    """
    
    fprefix = os.path.join(sc_path, 'GSE145926')
    
    info = pandas.read_excel(os.path.join(fprefix, 'patient_group.xlsx'), index_col=0)['Group']
    
    # the initial h5 files are downloaded from GEO
    file_list = glob(os.path.join(fprefix, '*.h5'))
    
    merge = []
    
    for f in file_list:
        title = os.path.basename(f).split('_')[1]
        print(title)
        
        data = get_matrix_from_h5(f)
       
        cat = info.loc[title]
        
        data.columns = cat + '_' + title + '.' + data.columns
        merge.append(data)
    
    data = pandas.concat(merge, axis=1, join='inner')
    data = data.loc[(data == 0).mean(axis=1) < 1]
    
    # down-sizing by 10 
    ratio = 1E5/data.sum()
    data *= ratio
    
    data = numpy.log2(data + 1)
    
    # add cell types
    ID_map = {
        'S1': 'severe_C143',
        'S2': 'severe_C145',
        'S3': 'severe_C146',
        'S4': 'severe_C148',
        'S5': 'severe_C149',
        'S6': 'severe_C152',
            
        'M1': 'mild_C141',
        'M2': 'mild_C142',
        'M3': 'mild_C144',
            
        'H1': 'healthy_C51',
        'H2': 'healthy_C52',
        'H3': 'healthy_C100',
    }
        
    # general cell map, these annotations are labeled by ourselves
    info = pandas.read_csv(os.path.join(fprefix, 'cluster_results.tsv'), sep='\t', index_col=0)
    info = info.loc[:, 'clusterID']
        
    cluster = pandas.read_csv(os.path.join(fprefix, 'cluster_annotation.txt'), sep='\t', index_col=0)
    cluster = cluster.loc[:, 'cellTypeSub']
        
    cell_type = info.apply(lambda v: cluster.loc[v])
    cell_type.index = [ID_map[v.split('_')[0]] + '.' + v.split('_')[1] for v in info.index]
        
    # T cell map
    info = pandas.read_csv(os.path.join(fprefix, 'cluster_T_results.tsv'), sep='\t', index_col=0)
    info = info.loc[:, 'clusterID']
            
    cluster = pandas.read_csv(os.path.join(fprefix, 'cluster_T_annotation.txt'), sep='\t', index_col=0)
    cluster = cluster.loc[:, 'cellTypeSub']
        
    Tcell_type = info.apply(lambda v: cluster.loc[v])
    Tcell_type.index = [ID_map[v.split('_')[0]] + '.' + v.split('_')[1] for v in info.index]
        
    cntmap = cell_type.loc[Tcell_type.index].value_counts()
    assert cntmap.shape[0] == 1
        
    cell_type.drop(Tcell_type.index, inplace=True)    
    
    # put in cell type to data columns
    data_1 = data.loc[:, cell_type.index]
    data_1.columns = cell_type.loc[data_1.columns] + '.' + data_1.columns
        
    data_2 = data.loc[:, Tcell_type.index]
    data_2.columns = Tcell_type.loc[data_2.columns] + '.' + data_2.columns
        
    data = pandas.concat([data_1, data_2], axis=1, join='inner')
    
    print(data.shape)
    data.sparse.to_dense().to_pickle(fprefix + '.pickle.gz', compression='gzip')
            


def load_mtx_to_dataframe(data, gene, barcode):
    data = io.mmread(data)
    data = pandas.DataFrame.sparse.from_spmatrix(data)
    
    gene = pandas.read_csv(gene, sep='\t', header=None).iloc[:,0]
    barcode = pandas.read_csv(barcode, sep='\t', header=None).iloc[:,0]
    
    data.index = gene
    data.columns = barcode
    
    return filter_matrix(data)



def load_COVID19_EGAS00001004571():
    fpath = os.path.join(sc_path, 'EGAS00001004571')
    
    # step 1: Seurat object to mtx matrix with normalized data
    # The initial rds file is downloaded from https://beta.fastgenomics.org/datasets/detail-dataset-952687f71ef34322a850553c4a24e82e#Files
    os.system(' '.join([
        os.path.join(src_path, 'convert_seurat.R'),
        os.path.join(fpath, 'seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds'),
        os.path.join(fpath, 'convert')
        ]))
    
    # step 2: sparse matrix to python pickle
    metadata = os.path.join(fpath, 'convert', 'metadata')
    metadata = pandas.read_csv(metadata, sep='\t', index_col=0, low_memory=False)
    
    postfix = re.compile('_[0-9]+$')
    
    cell_type = metadata['cluster_labels_res.0.4'].apply(lambda v: re.sub(postfix, '', v).replace('.', '-'))
    info = cell_type + '.' + metadata['group_per_sample'] + '_' + metadata['sampleID'].apply(lambda v: v.replace('.', '-'))
    
    # these data are seurat normalized already, no need to transform again
    data = load_mtx_to_dataframe(
        os.path.join(fpath, 'convert', 'matrix.mtx'),
        os.path.join(fpath, 'convert', 'genes.tsv'),
        os.path.join(fpath, 'convert', 'barcodes.tsv')
    )
    
    data.columns = info.loc[data.columns] + '.' + data.columns
    
    data.sparse.to_dense().to_pickle(os.path.join(sc_path, 'EGAS00001004571_cohort1.pickle.gz'))
    



def main():
    load_COVID19_GSE145926()
    load_COVID19_EGAS00001004571()
        
    return 0

if __name__ == '__main__': main()
