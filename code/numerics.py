#!/usr/bin/env python
import time
import mpathic.simulate_library
import mpathic.fast as fast
import mpathic.qc as qc
from mpathic.profile_mut import main as profile_mut
import numpy as np
from scipy.sparse import csr, csr_matrix, lil_matrix
import pdb
import sys

# Given a large list of sequences, produce a (sparse) mutation array
def dataset2mutarray(dataset_df, modeltype, chunksize=1000, rowsforwtcalc=100):

    # Determine the type of model and set seq2array function appropriately
    if modeltype=='MAT':
        seqs2array = mpathic.fast.seqs2array_for_matmodel
    elif modeltype=='NBR':
        seqs2array = mpathic.fast.seqs2array_for_nbrmodel
    else:
        raise SortSeqError('Unknown model type: %s'%modeltype)

    # Determine seqtype, etc.
    seqcol = 'seq'
    seqtype = qc.colname_to_seqtype_dict[seqcol]
    wtcol = qc.seqtype_to_wtcolname_dict[seqtype]

    # Compute the wt sequence
    rowsforwtcalc = min(rowsforwtcalc,dataset_df.shape[0])
    dataset_head_df = dataset_df.head(rowsforwtcalc)
    mut_df = profile_mut(dataset_head_df)
    wtseq = ''.join(list(mut_df[wtcol]))
    wtrow = seqs2array([wtseq], seq_type=seqtype).ravel().astype(bool)
    numfeatures = len(wtrow)

    # Process dataframe in chunks
    startrow = 0
    endrow = startrow+chunksize-1
    numrows = dataset_df.shape[0]

    # Fill in mutarray (a lil matrix) chunk by chunk
    mutarray_lil = lil_matrix((numrows,numfeatures),dtype=int)
    matrix_filled = False
    while not matrix_filled:

        if startrow >= numrows:
            matrix_filled = True
            continue
        elif endrow >= numrows:
            endrow = numrows-1
            matrix_filled = True

        # Compute seqarray
        seqlist = list(dataset_df[seqcol][startrow:(endrow+1)])
        seqarray = seqs2array(seqlist, seq_type=seqtype)

        # Remove wt entries
        tmp = seqarray.copy()
        tmp[:,wtrow] = 0

        # Store results from this chunk
        mutarray_lil[startrow:(endrow+1),:] = tmp

        # Increment rows
        startrow = endrow+1
        endrow = startrow + chunksize - 1

    # Convert to csr matrix
    mutarray_csr = mutarray_lil.tocsr()

    # Return vararray as well as binary representation of wt seq
    return mutarray_csr, wtrow


def eval_modelmatrix_on_mutarray(modelmatrix, mutarray, wtrow):

    # Compute constant contribution to model prediciton
    modelmatrix_vec = modelmatrix.ravel()
    const_val = np.dot(wtrow,modelmatrix_vec)

    # Prepare matrix for scanning mutarray
    tmp_matrix = modelmatrix.copy()
    indices = wtrow.reshape(modelmatrix.shape).astype(bool)
    wt_matrix_vals = tmp_matrix[indices]
    tmp_matrix -= wt_matrix_vals[:,np.newaxis]
    modelmatrix_for_mutarray = csr_matrix(np.matrix(tmp_matrix.ravel()).T)

    # Compute values
    mutarray_vals = mutarray*modelmatrix_for_mutarray
    vals = const_val + mutarray_vals.toarray().ravel()
    return vals


