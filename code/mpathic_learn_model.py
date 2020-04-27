#!/usr/bin/env python

'''A script which produces linear energy matrix models for a given data set.'''
from __future__ import division
#Our standard Modules
import argparse
import numpy as np
import scipy as sp
import sys
#Our miscellaneous functions
import pandas as pd
import utils
import pymc
import os
import gauge
import numerics


def MaximizeMI_memsaver(
        seq_mat,df,emat_0,wtrow,db=None,burnin=1000,iteration=30000,thin=10,
        runnum=0,verbose=False,temp=4200):
    '''This function will run a MCMC sampling with pymc'''
    #set n_seqs, this should be equal to the number of independent promoter variiants
    n_seqs = temp
    #set stochastic variables for MCMC
    @pymc.stochastic(observed=True,dtype=pd.DataFrame)
    def pymcdf(value=df):
        return 0
    @pymc.stochastic(dtype=float)
    def emat(p=pymcdf,value=emat_0):
        '''eval the current model, then return the log likelihood'''
        p['val'] = numerics.eval_modelmatrix_on_mutarray(np.transpose(value),seq_mat,wtrow)
        MI = EstimateMutualInfoforMImax.alt4(p.copy())  # New and improved
        return n_seqs*MI
    #save current run as an sqlite database
    if db:
        dbname = db + '_' + str(runnum) + '.sql'
        M = pymc.MCMC([pymcdf,emat],db='sqlite',dbname=dbname)
    else:
        M = pymc.MCMC([pymcdf,emat])
    M.use_step_method(stepper.GaugePreservingStepper,emat)
    if not verbose:
        M.sample = shutthefuckup(M.sample)
    M.sample(iteration,thin=thin)
    emat_mean = np.mean(M.trace('emat')[burnin:],axis=0)
    return emat_mean

def main(df,lm='IM',modeltype='MAT',LS_means_std=None,\
    db=None,iteration=30000,burnin=1000,thin=10,\
    runnum=0,initialize='LS',start=0,end=None,foreground=1,\
    background=0,alpha=0,pseudocounts=1,test=False,drop_library=False,\
    verbose=False,init_name=None,temp=4200):

    # Determine dictionary
    seq_cols = 'seqs'
    dicttype = 'dna'

    seq_dict,inv_dict = utils.choose_dict(dicttype,modeltype=modeltype)

    '''Check to make sure the chosen dictionary type correctly describes
         the sequences. An issue with this test is that if you have DNA sequence
         but choose a protein dictionary, you will still pass this test bc A,C,
         G,T are also valid amino acids'''
    #set name of sequences column based on type of sequence
    type_name_dict = {'dna':'seq','rna':'seq_rna','protein':'seq_pro'}
    seq_col_name = type_name_dict[dicttype]
    lin_seq_dict,lin_inv_dict = utils.choose_dict(dicttype,modeltype='MAT')
    par_seq_dict = {v:k for v,k in seq_dict.items() if k != (len(seq_dict)-1)}
    #drop any rows with ct = 0
    df = df[df.loc[:,'ct'] != 0]
    df.reset_index(drop=True,inplace=True)

    #If there are sequences of different lengths, then print error but continue
    if len(set(df[seq_col_name].apply(len))) > 1:
         sys.stderr.write('Lengths of all sequences are not the same!')
    #select target sequence region
    df.loc[:,seq_col_name] = df.loc[:,seq_col_name].str.slice(start,end)
    df = utils.collapse_further(df)
    col_headers = utils.get_column_headers(df)
    #make sure all counts are ints
    df[col_headers] = df[col_headers].astype(int)
    #create vector of column names
    val_cols = ['val_' + inv_dict[i] for i in range(len(seq_dict))]
    df.reset_index(inplace=True,drop=True)
    #Drop any sequences with incorrect length
    if not end:
        '''is no value for end of sequence was supplied, assume first seq is
            correct length'''
        seqL = len(df[seq_col_name][0]) - start
    else:
        seqL = end-start
    df = df[df[seq_col_name].apply(len) == (seqL)]
    df.reset_index(inplace=True,drop=True)
    #Do something different for each type of learning method (lm)

    seq_mat,wtrow = numerics.dataset2mutarray(df.copy(),modeltype)
    #this is also an MCMC routine, do the same as above.
    if initialize == 'rand':
        if modeltype == 'MAT':
            emat_0 = utils.RandEmat(len(df[seq_col_name][0]),len(seq_dict))
        elif modeltype == 'NBR':
            emat_0 = utils.RandEmat(len(df['seq'][0])-1,len(seq_dict))
    emat = MaximizeMI_memsaver(
            seq_mat,df.copy(),emat_0,wtrow,db=db,iteration=iteration,burnin=burnin,
            thin=thin,runnum=runnum,verbose=verbose,temp=temp)

    emat_typical = gauge.fix_matrix(np.transpose(emat))

    em = pd.DataFrame(emat_typical)
    em.columns = val_cols
    #add position column
    pos = pd.Series(range(start,start + len(df[seq_col_name][0])),name='pos')
    output_df = pd.concat([pos,em],axis=1)

    # Validate model and return
    output_df = qc.validate_model(output_df,fix=True)
    return output_df

#run full script
output_df = main()
