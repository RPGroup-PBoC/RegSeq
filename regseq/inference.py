import numpy as np
import pandas as pd

from sklearn import linear_model
from sklearn.preprocessing import StandardScaler

import scipy.sparse
import scipy as sp

import pdb

from .utils import choose_dict

def least_squares(raveledmat, batch):
    """Linear regression of effects on gene expression of mutation of the
    corresponding base.
    Uses sklearn to do this regression."""
    
    clf = linear_model.LinearRegression()
    clf.fit(raveledmat, batch)
    emat = clf.coef_
    return emat


def lin_reg(inputname, outputname, wildtypefile='../data/prior_designs/wtsequences.csv'):

    # Load data
    df = pd.read_csv(inputname)
    
    # Load wild type sequences
    genedf = pd.read_csv('../data/prior_designs/wtsequences.csv')
    
    # Gene name
    gene = inputname.split('/')[-1].split('_')[0]
    
    # Extract the wild type sequence of gene of interest
    wt = str(genedf.loc[genedf['name'] == gene,'geneseq'].tolist()[0])

    # Convert to list.
    wtlist = np.array(list(wt))

    # Barcode length
    taglength = 20

    # Total promoter length
    seqlength = len(df['seq'][0]) - taglength #160 bp
    
    #we create dictionaries that relate A,C,G,T to the number 1,2,3,4
    seq_dict,inv_dict = choose_dict('dna')

    # Initialize array to paramaterize sequence. Mutations are denoted by 1
    all_mutarr = np.zeros((len(df.index),seqlength))

    # Parameterize sequences
    for i,row in df.iterrows():
        s = np.array(list(row['seq']))
        # Clip off any sequence past the 160 bp mutated sequence length.
        s_clipped = s[:seqlength]
        # Find mutations
        all_mutarr[i,:seqlength] = (wtlist != s_clipped)
    # IUse the ratio of mRNA counts to DNA counts to regress against. Add a pseudocount of 1.
    thetarget = np.array((df['ct_0']+1)/(df['ct_1']+1))
    
    # Center the mean
    thetarget = thetarget - np.mean(thetarget)
    # Fit mutation effects using linear regression
    emat = least_squares(all_mutarr, thetarget)

    #output results
    np.savetxt(outputname,emat)
    
    
