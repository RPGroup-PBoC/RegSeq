from .utils import choose_dict, seq2mat
import os

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import pandas as pd
import numpy as np
import scipy as sp

import logomaker


def sliding_window(y,windowsize=3):
    """Average information values with neighbors."""
    
    if windowsize%2 != 1:
        raise RuntimeError("Window size has to be odd.")
        
    cut = int((windowsize - 1) / 2)
    out_vec = np.zeros(len(y)-2*cut)
    for i in range(cut, len(y)-cut):
        out_vec[i-cut] = np.sum(y[i-cut:i+cut+1])/windowsize
    return out_vec


def calcinfo(prob, bg):
    '''Calculate mutual information from a probability matrix.'''
    prob = prob + 1e-7
    bg = bg + 1e-7
    return np.sum(prob*np.log2(prob/bg))


def convert_to_df(em):
    outdf = pd.DataFrame()
    outdf['val_mut'] = em
    outdf['val_wt'] = 0
    outdf['pos'] = range(len(outdf.index))
    outdf = outdf[['pos','val_wt','val_mut']]
    return outdf


def effect_df_to_prob_df(effect_df, bg_df, beta):
    prob_df = effect_df.copy()
    vals = effect_df.values
    vals -= vals.min(axis=1)[:,np.newaxis]
    weights = np.exp(-beta*vals)*bg_df.values
    prob_df.loc[:,:] = weights/weights.sum(axis=1)[:,np.newaxis]
    return prob_df


def get_prob_df_info(prob_df,bg_df):
    prob_df = prob_df + 1e-7
    bg_df = bg_df + 1e-7
    return np.sum(prob_df.values*np.log2(prob_df.values/bg_df.values),axis=1)

#determine the relative probabilities of being mutated/being wt. In reality this
#is about 10 percent towards being mutated. However, to control for possible
#differing mutation rates, we will just arbitrarily set the ratio to be 50/50

def footprint(inarr, for_clip=None, seqlength=160, windowsize=3):
    """
    Compute information footprint from expression changes per position.
    
    Parameters
    ----------
    inarr : numpy array
        Change of expression per position.
    for_clip : str, default None
        if "clip", clip of last 21 bases
    seqlength : int, default 160
        Length of generated sequences without barcodes
    windowsize : int, default 3
        number of bases used to averagee footprint (has to be odd)
        
    Returns
    -------
    smoothinfo : numpy array
        mutual information at every base
    shiftcolors : numpy array
        colors for barplot
    signs : numpy array
        List of 1 and -1, saying if mutation at a base increased or decreased expression
    """
    

    background_array =pd.DataFrame([[.5,.5]])

    energy_df = convert_to_df(inarr)

    val_cols = ['val_wt','val_mut']
    energyarr = np.array(energy_df[val_cols]).T

    seq_dict,inv_dict = choose_dict('dna')

    #Create a matrix and use it to see if we need to invert all matrix values. We
    #will also use this to set color type of plot.
    y_sub = np.array(energyarr[1,:])

    if y_sub.sum() > 0:
        y_sub = y_sub*-1

    y_sub_smoothed = sliding_window(y_sub, windowsize=windowsize)
    abs_sub = np.abs(y_sub_smoothed)
    maxval = np.max(abs_sub)
    y_sub_normed = y_sub_smoothed/maxval/2 + 0.5
    colorinputs = np.zeros((len(y_sub_smoothed)))
    for i in range(seqlength - windowsize):
        if y_sub_smoothed[i] < 0:
            colorinputs[i] = 0.0
        else:
            colorinputs[i] = 1.0

    if for_clip == 'clip':
        total_length = len(energy_df.index)
        energy_df = energy_df.loc[:total_length - 21, :]

    energy_df = energy_df[['val_wt', 'val_mut']]

    energy_df_scaled = energy_df

    background_df = pd.DataFrame(np.tile(background_array,
                        (len(energy_df_scaled), 1)), columns=['val_wt', 'val_mut'])
    emat_min = -2
    emat_max = 2
    mid_val=0

    #calculate the probability matrices
    mutlogo = effect_df_to_prob_df(energy_df_scaled, background_df, 1)
    mut2logo = effect_df_to_prob_df(-1*energy_df_scaled, background_df, 1)
    mutarr = np.array(mutlogo).T
    mut2arr = np.array(mut2logo).T


    mutinfo = np.zeros(seqlength)
    mutbgprob = np.zeros(2)
    
    # store sign of expression change
    exp_change = np.ones(seqlength)
    
    for i in range(seqlength):
        mutbgprob[0] = .5
        mutbgprob[1] = 1-mutbgprob[0]
        
        mutinfo[i] = calcinfo(mutarr[:,i],mutbgprob)
        if y_sub[i] < 0:
            exp_change[i] *= -1
    
    cut = int((windowsize - 1)/2)
    smoothinfo = sliding_window(np.abs(mutinfo),windowsize=windowsize)
    shiftcolors = plt.cm.bwr(colorinputs)
    
    if cut != 0:
        return smoothinfo, shiftcolors, exp_change[cut:-cut]
    else:
        return smoothinfo, shiftcolors, exp_change


def get_info(df, gc=.508):
    '''Find the total information content of a binding site.'''
    

    background_array = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])
    
    #add in small value to make sure no probabilities are exactly zero.
    df = df + 1e-9
    
    return np.sum(df.values * np.log2(df.values/background_array))


def get_beta_for_effect_df(
    effect_df,
    target_info,
    min_beta=.001,
    max_beta=100,
    num_betas=1000):
    '''Find the appropriate scaling factor for displaying sequence
    logos. '''
    betas = np.exp(np.linspace(np.log(min_beta),np.log(max_beta),num_betas))
    infos = np.zeros(len(betas))
    for i, beta in enumerate(betas):
        prob_df = logomaker.transform_matrix(df=beta*effect_df,from_type='weight',to_type='probability')
        infos[i] = get_info(prob_df)
    i = np.argmin(np.abs(infos-target_info))
    beta = betas[i]
    return beta


def emat_to_information(
        file, 
        wildtype_file='../data/prior_designs/wtsequences.csv', 
        clip=False, 
        invert=False,
        old_format=False,
        gene=None
    ):
    """
    
    Parameters
    ----------
    file: str
        Energy matrix file
    wildtypefile: str, default
        File with wild type sequences
    clip : boolean, default False
        If True, clip off last 20 bases.
    invert : boolean, default True
        If True, flip signs for binding energy, if wildtype energy is positive.
    old_format : boolean, default False
        Determines if files are of old format
    gene : None
        For old file names, name of gene has to be given
    """
    
    #input the gene name, so we can get the wt sequence.
    if gene == None:
        gene = file.split("/")[-1].split("_")[0]
    elif type(gene) != str:
        raise RuntimeError("gene has to be of type string.")

    genedf = pd.read_csv(wildtype_file)
    am = str(genedf.loc[genedf['name'] == gene, 'geneseq'].tolist()[0])
    length_wt = len(am)
    
    # Load in the energy matrix

    energy_df = pd.read_csv(file, delim_whitespace=old_format)

    #convert to a numpy array
    val_cols = ['val_A','val_C','val_G','val_T']
    energyarr = np.array(energy_df[val_cols]).T
        
    #load in sequence dictionary and get matrix representation of wt seq.
    seq_dict,inv_dict = choose_dict('dna')
    wt_mat = seq2mat(am, seq_dict)

    #we now get a matrix that contains wt vs non-wt entries (averaged).
    wt_val = energyarr[:,:length_wt]*wt_mat
    submat = energyarr[:,:length_wt] - wt_val.sum(axis=0)
    y_sub = -1*submat.sum(axis=0)/3

    #make sure wild type energy is negative, if not flip the entries.
    if invert and y_sub.sum() > 0:
        y_sub = y_sub*-1

    return y_sub
    