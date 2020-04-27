
# coding: utf-8

# In[2]:

#import this just for backwards compatability with python2.7
from __future__ import division
import os

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import the project utils
import sys


# Import matplotlib stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Logo-generating module
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import sys
import scipy as sp




#define functions
def seq2mat(seq,seq_dict):
    '''Function that tells us how to convert a sequence into a matrix representation
    of shape (4,L)'''
    mat = sp.zeros((len(seq_dict),len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat

def choose_dict(dicttype,modeltype='MAT'):
    '''Creates a necessary tool to convert from bp to an index'''
    if dicttype == 'dna':
        seq_dict = {'A':0,'C':1,'G':2,'T':3}
        inv_dict = {0:'A',1:'C',2:'G',3:'T'}
    elif dicttype == 'rna':
        seq_dict = {'A':0,'C':1,'G':2,'U':3}
        inv_dict = {0:'A',1:'C',2:'G',3:'U'}
    elif dicttype == 'protein':
        seq_dict = {
            '*':0,'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,
            'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20}
        inv_dict = {v:k for k,v in seq_dict.items()}
    else:
        raise SortSeqError('Unkonwn dicttype: %s'%dicttype)
    return seq_dict,inv_dict

def sliding_window(y,windowsize=3):
    '''Function that averages the effect of bases in a 3 base pair window.
    This will be used on the information values at the end of the script'''
    out_vec = np.zeros((len(y)-windowsize))
    for i in range(len(y)-windowsize):
        out_vec[i] = np.sum(y[i:i+windowsize])/windowsize
    return out_vec

def calcinfo(prob,bg):
    '''Calculates mutual information from a probability matrix. We start
    from a unitless number representing an effective energy change
    upon mutation. We then convert to a probability matrix using a
    Boltzmann distribution. We can then use the below function to calculate
    the mutual information.'''
    prob = prob + 1e-7
    bg = bg + 1e-7
    return np.sum(prob*np.log2(prob/bg))

def convert_to_df(em):
    '''This function converts an numpy matrix format for mutation effect
    to a pandas data frame'''
    outdf = pd.DataFrame()
    outdf['val_mut'] = em
    outdf['val_wt'] = 0
    outdf['pos'] = range(len(outdf.index))
    outdf = outdf[['pos','val_wt','val_mut']]
    return outdf

def effect_df_to_prob_df(effect_df,bg_df,beta):
    '''This function converts effective energy matrices to a probability matrix.
    The probability matrix represents how likely that base would be to appear
    in a naturally evolved version of the promoter.'''
    prob_df = effect_df.copy()
    vals = effect_df.values
    vals -= vals.min(axis=1)[:,np.newaxis]
    weights = np.exp(-beta*vals)*bg_df.values
    prob_df.loc[:,:] = weights/weights.sum(axis=1)[:,np.newaxis]
    return prob_df

#determine the relative probabilities of being mutated/being wt. In reality this
#is about 10 percent towards being mutated. However, to control for possible
#differing mutation rates, we will just arbitrarily set the ratio to be 50/50

def main(inarr,for_clip=False,seqlength=160,for_invert=False):
    windowsize=3
    #base probabilites of being mutated vs. unmutated will be set to 50 percent.
    background_array =pd.DataFrame( [[.5,.5]])

    #convert the numpy matrix to a pandas dataframe
    energy_df = convert_to_df(inarr)

    #the columns in the pandas data frame corresponding to values.
    val_cols = ['val_wt','val_mut']
    
    #convert to a numpy matrix from the pandas dataframe
    energyarr = np.array(energy_df[val_cols]).T

    #get a sequence dictionary (for example A = 0, C = 1, etc)
    seq_dict,inv_dict = choose_dict('dna')


    #Create a matrix and use it to see if we need to invert all matrix values. We
    #will also use this to set color type of plot.
    y_sub = np.array(energyarr[1,:])

    #We want to fix the overall effect of a mutation to be negative.
    #so if the total effect is that mutation is positive, we will flip the
    #whole matrix.
    if y_sub.sum() > 0:
        y_sub = y_sub*-1

    #If we want to override this we can pass for_invert.
    if for_invert == 'invert':
        y_sub = y_sub*-1

    #we will smooth the matrix using the sliding window function.
    y_sub_smoothed = sliding_window(y_sub)
    
    #make sure everything is positive
    abs_sub = np.abs(y_sub_smoothed)

    #We will find the highest value of the array and use this value to set
    #the scale of the plots
    maxval = np.max(abs_sub)
    y_sub_normed = y_sub_smoothed/maxval/2 + 0.5

    #We need to color code the final plots depending on whether or not 
    #mutation tends to increase or decrease gene expression.
    colorinputs = np.zeros((seqlength))
    for i in range(seqlength-windowsize):
        if y_sub_smoothed[i] < 0:
            colorinputs[i] = 0.0
        else:
            colorinputs[i] = 1.0

    #If we wnat to remove the last 20 base pairs then we can pass in 'clip'
    if for_clip == 'clip':
        total_length = len(energy_df.index)
        energy_df = energy_df.loc[:total_length-21,:]

    #we only need to value columns, so we will remove the rest.
    energy_df_scaled = energy_df[['val_wt','val_mut']]

    #make an data frame out of the background values (50/50 probability for
    #mutated vs unmuated). It will be shaped like the input energy matrix (for
    #us that is generally 160 base pairs long.
    background_df = pd.DataFrame(pd.np.tile(background_array,
                        (len(energy_df_scaled), 1)), columns=['val_wt','val_mut'])

    emat_min = -2
    emat_max = 2
    mid_val=0

    #calculate the probability matrices. We will do separate ones for the bases
    #that have a negative effect on gene expression vs. a positive influence.

    mutlogo = effect_df_to_prob_df(energy_df_scaled,background_df,1)
    mut2logo = effect_df_to_prob_df(-1*energy_df_scaled,background_df,1)
    mutarr = np.array(mutlogo).T
    mut2arr = np.array(mut2logo).T

    #calculate info on each base in the sequence by looping. 
    mutinfo = np.zeros(seqlength)
    mutbgprob = np.zeros(2)
    for i in range(seqlength):
        mutbgprob[0] = .5
        mutbgprob[1] = 1-mutbgprob[0]
        if y_sub[i] > 0:
            mutinfo[i] = calcinfo(mutarr[:,i],mutbgprob)
        else:
            mutinfo[i] = -1*calcinfo(mut2arr[:,i],mutbgprob)

    #format output data frame, with a 'position' column and an info column.
    tempoutdf = pd.DataFrame()
    tempoutdf['pos'] = range(1,seqlength-windowsize+1)
    tempoutdf['info'] = sliding_window(mutinfo)

    #we will smooth the information footprint with the sliding window function
    #we first need to get the abolute value of the information footprint
    #because we used negative values to keep track of which bases had a negative
    #effect on gene expression.
    smoothinfo = sliding_window(np.abs(mutinfo),windowsize=windowsize)

    #generates colors for the bars in the output file. positive effects on mutation
    #get red bars and negative effects have blue bars. 
    shiftcolors = plt.cm.bwr(colorinputs)
    return np.abs(smoothinfo), shiftcolors
