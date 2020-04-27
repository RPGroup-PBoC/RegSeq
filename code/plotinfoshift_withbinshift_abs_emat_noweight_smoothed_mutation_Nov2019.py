
# coding: utf-8

# In[2]:

from __future__ import division
import os
import glob
import pickle
import re

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import the project utils
import sys

# Import matplotlib stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Logo-generating module
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from Bio import Seq
import sys
import scipy as sp
from scipy import stats
import scipy.stats

#we have 5 input arguments using sys.argv
#1) the input mutational model
#2) whether to clip off the back 20 base pairs
#3) whether to invert the mutational model
#4) file location to save information footprint eps file
#5) file location to save information footprint text file


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

    if modeltype == 'NBR' or modeltype == 'PAIR':
        seq_dict = {
            ''.join([inv_dict[i],inv_dict[z]]):i*len(seq_dict)+z 
            for i in range(len(seq_dict)) for z in range(len(seq_dict))}
        inv_dict = {seq_dict[i]:i for i in seq_dict.keys()}
    return seq_dict,inv_dict

def sliding_window(y,windowsize=3):
    '''We apply a 3 bp sliding window to results.'''
    out_vec = np.zeros_like(y)
    for i in range(len(y)-windowsize):
        out_vec[i] = np.sum(y[i:i+windowsize])/windowsize
    return out_vec

def calcinfo(prob,bg):
    '''Calculates mutual information from a probability matrix'''
    prob = prob + 1e-7
    bg = bg + 1e-7
    return np.sum(prob*np.log2(prob/bg))

def effect_df_to_prob_df(effect_df,bg_df,beta):
    prob_df = effect_df.copy()
    vals = effect_df.values
    vals -= vals.min(axis=1)[:,np.newaxis]
    weights = np.exp(-beta*vals)*bg_df.values
    prob_df.loc[:,:] = weights/weights.sum(axis=1)[:,np.newaxis]
    return prob_df

#determine the relative probabilities of being mutated/being wt. In reality this
#is about 10 percent towards being mutated. However, to control for possible
#differing mutation rates, we will just arbitrarily set the ratio to be 50/50

background_array =pd.DataFrame( [[.5,.5]])

energy_df = pd.read_csv(sys.argv[1],delim_whitespace=True)

val_cols = ['val_wt','val_mut']
energyarr = np.array(energy_df[val_cols]).T

seq_dict,inv_dict = choose_dict('dna')


#Create a matrix and use it to see if we need to invert all matrix values. We
#will also use this to set color type of plot.
y_sub = np.array(energyarr[1,:])

if y_sub.sum() > 0:
    y_sub = y_sub*-1

if sys.argv[3] == 'invert':
    y_sub = y_sub*-1
y_sub_smoothed = sliding_window(y_sub)
abs_sub = np.abs(y_sub_smoothed)
maxval = np.max(abs_sub)
y_sub_normed = y_sub_smoothed/maxval/2 + 0.5



colorinputs = np.zeros((160))
for i in range(160):
    if y_sub_smoothed[i] < 0:
        colorinputs[i] = 0.0
    else:
        colorinputs[i] = 1.0


if sys.argv[2] == 'clip':
    total_length = len(energy_df.index)
    
    energy_df = energy_df.loc[:total_length-21,:]



energy_df = energy_df[['val_wt','val_mut']]

energy_df_scaled = energy_df

background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['val_wt','val_mut'])
emat_min = -2
emat_max = 2
mid_val=0

#calculate the probability matrices

mutlogo = effect_df_to_prob_df(energy_df_scaled,background_df,1)
mut2logo = effect_df_to_prob_df(-1*energy_df_scaled,background_df,1)
mutarr = np.array(mutlogo).T
mut2arr = np.array(mut2logo).T


#calculate info

mutinfo = np.zeros(160)
mutbgprob = np.zeros(2)
for i in range(160):
    mutbgprob[0] = .5
    mutbgprob[1] = 1-mutbgprob[0]
    if y_sub[i] > 0:
        mutinfo[i] = calcinfo(mutarr[:,i],mutbgprob)
    else:
        mutinfo[i] = -1*calcinfo(mut2arr[:,i],mutbgprob)


tempoutdf = pd.DataFrame()
tempoutdf['pos'] = range(160)
tempoutdf['info'] = sliding_window(mutinfo)
pd.set_option('max_colwidth',int(1e8))
    
tempoutdf.to_string(
        open(sys.argv[5],'w'), index=False,col_space=10)

smoothinfo = sliding_window(np.abs(mutinfo))

matplotlib.rc('xtick', labelsize=40) 
matplotlib.rc('ytick', labelsize=40) 
shiftcolors = plt.cm.bwr(colorinputs)
fig,ax = plt.subplots(figsize=(50,10))
ax.bar(range(160),np.abs(smoothinfo),color=shiftcolors)
plt.savefig(sys.argv[4],format='eps')

