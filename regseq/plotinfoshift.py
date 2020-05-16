
from __future__ import division
import os
import glob
import pickle
import re
import sys

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import matplotlib stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from Bio import Seq
import sys
import scipy as sp
from scipy import stats
import scipy.stats

#E coli gc content for background frequencies.
gc = .508
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])


#input the gene name, so we can get the wt sequence.
genelabel = sys.argv[1]
genedf = pd.io.parsers.read_csv('/home/bill/next100genes/compedgenesv4.csv')
am = str(genedf.loc[genedf['name'] == genelabel,'geneseq'].tolist()[0])
length_wt = len(am)

#load in the energy matrix to use for generating an information footprint
energy_df = pd.read_csv(sys.argv[2],delim_whitespace=True)

#convert to a numpy array
val_cols = ['val_A','val_C','val_G','val_T']
energyarr = np.array(energy_df[val_cols]).T

def seq2mat(seq,seq_dict):
    '''convert a sequence to a numeric representation'''
    mat = sp.zeros((len(seq_dict),len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat

def choose_dict(dicttype,modeltype='MAT'):
    '''create a dictionary that matches the bases A C G T to the number 0 1 2 3'''
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

def calcinfo(prob,bg):
    prob = prob + 1e-7
    bg = bg + 1e-7
    return np.sum(prob*np.log2(prob/bg))

def sliding_window(y,windowsize=3):
    '''average mutual information footprint with neighbors'''
    out_vec = np.zeros_like(y)
    for i in range(len(y)-windowsize+1):
        out_vec[i] = np.sum(y[i:i+windowsize])/windowsize
    return out_vec

def effect_df_to_prob_df(effect_df,bg_df,beta):
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

#load in sequence dictionary and get matrix representation of wt seq.
seq_dict,inv_dict = choose_dict('dna')
wt_mat = seq2mat(am,seq_dict)

#we now get a matrix that contains wt vs non-wt entries (averaged).
wt_val = energyarr[:,:length_wt]*wt_mat
submat = energyarr[:,:length_wt] - wt_val.sum(axis=0)
y_sub = -1*submat.sum(axis=0)/3

#make sure wild type energy is negative, if not flip the entries.
if y_sub.sum() > 0:
    y_sub = y_sub*-1

#in case you want to manually override this, you can pass this 
if sys.argv[4] == 'invert':
    y_sub = y_sub*-1


#smooth the y_sub array by averaging with neighbors
y_sub_smoothed = sliding_window(y_sub)

#we color code output files to show repressor-like or activator-like.
colorinputs = np.zeros((length_wt))
for i in range(length_wt):
    if y_sub_smoothed[i] < 0:
        colorinputs[i] = 0.0
    else:
        colorinputs[i] = 1.0

#rename columns
energy_df = energy_df.rename(columns={'val_A':'A','val_C':'C','val_G':'G','val_T':'T'})


#if designated clip off last 20 bases.
if sys.argv[3] == 'clip':
    total_length = len(energy_df.index)
    
    energy_df = energy_df.loc[:total_length-21,:]

#only take sequence columns
energy_df = energy_df[['A','C','G','T']]
energy_df_scaled = energy_df

background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])
emat_min = -2
emat_max = 2
mid_val=0


#get probability array
prob = effect_df_to_prob_df(energy_df_scaled,background_df,1)
prob2 = effect_df_to_prob_df(-1*energy_df_scaled,background_df,1)
probarr = np.array(prob).T
prob2arr = np.array(prob2).T

nonwt_mat = -1*(wt_mat - 1)
mutprob = np.zeros((2,length_wt))
mut2prob = np.zeros((2,length_wt))
mutprob[0,:] = np.sum(wt_mat*probarr,axis=0)
mut2prob[0,:] = np.sum(wt_mat*prob2arr,axis=0)
mutprob[1,:] = np.sum(nonwt_mat*probarr,axis=0)
mut2prob[1,:] = np.sum(nonwt_mat*prob2arr,axis=0)


#calculate information
mutinfo = np.zeros(length_wt)
mutbgprob = np.zeros(2)
for i in range(length_wt):
    if am[i] in ['C','G']:
        mutbgprob[0] = .508/2
    else:
        mutbgprob[0] = .492/2
    mutbgprob[1] = 1-mutbgprob[0]
    if y_sub[i] > 0:
        mutinfo[i] = calcinfo(mutprob[:,i],mutbgprob)
    else:
        mutinfo[i] = -1*calcinfo(mut2prob[:,i],mutbgprob)

#create output dataframe
tempoutdf = pd.DataFrame()
tempoutdf['pos'] = range(1,length_wt-1)
tempoutdf['info'] = sliding_window(np.abs(mutinfo))[:-2]
pd.set_option('max_colwidth',int(1e8))
    
tempoutdf.to_string(
        open(sys.argv[6],'w'), index=False,col_space=10)



smoothinfo = sliding_window(np.abs(mutinfo))

matplotlib.rc('xtick', labelsize=40) 
matplotlib.rc('ytick', labelsize=40) 
shiftcolors = plt.cm.bwr(colorinputs)
fig,ax = plt.subplots(figsize=(50,10))
ax.bar(range(length_wt),np.abs(smoothinfo),color=shiftcolors)
plt.savefig(sys.argv[5],format='eps')


