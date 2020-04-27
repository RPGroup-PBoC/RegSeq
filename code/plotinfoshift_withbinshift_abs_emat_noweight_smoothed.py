
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
sys.path.insert(0, '/media/bill/New_Volume/20171214RNAseq/plot/')
import NB_sortseq_utils as utils

# Import matplotlib stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Logo-generating module
import anylogo
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from Bio import Seq
import sys
import scipy as sp
from scipy import stats
import scipy.stats

gc = .508
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])


genelabel = sys.argv[1]

# In[3]:

genedf = pd.io.parsers.read_csv('/home/bill/next100genes/compedgenesv4.csv')

am = str(genedf.loc[genedf['name'] == genelabel,'geneseq'].tolist()[0])
#genedf = pd.io.parsers.read_csv('/home/bill/Downloads/gene_seqsv3(1)',delim_whitespace=True)
#print(genedf)
#am = str(genedf.loc[genedf['gene'] == genelabel,'seq'].tolist()[0])


#df1 = pd.io.parsers.read_csv('5760' + genelabel + '_goodlength_df',delim_whitespace=True)


# In[4]:

#df2 = pd.io.parsers.read_csv('../all_5759/5759' + genelabel + '_goodlength_df',delim_whitespace=True)


# In[5]:


#fr = np.load(sys.argv[2])

# In[59]:

#bad_shared2 = pd.io.parsers.read_csv('../all_5859/5859' + genelabel + '_goodlength_df_badtag_shared',delim_whitespace=True)


# In[60]:






# In[180]:

#now create logo using anylogo

energy_df = pd.read_csv(sys.argv[3],delim_whitespace=True)

val_cols = ['val_A','val_C','val_G','val_T']
energyarr = np.array(energy_df[val_cols]).T

def seq2mat(seq,seq_dict):
    
    mat = sp.zeros((len(seq_dict),len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat

def choose_dict(dicttype,modeltype='MAT'):
    
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

seq_dict,inv_dict = choose_dict('dna')

wt_mat = seq2mat(am,seq_dict)


def sliding_window(y,windowsize=3):
    out_vec = np.zeros_like(y)
    for i in range(len(y)-windowsize):
        out_vec[i] = np.sum(y[i:i+windowsize])/windowsize
    return out_vec
# In[173]:


wt_val = energyarr[:,:160]*wt_mat
submat = energyarr[:,:160] - wt_val.sum(axis=0)
y_sub = -1*submat.sum(axis=0)/3

if y_sub.sum() > 0:
    y_sub = y_sub*-1

if sys.argv[5] == 'invert':
    #energy_df = energy_df*-1
    y_sub = y_sub*-1
y_sub_smoothed = sliding_window(y_sub)
abs_sub = np.abs(y_sub_smoothed)
maxval = np.max(abs_sub)
y_sub_normed = y_sub_smoothed/maxval/2 + 0.5

thestats = np.zeros(160)
for i in range(160):
    q = submat[:,i]
    outq = q[np.nonzero(q)[0]]
    thestats[i] = scipy.stats.ttest_1samp(outq, 0)[1]
#if fr[0] > fr[499]:
#    y_sub_normed = y_sub_normed*-1

colorinputs = np.zeros((160))
for i in range(160):
    if y_sub_smoothed[i] < 0:
        colorinputs[i] = 0.0
    else:
        colorinputs[i] = 1.0


# energy_df = pd.read_csv(datadir2 + '20160710_purT_MG1655_M9glucose_na_mut1_4bins_pymc_rnap_MCMC_114_thermo.csv')
energy_df = energy_df.rename(columns={'val_A':'A','val_C':'C','val_G':'G','val_T':'T'})


#if designated clip off tag
if sys.argv[4] == 'clip':
    total_length = len(energy_df.index)
    
    energy_df = energy_df.loc[:total_length-21,:]

#invert if necessary

# energy_df = pd.read_csv(datadir2 + '20160824_purT_MG1655deltapurR_M9glucose_adenine_mut1_4bins_RNAP_emat_mean.csv')

# energy_df = pd.read_csv(datadir2 + '20160710_purT_MG1655_M9glucose_na_mut1_4bins_pymc_rnap_MCMC_114_thermo.csv')

energy_df = energy_df[['A','C','G','T']]
seq = 'AAAGACACACGCAAACGTTTTCGTTTATACTG'

#energy_df_scaled = -utils.estimate_scalefactor(np.array(energy_df).T)*energy_df.copy()
energy_df_scaled = energy_df
#if sys.argv[4]:
#    energy_df_scaled = energy_df_scaled/float(sys.argv[4])
# # scale factor for purR determined by MCMC
# scalefactor = 8.23

# for plotting logo and emat
#scalefactor = -10

#energy_df_scaled = energy_df.copy()
#energy_df_scaled = scalefactor * energy_df_scaled[['A','C','G','T']]
# emat_scaled = scalefactor * utils.zero_matrix_WT(np.array(energy_df.T), seq)
# create background nucleotide frequencies dataframe
background_df = pd.DataFrame(pd.np.tile(background_array,
                    (len(energy_df_scaled), 1)), columns=['A','C','G','T'])
# Set color scale - I want the colorbar to be symmetric
# emat_min=emat_scaled.min()
# emat_max=-emat_scaled.min()
# mid_val=0.0
emat_min = -2
emat_max = 2
mid_val=0


#relative_scale=1.5
#relative_spacing=.65
#emat_ymin = -2 * (relative_scale + relative_spacing)
#emat_ymax = -2 * relative_spacing
#yticks = np.linspace(emat_ymin, emat_ymax, 9)[[1, 3, 5, 7]]
#yticklabels = list('TGCA')

logo = anylogo.effect_df_to_prob_df(energy_df_scaled,background_df,1)
logo2 = anylogo.effect_df_to_prob_df(-1*energy_df_scaled,background_df,1)
logoarr = np.array(logo).T
logo2arr = np.array(logo).T

#collapse logos into wt and nonwt

nonwt_mat = -1*(wt_mat - 1)
mutlogo = np.zeros((2,160))
mut2logo = np.zeros((2,160))
mutlogo[0,:] = np.sum(wt_mat*logoarr,axis=0)
mut2logo[0,:] = np.sum(wt_mat*logo2arr,axis=0)
mutlogo[1,:] = np.sum(nonwt_mat*logoarr,axis=0)
mut2logo[1,:] = np.sum(nonwt_mat*logo2arr,axis=0)

def calcinfo(prob,bg):
    prob = prob + 1e-7
    bg = bg + 1e-7
    #return np.sum(prob_df.values*np.log2(prob_df.values),axis=1) - np.sum(bg_df.values*np.log2(bg_df.values),axis=1)
    return np.sum(prob*np.log2(prob/bg))

mutinfo = np.zeros(160)
mutbgprob = np.zeros(2)
for i in range(160):
    if am[i] in ['C','G']:
        mutbgprob[0] = .508/2
    else:
        mutbgprob[0] = .492/2
    mutbgprob[1] = 1-mutbgprob[0]
    if y_sub[i] > 0:
        mutinfo[i] = calcinfo(mutlogo[:,i],mutbgprob)
    else:
        mutinfo[i] = -1*calcinfo(mut2logo[:,i],mutbgprob)
tempoutdf = pd.DataFrame()
tempoutdf['pos'] = range(160)
tempoutdf['info'] = sliding_window(mutinfo)
pd.set_option('max_colwidth',int(1e8))
    
tempoutdf.to_string(
        open(sys.argv[7],'w'), index=False,col_space=10)

def get_prob_df_info(prob_df,bg_df):
    prob_df = prob_df + 1e-7
    bg_df = bg_df + 1e-7
    #return np.sum(prob_df.values*np.log2(prob_df.values),axis=1) - np.sum(bg_df.values*np.log2(bg_df.values),axis=1)
    return np.sum(prob_df.values*np.log2(prob_df.values/bg_df.values),axis=1)

infos = get_prob_df_info(logo,background_df)
infos2 = get_prob_df_info(logo2,background_df)
out_infos = np.zeros(160)
for i in range(160):
    if y_sub[i] > 0:
        out_infos[i] = infos[i]
    else:
        out_infos[i] = infos2[i]*-1


#negindexes = np.nonzero(mutinfo<0)[0]
#posindexes = np.nonzero(mutinfo>0)[0]
smoothinfo = sliding_window(np.abs(mutinfo))
#smoothinfopos = smoothinfo.copy()
#smoothinfopos[negindexes] = 0
#smoothinfoneg = smoothinfo.copy()
#smoothinfoneg[posindexes] = 0

#mutinfopos = mutinfo.copy()
#mutinfopos[negindexes] = 0
#mutinfoneg = mutinfo.copy()
#mutinfoneg[posindexes] = 0
#fig,ax = plt.subplots()

#ax.bar(range(160),smoothinfoneg,color='r')
#ax.bar(range(160),smoothinfopos)
#plt.show()
#ax.bar(range(160),out_infos)
#plt.show()

#fig,ax = plt.subplots()
#ax.bar(range(160),sliding_window(mutinfo))
#plt.show()

matplotlib.rc('xtick', labelsize=40) 
matplotlib.rc('ytick', labelsize=40) 
shiftcolors = plt.cm.bwr(colorinputs)
fig,ax = plt.subplots(figsize=(50,10))
ax.bar(range(160),np.abs(smoothinfo),color=shiftcolors)
#ax.bar(range(160),np.abs(smoothinfopos))
plt.savefig(sys.argv[6],format='eps')
plt.show()
'''

def seq2mat(seq,seq_dict):
    
    mat = sp.zeros((len(seq_dict),len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat

def choose_dict(dicttype,modeltype='MAT'):
    
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

seq_dict,inv_dict = choose_dict('dna')

wt_mat = seq2mat(am,seq_dict).T

print wt_mat

# In[173]:


wt_val = np.array(energy_df_scaled)*wt_mat
print wt_val

# In[ ]:




# In[174]:

submat = energy_df_scaled.T - wt_val.sum(axis=1).T




# In[176]:

y_sub = submat.sum(axis=0)/3


ax.bar(range(160),y_sub)
'''
plt.show()

