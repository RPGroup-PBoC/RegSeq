import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp

from deprecated import deprecated

@deprecated(version='0.0.2', reason="This module is moved to utils.py")
def seq2mat(seq,seq_dict):

    mat = sp.zeros((len(seq_dict),len(seq)),dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp],i] = 1
    return mat

@deprecated(version='0.0.2', reason="This module is moved to utils.py")
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
        raise SortSeqError('Unknown dicttype: %s'%dicttype)

    if modeltype == 'NBR' or modeltype == 'PAIR':
        seq_dict = {
            ''.join([inv_dict[i],inv_dict[z]]):i*len(seq_dict)+z 
            for i in range(len(seq_dict)) for z in range(len(seq_dict))}
        inv_dict = {seq_dict[i]:i for i in seq_dict.keys()}
    return seq_dict,inv_dict

@deprecated(version='0.0.2', reason="This module is moved to utils.py")
def sliding_window(y,windowsize=3):
        out_vec = np.zeros_like(y)
        for i in range(len(y)-windowsize):
            out_vec[i] = np.sum(y[i:i+windowsize])/windowsize
        return out_vec

