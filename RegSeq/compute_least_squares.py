#Import basic stuff
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import linear_model

#import the custom analysis software
import scipy as sp
import plot_informationfootprint as pli
import seaborn as sns
import sys

#The input file will be input via sys.argv. The arguments already
#1: input file name
#2: output information footprint file name

inputname = sys.argv[1]
outputname = sys.argv[2]
gene = sys.argv[3]

def Compute_Least_Squares(raveledmat,batch):
    '''this is a simple linear regression function that will return the coefficients from the regression.
    In this case each coefficient represents the effect on gene expression of mutation of the corresponding base.
    We use sklearn to do this regression.'''
    clf = linear_model.LinearRegression()
    clf.fit(raveledmat,batch)
    emat = clf.coef_
    return emat

#We handle most data using the Pandas package, we will load in the target data set now.
df = pd.io.parsers.read_csv(inputname,delim_whitespace=True)

#we have a file with all wild type sequences for our genes. We load it in now
genedf = pd.io.parsers.read_csv('../data/test_data/wtsequences.csv')

#we extract the wild type sequence for aphA from that file.
wt = str(genedf.loc[genedf['name'] == gene,'geneseq'].tolist()[0])

#we convert the wild type sequence to a list.
wtlist = np.array(list(wt))

#some basic parameters of our sequences

#barcode length
taglength = 20
#total promoter length
seqlength = len(df['seq'][0]) - taglength #160 bp
#we create dictionaries that relate A,C,G,T to the number 1,2,3,4
seq_dict,inv_dict = pli.choose_dict('dna')

'''we initialize our array where we parameterize the sequence. There is one entry per base pair which
is equal to 1 for mutated, or 0 for wild type.'''
all_mutarr = np.zeros((len(df.index),seqlength))

#We will now parameterize our sequences
for i,row in df.iterrows():
    s = np.array(list(row['seq']))
    #clip off any sequence past the 160 bp mutated sequence length.
    s_clipped = s[:seqlength]
    #determine which bases are mutated
    all_mutarr[i,:seqlength] = (wtlist != s_clipped)

#We will use the ratio of mRNA counts to DNA counts to regress against. We add a pseudocount of 1.
thetarget = np.array((df['ct_0']+1)/(df['ct_1']+1))
#We will center the mean to be 0.
thetarget = thetarget - np.mean(thetarget)

#We then need to fit the effect of mutation from the data.
#For illustration purposes we now show how this can be done using linear regression.

emat = Compute_Least_Squares(all_mutarr,thetarget)

#output results
np.savetxt(outputname,emat)
