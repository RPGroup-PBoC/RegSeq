import pandas as pd
import os
import matplotlib.pyplot as plt
import sys


#This is a helper script that runs along with matchdatasets_all.py

#The inputs here are the fastq file name for RNA and plasmid
name = sys.argv[1]
nameplas = sys.argv[2]

#read in file names
df = pd.io.parsers.read_csv(name,header=None)

#clip out only the sequences (not quality scores)
df = df.loc[1::4]

#there are 18 specific barcodes in our sequences, with different tag labels.
groupnum = int(name[-3:])

#get tags
tags = df[0].str.slice(-41,-21)

#count freuqncies of tags
tagcounts = tags.value_counts()

#do same things for plasmids
dfplas = pd.io.parsers.read_csv(nameplas,header=None)

dfplas = dfplas.loc[1::4]

groupnum = int(name[-3:])

tagsplas = dfplas[0].str.slice(-41,-21)

tagcountsplas = tagsplas.value_counts()


#load in the file that translates between barcode and the group number
genecodes = pd.io.parsers.read_csv('genetogroupnum')


#if we only want to do a partial match of the 18 groups, we check that here.
genestodo = list(genecodes.loc[genecodes['pnum'] == groupnum,'genename'])


def comb_tag(s):
    '''For ease of fitting, we append the tag sequence to the end of the dataset'''
    return s['seq'] + s['tag']

def format_string(x):
    '''This function is crucial for getting out ability to save the datasets with 
    proper float formatting'''
    return '%10.6f' %x
pd.set_option('max_colwidth',int(1e8))


for gene in genestodo:
    #load in the file that allow you to match tag sequence to promoter sequence
    # we have previously sequenced the total library.
    tagkeyname = '../107gene/' + gene + 'gooddataset_nslater_with_large'
    tagkey = pd.io.parsers.read_csv(tagkeyname,delim_whitespace=True)
    tagkey = tagkey.set_index('tag')

    #get counts for both mRNA and plasmid counts.
    tempdf = tagkey.loc[tagcounts.copy().index]

    tempdf['ct_1'] = tagcounts.copy()
    tempdf = tempdf.dropna()

    c = tagkey.loc[tagcountsplas.copy().index]
    c['ct_0'] = tagcountsplas.copy()
    c = c.dropna()

    outdf = pd.concat([tempdf,c],axis=0)

    #if there arent any counts then there will be a NaN in the dataset, we need to replace this with 0
    outdf = outdf.fillna(0)



    outdf = outdf[['ct_0','ct_1','seq']]


    outdf['ct'] = outdf[['ct_0','ct_1']].sum(axis=1)


    outdf = outdf.reset_index()

    outdf = outdf.rename(columns={'index':'tag'})

    outdf['seqall'] = outdf.apply(comb_tag,axis=1)

    outdf2 = outdf[['ct','ct_0','ct_1','seqall']]

    outdf2 = outdf2.rename(columns={'seqall':'seq'})

    temp = outdf2.groupby(by='seq').sum()

    temp = temp.reset_index()

    outdf3 = temp[['ct','ct_0','ct_1','seq']]
    #output the dataset
    outdf3.to_string(open(gene + sys.argv[3] + 'dataset_alldone_with_large','w'), index=False,col_space=10,float_format=format_string,justify='left')

