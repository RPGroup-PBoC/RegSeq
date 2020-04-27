import pandas as pd
import os
import matplotlib.pyplot as plt
import sys

#This script will process two input sequencing files for mRNA
#and DNA into a data set for all genes 

# The four inline arguments passed to this script are
# 1: file name of mRNA sequencing .fastq file.
# 2: file name of DNA sequencing .fastq file.
# 3: output name prefix. The output name for each gene file will be given
# as gene_name + this_input + 'dataset'
# 4  group number. Genes are separated into 18 groups, labeled 101 to 118.
# Only those genes in the given group number will have their datasets generated
# by this script. The associated group number to gene association is given
# by the genetogroupnum file.

name = sys.argv[1]
nameplas = sys.argv[2]

groupnum = int(sys.argv[5])

#define needed functions

def comb_tag(s):
    '''function to combine mutated sequence with barcode'''
    return s['seq'] + s['tag']

def format_string(x):
    'basic function to format output string'
    return '%10.6f' %x

#set no maximum length on output column size.
pd.set_option('max_colwidth',int(1e8))


#load in dataframe version of mRNA sequences.
df = pd.io.parsers.read_csv(name,header=None)

#extract the sequences from the fastq format.
df = df.loc[1::4]

#we will select out the barcodes from each sequence. They will be located
#from -41 to -21 bases from the end of the sequence.

tags = df[0].str.slice(-41,-21)

#we will get the numbers of each barcode.

tagcounts = tags.value_counts()

#We will now preform an identical procedure for the DNA sequences.

dfplas = pd.io.parsers.read_csv(nameplas,header=None)

dfplas = dfplas.loc[1::4]

tagsplas = dfplas[0].str.slice(-41,-21)

tagcountsplas = tagsplas.value_counts()


#we will get the genes for the associated group number. This is generally 6
#genes.

#load in key for group number for each gene

genecodes = pd.io.parsers.read_csv('../data/test_data/genetogroupnum')

#use group code to find the genes we need to make datasets for.

genestodo = list(genecodes.loc[genecodes['pnum'] == groupnum,'genename'])


#loop through the genes for the appropriate group number.
for gene in genestodo:
    #load in the file that relates barcode to mutated sequence.
    tagkeyname = sys.argv[3]
    tagkey = pd.io.parsers.read_csv(tagkeyname,delim_whitespace=True)
    #reset the barcode to be the pandas index.
    tagkey = tagkey.set_index('tag')

    #make a dataframe that has the number of counts for the assocated barcode
    #for the mRNA sequencing. After this we will do the same for the DNA plasmid
    #sequencing then combine them.

    #make dataframe with mutated sequence and barcode

    tempdf = tagkey.reindex(tagcounts.copy().index)

    #assign sequencing counts based on mRNA sequencing file.

    tempdf['ct_1'] = tagcounts.copy()
    tempdf = tempdf.dropna()

    #we now do the same thing for the DNA plasmid sequencing.

    c = tagkey.reindex(tagcountsplas.copy().index)
    c['ct_0'] = tagcountsplas.copy()
    c = c.dropna()

    #combine the dataframes to get sequencing counts for mRNA and DNA

    outdf = pd.concat([tempdf,c],axis=0,sort=True)

    #any entries with no counts for mRNA or DNA will be automatically filled with
    #np.NaN, we want to replace these with 0 counts.

    outdf = outdf.fillna(0)
    
    #remove unnecessary columns.
    outdf = outdf[['ct_0','ct_1','seq']]

    #add in a total counts column.
    outdf['ct'] = outdf[['ct_0','ct_1']].sum(axis=1)

    #we want a column for the barcode. barcode is currently the index. we
    #need to reset this.

    outdf = outdf.reset_index()

    #rename barcode column.

    outdf = outdf.rename(columns={'index':'tag'})

    #combine the mutated sequence and barcode into one sequence that is 180 bp
    #long.

    outdf['seqall'] = outdf.apply(comb_tag,axis=1)

    #reorder columns.

    outdf2 = outdf[['ct','ct_0','ct_1','seqall']]

    #rename columns.

    outdf2 = outdf2.rename(columns={'seqall':'seq'})

    temp = outdf2.groupby(by='seq').sum()

    temp = temp.reset_index()

    outdf3 = temp[['ct','ct_0','ct_1','seq']]

    #output the final dataset for associated gene.
    outdf3.to_string(open(sys.argv[4],'w'), index=False,col_space=10,float_format=format_string,justify='left')


