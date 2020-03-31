import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

#we use inline arguments using sys.argv to input file names for processing.

#The 1 input argument is the file name for the fastq sequencing.
#each sequence will include the 160 bp of mutated sequence, with the associated
#barcode (20bp) at the end of the sequence

#define sequence parameters These might need to be changed depending on your
#sequencing set up.

input_file_name = sys.argv[1]
barcode_length = 20
sequence_length = 160

#There will be some number of base pairs after the barcode in each sequence

#For our sequencing set up we have either 24 bp or 20 bp at the end of the
#sequences. We will set the correct parameters for this now.

#if total sequence length is 299 bp. then the toal trailing sequence length is
#24 bp, otherwise if the total length is 295 hen trailing sequence length is
#20 bp.
longth_trailing_length = 24
short_trailing_length = 20

#all sequences have 20 base pairs at the start of the sequence.

starting_length = 20


#argument 1 = file name for fastq file for sequencing
df = pd.io.parsers.read_csv(input_file_name,delim_whitespace=True,header=None)


#as it is a fastq file, the actual sequences will occur every 4 lines (starting
#with line index 1). The
#sequences will already have been filtered for quality score (phred > 20). 

#select only rows with sequences.
df = df.loc[1::4,0]


#we need to filter out incorrect sequence lengths to make sure they don't
#have improper insertions or deletions

#find all lengths
lengths = df.apply(len)

#find the correct length by finding the most common length
lengthsmax = list(lengths.value_counts().index)[0]
print('optimal length is ' + str(lengthsmax))

#find all sequences with incorrect length
badlength = (df.apply(len) != lengthsmax)
baddf = df.loc[badlength]

#find all sequences with correct length.
goodlength = (df.apply(len) == lengthsmax)
df = df.loc[goodlength]


#We will also be discluding those sequences with sequences with undetermined
#bases (shown as N's).

def check_N(s):
    if 'N' in s:
        return False
    else:
        return True

if lengthsmax == 295:
    sliceddf = df.str.slice(starting_length,-short_trailing_length)
elif lengthsmax == 299:
    sliceddf = df.str.slice(starting_length,-long_trailing_length)
else:
    print('Sequences not either 299 or 295 bp')
    raise

def stitch(s):
    #this function will combine the mutated sequence with barcode.
    return s['seq'] + s['tag']

#We will now construct a data frame with mutated sequence and barcode
tempdf = pd.DataFrame()
tempdf['seq'] =sliceddf.str.slice(0,sequence_length)
tempdf['tag'] =sliceddf.str.slice(-barcode_length,)

#stitch them together into one sequence (with barcode at end).
tempstitched = tempdf.apply(stitch,axis=1)

#remove all sequences with N's in them.
noN = tempstitched.apply(check_N)
tempstitched = tempstitched.loc[noN]

#these are the bad sequences.
badnoN = baddf.apply(check_N)
baddf = baddf.loc[badnoN]

#To determine which barcode is linked to which sequence we need to have only
#sequence matched to each sequence.


#first we check to see how many time each sequence/barcode combination is
#sequenced.
q2 = tempstitched.value_counts()

#we then drop those that only are sequenced 1 time. This will help to cut
#out some sequencing errors.
bigdf = q2.loc[q2 > 1]


#we will now return the data frame to have separate columns for sequence
#and barcode.
outbigdf = pd.DataFrame()

tempbigdf = pd.DataFrame(bigdf.index)

outbigdf['seq'] = tempbigdf[0].str.slice(0,sequence_length)
outbigdf['tag'] = tempbigdf[0].str.slice(-barcode_length,)

#we will now check to make sure that each barcode does not map incorrectly to
#multiple sequences

#count number of sequences associated with a given tag

q3 = outbigdf['tag'].value_counts()

temp = bigdf.reset_index()
outbigdf['ct'] = temp[0]

#this shows the number of sequences for a given tag
oq = outbigdf['tag'].value_counts()


#we will find those with exactly 1 sequence assocated with them.

goodtags = oq.loc[oq == 1].index

toomany = (oq > 1)

outbigdfindexed = outbigdf.set_index('tag')

toobigdf = outbigdfindexed.loc[toomany]

toobigdf = toobigdf.reset_index()

#because even for low sequencing error rates you could have a single error
#that could make it appear that a barcode maps to multiple sequences even
#if it doesn't. As a result, if the most common sequence is 99% of all sequences
#for a given barcode we will include it in the final mapping.

def findshare(s):
    'Finds what percent of all sequences are composed of a single sequence'
    ourmax = s['ct'].max()
    oursum = s['ct'].sum()
    return ourmax/oursum

def findmaxseq(s):
    'find most common sequence'
    ourmax = s['ct'].idxmax()
    maxseq = s.loc[ourmax,'seq']
    return maxseq

#find sequences which have 99% of sequences mapping a single barcode to a
#single sequence.

shares = toobigdf[['tag','seq','ct']].groupby(by='tag').apply(findshare)

sharesseq = toobigdf[['tag','seq','ct']].groupby(by='tag').apply(findmaxseq)

bigshares = sharesseq.loc[shares>.99]

outbigshares = pd.DataFrame()

outbigshares['seq'] = sharesseq.loc[shares>.99]

outbigshares['ct'] = shares.loc[shares>.99]

outdf5 = outbigshares.reset_index()


print('number of good sequencing counts ' + str(outdf5['ct'].sum()))




#we now want to automatically identify and split the identified sequences
#into genes.
genedf = pd.io.parsers.read_csv('/home/bill/next100genes/compedgenesv4.csv')



def check_all_muts(s,seq):
    '''This functions will check how many differences there are between
    a wild type sequence and a given sequence.'''
    temparr = np.array(list(s))
    return np.sum(temparr==seq)


def findgene(s):
    '''this function will return the genename from a set of sequences'''
    tempdf = genedf['geneseq'].apply(check_all_muts,args=(np.array(list(s)),))
    nmut = np.max(tempdf)
    maxind = np.argmax(tempdf)
    genename = genedf.loc[maxind,'name']
    return genename


def findgene_nmut(s):
    '''this function returns the number of mutations from a wild type sequence
    '''
    tempdf = genedf['geneseq'].apply(check_all_muts,args=(np.array(list(s)),))
    nmut = np.max(tempdf)
    maxind = np.argmax(tempdf)
    genename = genedf.loc[maxind,'name']
    return nmut


outdf5['gene'] = outdf5['seq'].apply(findgene)


outdf5['nmut'] = outdf5['seq'].apply(findgene_nmut)


#we will now filter overly mutated sequences, if a sequence has an over 50%
#mutation rate, we will remove it.
goodmut = outdf5['nmut'] > 80
outdf6 = outdf5.loc[goodmut]


def format_string(x):
    '''We use this to format output string for saving'''
    return '%10.6f' %x

pd.set_option('max_colwidth',int(1e8))

#loop through all genes in the experiment, outputing a key between barcode
#and sequence for each gene.
for gene in outdf6['gene'].value_counts().index:
    outdfgene = outdf6.loc[outdf6['gene'] == gene]
    print(len(outdfgene['seq'].value_counts().index))
    
    outdfgene.to_string(open(str(gene) + 'gooddataset_nslater','w'), index=False,col_space=10,float_format=format_string,justify='left')