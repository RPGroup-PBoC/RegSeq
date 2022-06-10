import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from .seq_utils import stitch, format_string


def check_length(
    input_file_name, 
    barcode_length=20, 
    sequence_length=160, 
    optimal_lengths= np.array([299, 295]), 
    trailing_lengths= np.array([24, 20]),
    starting_length=20
):
    """
    Check length of sequences.
    
    Return sequences that have 24 or 20 trailing bp, and have 20 bp at the end. Sequences
    with varying length have insertions and deletions and are filtered out.
    Paramters:
    ----------
    input_file_name : str
        file name for the fastq sequencing
    barcode_length : int, default 20
    sequence_length : int, default 160
    optimal_lengths : array-like, default np.array([299, 295])
        Numpy array containing possible lengths of sequences
    trailing_lengths : array-like, default np.array([24, 20])
        Numpy array containing possible number of trailing bp after barcode
    starting_length : int, default 20
        Number of bp preceeding the sequence of interest
    
    Return
    ------
    sliceddf : Pandas DataFrame
        DataFrame containing only the mutated sequence and the barcode
    
    """
    if type(trailing_lengths) == list:
        trailing_lengths = np.array(trailing_lengths)
    elif type(trailing_lengths) == np.ndarray:
        pass
    else:
        raise RuntimeError("`trailing_lengths` has to be a numpy array or list.")
        
    if type(optimal_lengths) == list:
        optimal_lengths = np.array(optimal_lengths)
    elif type(optimal_lengths) == np.ndarray:
        pass
    else:
        raise RuntimeError("`optimal_lengths` has to be a numpy array or list.")
        
    # Load data
    df = pd.io.parsers.read_csv(input_file_name,delim_whitespace=True,header=None)

    #Select only rows with sequences (fastq format).
    df = df.loc[1::4,0]

    # Find all lengths
    lengths = df.apply(len)

    # Find the correct length by finding the most common length
    lengthsmax = list(lengths.value_counts().index)[0]
    print('optimal length is ' + str(lengthsmax))


    # Find all sequences with correct length.
    goodlength = (df.apply(len) == lengthsmax)
    df = df.loc[goodlength]
    ind = np.where(optimal_lengths ==lengthsmax)[0]
    if len(ind) == 0:
        raise ValueError('Sequence length not in the list of required lengths.')
    else:
        sliceddf = df.str.slice(starting_length,-trailing_lengths[ind[0]])
        
    return sliceddf


def check_N(s):
    "Exclude sequences with undetermined bases."
    if 'N' in s:
        return False
    else:
        return True


def stitch_barcode_sequence(df, sequence_length=160, barcode_length=20):
    """Construct a data frame with mutated sequence and barcode"""
    #We will now construct a data frame with mutated sequence and barcode
    tempdf = pd.DataFrame()
    tempdf['seq'] =df.str.slice(0,sequence_length)
    tempdf['tag'] =df.str.slice(-barcode_length,)

    # Stitch them together into one sequence (with barcode at end).
    tempstitched = tempdf.apply(stitch,axis=1)

    # Remove all sequences with N's in them.
    noN = tempstitched.apply(check_N)
    tempstitched = tempstitched.loc[noN]
    return tempstitched


def check_barcode_uniqueness(df, sequence_length=160, barcode_length=20):
    """Check unique mapping of barcodes and sequences.
    
    Parameters
    ----------
    df : Pandas.Series
        Series of sequences of length `sequence_length`, stitched together
        with barcode of length `barcode_length.
    sequence_length : int, default 160
    barcode_length : int, default 20
    
    Returns
    -------
    goodtags : pd.Series
        Tags that are related to unique sequence
    counts : pd.Series
        Number of sequences related to a barcode
    output_df
    
    """
    # Number of barcode/sequence combinations
    q = df.value_counts()

    #Drop those that only are sequenced 1 time. Helps to cut
    #out some sequencing errors.
    #temp_df = q.loc[q > 1]
    count_df = q

    # Seperate sequence and barcode
    output_df = pd.DataFrame()

    temp_df = pd.DataFrame(count_df.index)
    output_df['seq'] = temp_df[0].str.slice(0,sequence_length)
    output_df['tag'] = temp_df[0].str.slice(-barcode_length,)

    count_df = count_df.reset_index()
    output_df['count'] = count_df[0]

    # Number of sequences for a given tag
    counts = output_df['tag'].value_counts()

    # Find unique tag/sequence 
    goodtags = counts.loc[counts == 1].index
    
    return goodtags, counts, output_df


def check_rare_barcode_errors(goodtags, counts, df):
    """Even for low sequencing error rates you could have a single error
    that could make it appear that a barcode maps to multiple sequences even
    if it doesn't. As a result, if the most common sequence is 99% of all sequences
    for a given barcode we will include it in the final mapping."""
    exceed = (counts > 1)

    indexed_df = df.set_index('tag')

    gooddf = indexed_df.loc[counts == 1]

    gooddf = gooddf.reset_index()

    toobigdf = indexed_df.loc[exceed]

    toobigdf = toobigdf.reset_index()

    shares = toobigdf[['tag','seq','count']].groupby(by='tag').apply(findshare)

    sharesseq = toobigdf[['tag','seq','count']].groupby(by='tag').apply(findmaxseq)
    

    if len(toobigdf.index) != 0:
        bigshares = sharesseq.loc[shares>.99]

        outbigshares = pd.DataFrame()

        outbigshares['seq'] = sharesseq.loc[shares>.99]

        outbigshares['count'] = shares.loc[shares>.99]

        output_df = outbigshares.reset_index()
    else:
        output_df = pd.DataFrame()

    output_df = pd.concat([gooddf, output_df])
    print('number of good sequencing counts ' + str(output_df['count'].sum()))
    return output_df
  

def detect_genes(df, wildtypefile):
    """
    Find genes that sequences belong to.
    
    Output is written into file.
    
    Parameters
    ----------
    df : Pandas.Dataframe
        DataFrame containg barcodes and sequences that have been
        checked to unique mapping.
    wildtypefile : str
        Path to file containing genes and sequences of wildtype
        
    Returns
    -------
    temp_df : Pandas DataFrame
        DataFrame containing found sequences for each gene and number of mutations
    """
    
    wildtype_df = pd.read_csv(wildtypefile)
    df['gene'] = df['seq'].apply(findgene, args=(wildtype_df,))
    df['nmut'] = df['seq'].apply(findgene_nmut, args=(wildtype_df,))


    # Filter overly mutated sequences, if >50% is mutated
    goodmut = df['nmut'] < 80
    temp_df = df.loc[goodmut]

    pd.set_option('max_colwidth',int(1e8))
    return temp_df
    
    
    
def key_barcode_sequence(data_file, output_path, wildtypefile='../data/prior_designs/wtsequences.csv', genes=None): 
    """
    Go through functions to create unique map of barcode to sequence and gene in wiltype.
    
    The sequences are checked for correct lengths, to exlude insertion and deletion events.
    Then, created sequences and barcodes are extracted (removing overhangs) and unique barcode/
    sequence maps are found. Possible sequencing errors that lead to false negatives in uniqueness
    are considered. Sequences are compared to gene sequences in wildtype.
    
    Parameters
    ----------
    data_file : str
        Path to file containing sequencing data.
    output_path : str
        Path to folder where results are stored.
    wildtypefile : str
        Path for file containing wild type genetic sequences.
    genes : List, default None
        List of genes for which mapping is returned. If None, all maps are returned.
    -------
    """
    if genes != None:
        if type(gene) != list:
            raise RuntimeError("Type of `genes` has to be list.")
    
    # Find sequences with correct length
    correct_seq = check_length(data_file)
    
    # Extract sequence and barcode
    stitched = stitch_barcode_sequence(correct_seq)
    
    # Find unique barcode/sequence relations
    barcodes, counts, seq_tag_df = check_barcode_uniqueness(stitched)
    
    # Check for possible sequencing errors in mapping
    seq_tag_df = check_rare_barcode_errors(barcodes, counts, seq_tag_df)
    
    # Find gene relating to sequence and store result
    df = detect_genes(seq_tag_df, wildtypefile)
    if genes == None:
        for gene in df["gene"].unique():
            genedf = df.loc[df["gene"] ==gene]
            genedf.drop(['gene'], axis=1).to_csv(output_path + gene + "_barcode_key.csv", index=False)
    else:
        for gene in genes:
            genedf = df.loc[df["gene"] ==gene]
            genedf.drop(['gene'], axis=1).to_csv(output_path + gene + "_barcode_key.csv", index=False)
    
    
def findshare(df):
    """Finds percentage of single counted sequences."""
    ourmax = df['count'].max()
    oursum = df['count'].sum()
    return ourmax/oursum


def findmaxseq(df):
    """Find most common sequence.
    
    Parameters
    -----------
    df : Pandas.DataFrame
        Dataframe with colums `seq`: sequences, and `count`: sequence abundance.
    
    Returns
    -------
    maxseq : 
    """
    ourmax = df['count'].idxmax()
    maxseq = df.loc[ourmax, 'seq']
    return maxseq


def check_all_muts(s, seq):
    """Compare two sequences for identical bases."""
    temparr = np.array(list(s))
    return np.sum(temparr!=seq)


def findgene(s, df):
    """Return the gene name from a set of sequences.
    
    Looks for wildtype gene in given table whose sequence is most similar to
    provided sequence.
    
    Parameters
    ----------
    s : str
        Sequence
    df : Pandas.Dataframe
        Table of genes and wiltype sequences
    
    Returns
    -------
    genename : str
        Name of gene 
    """

    tempdf = df['geneseq'].apply(check_all_muts, args=(np.array(list(s)),))
    maxind = np.argmin(np.array(tempdf))
    genename = df.loc[maxind, 'name']
    return genename


def findgene_nmut(s, df):
    """Return the number equal bp of sequence compared to wt.
    
    Looks for wildtype gene in given table whose sequence is most similar to
    provided sequence. Then computes number of identical bases. 
    
    Parameters
    ----------
    s : str
        Sequence
    df : Pandas.Dataframe
        Table of genes and wiltype sequences
    
    Returns
    -------
    nmut : int
        Number of mutations in sequence
    """

    tempdf = df['geneseq'].apply(check_all_muts,args=(np.array(list(s)),))
    nmut = np.min(tempdf)
    return nmut

