import pandas as pd
import os
import matplotlib.pyplot as plt
from .seq_utils import stitch, format_string


def read_data(file, barcode_length=20, trailing_sequence_length=21):
    """Load fastq file and extract sequence and barcode.
    
    Parameters
    ----------
    file : str
        Path to fastq file
    barcode_length : int, default 20
        Length of barcode sequence
    trailing_sequence_length : int, default 21
        Number of bp trailing barcode
        
    Returns
    -------
    tagcounts : Pandas.Series
        Panda Series contaning barcodes as index and how often they are counted
    """
    # No maximum length for columns
    pd.set_option('max_colwidth',int(1e8))
    
    # Load sequences
    df = pd.read_csv(file, header=None)
    
    # Extract the sequences from fastq format
    df = df.loc[1::4]
    
    # Extract barcodes
    tags = df[0].str.slice(
        -trailing_sequence_length - barcode_length,
        -trailing_sequence_length
    )
    
    # Count barcodes
    tagcounts = tags.value_counts()
    return tagcounts


def combine_counts(
    mRNA_file,
    DNA_file,
    tag_key_file,
    output_file
    ):
    """Compute counts for sequences from mRNA and DNA.
    
    Parameters
    ----------
    mRNA_file : str
        Path of file for mRNA sequencing
    DNA_file : str
        Path of file for DNA sequencing
    tag_key_file : str
        Path of file for barcode/sequence mapping
    output_file : str
        Path of file constructed for output
    """
    
    # Load barcode key
    tagkey = pd.read_csv(tag_key_file)
    tagkey = tagkey.set_index('tag')
    
    # Load mRNA data
    mRNA_counts = read_data(mRNA_file)
    
    # Make dataframe with mutated sequence and barcode
    tempdf = tagkey.reindex(mRNA_counts.copy().index)

    # Assign sequencing counts based on mRNA sequencing file
    tempdf['ct_1'] = mRNA_counts.copy()
    tempdf = tempdf.dropna()

    # Load DNA plasmid data
    DNA_counts = read_data(DNA_file)
    
    # Make dataframe with mutated sequence and barcode
    c = tagkey.reindex(DNA_counts.copy().index)
    
    # Assign sequencing counts based on DNA plasmid sequencing file
    c['ct_0'] = DNA_counts.copy()
    c = c.dropna()
    
    # Concatenate DataFrames
    output_df = pd.concat([tempdf, c], axis=0, sort=True)
    
    # Assign zero counts to NaN
    output_df = output_df.fillna(0)
    
    # Remove unnecessary columns
    output_df = output_df[['ct_0','ct_1','gene','seq']]

    # Total counts column
    output_df['ct'] = output_df[['ct_0','ct_1']].sum(axis=1)
    output_df = output_df.reset_index()
    
    # Rename barcode column
    output_df = output_df.rename(columns={'index':'tag'})
    
    # Combine sequence and barcode
    output_df['seq'] = output_df.apply(stitch, axis=1)

    output_df = output_df.groupby(by=['seq',"gene"]).mean()

    output_df = output_df.reset_index()
    output_df = output_df[['ct','ct_0','ct_1','gene','seq']]

    #Output the final dataset for associated gene.
    output_df.to_csv(output_file, index=False)
    
    # Old output format, might be needed again
    #output_df.to_string(
    #    open(output_file,'w'), 
    #    index=False,
    #    col_space=10,
    #    float_format=format_string,
    #    justify='left'
    #)