import pymc
import pandas as pd
import scipy as sp
import numpy as np
import Bio.Seq as Seq
import mpathic.simulate_library as simulate_library
from mpathic.utils import collapse_further
import mpathic.profile_mut as profile_mut
import mpathic.profile_freq as profile_freq
from scipy.special import erfc
import logomaker
import matplotlib.pyplot as plt
import seaborn as sns


def seq2mat(seq,seq_dict):

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


def seq2array_for_matmodel(dataset_df, chunksize=1000):
    """Compute the wt sequence. """
    #use first sequence to figure out how many basepairs (and so number of
    #features to fit)
    seq_dict,inv_dict = choose_dict('dna')
    n_bases = len(dataset_df.loc[0,'seq'])
    numfeatures = len(dataset_df.loc[0,'seq'])*len(seq_dict)
    

    # Process dataframe in chunks
    startrow = 0
    endrow = startrow+chunksize-1
    numrows = dataset_df.shape[0]

    temp_mat = np.zeros((numrows,numfeatures))
    for i in range(numrows):
        temp_seq = seq2mat(dataset_df.loc[i,'seq'],seq_dict)
        temp_mat[i,:] = temp_seq.ravel(order='F')
    
    #now we are going to calculate the wild type sequence
    #first get counts 
    summed_mat = temp_mat.sum(axis=0)
    #now we find the highest number number of bases.
    wt_mat = np.zeros((1,numfeatures))
    for i in range(n_bases):
        max_index = np.argmax(summed_mat[i*len(seq_dict):(i+1)*len(seq_dict)])[0]
        #find the corresponding base
        max_base = seq_dict_inv[max_index]
        wt_mat[i*4 + max_index] = 1

    return scipy.sparse.csr(temp_mat),wt_mat

def calc_test_stat(allratios, r1, r0):
    """
    """
    return (allratios - r0) / (r1 - r0)


def cox_mann_p_values(files, output_file='test_pval.txt'):
    """Compute p-values following Cox and Mann, Nature Biotechnology 2008."""
    for z, name in enumerate(files):
        noleader = name.split('/')[-1]
        # load in file with proteins and enrichments
        indf = pd.read_csv(name, delimiter=",")

        # get the correct column name that contains the heavy to light ratio.
        indf_ratio_col = "Ratio H/L normalized"
        row = indf.loc[0,:]
        name_protein = row['Protein names']
        
        # Get the column of all enrichment ratios.
        q = indf[indf_ratio_col]

        # Drop any ratios equal to zero.
        q = q[q > 1e-8]
        
        # log transform ratios 
        allratios = np.log(np.array(list(q)))
        
        # Estimate the S.D.
        [rlow,r0,r1] = np.percentile(allratios,[15.87,50,84.13])
        
        # Compute the test stat z 
        test_stats = calc_test_stat(allratios,r0=r0,r1=r1)
        # Calculate the p_values f
        p = .5*erfc(test_stats/np.sqrt(2))
        # Check lowest p-value to see if the most enriched protein is an outlier
        pval = p.min()
        # Multiply by the number of enrichment ratios to correct for multiple hypothesis testing.
        pval = pval*len(p)
        # Write results to file.
        if z == 0:
            with open(output_file,'w') as f:
                f.write(noleader + ',p_val' +  '\n')
                f.write(name_protein + ',' + str(p.min()*len(allratios)) + '\n')
        else:
            with open(output_file,'a') as f:
                f.write(noleader + ',p_val' +  '\n')
                f.write(name_protein + ',' + str(p.min()*len(allratios)) + '\n')
        
    

#set the plotting style
def pboc_style_mpl():
    """
    Formats matplotlib plotting enviroment to the style used in
    Physical Biology of the Cell, 2nd edition.
    """
    rc = {'lines.linewidth': 1.25,
          'axes.labelsize': 8,
          'axes.titlesize': 9,
          'axes.facecolor': '#E3DCD0',
          'xtick.labelsize': 7,
          'ytick.labelsize': 7,
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.color': '#ffffff',
          'legend.fontsize': 8,
          'figure.dpi': 300,
          'savefig.dpi': 300}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('xtick.major', pad=-1)
    plt.rc('ytick.major', pad=-1)
    plt.rc('mathtext', fontset='stixsans', sf='sansserif')
    plt.rc('figure', figsize=[3.5, 2.5])
    plt.rc('svg', fonttype='none')
    plt.rc('legend', title_fontsize='8', frameon=True,
           facecolor='#E3DCD0', framealpha=1)
    sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)
