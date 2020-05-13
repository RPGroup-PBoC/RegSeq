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
    '''Get numbering dictionary for either dna,rna, or proteins'''
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

    if modeltype == 'NBR' or modeltype == 'PAIR' or "Neighbor":
        seq_dict = {
            ''.join([inv_dict[i],inv_dict[z]]):i*len(seq_dict)+z 
            for i in range(len(seq_dict)) for z in range(len(seq_dict))}
        inv_dict = {seq_dict[i]:i for i in seq_dict.keys()}
    return seq_dict,inv_dict


def sliding_window(y,windowsize=3):
        out_vec = np.zeros_like(y)
        for i in range(len(y)-windowsize):
            out_vec[i] = np.sum(y[i:i+windowsize])/windowsize
        return out_vec


class GaugePreservingStepper(pymc.Metropolis):
    """Perform monte carlo steps that preserve the following choise of gauge:
    sum of elements in each column = 0, overall matrix norm = 1."""
    def __init__(self,stochastic):
        pymc.Metropolis.__init__(self,stochastic)

    def propose(self):
        # to avoid overwriting weirdness
        
        emat_temp = self.stochastic.value.copy()
        # number of columns in matrix
        num_col = emat_temp.shape[1]
        # first choice a random 4L dimensional vector
        r = sp.random.standard_normal(emat_temp.shape)
        # dot product of proposed direction with current vector
        lambda_0 = sp.sum(emat_temp*r)
        # dot products of each column with basis vectors
        lambda_vec = 0.5*sp.sum(r,axis=0)

        s = sp.zeros_like(r)
        # s is a vector based on r, but along the surface of the hypersphere
        for j in range(emat_temp.shape[1]):
            s[:,j] = r[:,j] - lambda_0*emat_temp[:,j] - lambda_vec[j]*(0.5*sp.ones(emat_temp.shape[0]))
        dx = self.adaptive_scale_factor*s/sp.sqrt(sp.sum(s*s))
        self.stochastic.value = (emat_temp+dx)/sp.sqrt(sp.sum((emat_temp+dx)**2))

def get_beta_for_effect_df(effect_df,target_info,\
    min_beta=.001,max_beta=100,num_betas=1000):
    '''This function finds the appropriate scaling factor for displaying sequence
    logos. From empirical results, most binding sites will ahve approximately
    1 bit per base pair of information'''
    betas = np.exp(np.linspace(np.log(min_beta),np.log(max_beta),num_betas))
    infos = np.zeros(len(betas))
    for i, beta in enumerate(betas):
        prob_df = logomaker.transform_matrix(df=beta*effect_df,from_type='weight',to_type='probability')
        infos[i] = get_info(prob_df)
    i = np.argmin(np.abs(infos-target_info))
    beta = betas[i]
    return beta


def get_info(df):
    '''This function finds the total information content of a binding site'''
    #define the background probabilities for E. coli bases.
    gc = .508
    background_array =np.array([(1-gc)/2,gc/2,gc/2,(1-gc)/2])
    #add in small value to make sure no probabilities are exactly zero.
    df = df + 1e-7
    
    return np.sum(df.values*np.log2(df.values/background_array))

def alt4(df, coarse_graining_level = 0.01,return_freg=False):
    '''
    MI ESTIMATOR EDITED BY JBK 
    Used when lm=memsaver 
    REQUIRES TESTING AND PROFILING.
    '''
    n_groups=500
    n_seqs = len(df.index)
    binheaders = utils.get_column_headers(df)
    n_batches = len(binheaders)
    cts_grouped = sp.zeros([n_groups,n_batches])
    group_num = 0
    frac_empty = 1.0
    
    #copy dataframe
    tmp_df = df.copy(binheaders+['val'])

    # Speed computation by coarse-graining model predictions
    if coarse_graining_level:
        assert type(coarse_graining_level)==float
        assert coarse_graining_level > 0
        vals = tmp_df['val'].values
        scale = np.std(vals)
        coarse_vals = np.floor((vals/scale)/coarse_graining_level)
        tmp_df['val'] = coarse_vals
        grouped = tmp_df.groupby('val')
        grouped_tmp_df = grouped.aggregate(np.sum)
        grouped_tmp_df.sort_index(inplace=True)
    else:
        grouped_tmp_df = tmp_df
        grouped_tmp_df.sort_values(by='val',inplace=True)
    # Get ct_xxx columns
    ct_df = grouped_tmp_df[binheaders].astype(float)
    cts_per_group = ct_df.sum(axis=0).sum()/n_groups
    # Histogram counts in groups. This is a bit tricky
    group_vec = np.zeros(n_batches)
    for i,row in ct_df.iterrows():
        row_ct_tot = row.sum()
        row_ct_vec = row.values
        row_frac_vec = row_ct_vec/row_ct_tot 

        while row_ct_tot >= cts_per_group*frac_empty:
            group_vec = group_vec + row_frac_vec*(cts_per_group*frac_empty)
            row_ct_tot -= cts_per_group*frac_empty

            # Only do once per group_num
            cts_grouped[group_num,:] = group_vec.copy() 
            # Reset for new group_num
            group_num += 1
            frac_empty = 1.0
            group_vec[:] = 0.0
        group_vec += row_frac_vec*row_ct_tot
        
        frac_empty -= row_ct_tot/cts_per_group
    if group_num == n_groups-1:
        cts_grouped[group_num,:] = group_vec.copy()
    elif group_num == n_groups:
        pass
    else:
        raise TypeError(\
            'group_num=%d does not match n_groups=%s'%(group_num,n_groups))
    # Smooth empirical distribution with gaussian KDE
    f_reg = scipy.ndimage.gaussian_filter1d(cts_grouped,0.04*n_groups,axis=0)
    if return_freg:
    # Return mutual information
    	return info.mutualinfo(f_reg),f_reg
    else:
        return info.mutualinfo(f_reg)

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
        
    
def dataset2mutarray(dataset_df, modeltype, chunksize=1000, rowsforwtcalc=100):
    """Produce sparse mutation array for large list of sequences."""
    # Determine the type of model and set seq2array function appropriately
    if modeltype=='MAT':
        seqs2array = seqs2array_for_matmodel
    else:
        raise SortSeqError('Unknown model type: %s'%modeltype)

    # Determine seqtype, etc.
    seqcol = qc.get_cols_from_df(dataset_df,'seqs')[0]
    seqtype = qc.colname_to_seqtype_dict[seqcol]
    wtcol = qc.seqtype_to_wtcolname_dict[seqtype]

    # Compute the wt sequence
    rowsforwtcalc = min(rowsforwtcalc,dataset_df.shape[0])
    dataset_head_df = dataset_df.head(rowsforwtcalc)
    mut_df = profile_mut(dataset_head_df)
    wtseq = ''.join(list(mut_df[wtcol]))
    wtrow = seqs2array([wtseq], seq_type=seqtype).ravel().astype(bool)
    numfeatures = len(wtrow)

    # Process dataframe in chunks
    startrow = 0
    endrow = startrow+chunksize-1
    numrows = dataset_df.shape[0]

    # Fill in mutarray (a lil matrix) chunk by chunk
    mutarray_lil = lil_matrix((numrows,numfeatures),dtype=int)
    matrix_filled = False
    while not matrix_filled:

        if startrow >= numrows:
            matrix_filled = True
            continue
        elif endrow >= numrows:
            endrow = numrows-1
            matrix_filled = True


        # Compute seqarray
        seqlist = list(dataset_df[seqcol][startrow:(endrow+1)])
        seqarray = seqs2array(seqlist, seq_type=seqtype)

        # Remove wt entries
        tmp = seqarray.copy()
        tmp[:,wtrow] = 0

        # Store results from this chunk
        mutarray_lil[startrow:(endrow+1),:] = tmp

        # Increment rows
        startrow = endrow+1
        endrow = startrow + chunksize - 1

    # Convert to csr matrix
    mutarray_csr = mutarray_lil.tocsr()

    # Return vararray as well as binary representation of wt seq
    return mutarray_csr, wtrow


def find_seq(s, genome):
    '''
    Find sequence 115 bp upstream and 45 bp downstream of TSS. If gene is reverse 
    transcribed, return reverse complement of sequence.
    
    Parameters
    ----------
    s : array-like
        Name, TSS and transcription direction of a gene
    genome : string
        Reference genome
        
    Returns
    -------
    output : string
        Genetic sequence 115 bp upstream and 45 bp downstream of TSS
    '''
    ss = int(s[1])
    if s[2] == 'rev':
        
        gene = genome[ss-45:ss+115]
        tempgene = Seq.Seq(gene)
        outgene = str(tempgene.reverse_complement())
    elif s[2] == 'fwd':
        outgene = genome[ss-115:ss+45]
        
    return outgene


def mut_seq_random(s, fwdprimer, revprimer, nseqs):
    """Mutate sequence at given rate and add primers.
    
    Parameters
    ----------
    s : array-like
        Name, TSS and transcription direction of a gene
    fwdprimer : string
        Forward primer
    revprimer : string
        Reverse primer
    nseqs : Int
        Number of generated sequences
        
    Returns
    -------
    outdf : Pandas.DataFrame
        Dataframe of mutated sequences
    """
    
    # mutate interal 140 bp.
    firstbase = 10
    finalbase = 150
    
    # use mpathic function to generate randomly mutated library
    seqs_df = simulate_library.main(wtseq=s[firstbase:finalbase], numseq=nseqs)
    # only take independent sequences
    seqs_df = collapse_further(seqs_df)
    seqs = seqs_df['seq']
    # remove wild type sequences
    outseqs = [fwdprimer + s + revprimer]  + [fwdprimer + s[:firstbase] + str(x) + s[finalbase:] + revprimer  for x in seqs]
    outdf = pd.DataFrame()
    outdf['seq'] = outseqs
    return outdf


def get_primers(fwd, rev, start, end):
    """Return primers from Bio.Seq iterator objects.
    
    Parameters
    ----------
    fwd : Bio.Seq iterator
        Library for forward primers
    rev : Bio.Seq iterator
        Library for reverse primers
    start : Int
        Index of first primers to be read
    end : Int
        Index of last primers to be read
        
    Returns
    -------
    fwd_list : List
        List of forward primers
    rev_list : List
        List of reverse primers
    """
    
    
    fwd_list = []
    rev_list = []
    records = list(enumerate(zip(fwd, rev)))[start:end+1]
    for i,(fwd_rec, rev_rec) in records:
        fwd_prime = str(fwd_rec.seq)
        rev_prime = str(rev_rec.seq)
        fwd_list.append(fwd_prime)
        rev_list.append(rev_prime)
        
    return fwd_list, rev_list
        
    
    
def mutation_sequences(genes, fwdprimes, revprimes, nseqs):
    """
    Get mutated sequences for a set of genes and primers.
    
    Parameters
    ----------
    genes : Pandas Dataframe
        Gene names and sequences
    fwdprimes : array-like
        List of forward primers
    revprimes : array-like
        List of reverse primers
    nseqs : int
        Number of mutated sequences
        
    Returns
    -------
    allseqs : Pandas Dataframe
        Dataframe of mutated sequences
    """
    
    allseqs = pd.DataFrame()
    primer_df = pd.DataFrame()
    for i,row in genes.iterrows():
        
        primernum = int(np.floor(i/5))
        # get fwd and rev primer.
        thefwdprimer = fwdprimes[primernum]
        therevprimer = revprimes[primernum]
        # mutate the sequence
        tempdf =  mut_seq_random(row['geneseq'], thefwdprimer, therevprimer, nseqs)
        # set which gene the group of sequences comes from.
        tempdf['gene'] = row['name']
        #we will build up the final dataframe of mutated sequences from the individual
        #dataframes of sequences. 
        allseqs = pd.concat([allseqs,tempdf],axis=0)
        primer_df = pd.concat(
            [primer_df,
             pd.DataFrame([[row['name'], row['geneseq'], thefwdprimer, therevprimer]], columns=['name', 'geneseq', "fwd_primer", "rev_primer"])],
          ignore_index=True
        )
    return allseqs, primer_df


def check_mutation_rate(
    genes, 
    sequences, 
    buffer=10, 
    fix_ex_rate=False, 
    fix_rep=False,
    fix_low_base=False,
    primer_df=pd.DataFrame()
):
    """
    Check mutation rate of generated sequences. Sequences can be generated again,
    if certain tests are not passed.
    
    Parameters
    ----------
    genes : Pandas.DataFrame
        DataFrame of gene names and sequences
    sequences : Pandas.DataFrame
        DataFrame of mutated sequences for each gene
    buffer : int
        Number of bases as buffer when checking non-primer and non-wildtype bases
    fix_ex_rate : boolean, default False
        States if sequences should be re-generated if mutation rate is far away from
        expectation.
    fix_rep : boolean, default False
        States if sequences should be re-generated if sequences are repeated.
    fix_low_base : boolean, default False
        States if sequences should be re-generated if there is a base with outstanding
        low mutation rate.
    primer_df : Pandas.DataFrame
        DataFrame containing information about wild type sequence and primer pairs, 
        needed to generate new sequences if any of the fix* arguments is True.
    
    Returns
    -------
    sequences : Pandas.DataFrame
        DataFrame of mutated sequences for each gene. Possibly with new sequences, if
        input sequences did not match criteria.
    """
    
    # Check that primer_df is given i
    if np.array([fix_ex_rate, fix_rep, fix_low_base]).any():
        if primer_df.empty:
            raise Exception('If sequences are supposed to be fixed, primers have to be given.')

    
    
    for i,row in genes.iterrows():
        gene = row['name']
        #dnaE gene is bugged out right now, just skip it for the moment.
        if gene == 'dnaE':
            continue

        # get sequences from target gene
        partialdf = _make_partial(gene, sequences)
        if np.array([fix_ex_rate, fix_rep, fix_low_base]).any():
            if np.array([fix_ex_rate * _check_mut_rate(partialdf, buffer),
                fix_rep * _check_repeated(partialdf),
                fix_low_base * _check_bases(partialdf, buffer)]).any() :
                fixed = False
                max_it = 10
                it = 0
                fwd = primer_df.loc[primer_df["name"] == gene, "fwd_primer"].values
                rev = primer_df.loc[primer_df["name"] == gene, "rev_primer"].values
                while fixed==False and it<max_it:
                    gene_df = primer_df.loc[primer_df["name"] == gene]
                    gene_df.index = [0]
                    temp, _ = mutation_sequences(gene_df, fwd, rev, len(partialdf))
                    temp = _make_partial(gene, temp)
                    if not np.array(
                        [fix_ex_rate * _check_mut_rate(temp, buffer),
                        fix_rep * _check_repeated(temp),
                        fix_low_base * _check_bases(temp, buffer)
                        ]).any() :
                        fixed = True
                    it += 1
                    
                sequences.loc[sequences['gene'] == gene,['seq']] = temp['seq']
                if it == max_it:
                    if _check_mut_rate(temp, buffer):
                        print('Bad mutation rate in gene {}!'.format(gene))

                    if _check_repeated(temp):
                        print('Repeated sequence for gene {}'.format(gene))

                    if _check_bases(temp, buffer):
                        print("Base with low frequency in gene {}!".format(gene))
        else:

            if _check_mut_rate(partialdf, buffer):
                print('Bad mutation rate in gene {}!'.format(gene))
            
            if _check_repeated(partialdf):
                print('Repeated sequence for gene {}'.format(gene))
            
            if _check_bases(partialdf, buffer):
                print("Base with low frequency in gene {}!".format(gene))
        print("Gene {} done.".format(gene))
                    
    return sequences
                    
                    

            
            
            
def calc_test_stat(allratios, r1, r0):
    return (allratios - r0) / (r1 - r0)


def cox_mann_p_values(files, output_file='test_pval.txt'):
    """Compute p-values following Cox and Mann, Nature Biotechnology 2008."""
    for z, name in enumerate(files):
        # load in file with proteins and enrichments
        indf = pd.read_csv(name, delimiter="\t")
        # get the correct column name that contains the heavy to light ratio.
        indf_ratio_col = "Ratio H/L normalized"
        row = indf.loc[0,:]
        name = row['Protein names']
        
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
        with open(output_file,'a') as f:
            f.write("Gene: {}, {}: {} \n ".format(name, indf_ratio_col, str(p.min()*len(allratios))))
            


            
# Helper functions
def _check_mut_rate(partialdf, buffer):
    mut = profile_mut.main(partialdf)

    #Check mutation rate of non-primer and non-wildtype bases
    relevantmuts = mut.loc[buffer + 20:179-buffer*2,'mut']
    return ((relevantmuts > .14).any() or (relevantmuts < .07).any())
    
    
def _check_repeated(partialdf):
    # Check for repeated sequences
    return len(partialdf['seq'].index) != len(partialdf['seq'].unique())
     
        
def _check_bases(partialdf, buffer):
    # Check frequencies of each base in the sequence range
    freqmat = profile_freq.main(partialdf)
    relevantfreq = freqmat.loc[20 + buffer:179-buffer*2,:]

    # Check mutation rates of all bases
    freqmin = relevantfreq[['freq_A','freq_C','freq_G','freq_T']].min(axis=1)
    relevantfreqmat = np.array(relevantfreq[['freq_A','freq_C','freq_G','freq_T']])
    return (freqmin < .02).any()
    
    
def _make_partial(gene, sequences):
    # get sequences from target gene
    partialdf = sequences.loc[sequences['gene'] == gene,['seq']]

    # Assign test count
    partialdf['ct'] = 1
    partialdf = partialdf[['ct','seq']]
    partialdf = partialdf.reset_index(drop=True)
    return partialdf


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
pboc_style_mpl()