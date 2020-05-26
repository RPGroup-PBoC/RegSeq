import pandas as pd
import Bio.Seq as Seq
import Bio.SeqIO
import numpy as np
import mpathic.simulate_library as simulate_library
from mpathic.utils import collapse_further
import mpathic.profile_mut as profile_mut
import mpathic.profile_freq as profile_freq



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


def get_wt_seqs(file, output, wild_type_genome='../data/wild_type/sequencev3.fasta'):
    """
    Read table and find wildtype sequences for given genes 
    115 bp upstream and 45 bp downstream of TSS.
    
    Parameters
    ----------
    file : str
        Path to file used as input
    output : str
        Path to store dataframe
    wild_type_genome : str, default '../data/wild_type/sequencev3.fasta'
        Path of file containing wildtype genome
    """
    for record in Bio.SeqIO.parse(wild_type_genome, "fasta"):
        genome = str(record.seq)
        
    totdf = pd.read_csv(file)
    totdf = totdf.loc[:,['name','start_site','rev']]
    totdf['geneseq'] = totdf.apply(find_seq, axis=1, args=(genome,))
    totdf.dropna()
    totdf.to_csv(output, index=False)

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
    
    # mutate interal all bp. Can be changed if o
    firstbase = 0
    finalbase = 160
    
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
    buffer=0,
    max_it=10,
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
    max_int : int
        Number of maximal interations.
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
            


def gen_mutated_seq(
    file,
    output,
    forward_primers='../data/primers/forward_finalprimers.fasta',
    reverse_primers='../data/primers/reverse_finalprimers.fasta',
    norder=150000,
    fix_ex_rate=False, 
    fix_rep=False,
    fix_low_base=False,
):
    """
    Generate mutated sequences from given sequences.
    
    Parameters
    ----------
    file : str
        Path to file with wild type sequences
    output : str
        Path to file where results are stored
    forward_primers : str, default '../data/primers/forward_finalprimers.fasta'
        Path to file containing list of forward primers.
    reverse_primers : str, default '../data/primers/reverse_finalprimers.fasta'  
        Path to file containing list of reverse primers.
    norder : int, default 150000
        Number of total sequences
    fix_ex_rate : boolean, default False
        If True, generate new sequences if genes have expectional mutation rates.
    fix_rep : boolean, default False
        If True, generate new sequences if sequences are not unique.
    fix_low_base : boolean, default False
        If True, generate new sequences if some bases are underrepresented.
    """
    # Read files
    df = pd.read_csv(file)
    kosprimefwd = Bio.SeqIO.parse(forward_primers,'fasta')
    kosprimerev = Bio.SeqIO.parse(reverse_primers,'fasta')
    ngenes = len(df.index)
    
    # Get primer pairs
    num_groups = int(np.ceil(ngenes/5))
    fwdprimes, revprimes = get_primers(kosprimefwd, kosprimerev, 100, 100+num_groups-1)
    nseqs = int(np.floor(norder/ngenes))
    
    # Generate sequences
    allseqs, primer_df = mutation_sequences(df, fwdprimes, revprimes, nseqs)
    
    # Check criteria
    allseqs = check_mutation_rate(
        df, 
        allseqs, 
        fix_ex_rate=fix_ex_rate, 
        fix_rep=fix_rep,
        fix_low_base=fix_low_base, 
        primer_df=primer_df
    )
    
    # Store results
    pd.set_option('max_colwidth',int(1e8))
    allseqs.to_csv(output, index=False)
    
    
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