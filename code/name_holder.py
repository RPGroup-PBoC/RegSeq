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
    seqs_df = utils.collapse_further(seqs_df)
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
    records = list(enumerate(zip(kosprimefwd, kosprimerev)))[start:end+1]
    for i,(fwd_rec, rev_rec) in records:
        fwd_prime = str(fwd_rec.seq)
        rev_prime = str(rev_rec.seq)
        fwd_list.append(fwd_prime)
        rev_list.append(rev_prime)
        
    return fwd_list, rev_list
        