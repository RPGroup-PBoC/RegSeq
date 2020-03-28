import pymc
import scipy as sp
import numpy as np
import 

def choose_dict(dicttype,modeltype='LinearEmat'):
    '''Get numbering dictionary for either dna,rna, or proteins'''
    if dicttype == 'dna':
        seq_dict = {'A':0,'C':1,'G':2,'T':3}
        inv_dict = {0:'A',1:'C',2:'G',3:'T'}
    if dicttype == 'rna':
        seq_dict = {'A':0,'C':1,'G':2,'U':3}
        inv_dict = {0:'A',1:'C',2:'G',3:'U'}
    if dicttype == 'protein':
        seq_dict = {
            '*':0,'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,
            'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20}
        inv_dict = {v:k for k,v in seq_dict.items()}
    if modeltype == 'Neighbor':
        seq_dict = {
            ''.join([inv_dict[i],inv_dict[z]]):i*len(seq_dict)+z 
            for i in range(len(seq_dict)) for z in range(len(seq_dict))}
        inv_dict = {seq_dict[i]:i for i in seq_dict.keys()}
    return seq_dict,inv_dict


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
    # Compute the wt sequence
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

    
        
    


# Given a large list of sequences, produce a (sparse) mutation array
def dataset2mutarray(dataset_df, modeltype, chunksize=1000, rowsforwtcalc=100):

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
