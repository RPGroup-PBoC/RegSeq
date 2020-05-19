import numpy as np
import pandas as pd

from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
import scipy.sparse
import scipy as sp
import pdb
import pymc

from .utils import choose_dict

def least_squares(raveledmat, batch):
    """Linear regression of effects on gene expression of mutation of the
    corresponding base.
    Uses sklearn to do this regression."""
    
    clf = linear_model.LinearRegression()
    clf.fit(raveledmat, batch)
    emat = clf.coef_
    return emat


def lin_reg(inputname, outputname, wildtypefile='../data/prior_designs/wtsequences.csv'):

    # Load data
    df = pd.read_csv(inputname)
    
    # Load wild type sequences
    genedf = pd.read_csv('../data/prior_designs/wtsequences.csv')
    
    # Gene name
    gene = inputname.split('/')[-1].split('_')[0]
    
    # Extract the wild type sequence of gene of interest
    wt = str(genedf.loc[genedf['name'] == gene,'geneseq'].tolist()[0])

    # Convert to list.
    wtlist = np.array(list(wt))

    # Barcode length
    taglength = 20

    # Total promoter length
    seqlength = len(df['seq'][0]) - taglength #160 bp
    
    #we create dictionaries that relate A,C,G,T to the number 1,2,3,4
    seq_dict,inv_dict = choose_dict('dna')

    # Initialize array to paramaterize sequence. Mutations are denoted by 1
    all_mutarr = np.zeros((len(df.index),seqlength))

    # Parameterize sequences
    for i,row in df.iterrows():
        s = np.array(list(row['seq']))
        # Clip off any sequence past the 160 bp mutated sequence length.
        s_clipped = s[:seqlength]
        # Find mutations
        all_mutarr[i,:seqlength] = (wtlist != s_clipped)
    # IUse the ratio of mRNA counts to DNA counts to regress against. Add a pseudocount of 1.
    thetarget = np.array((df['ct_0']+1)/(df['ct_1']+1))
    
    # Center the mean
    thetarget = thetarget - np.mean(thetarget)
    # Fit mutation effects using linear regression
    emat = least_squares(all_mutarr, thetarget)

    #output results
    np.savetxt(outputname,emat)
    

class GaugePreservingStepper(pymc.Metropolis):
    """Perform monte carlo steps that preserve the following choise of gauge:
    sum of elements in each column = 0, overall matrix norm = 1."""
    def __init__(self,stochastic):
        pymc.Metropolis.__init__(self,stochastic)

    def propose(self):
        # to avoid overwriting weirdness
        
        emat_temp = self.stochastic.value.copy()
        emat_temp1 = emat_temp[:2, :160].copy()
        emat_temp2 = emat_temp[:, 160:].copy()
        
        # Number of columns in matrix
        num_col1 = emat_temp1.shape[1]
        
        # First choice a random 4L dimensional vector
        r = sp.random.standard_normal(emat_temp1.shape)
        
        # Dot product of proposed direction with current vector
        lambda_0 = sp.sum(emat_temp1 * r)
        
        # Dot products of each column with basis vectors
        lambda_vec = 0.5 * sp.sum(r, axis=0) 
        s = sp.zeros_like(r)
        
        # s is a vector based on r, but along the surface of the hypersphere
        for j in range(emat_temp1.shape[1]):
            s[:,j] = r[:,j] - lambda_0*emat_temp1[:,j] - lambda_vec[j] * \
            (0.5 * sp.ones(emat_temp1.shape[0]))
        
        dx = self.adaptive_scale_factor * s / sp.sqrt(sp.sum(s * s))
        emat_temp1_new = (emat_temp1 + dx)

        num_col = emat_temp2.shape[1]
        
        # First choice a random 4L dimensional vector
        r = sp.random.standard_normal(emat_temp2.shape)
        
        # Dot product of proposed direction with current vector
        lambda_0 = sp.sum(emat_temp2*r)
        
        # Dot products of each column with basis vectors
        lambda_vec = 0.5*sp.sum(r,axis=0)

        s = sp.zeros_like(r)
        
        # s is a vector based on r, but along the surface of the hypersphere
        for j in range(emat_temp2.shape[1]):
            s[:,j] = r[:,j] - lambda_0*emat_temp2[:,j] - lambda_vec[j] * \
            (0.5 * sp.ones(emat_temp2.shape[0]))
            
        dx = self.adaptive_scale_factor*s/sp.sqrt(sp.sum(s*s))
        emat_temp2_new = (emat_temp2+dx)
        emat_temp[:2,:160] = emat_temp1_new
        emat_temp[:,160:] = emat_temp2_new
        self.stochastic.value = (emat_temp) / sp.sqrt(sp.sum((emat_temp1_new)**2) \
                                                    +sp.sum((emat_temp2_new)**2))
        

def MaximizeMI_memsaver(
    seq_mat,
    df,
    emat_0,
    wtrow,
    db=None,
    burnin=1000,
    iteration=30000,
    thin=10,
    runnum=0,
    verbose=False,
    temp=4200
    ):
    """Performs MCMC MI maximzation in the case where lm = memsaver
    
    
    Parameters
    ----------
    seq_mat : 
    
    df : 
    
    emat_0 : 
    
    wtrow : 
    
    db : 
    
    burnin : 
    
    iteration : 
    
    thin : 
    
    runnum : 
    
    verbose : 
    
    temp : 
    
    Returns
    -------
    """
    
    # Number of sequences in the dataset.
    n_seqs = len(df.index)
    # Model for MCMC inference with pymc2.
    @pymc.stochastic(observed=True, dtype=pd.DataFrame)
    def pymcdf(value=df):
        return 0
    
    @pymc.stochastic(dtype=float)
    def emat(p=pymcdf, value=emat_0):
        #convert a L_bases array to a 2xL array for gauge fixing function
        len_seq = len(df.loc[0,'seq'])
        len_barcode = 20
        len_outputseq = len_seq - len_barcode 
        total_params = len_outputseq + len_barcode*4
        temp_val = np.zeros(total_params)
        temp_val[:len_outputseq] = value[1,:len_outputseq] - value[0,:len_outputseq]
        temp_val[len_outputseq:] = value[:,len_outputseq:].ravel(order='F')
        #Evaluate each sequence in the dataset with the temporary model.         
        p['val'] = seq_mat*temp_val      
        #evaluate mutual information, likelihood is the log             
        MI = EstimateMutualInfoforMImax.alt4(p.copy())  # New and improved
        return n_seqs*MI
    if db:
        dbname = db + '_' + str(runnum) + '.sql'
        M = pymc.MCMC([pymcdf,emat],db='sqlite',dbname=dbname)
    else:
        M = pymc.MCMC([pymcdf,emat])
    M.use_step_method(GaugePreservingStepper,emat)

    if not verbose:
        M.sample = shutthefuckup(M.sample)

    #M.sample(iteration,thin=thin,tune_interval=20000)
    M.sample(iteration,thin=thin)
    emat_mean = np.mean(M.trace('emat')[burnin:],axis=0)
    return emat_mean


def max_inf_mcmc(
    data_set_file, 
    output_db, 
    output_mean,
    burnin=1000,
    iteration=30000,
    thin=10,
    runnum=0,
    verbose=False,
    temp=4200
    ):
    """ Perform MCMC.
    
    Parameters
    ----------
    data_set_file : str
        Path to file containing data set formatted `match_data` module.
    output_db : str
        Path to file for output of database of MCMC
    output_mean : str
        Path to file for output of mean
    Returns
    -------
    
    """
    # Load data set
    df = pd.io.parsers.read_csv(data_set_file,delim_whitespace=True)#, delim_whitespace=True)

    temp_seq_mat, wtrow = numerics.dataset2mutarray(df.copy(),'MAT')
    temp_seq_mat2 = temp_seq_mat.toarray()

    #we need a parameter for the effect of mutation at each base pair.
    #we remove 20 base pairs to remove the barcode sequence on the end of the sequence.
    len_seq = len(df.loc[0,'seq'])
    len_outputseq = len_seq - 20
    len_barcode = 20
    # 4 parameters for the barcode
    total_params = len_outputseq + len_barcode * 4
    seq_mat = np.zeros((temp_seq_mat2.shape[0],total_params))
    for i in range(len_outputseq):
        seq_mat[:, i] = np.sum(temp_seq_mat2[:, i * 4:( i * 4 + 4)], axis=1)
        
    seq_mat[:, len_outputseq:] = temp_seq_mat2[:, -len_barcode * 4:]
    seq_mat = scipy.sparse.csr_matrix(seq_mat)
    emat_0 = np.zeros((4, len_seq))
    emat_0[:2, :len_outputseq] = mpathic.utils.RandEmat(len_outputseq, 2)
    emat_0[:, len_outputseq:] = mpathic.utils.RandEmat(20, 4)
    emat = MaximizeMI_memsaver(
        seq_mat,
        df.copy(),
        emat_0,
        wtrow,
        db=output_db,
        iteration=iteration,
        burnin=burnin,
        thin=thin,
        runnum=runnum,
        verbose=verbose,
        temp=temp
    )

    np.savetxt(output_mean, emat)
