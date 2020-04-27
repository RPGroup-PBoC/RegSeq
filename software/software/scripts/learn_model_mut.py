from __future__ import division
#Our standard Modules
import argparse
import numpy as np
import scipy as sp
import sys
#Our miscellaneous functions
import pandas as pd
import mpathic.utils as utils
from sklearn import linear_model
import mpathic.EstimateMutualInfoforMImax as EstimateMutualInfoforMImax
import pymc
import mpathic.stepper as stepper
import os
from mpathic import SortSeqError
import mpathic.io as io
import mpathic.gauge as gauge
import mpathic.qc as qc
import pdb
from mpathic import shutthefuckup
import mpathic.numerics as numerics
from sklearn.preprocessing import StandardScaler
import scipy.sparse   

#this function is used to infer the effect of mutation on gene expression for
#each base pair using MCMC. This is crucial for converting to an information 
#footprint.

class GaugePreservingStepper(pymc.Metropolis):
    """Perform monte carlo steps that preserve the following choise of gauge:
    sum of elements in each column = 0, overall matrix norm = 1."""
    def __init__(self,stochastic):
        pymc.Metropolis.__init__(self,stochastic)

    def propose(self):
        # to avoid overwriting weirdness
        
        emat_temp = self.stochastic.value.copy()
        emat_temp1 = emat_temp[:2,:160].copy()
        emat_temp2 = emat_temp[:,160:].copy()
        # number of columns in matrix
        num_col1 = emat_temp1.shape[1]
        # first choice a random 4L dimensional vector
        r = sp.random.standard_normal(emat_temp1.shape)
        # dot product of proposed direction with current vector
        lambda_0 = sp.sum(emat_temp1*r)
        # dot products of each column with basis vectors
        lambda_vec = 0.5*sp.sum(r,axis=0)

        s = sp.zeros_like(r)
        # s is a vector based on r, but along the surface of the hypersphere
        for j in range(emat_temp1.shape[1]):
            s[:,j] = r[:,j] - lambda_0*emat_temp1[:,j] - lambda_vec[j]*(0.5*sp.ones(emat_temp1.shape[0]))
        dx = self.adaptive_scale_factor*s/sp.sqrt(sp.sum(s*s))
        emat_temp1_new = (emat_temp1+dx)


        num_col = emat_temp2.shape[1]
        # first choice a random 4L dimensional vector
        r = sp.random.standard_normal(emat_temp2.shape)
        # dot product of proposed direction with current vector
        lambda_0 = sp.sum(emat_temp2*r)
        # dot products of each column with basis vectors
        lambda_vec = 0.5*sp.sum(r,axis=0)

        s = sp.zeros_like(r)
        # s is a vector based on r, but along the surface of the hypersphere
        for j in range(emat_temp2.shape[1]):
            s[:,j] = r[:,j] - lambda_0*emat_temp2[:,j] - lambda_vec[j]*(0.5*sp.ones(emat_temp2.shape[0]))
        dx = self.adaptive_scale_factor*s/sp.sqrt(sp.sum(s*s))
        emat_temp2_new = (emat_temp2+dx)
        emat_temp[:2,:160] = emat_temp1_new
        emat_temp[:,160:] = emat_temp2_new
        self.stochastic.value = (emat_temp)/sp.sqrt(sp.sum((emat_temp1_new)**2) + sp.sum((emat_temp2_new)**2))


def MaximizeMI_memsaver(
        seq_mat,df,emat_0,wtrow,db=None,burnin=1000,iteration=30000,thin=10,
        runnum=0,verbose=False,temp=4200):
    '''Performs MCMC MI maximzation in the case where lm = memsaver'''    
    
    #sets the number of sequences in the dataset.
    n_seqs = len(df.index)
    #set up the model for MCMC inference with pymc2.
    @pymc.stochastic(observed=True,dtype=pd.DataFrame)
    def pymcdf(value=df):
        return 0
    @pymc.stochastic(dtype=float)
    def emat(p=pymcdf,value=emat_0):
        #on each iteration we will need to convert a L_bases array to a 
        #2xL array so we can use the gauge fixing function
        len_seq = len(df.loc[0,'seq'])
        len_barcode = 20
        len_outputseq = len_seq - len_barcode 
        total_params = len_outputseq + len_barcode*4
        temp_val = np.zeros(total_params)
        temp_val[:len_outputseq] = value[1,:len_outputseq] - value[0,:len_outputseq]
        temp_val[len_outputseq:] = value[:,len_outputseq:].ravel(order='F')
        #we now evaluate each sequence in the dataset with the temporary model.         
        p['val'] = seq_mat*temp_val      
        #evaluate mutual information likelihood is the log of this.             
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

df = pd.io.parsers.read_csv(sys.argv[1],delim_whitespace=True)

temp_seq_mat,wtrow = numerics.dataset2mutarray(df.copy(),'MAT')
temp_seq_mat2 = temp_seq_mat.toarray()


seq_mat = np.zeros((temp_seq_mat2.shape[0],total_params))
for i in range(len_outputseq):
    seq_mat[:,i] = np.sum(temp_seq_mat2[:,i*4:(i*4+4)],axis=1)
seq_mat[:,len_outputseq:] = temp_seq_mat2[:,-80:]
seq_mat = scipy.sparse.csr_matrix(seq_mat)
emat_0 = np.zeros((4,len_seq))
emat_0[:2,:len_outputseq] = utils.RandEmat(len_outputseq,2)
emat_0[:,len_outputseq:] = utils.RandEmat(20,4)
emat = MaximizeMI_memsaver(
                seq_mat,df.copy(),emat_0,wtrow,db=sys.argv[2],iteration=600000,burnin=1000,
                thin=60,runnum=0,verbose=True,temp=4200)

np.savetxt(sys.argv[3],emat)

