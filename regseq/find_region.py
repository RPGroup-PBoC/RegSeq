import pandas as pd
import scipy as sp
import numpy as np
from regseq import information

#we define some column names. The only purpose of that these are the columns that are used in
#the any energy matrices.
val_cols = ['val_A','val_C','val_G','val_T']

def find_edges(em, start, end, thresh=0.00025):
    newstart = np.nan
    newend = np.nan
    em_abs = np.abs(em)
    for s in range(start,end+15):
        if em_abs[s] > thresh:
            newstart = s
            break
    for s in range(end+14, end-1, -1):
        if em_abs[s] > thresh:
            newend = s
            break
        elif s == end:
            newend = s
    return newstart, newend



def select_region(temp_significant, gene, growth, pos, windowsize=15):
    """Find all transcription factor binding sites in a sequence given thresholded 
    information about expression shift for each base in case of mutation.
    
    If there is only a 1-4 base pair break in significance, the part
    is likely still part of one binding site. If none of those bases are signifanct, the 
    binding site has ended. 
    If the binding site is an activator, and none of the next four bases are also activator like
    we end the binding site (if instead they are repressor like, in other words if
    significant is > 0, then we can still assume the binding site has ended).
    
    Parameters
    ----------
    temp_significant : array-like
        List if significance for each base in the sequence. -1 if expression is
        significantly reduced, 1 if expression is significantly increased, and 0 if no
        significance for base
    gene : str
        Name of the gene
    growth : str
        Name of the growth condition
    pos : int
        Positions of the first base
    windowsize : Int, default 15
        Minimal size of binding site
        
    Returns
    -------
    outdf : Pandas DataFrame
    """
    # Shows whether or not the current base is part of a binding site or if its starting a new one.
    ongoing = 0
    
    info_length = len(temp_significant)
    outdf = pd.DataFrame(columns=['gene', 'growth', 'feat_num', 'start', 'end', 'type'])
    
    TF_type = 'None'
    num_feat = 0
    counter = 0
    
    for i in range(0, info_length):

        # Look for new binding site if not in one currently
        if not ongoing and i < info_length - windowsize:
            # Log new repressor binding site (expression is increased by mutation)
            if (temp_significant[i - 1] == 0 or temp_significant[i - 1] == -1 or i == 0) and temp_significant[i] == 1:
                start = i
                ongoing = 1
                TF_type = 'rep'
                
            # Log new activator binding site (expression is decreased by mutation)
            elif (temp_significant[i-1] == 0 or temp_significant[i-1] == 1 or i == 0) and temp_significant[i] == -1:
                start = i
                ongoing = 1
                TF_type = 'act'
                
        elif ongoing:
            # Compute significance of next 4 bases
            future_sum = temp_significant[i:i+4].sum()
            if (temp_significant[i] == 0):
                # Activator binding site ends if next bases are not activator like (negative significance)
                if future_sum > -.5 and (TF_type == 'act'):
                    end = i
                    
                    # Store information about binding site in data frame
                    outdf.loc[counter,['gene', 'growth', 'feat_num', 'start', 'end', 'type']] = \
                        [gene, growth, num_feat, start, end, TF_type]
                    
                    # Reset variables and increase counters
                    ongoing = 0
                    num_feat = num_feat + 1
                    TF_type = 'None'
                    counter = counter + 1
                    
                # Repressor binding site ends if next bases are not repressor like (positive significance)
                elif future_sum < .5 and (TF_type == 'rep'):
                    end = i
                    
                    # Store information about binding site in data frame
                    outdf.loc[counter, ['gene', 'growth', 'feat_num', 'start', 'end', 'type']] =\
                        [gene, growth, num_feat, start, end, TF_type]
                    
                    # Reset variables and increase counters
                    ongoing = 0
                    num_feat = num_feat + 1
                    TF_type = 'None'
                    counter = counter + 1
                else:
                    pass
            
            # Possibly adjacent binding sites 
            elif (temp_significant[i - 1] == 1 and temp_significant[i] == -1 and future_sum < .5):
                end = i
                # Store information about binding site in data frame
                outdf.loc[counter,['gene', 'growth', 'feat_num', 'start', 'end', 'type']] =\
                    [gene, growth, num_feat, start, end, TF_type]
                
                # Start new binding site
                start = i
                num_feat = num_feat + 1
                TF_type = 'act'
                counter = counter + 1
            elif (temp_significant[i-1] == -1 and temp_significant[i] == 1 and future_sum > -.5):
                end = i
                # Store information about binding site in data frame
                outdf.loc[counter,['gene','growth','feat_num','start','end','type']] =\
                    [gene, growth, num_feat, start, end, TF_type]
                
                # Start new binding site
                start = i
                num_feat = num_feat + 1
                TF_type = 'rep'
                counter = counter + 1

    return outdf


def do_sum2(s, windowsize=15):
    """Compute the sum of expression shifts of all sub-sequences of a given length
    
    Parameters
    ----------
    s : array-like
        List of expression shifts for single mutations.
    windowsize : Int, default 15
        Length of sequences to be summed over.
        
    Returns
    -------
    outarr : array-like
        Array of sums of expression shifts from one base of the next number of bases
        defined by windowsize.
    """
    
    info_length = len(s)
    outarr = np.zeros(info_length - windowsize)
    for i in range(info_length-windowsize):
        outarr[i] = s[i:(i+windowsize)].sum()
    return outarr


def find_region(file, gene, growth, windowsize=15, thresh=0.00025, old_format=False):
    """Find activator and repressor binding sites in sequence.
    
    Parameters
    ----------
    
    Returns
    -------
    
    """
    infofootprint = information.emat_to_information(
        file, 
        old_format=old_format,
        gene=gene
    )
    info, _, signs = information.footprint(infofootprint)
    pos = np.arange(-114, 44)
    counter = 0
    outdf = pd.DataFrame(columns=['gene', 'growth', 'feat_num', 'start', 'end', 'type'])
    info_length = len(info)
    
    em = info
    em_noabs = info * signs

    em_noabs = em_noabs*-1

    #sum into groupings of 15 base pairs so that we can see if large regions are statistically
    #signficant for expression.
    summedarr2 = do_sum2(em)
    summedarr_noabs = do_sum2(em_noabs)
    for q in range(len(summedarr_noabs)):
        if summedarr_noabs[q] > 0:
            summedarr_noabs[q] = 1
        else:
            summedarr_noabs[q] = -1

    #initialize array where we will store whether or not the outcome is signficant.
    is_significant = np.zeros(info_length - windowsize)
    for i in range(info_length - windowsize):
        #make a 99.5 percent confidence interval
        #to do this we will check from .5 to .95 because we want to know if 0
        #is in the range (.5% to 100%) or (0 to 99.5%) range, depending on
        #whether the shift is positive or negative.
        #is_significant[q,i] = summedarr2[i] > thresh*meanval*windowsize
        is_significant[i] = summedarr2[i] > thresh * windowsize
        is_significant[i] = is_significant[i]*summedarr_noabs[i]
        #if zero is in interval, the base is not signficant, otherwise it is.

    outdf_temp = select_region(is_significant, gene, growth, pos)
    for i,row in outdf_temp.iterrows():
        start = row['start']
        end = row['end']
        newstart,newend = find_edges(em,start,end)
        outdf.loc[counter, ['gene', 'growth', 'feat_num', 'start', 'end', 'type']] =\
            [row['gene'], growth, row['feat_num'], newstart, newend, row['type']]
        counter = counter + 1
    output_merged = merge_growths(outdf, windowsize, info_length=info_length)
    return output_merged





def merge_growths(df, windowsize, info_length=160):
    allgenes = list(set(df['gene']))
    processed_df = pd.DataFrame(columns=['gene', 'feat_num', 'start', 'end', 'type'])
    counter = 0
    for gene in allgenes:
        tempdf = df.loc[df['gene'] == gene]
        all_acts_df = tempdf.loc[tempdf['type'] == 'act']
        all_rep_df = tempdf.loc[tempdf['type'] == 'rep']
        act_array = np.zeros(info_length)
        rep_array = np.zeros(info_length)
        for i, row in all_acts_df.iterrows():
            act_array[int(row['start']):int(row['end'])] = 1
        for i, row in all_rep_df.iterrows():
            rep_array[int(row['start']):int(row['end'])] = 1

        outdf_temp_act = select_region(act_array, gene, pos=0, growth='combined')
        outdf_temp_rep = select_region(rep_array, gene, pos=0, growth='combined')
        if len(all_acts_df.index) != 0:
            for i,row in outdf_temp_act.iterrows():
                start = row['start']
                end = row['end']
                processed_df.loc[counter,['gene', 'feat_num', 'start', 'end', 'type']] = [row['gene'], row['feat_num'], start, end, 'act']
                counter = counter + 1
        if len(all_rep_df.index) != 0:
            for i,row in outdf_temp_rep.iterrows():
                start = row['start']
                end = row['end']
                processed_df.loc[counter,['gene','feat_num', 'start', 'end', 'type']] = [row['gene'], row['feat_num'], start, end, 'rep']
                counter = counter + 1

    return processed_df
