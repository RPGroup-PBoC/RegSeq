import pandas as pd
import scipy as sp
import numpy as np

#we define some column names. The only purpose of that these are the columns that are used in
#the any energy matrices.
val_cols = ['val_A','val_C','val_G','val_T']

def find_edges(em,start,end,thresh=0.00025):
    newstart = np.nan
    newend=np.nan
    em_abs = np.abs(em)
    for s in range(start,end+15):
        if em_abs[s] > thresh:
            newstart = s
            break
    for s in range(end+14,end-1,-1):
        if em_abs[s] > thresh:
            newend = s
            break
        elif s == end:
            newend = s
    return newstart, newend

def select_region(temp_significant,gene,growth,thresh=0.00025):
    #initialize the variable 'ongoing'. this shows whether or not the current base is part of a binding site
    #or if its starting a new one.
    info_length = len(temp_significant)
    outdf = pd.DataFrame(columns=['gene','growth','feat_num','start','end','type'])
    ongoing = 0
    TF_type = 'None'
    num_feat = 0
    counter = 0
    #loop through only 145 base pairs (we don't go to 160 becuase the regions are 15 base pairs long.)
    for i in range(0,info_length-windowsize):

        #if we are not currently part of a binding site then we do this...
        if not ongoing:
            #if we have a new significant base pair we will start to log a new binding site.
            if (temp_significant[i-1] == 0 or temp_significant[i-1] == -1 or i == 0) and temp_significant[i] == 1:
                start = i
                ongoing = 1
                #we checked whether temp_signficant was '1' which means the effect of mutation on expression
                #is positive. This would indicate that the transcription factor type is 'repressor'
                TF_type = 'rep'
            elif (temp_significant[i-1] == 0 or temp_significant[i-1] == 1 or i == 0) and temp_significant[i] == -1:
                start = i
                ongoing = 1
                #we checked whether temp_signficant was '-1' which means the effect of mutation on expression
                #is negative. This would indicate that the transcription factor type is 'activator'
                TF_type = 'act'
        elif ongoing:
            #if we are currently within an ongoing binding site we need to see whether or not the binding site
            #has ended at the current base. To do that we first check if the new base is not significant.
            future_sum = temp_significant[i:i+4].sum()
            if (temp_significant[i] == 0):
                #next, if there is only a 1-4 base pair break in which bases are significant, the whole
                #thing is probably still part of one binding site. So we check whether or not the next 4
                #base pairs are not significant, if they are not we declare the binding site ended.
                #if the binding site is an activator, and none of the next four bases are also activator like
                #we end the binding site (if instead they are repressor like, in other words if
                #significant is > 0, then we can still assume the binding site has ended)
                if future_sum > -.5 and (TF_type == 'act'):
                    end = i
                    ongoing = 0

                    #now that the current binding site has ended we will update the list of binding sites.
                    outdf.loc[counter,['gene','growth','feat_num','start','end','type']] = [gene,growth,num_feat,start,end,TF_type]
                    num_feat = num_feat + 1
                    TF_type = 'None'
                    counter = counter + 1
                #now do the same in the case that the current binding site is a repressor.
                elif future_sum < .5 and (TF_type == 'rep'):
                    end = i
                    ongoing = 0
                    outdf.loc[counter,['gene','growth','feat_num','start','end','type']] = [gene,growth,num_feat,start,end,TF_type]
                    num_feat = num_feat + 1
                    TF_type = 'None'
                    counter = counter + 1
                else:
                    pass
            elif (temp_significant[i-1] == 1 and temp_significant[i] == -1 and future_sum < .5):
                end = i
                outdf.loc[counter,['gene','growth','feat_num','start','end','type']] = [gene,growth,num_feat,start,end,TF_type]
                start = i
                num_feat = num_feat + 1
                TF_type = 'act'
                counter = counter + 1
            elif (temp_significant[i-1] == -1 and temp_significant[i] == 1 and future_sum > -.5):
                end = i
                outdf.loc[counter,['gene','growth','feat_num','start','end','type']] = [gene,growth,num_feat,start,end,TF_type]
                start = i
                num_feat = num_feat + 1
                TF_type = 'rep'
                counter = counter + 1

    return outdf

windowsize = 15
def do_sum2(s):
    '''this function does a summation 15 base pairs from the experession shifts models.
    We will be seeing if the summation of 15 consecutive base pairs are significant for gene
    expression (99% confidence interval).'''
    info_length = len(s)
    outarr = np.zeros(info_length - windowsize)
    for i in range(info_length-windowsize):
        outarr[i] = s[i:(i+windowsize)].sum()
    return outarr

def merge_growths(df,info_length=160):
    allgenes = list(set(df['gene']))
    processed_df = pd.DataFrame(columns=['gene','feat_num','start','end','type'])
    counter = 0
    for gene in allgenes:
        tempdf = df.loc[df['gene'] == gene]
        all_acts_df = tempdf.loc[tempdf['type']=='act']
        all_rep_df = tempdf.loc[tempdf['type']=='rep']
        act_array = np.zeros(info_length)
        rep_array = np.zeros(info_length)
        for i,row in all_acts_df.iterrows():
            act_array[int(row['start']):int(row['end'])] = 1
        for i,row in all_rep_df.iterrows():
            rep_array[int(row['start']):int(row['end'])] = 1

        outdf_temp_act = select_region(act_array,gene,growth='combined')
        outdf_temp_rep = select_region(rep_array,gene,growth='combined')
        if len(all_acts_df.index) != 0:
            for i,row in outdf_temp_act.iterrows():
                start = row['start']
                end = row['end']
                processed_df.loc[counter,['gene','feat_num','start','end','type']] = [row['gene'],row['feat_num'],start,end,'act']
                counter = counter + 1
        if len(all_rep_df.index) != 0:
            for i,row in outdf_temp_rep.iterrows():
                start = row['start']
                end = row['end']
                processed_df.loc[counter,['gene','feat_num','start','end','type']] = [row['gene'],row['feat_num'],start,end,'rep']
                counter = counter + 1

    return processed_df




def find_region(df,gene,growth):
        #information threshhold
        thresh = 0.00025
        counter = 0
        outdf = pd.DataFrame(columns=['gene','growth','feat_num','start','end','type'])
        info_length = len(df.index)
        #we will use pymc to load in the database of all MCMC steps.
        em = np.abs(np.array(list(df['info'])))
        em_noabs = np.array(list(df['info']))
        #meanval = np.average(em)
        #meanval_noabs = np.average(em_noabs)
        #if meanval_noabs > 0:
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
        is_significant = np.zeros(info_length-windowsize)
        for i in range(info_length-windowsize):
            #make a 99.5 percent confidence interval
            #to do this we will check from .5 to .95 because we want to know if 0
            #is in the range (.5% to 100%) or (0 to 99.5%) range, depending on
            #whether the shift is positive or negative.
            #is_significant[q,i] = summedarr2[i] > thresh*meanval*windowsize
            is_significant[i] = summedarr2[i] > thresh*windowsize
            is_significant[i] = is_significant[i]*summedarr_noabs[i]
            #if zero is in interval, the base is not signficant, otherwise it is.

        #both plot and save signficance results. signficant sites will end 15 base pairs after the last
        #shown signficant base because we average over 15 bases.
        #plt.plot(range(145),is_significant)
        #plt.show()
        outdf_temp = select_region(is_significant,gene,growth,thresh=thresh)
        for i,row in outdf_temp.iterrows():
            start = row['start']
            end = row['end']
            newstart,newend = find_edges(em,start,end)
            outdf.loc[counter,['gene','growth','feat_num','start','end','type']] = [row['gene'],growth,row['feat_num'],newstart,newend,row['type']]
            counter = counter + 1
        output_merged = merge_growths(outdf,info_length=info_length)
        return output_merged
