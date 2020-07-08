import pandas as pd
import scipy as sp
import numpy as np
from regseq import information

#we define some column names. The only purpose of that these are the columns that are used in
#the any energy matrices.
val_cols = ['val_A','val_C','val_G','val_T']

def find_edges(em, start, end, thresh=0.00020):
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




def select_region(temp_significant, gene, growth, windowsize=15):
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
    if ongoing:
        if TF_type == 'rep':
            end = i
            ongoing = 0
            outdf.loc[counter,['gene','growth','feat_num','start','end','type']] =\
                [gene, growth, num_feat, start, end, TF_type]
            num_feat = num_feat + 1
            TF_type = 'None'
            counter = counter + 1
        elif TF_type == 'act':
            end = i
            ongoing = 0
            #now that the current binding site has ended we will update the list of binding sites.
            outdf.loc[counter,['gene','growth','feat_num','start','end','type']] =\
                [gene, growth, num_feat, start, end, TF_type]
            num_feat = num_feat + 1
            TF_type = 'None'
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


def find_region(
    file, 
    gene, 
    growth, 
    windowsize=15, 
    thresh=0.00025, 
    old_format=False, 
    wildtype_file='../data/prior_designs/wtsequences.csv'
):
    """Find activator and repressor binding sites in sequence.

    Parameters
    ----------
    file : str
        path to file containing energy matrix
    gene : str,
        name of observed gene
    growth : str
        name of growth condition
    windowsize : int, default 15
        size of window to average over to find regions
    thresh : float, default 0.00025
        emperically determined threshold to identify significant sites
    old_format : boolean, default False
        if True, loads file with argument `delim_whitespace`
    wildtype_file : str, default '../data/prior_designs/wtsequences.csv'
        location of file containing wildtype promoter sequences
        
    Returns
    -------
    output_merged : Pandas DataFrame
        Dataframe with start and end of significant sites, and if activator or repressor
    """
    infofootprint = information.emat_to_information(
        file,
        old_format=old_format,
        gene=gene,
        wildtype_file=wildtype_file
    )

    genedf = pd.read_csv(wildtype_file)
    
    # Compute unsmoothed information footprint
    info, _, signs = information.footprint(infofootprint, windowsize=1)
    counter = 0
    outdf = pd.DataFrame(columns=['gene', 'growth', 'feat_num', 'start', 'end', 'type'])
    info_length = len(info)

    em = info
    em_noabs = info * signs


    pos_mat = np.zeros(160)
    neg_mat = np.zeros(160)
    for q in range(len(em_noabs)):
        if em_noabs[q] > 0:
            pos_mat[q] = np.abs(em_noabs[q])
        else:
            neg_mat[q] = np.abs(em_noabs[q])
    #sum into groupings of 15 base pairs so that we can see if large regions are statistically
    #signficant for expression.
    summedarr2 = do_sum2(pos_mat)
    summedarr2_neg = do_sum2(neg_mat)
    # initialize array where we will store whether or not the outcome is signficant.
    is_significant = np.zeros(info_length-windowsize)
    is_significant_neg = np.zeros(info_length-windowsize) 
    for i in range(info_length - windowsize):
        is_significant[i] = summedarr2[i] > thresh*windowsize
        is_significant_neg[i] = summedarr2_neg[i] > thresh*windowsize
        is_significant_neg[i] = is_significant_neg[i]*-1

    outdf_temp = select_region(is_significant,gene,growth)
    outdf_temp2 = select_region(is_significant_neg,gene,growth)
    for i,row in outdf_temp.iterrows():
        start = row['start']
        end = row['end']
        newstart,newend = find_edges(pos_mat, start, end)
        outdf.loc[counter,['gene','growth','feat_num','start','end','type']] =\
            [row['gene'],growth,row['feat_num'],newstart,newend,row['type']]
        counter = counter + 1
    for i,row in outdf_temp2.iterrows():
        start = row['start']
        end = row['end']
        newstart,newend = find_edges(neg_mat, start, end)
        outdf.loc[counter,['gene','growth','feat_num','start','end','type']] =\
            [row['gene'],growth,row['feat_num'],newstart,newend,row['type']]
        counter = counter + 1

    output_merged = merge_growths(outdf, windowsize, info_length=info_length)
    output_merged["start"] -= 115
    output_merged["end"] -= 115
    return output_merged





def merge_growths(df, windowsize, info_length=160):
    """Merge close binding sites.
    
    Parameters
    ----------
    df : Pandas DataFrame
        Dataframe of observed significant sites
    windowsize : int
        size of minimal number of bases between sites to be distinguished\
        
    Returns
    -------
    df : Pandas DataFrame
        Dataframe of merged sites
    """
    
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

        outdf_temp_act = select_region(act_array, gene, growth='combined')
        outdf_temp_rep = select_region(rep_array, gene, growth='combined')
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



def RNAP_sites(files, growth_conditions, alldf_RNAP, thresh=0.00025, windowsize=15):
    """Find RNAP binding sites in information footprints.
    
    Parameters
    ----------
    - files : list of files containing information footprints
    - growth_conditions : list of used growth conditions
    - alldf_RNAP: List of known RNAP binding sites for genes
    - thresh : Threshold used to identify significant sites
    - windowsize : number of bp to average over when smooting the information footprint
    """
    counter = 0
    outdf = pd.DataFrame(columns=["gene", "growth", "feat_num", "start", "end", "type"])
    for name in files:
        try:
            # Get gene name
            noleader = name.split("/")[-1]

            for x in growth_conditions:
                if (
                    x in noleader
                    and "500cAMP" not in x
                    and noleader.split("dataset")[0] != ""
                ):
                    growth = x
                    if (
                        growth == "acetate"
                        or growth == "42"
                        or growth == "fructose"
                        or growth == "_comb"
                    ):
                        gene = noleader.split("dataset")[0]
                    else:
                        gene = noleader.split(growth)[0]
            gene = gene.split("_")[-1]

            # load in information footprint
            df = pd.read_csv(name, delim_whitespace=True)
            length_info = len(df["info"])
            em = np.abs(np.array(list(df["info"])))
            em_noabs = np.array(list(df["info"]))

            # invert information values.
            em_noabs = em_noabs * -1

            # split the information footprint into two groups of repressor like and activator like bases.
            pos_mat = np.zeros(length_info)
            neg_mat = np.zeros(length_info)
            for q in range(len(em_noabs)):
                if em_noabs[q] > 0:
                    pos_mat[q] = np.abs(em_noabs[q])
                else:
                    neg_mat[q] = np.abs(em_noabs[q])

            # sum over 15 base pair windows
            summedarr2 = do_sum2(pos_mat)
            summedarr2_neg = do_sum2(neg_mat)

            # initialize array where we will store whether or not the average outcome beats threshold.
            is_thresh = np.zeros(length_info - windowsize)
            is_thresh_neg = np.zeros(length_info - windowsize)
            for i in range(length_info - windowsize):
                is_thresh[i] = summedarr2[i] > thresh * windowsize
                is_thresh_neg[i] = summedarr2_neg[i] > thresh * windowsize
                is_thresh_neg[i] = is_thresh_neg[i] * -1

            outdf_temp = select_region(is_thresh, gene, growth)
            outdf_temp2 = select_region(is_thresh_neg, gene, growth)
            for i, row in outdf_temp.iterrows():
                start = row["start"]
                end = row["end"]
                newstart, newend = find_edges(pos_mat, start, end)
                outdf.loc[
                    counter, ["gene", "growth", "feat_num", "start", "end", "type"]
                ] = [
                    row["gene"],
                    growth,
                    row["feat_num"],
                    newstart,
                    newend,
                    row["type"],
                ]
                counter = counter + 1
            for i, row in outdf_temp2.iterrows():
                start = row["start"]
                end = row["end"]
                newstart, newend = find_edges(neg_mat, start, end)
                outdf.loc[
                    counter, ["gene", "growth", "feat_num", "start", "end", "type"]
                ] = [
                    row["gene"],
                    growth,
                    row["feat_num"],
                    newstart,
                    newend,
                    row["type"],
                ]
                counter = counter + 1
        except Exception as e:
            print(e)
            
        
    outdf["RNAP"] = outdf.apply(find_RNAP, axis=1, args=(alldf_RNAP,))
    outdf.to_csv("all_features_auto20_for_check_split_" + str(thresh) + ".csv", index=False)
    z = merge_growths(outdf, 15)
    z["contains RNAP"] = z.apply(find_RNAP, axis=1, args=(alldf_RNAP,))
    z.to_csv(
        "all_features_auto_merged_20_scaled_for_RNAP_check_" + str(thresh) + ".csv",
        index=False
    )
    
    
def concat_RNAP_sites(df):
    """Concatenate all binding sites from DataFrame into list."""
    outsites = []
    for i, row in df.iterrows():
        q = row["sites"]
        q2 = q.split()
        q3 = [int(x) for x in q2]
        outsites = outsites + q3
    return list(set(outsites))


def find_RNAP(df, alldf_RNAP):
    """Identify if a given putative regulatory region contains an RNAP binding site.
    To do this we cross reference the location of the putative binding site with a list of identified
    RNAP sites."""
    gene = df["gene"]
    try:
        RNAP_sites = concat_RNAP_sites(alldf_RNAP[alldf_RNAP["gene"] == gene])
    except:
        RNAP_sites = []
    is_RNAP = False
    start = df["start"]
    end = df["end"]
    for RNAP in RNAP_sites:
        """Our RNAP sites are listed based on the location of the RNAP minus 35. We check whether that
        or the minus 10 region (which is why we have to possible regions), exist in the region."""
        if ((start - 2 < RNAP) and (end > RNAP + 5)) or (
            (start < RNAP + 29) and (end > RNAP + 21)
        ):
            is_RNAP = True
    return is_RNAP