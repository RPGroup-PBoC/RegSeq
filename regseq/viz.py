from regseq import information
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import pandas as pd
import seaborn as sns


def footprint(
    matrix, 
    output_file=None, 
    wildtype_file='../data/prior_designs/wtsequences.csv', 
    gene=None, 
    show_real_pos=False,
    windowsize=3
):
    """ Plot information footprint.
    
    Footprint is smoothed with a window of size windowsize. Bars are colored by
    increased or decreased expression.
    
    Parameters
    ----------
    matrix : str
        File path for energy matrix
    output_file : 
        File path to store plot
    wildtype_file : str, default None
        Path to file containing information about gene TSS and direction of transcription
    gene : str, default None
        Gene name in to find TSS and direction of transcription in gene file.
    show_real_pos : boolean, default False
        If True, positions will be given as in the genome. If False, positions are shown
        relative to TSS.
    windowsize: Int, default 3
        Size of sliding window used to average mutual information
    
    Returns
    -------
    ax : matplotlib.pyplot object
        Information Footprint
    
    """
   
    if gene == None:
        raise RuntimeError("Give gene name to find TSS and direction of trancription.")
    
    genedf = pd.read_csv(wildtype_file)
    
    direction = 0 if genedf.loc[genedf["name"] == gene, "rev"].values == "fwd" else 1
    
    cut = int((windowsize - 1) / 2)
    if direction == 0:
        #x = np.arange(-45 + cut, 115 - cut)
        # We always plot upstream to downstream
        x = np.arange(-115 + cut, 45 - cut)
    else:
        x = np.arange(-115 + cut, 44 - cut)
    
    if show_real_pos:
        TSS = int(genedf.loc[genedf["name"] == gene, "start_site"].values)
        x = x + TSS
        
    try :
        emat = np.loadtxt(matrix)
    except ValueError:
        emat = np.loadtxt(matrix, skiprows=1)[:, 1]
        x = np.arange(0, len(emat)-2*cut)

    smoothinfo, shiftcolors = information.footprint(emat, windowsize=windowsize)[:2]
    fig,ax = plt.subplots(figsize=(10,2))
    ax.set_ylabel('Information (bits)',fontname='DejaVu Sans',fontsize=12)
    ax.set_xlabel('position',fontname='DejaVu Sans',fontsize=12)
    ax.bar(x,np.abs(smoothinfo),color=shiftcolors)
    if not output_file == None:
        plt.savefig(output_file)
    return ax


def footprint_from_emat(
    file, 
    output_file=None, 
    old_format=False, 
    wildtype_file='../data/prior_designs/wtsequences.csv', 
    gene=None, 
    show_real_pos=False,
    windowsize=3
):
    """ Plot information footprint.
    
    Footprint is smoothed with a window of size=3. Bars are colored by
    increased or decreased expression.
    
    Parameters
    ----------
    matrix : str
        File path for energy matrix
    output_file : str
        If not None, path where figure is saved as pdf.
    old_format : boolean, default False
        Determines if file is loaded from old format.
    wildtype_file : str, default None
        Path to file containing information about gene TSS and direction of transcription
    gene : str, default None
        Name of gene has to be given in case of old file format
    show_real_pos : boolean, default False
        If True, positions will be given as in the genome. If False, positions are shown
        relative to TSS.
    windowsize: Int, default 3
        Size of sliding window used to average mutual information
        
    Returns
    -------
    plt : matplotlib.pyplot object
        Information Footprint
    
    """
    if gene == None:
        raise RuntimeError("Give gene name to find TSS and direction of trancription.")
    
    genedf = pd.read_csv(wildtype_file)
    
    
    direction = 0 if genedf.loc[genedf["name"] == gene, "rev"].values == "fwd" else 1
    
    cut = int((windowsize - 1) / 2)
    
    if direction == 0:
        #x = np.arange(-45 + cut, 115 - cut)
        # We always plot upstream to downstream
        x = np.arange(-115 + cut, 45 - cut)
    else:
        x = np.arange(-115 + cut, 45 - cut)
    
    if show_real_pos:
        TSS = int(genedf.loc[genedf["name"] == gene, "start_site"].values)
        x = x + TSS
        
    arr = information.emat_to_information(file, old_format=old_format, gene=gene)
    smoothinfo, shiftcolors = information.footprint(arr, windowsize=windowsize)[:2]
    fig,ax = plt.subplots(figsize=(10,2))
    ax.set_ylabel('Information (bits)',fontname='DejaVu Sans',fontsize=12)
    ax.set_xlabel('position',fontname='DejaVu Sans',fontsize=12)
    ax.bar(x, np.abs(smoothinfo), color=shiftcolors)
    if not output_file == None:
        plt.savefig(output_file,format='pdf')
    return plt


def logo(
    file, 
    limit=(),
    output_file=None, 
    old_format=False,
    min_beta=.001, 
    max_beta=100., 
    num_betas=1000
   ):
    """Plot sequence logos using logomaker.
    
    Parameters
    ----------
    file : str
        path to file containing energy matrix
    limit : Tuple, default ()
        first and last base of sequence that is converted into logo
    output_file : str, default None
        path where plot is saved to, if not None
    old_format : boolean, default False
        If True, file is loaded with extra argument "delim_whitespace"
    min_beta : float, default 0.001
        minimal scaling factor
    max_beta : float, default 100
        maximal scaling factor
    max_beta : int, default 1000
        number of tested scaling factors
    Returns
    -------
    binding_logo : logomaker.src.Logo.Logo
    """
    
    # Load in a binding site matrix.
    arraydf = pd.read_csv(file, index_col="pos", delim_whitespace=old_format)
    
    # Rename columns to be useable by the logomaker package
    arraydf = arraydf.rename(columns={'val_A':'A','val_C':'C','val_G':'G','val_T':'T'})
    
    if len(limit) != 0:
        if len(limit) != 2:
            raise RuntimeError("limit must have length 2.")
        else:
            arraydf = arraydf.iloc[limit[0]:limit[1]+1]
    # finding scaling factor
    target_info = len(arraydf.index)
    beta = information.get_beta_for_effect_df(
        arraydf,
        target_info,
        min_beta=min_beta,
        max_beta=max_beta,
        num_betas=num_betas
    )

    # use logomaker to convert energy matrix to information matrix
    binding_info = logomaker.transform_matrix(df=beta*arraydf,from_type='weight',to_type='information')
    binding_logo = logomaker.Logo(binding_info,
                             #font_name='Stencil Std',
                             vpad=.1,
                             width=.8)

    # style using Logo methods
    binding_logo.style_spines(visible=False)
    binding_logo.style_spines(spines=['left', 'bottom'], visible=True)
    binding_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    # style using Axes methods
    binding_logo.ax.set_ylabel("Information (bits)", labelpad=-1)
    binding_logo.ax.xaxis.set_ticks_position('none')
    binding_logo.ax.xaxis.set_tick_params(pad=-1)
    binding_logo.ax.grid(False)
    binding_logo.ax.set_xticklabels(np.arange(limit[0], limit[1]+1))
    
    if output_file != None:
        plt.savefig(output_file)
        
    return binding_logo


def energy_matrix(file, limit=(), output_file=None, old_format=False):
    """ Plot energy matrix.
    
    Parameters
    ----------
    file : str
        File for energy matrix.
    limit : Tuple, default ()
        range of positions plotted
    output_file : str, default None
        path where plot is saved to, if not None
        
    Returns 
    -------
    ax : atplotlib.axes._subplots.AxesSubplot
        Plotted energy matrix
    """
    
     # Load in a binding site matrix.
    arraydf = pd.read_csv(file, index_col="pos", delim_whitespace=old_format)
    
    # Apply Limit if given
    if len(limit) != 0:
        if len(limit) != 2:
            raise RuntimeError("limit must have length 2.")
        else:
            arraydf = arraydf.iloc[limit[0]:limit[1]+1]
            
    # Convert to numpy array
    temparr = np.array(arraydf[['val_A','val_C','val_G','val_T']])
    
    # find maximum absolute to center colorbar
    maximum = np.max(np.abs(temparr))
    #now plot using matplotlib
    fig,ax = plt.subplots(figsize=((10,2)))
    plt.imshow(
        temparr.T,
        aspect='auto', 
        interpolation='nearest',
        cmap='coolwarm',
        vmin=-maximum,
        vmax=maximum
    )
    plt.colorbar()
    plt.xlabel('Position')
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(['A','C','G','T'])
    ax.set_xticks(np.arange(0, len(temparr), step=5))
    if len(limit) != 0:
        ax.set_xticklabels(np.arange(limit[0], limit[1]+1, step=5))
    ax.grid(False)

    for item in [fig, ax]:
        item.patch.set_visible(True)
    
    if output_file != None:
        plt.savefig(output_file)
    
    return ax


def mass_spec(file, good_column="Ratio H/L normalized", output_file=None, filtered=True):
    """Jitter plot for protein enrichment.
    
    Parameters
    ---------
    file : str
        path to mass spec file
    good_column : str, default "Ratio H/L normalized"
        plotted column
    output_file : str, default None
        path where plot is saved to, if not None
    filtered : boolean, default True
        If False, check which proteins bind DNA
        
    Returns
    -------
    ax : pyplot axes
        Plot
    """
    
    def check_DNA(s):
        '''Return only proteins which have DNA binding activity.'''
        with open('../data/massspec/DNAbinding_genenames.txt') as f:
            genenames = f.read()
            genenames = genenames.split(',\n')
        if s in genenames:
            return True
        else:
            return False
    
    
    df = pd.io.parsers.read_csv(file, sep=',')
    enrichment = df[['Protein names',good_column]]
    
    enrichment2 = enrichment.dropna()
    enrichment2 = enrichment2.sort_values(by=good_column,ascending=False)
    
   
    fig,ax = plt.subplots()
    ax.set_xlabel('')
    ax.set_xticks([])
    ax.set_ylabel('enrichment')
    if filtered:
        sns.stripplot(data=list(enrichment2[good_column]),jitter=True,size=12)
    else:
        goodrows = enrichment2['Protein names'].apply(check_DNA)
        enrichment3 = enrichment2[goodrows]
        sns.stripplot(data=list(enrichment2[good_column]),jitter=True,size=12)
        
    ax.set_xticks([])
    
    if output_file != None:
        plt.savefig(output_file)
        
    return ax
