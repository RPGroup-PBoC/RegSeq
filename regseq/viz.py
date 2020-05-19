from regseq import information
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import pandas as pd
import seaborn as sns


def footprint(matrix, output_file=None):
    """ Plot information footprint.
    
    Footprint is smoothed with a window of size=3. Bars are colored by
    increased or decreased expression.
    
    Parameters
    ----------
    matrix : str
        File path for energy matrix
    output_file : 
        File path to store plot

    Returns
    -------
    plt : matplotlib.pyplot object
        Information Footprint
    
    """
    
    emat = np.loadtxt(matrix)
    smoothinfo, shiftcolors = information.footprint(emat)
    fig,ax = plt.subplots(figsize=(10,2))
    ax.set_ylabel('Information (bits)',fontname='DejaVu Sans',fontsize=12)
    ax.set_xlabel('position',fontname='DejaVu Sans',fontsize=12)
    ax.bar(range(-114,43),np.abs(smoothinfo),color=shiftcolors)
    if not output_file == None:
        plt.savefig(output_file,format='pdf')
    return plt


def footprint_from_emat(file, output_file=None):
    """ Plot information footprint.
    
    Footprint is smoothed with a window of size=3. Bars are colored by
    increased or decreased expression.
    
    Parameters
    ----------
    matrix : str
        File path for energy matrix
    output_file : str
        If not None, path where figure is saved as pdf.

    Returns
    -------
    plt : matplotlib.pyplot object
        Information Footprint
    
    """
    
    arr = information.emat_to_information(file)
    smoothinfo, shiftcolors = information.footprint(arr)
    fig,ax = plt.subplots(figsize=(10,2))
    ax.set_ylabel('Information (bits)',fontname='DejaVu Sans',fontsize=12)
    ax.set_xlabel('position',fontname='DejaVu Sans',fontsize=12)
    ax.bar(range(-114,43),np.abs(smoothinfo),color=shiftcolors)
    if not output_file == None:
        plt.savefig(output_file,format='pdf')
    return plt


def logo(file, limit=(), min_beta=.001, max_beta=100, num_betas=1000, output_file=None,save=False):
    """ Plot logo using Logomaker package.
    
    Parameters
    ----------
    file : str
        file containing energy matrix
    limit : Tuple, default ()
        range of positions plotted
    out
    
    """
    
    # Load in a binding site matrix.
    arraydf = pd.read_csv(file, index_col="pos")
    
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
    
    if save:
        plt.savefig("../figures/"+file.split("/")[-1].split(".")[0]+'_logo.pdf',format='pdf')
        
    return binding_logo


def energy_matrix(file, limit=(), save=False):
    """ Plot energy matrix.
    
    Parameters
    ----------
    file : str
        File for energy matrix.
    limit : Tuple, default ()
        range of positions plotted
    save : boolean, default False
        If True, save result with "_array.pdf" attached in RegSeq/figures/.
        
    Returns 
    -------
    ax : atplotlib.axes._subplots.AxesSubplot
        Plotted energy matrix
    """
    
     # Load in a binding site matrix.
    tempdf = pd.read_csv(file, index_col="pos")
    
    # Apply Limit if given
    if len(limit) != 0:
        if len(limit) != 2:
            raise RuntimeError("limit must have length 2.")
        else:
            arraydf = arraydf.iloc[limit[0]:limit[1]+1]
            
    # Convert to numpy array
    temparr = np.array(tempdf[['val_A','val_C','val_G','val_T']])
    
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
    ax.set_yticklabels(['','A','C','G','T'])
    ax.grid(False)
    for item in [fig, ax]:
        item.patch.set_visible(True)
    
    if save:
        plt.savefig("../figures/"+file.split("/")[-1].split(".")[0]+'_array.pdf',format='pdf')
    
    return ax


def mass_spec(file, inname="A5", save=False):
    def check_DNA(s):
        '''Return only proteins which have DNA binding activity.'''
        with open('../data/massspec/DNAbinding_genenames.txt') as f:
            genenames = f.read()
            genenames = genenames.split(',\n')
        if s in genenames:
            return True
        else:
            return False
    
    
    df = pd.io.parsers.read_csv(file, sep='\t')
    # Column name that contains the normal
    good_column = 'Ratio H/L normalized ' + inname
    # Extract the only necessary columns, protein name and normalized ratio
    enrichment = df[['Protein names',good_column]]
    
    enrichment2 = enrichment.dropna()
    enrichment2 = enrichment2.sort_values(by=good_column,ascending=False)
    
    goodrows = enrichment2['Protein names'].apply(check_DNA)
    enrichment3 = enrichment2[goodrows]
    
    fig,ax = plt.subplots()
    ax.set_xlabel('')
    ax.set_xticks([])
    ax.set_ylabel('enrichment')
    sns.stripplot(data=list(enrichment3[good_column]),jitter=True,size=12)
    ax.set_xticks([])
    if save == True:
        plt.savefig(inname + '_output.eps', format='eps')
        
    return ax
