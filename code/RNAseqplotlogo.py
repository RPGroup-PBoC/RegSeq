import os
import glob

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import the project utils
import sys

# Import matplotlib stuff for plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Logo-generating module
import anylogo

gc = .5
background_array =pd.DataFrame( [[(1-gc)/2,gc/2,gc/2,(1-gc)/2]])

def plot_logo(energy_df,clip=False,invert=False,scaledown=None):

    energy_df = energy_df.rename(columns={'val_A':'A','val_C':'C','val_G':'G','val_T':'T'})
    #if designated clip off tag
    if clip == 'clip':
        total_length = len(energy_df.index)
        energy_df = energy_df.loc[:total_length-21,:]

    #invert if necessary
    if invert == 'invert':
        energy_df = energy_df*-1

    energy_df = energy_df[['A','C','G','T']]

    #no scaling factor
    energy_df_scaled = energy_df
    if scaledown:
        energy_df_scaled = energy_df_scaled/float(scaledown)

    background_df = pd.DataFrame(pd.np.tile(background_array,
                        (len(energy_df_scaled), 1)), columns=['A','C','G','T'])

    emat_min = -2
    emat_max = 2
    mid_val=0

    fig,ax = plt.subplots()

    logo = anylogo.effect_df_to_prob_df(energy_df_scaled,background_df,1)
    pd.set_option('max_colwidth',int(1e8)) #makes sure seq columns aren't shortened

    anylogo.draw(ax, effect_df=energy_df_scaled, logo_type='information',
                 background = background_df,
                 use_transparency=False,find_beta=False)
    ax.set_ylim([0,2])
    return fig
