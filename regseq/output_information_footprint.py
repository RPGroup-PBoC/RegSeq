#Import basic stuff
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import linear_model

#import the custom analysis software
import scipy as sp
import plot_informationfootprint as pli
import seaborn as sns
import sys

#The input file will be input via sys.argv. The arguments already
#1: input inferred matrix
#2: output information footprint file name

inputname = sys.argv[1]

emat = np.loadtxt(inputname)

'''The conversion script returns the information footprint smoothed with a window size = 3
(the values are averaged with their neighbors), and also a variable called shiftcolors. Shiftcolors is
a colorcoding scheme for showing whether mutation tends to increase gene expression (repressor like) or
decrease gene expression (activator like).'''
smoothinfo, shiftcolors = pli.main(emat)

fig,ax = plt.subplots(figsize=(10,2))
ax.set_ylabel('Information (bits)',fontname='DejaVu Sans',fontsize=12)
ax.set_xlabel('position',fontname='DejaVu Sans',fontsize=12)
ax.bar(range(-114,43),np.abs(smoothinfo),color=shiftcolors)
plt.savefig(sys.argv[2],format='pdf')
