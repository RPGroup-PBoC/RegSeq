import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import panel as pn
import glob

allnames = glob.glob('*infoshift*')

possible_growths = ['metal','LB','xylara','xanth','SS','arabinose','heat']
growths = []
thegenes = []
#output_dict
for name in allnames:
    for x in possible_growths:
        if x in name:
            growths.append(x)
            thegenes.append(name.split(x)[0])
            #infodf = pd.io.parsers.read_csv(name,delim_whitespace=True)
            #infoarr = np.array(df.loc[:,'info'])
growthsarr = np.array(growths)
allnamesarr = np.array(allnames)
geneset = set(thegenes)
outdict = {}
for gene in geneset:
    indexes = np.nonzero(np.array(thegenes) == gene)[0]
    genegrowths = growthsarr[indexes]
    genenames = allnamesarr[indexes]
    tempdict = {}
    for i,name in enumerate(genenames):
        infodf = pd.io.parsers.read_csv(name,delim_whitespace=True)
        infoarr = np.array(infodf.loc[:,'info'])
        tempdict[genegrowths[i]] = infoarr
    outdict[gene] = tempdict    

select_gene = pn.widgets.Select(name='Select gene', options=list(outdict.keys()))
select_growth = pn.widgets.Select(name='Select growth', options=list(outdict[select_gene.value].keys()))

def plot_infoshift(gene='aphA',growth='arabinose'):
    plt.clf()
    plt.cla()
    s = outdict[gene][growth]
    
    colorinputs = np.zeros((160))
    for i in range(160):
        if s[i] < 0:
            colorinputs[i] = 0.0
        else:
            colorinputs[i] = 1.0
    fig,ax = plt.subplots(figsize=(10,2))
    shiftcolors = plt.cm.bwr(colorinputs)
    s_abs = np.abs(s)
    ax.bar(range(160),np.abs(s),color=shiftcolors)
    #fig = plt.bar(range(160),np.abs(s),color=shiftcolors)
    plt.close(fig)
    return fig

m = pn.interact(plot_infoshift,gene=select_gene,growth=select_growth)
mout = pn.Row(pn.Column(m[0]), pn.Column(m[1])).servable()
mout.save('/home/bill/test.html')

