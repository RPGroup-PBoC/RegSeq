import sys
import glob
import os

#This script allows you to format the inputs for matchdatasets-all.py
rnanamesfwd = ['98']
plasnamesfwd = ['95']
growthname = ['glpR']
#choose which groups to fit, all 18 groups range from 101 to 118. 
allgroups = [111,113]


for i in range(len(rnanamesfwd)):
    for z in allgroups:
        '''these are the names of the input files. you may have to change this
        if you decide on a different naming scheme'''
        rnainput = 'BI' + rnanamesfwd[i] + '_' + str(z)
        plasinput = 'BI' + plasnamesfwd[i] + '_' + str(z)
        os.system('''python matchdatasets-all_glpR.py %s %s %s''' %(rnainput,plasinput,growthname[i]))
 
