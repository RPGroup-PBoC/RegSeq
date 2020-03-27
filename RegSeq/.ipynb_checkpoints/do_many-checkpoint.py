import os
import glob
ournames = ['tarphoPdataset_alldone_with_large']
os.system('/home/ubuntu/anaconda2/bin/mpathic learn_model -i /home/ubuntu/tarphoPdataset_alldone_with_large -o /home/ubuntu/results/tarphoPdataset_alldone_with_largeMCMC194 -lm IM --iteration 300000 --thin 60 -db /home/ubuntu/results/tarphoPdataset_alldone_with_largeMCMC194_0 --initialize rand') 
os.system('/home/ubuntu/anaconda2/bin/mpathic learn_model -i /home/ubuntu/tarphoPdataset_alldone_with_large -o /home/ubuntu/results/tarphoPdataset_alldone_with_largeMCMC195 -lm IM --iteration 300000 --thin 60 -db /home/ubuntu/results/tarphoPdataset_alldone_with_largeMCMC195_0 --initialize rand') 
