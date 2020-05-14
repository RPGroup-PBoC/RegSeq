import os

#first we will generate a key that matches barcode to mutated sequences

#take original .fastq file from sequencing constructs and construct
#the barcode to sequence key. This will operate on a subset of the whole
#sequencing dataset we originally used to avoid loading an extremely large
#file to github

#bdcR_matching_sequencing.fastq is the input file name for the bdcR gene
#bdcR_barcode_key is the output file name for the barcode matching.

os.system('''python create_key_to_match_sequence_to_barcode.py
             ../data/sequencing_data/mappingseqs.fastq 
             ../data/test_data/bdcR_barcode_key''')

#Next we create a dataset for the bdcR gene from sequencing the mRNAs and
#the DNA sequences in the construct libraryself.

#We use the fastx toolkit to do quality score filtering and barcode splitting

#More information on the fastx toolkit can be found on their online documentation
#an example command for filtering is
#./fastq_quality_filter -i infile.fastq -o filteredfile.fastq -Q33
#and for barcode splitting
#cat filteredfile.fastq | ./fastx_barcode_splitter.pl --bcfile bcfile.txt --bol --exact
# --prefix BI_split

#inline arguments:
#1: BI94_102_mRNA is the file name for the mRNA sequencing .fastq file.
#2: BI95_102_DNA is the file name for the DNA sequencing .fastq file.
#3: mapping key file
#4: Anaero is the growth condition of the experiment. This is used in constructing
#the output file name.
#5: the target gene name.
#group numbers for each gene can be found in the file genetogroupnum

os.system('''python matchdatasets.py ../data/sequencing_data/BI94_102_mRNA 
             ../data/sequencing_data/BI95_102_DNA ../data/test_data/bdcR_barcode_key ../data/sequencing_data/bdcRAnaerodataset bdcR''')

#The output file name will be bdcRAnaerodataset

# Now we can fit a model from which we can extract information footprints.

#The inline arguments are
#1: bdcRAnaerodataset is the file name for the input matchdataset.
#2: bdcRAnaerodataset_db is the output name for the database file for MCMC.
#3: bdcRAnaero_MCMC_mut is the output name for the mean file

#This MCMC command can take significant computational time. We typically use
#amazon AWS to run these longer computations. This command will be commented
#out for now, and replaced with a quick but less accurate least squares version
# of the inference.

#original command
os.system('''python learn_model_mut.py ../data/sequencing_data/bdcRAnaerodataset ../data/sequencing_data/bdcRAnaerodataset_db ../data/sequencing_data/bdcRAnaero_MCMC_mut''')

"""
#Least squares command
os.system('''python compute_least_squares.py ../data/sequencing_data/bdcRAnaerodataset
             ../data/test_data/bdcRAnaero_LS_mut bdcR''')

#Now we can create an information footprint from the previous models
os.system('''python output_information_footprint.py
             ../data/test_data/bdcRAnaero_LS_mut 
             ../data/test_data/bdcRAnaero_informationfootprint.pdf''')

#The following commands are used to fit energy matrices.... to do
"""