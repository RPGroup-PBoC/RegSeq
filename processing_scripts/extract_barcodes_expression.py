import pandas as pd
import git


# Barcode identification
barcode_dict = {
    "AGAG":"LB_RNA",
    "CAAG":"LB_DNA",
    "TCTA":"heatshock_RNA",
    "ATGC":"heatshock_DNA"
}

# Arrays to store barcodes
barcode1_arr = []
barcode2_arr = []
barcode_promoter_arr = []

# Find project parental directory
repo = git.Repo('./', search_parent_directories=True)
homedir = repo.working_dir

filename = 'LB_heatshock_w_barcode.txt'

with open(f'{homedir}/data/sequencing_data/processed_sequencing/{filename}') as f:
    while True:
        line = f.readline()
        if not line: 
            break
        ind = line.find("TATTAGGCTTCTCCTCAGCG")

        if ind > 0 :
            if (ind > 3) and (ind < 76):
                # Identify barcodes
                barcode1 = line[ind-4:ind]

                # Map barcode to condition
                if barcode1 in barcode_dict.keys():
                    barcode1_arr.append(barcode_dict[barcode1])
                else:
                    barcode1_arr.append("None")
                
            else:
                barcode1_arr.append("None")
                barcode2_arr.append("None")
            
            ind_prom_bc = line.find("TTTTACATGACTGACTGA")
            if ind_prom_bc > 0:
                if ind_prom_bc < 80-18:
                    barcode_promoter_arr.append(line[ind_prom_bc+18:ind_prom_bc+38])
                else:
                    barcode_promoter_arr.append("None")
                
            else:
                barcode_promoter_arr.append("None")


pd.DataFrame(data={
        "Barcode1":barcode1_arr,
        "Barcode_promoter":barcode_promoter_arr
    }).to_csv(
    f"{homedir}/data/sequencing_data/processed_sequencing/LB_heatshock_extracted.txt",
    index=False,
    header=False
)
