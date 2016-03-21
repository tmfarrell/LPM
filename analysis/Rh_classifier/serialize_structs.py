'''
Rh_classifier-serialize_structs.py

Serializes Rh and Sample data structs for Rh classifier project.

Tim Farrell, tmf@bu.edu
20160131
'''
import pickle 
import pandas as pd
import numpy as np

datadir = '/Volumes/Data/work/projects/lpm/data/'

## Build Rh data struct
# Of note, the provided gene positions are different from those cited by NCBI
# for assembly GRCH37.p13: RHCE: (25687853, 25747363), RHD: (25598884, 25656936)
Rh = {'RHCE': {'whole': {1: (25688739, 25757662)}},\
      'RHD':  {'whole': {1: (25588680, 25656935)}}}

# get exon data
rh_exons = pd.read_table(datadir + 'Rh_classifier/Rh_exons_grch37.txt',\
                    names=['Genomic Region','Source','Type','Start',\
                           'End','.','Strand','..','Notes'])
# preprocess it
rh_exons['Gene'] = rh_exons['Notes'].apply(\
        lambda n: n[n.find('GeneID:') + len('GeneID:'):\
                    n.find('GeneID:') + len('GeneID:') + 4])
rh_exons['Gene'] = rh_exons['Gene'].apply(lambda g: 'RHD' if g == '6007' else 'RHCE')
for col in ['.','..','Notes','Genomic Region','Source','Type','Strand']:
    del(rh_exons[col]) 
rh_exons = rh_exons.drop_duplicates(subset=['Start', 'Gene'])\
                   .drop_duplicates(subset=['End', 'Gene'])
rh_exons = rh_exons.sort(columns=['Start','End']).reset_index()

# index it and add to Rh data struct
# add the putative snv coordinates as well 
for gene in Rh:
    Rh[gene]['snps'] = list(map(int, list(pd.read_csv(\
                                datadir + 'Rh_classifier/'+ gene +'.csv').loc[0,:][2:])))
    if gene == 'RHCE':
        # since RHCE is backwards (i.e. minus strand) need to index exons differently
        rh_exons.loc[rh_exons['Gene'] == gene,'Exon #'] = \
                        pd.Series(range(10,0,-1), index=range(10,20)) 
        # and count from end in making snv coordinates absolute
        Rh[gene]['snps'] = {n: (Rh[gene]['whole'][1][1] - c, Rh[gene]['whole'][1][1] - c + 1)\
                                    for n, c in enumerate(Rh[gene]['snps'])}
    
    else:
        rh_exons.loc[rh_exons['Gene'] == gene,'Exon #'] = \
                        pd.Series(range(1,11), index=range(0,10))
        Rh[gene]['snps'] = {n: (Rh[gene]['whole'][1][0] + c, Rh[gene]['whole'][1][0] + c + 1)\
                                    for n, c in enumerate(Rh[gene]['snps'])}
            
    Rh[gene]['exome'] = {int(v[-1]): (v[1], v[2])\
                         for v in rh_exons[rh_exons['Gene'] == gene].values}



## Build Samples data struct
# process samples and their phenotype data
f = open(datadir + 'ngs_rhd_rhce/RHD_RHCE_Serology.csv','r')
serology = []
for line in f.readlines()[1:]:
    fields = line.strip().split(',')
    # ensure has readings for all antigens
    if fields[2:] != ['','','','']:
        # convert serological reading into integer. the mapping is as follows:
        #  n+ |-> n
        #  +  |-> 1
        #  -  |-> 0
        serology.append([fields[0]] + \
            list(map(lambda x: 0 if '-' in x else (int(x[0]) if x != '+' else 1),\
                    fields[1:])))
f.close()
# collect in dict 
Samples = {}
for s in serology:
    Samples[s[0]] = {'phenotypes': list(zip(['D','C','c','E','e'], s[1:]))}


## Save both
pickle.dump(Rh, open('Rh.p', 'wb'))
pickle.dump(Samples, open('Samples.p', 'wb'))