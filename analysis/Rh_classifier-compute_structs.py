'''
Rh_classifier-compute_structs.py

Computes Rh and Sample data structs for  Rh classifier development project.

Tim Farrell, tmf@bu.edu
20160131
'''
import pysam
import pandas as pd
import numpy as np

##
## Build Rh data struct
##
# of note, the provided gene positions are different from those cited by NCBI
# for assembly GRCH37.p13: RHCE: (25687853, 25747363), RHD: (25598884, 25656936)
Rh = {'RHCE': {'pos': (25688739, 25757662)}, \
      'RHD':  {'pos': (25588680, 25656935)}}

# get exon data
rh_exons = pd.read_table('./variant_data/Rh_exons_grch37.txt',\
                    names=['Genomic Region','Source','Type','Start',\
                           'End','.','Strand','..','Notes'])
# preprocess
rh_exons['Gene'] = rh_exons['Notes'].apply(\
        lambda n: n[n.find('GeneID:') + len('GeneID:'):\
                    n.find('GeneID:') + len('GeneID:') + 4])

rh_exons['Gene'] = rh_exons['Gene'].apply(lambda g:\
                                            'RHD' if g == '6007' else 'RHCE')

for col in ['.','..','Notes','Genomic Region','Source','Type','Strand']:
    del(rh_exons[col]);

rh_exons = rh_exons.drop_duplicates(subset=['Start', 'Gene'])\
                   .drop_duplicates(subset=['End', 'Gene'])

rh_exons = rh_exons.sort(columns=['Start','End']).reset_index()

# index exons and add to Rh data struct
for gene in Rh:
    if gene == 'RHCE':
        rh_exons.loc[rh_exons['Gene'] == gene,'Exon #'] = \
                pd.Series(range(10,0,-1), index=range(10,20)) # since minus strand
    else:
        rh_exons.loc[rh_exons['Gene'] == gene,'Exon #'] = \
                pd.Series(range(1,11), index=range(0,10))

    Rh[gene]['exons'] = {str(int(v[-1])): (v[1],v[2]) for v in \
                            rh_exons[rh_exons['Gene'] == gene].values}


##
## Build Samples data struct
##
# process serology data
f = open('/Volumes/Data/LPM/data/RHD_RHCE_Serology.csv','r')
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
            map(lambda x: 0 if '-' in x else (int(x[0]) if x != '+' else 1),\
                    fields[1:]))
f.close()

# add to samples struct
Samples = {}
for s in serology:
    Samples[s[0]] = {'phenotype': zip(['D','C','c','E','e'], s[1:])}

## Compute PFMs for each sample for each gene's exons
bases_to_idx = {'A':0, 'T':1, 'C':2, 'G':3, 'N':4, 'I':5}
idx_to_bases = {v: k for k, v in bases_to_idx.items()}
for sample in Samples:
    alignment = pysam.AlignmentFile('/Volumes/Data/LPM/data/' + sample +\
                                    '-extract_for_RHD_RHCE.bam', 'rb')
    Samples[sample]['PFM'] = {}
    #Samples[sample]['FASTA'] = {}

    for gene in Rh:
        Samples[sample]['PFM'][gene] = np.array([[]])

        for exon in Rh[gene]['exons']:
            #fasta_file = gene + '_exon' + exon + '.fasta'
            #fasta = open(fasta_file, 'w')

            exon_start = Rh[gene]['exons'][exon][0]
            exon_end   = Rh[gene]['exons'][exon][1]

            pfm = np.zeros((6, exon_end - exon_start + 1))
            for n, read in enumerate(alignment.fetch(reference='1',\
                                        start=exon_start, end=exon_end)):
                # (a) build exon specific PFM
                for aP in read.get_aligned_pairs():
                    ref_coor = aP[1]
                    if ref_coor != None \
                    and ref_coor < exon_end and ref_coor > exon_start:
                        base = read.query_sequence[aP[0]] if aP[0] != None else 'I'
                        pfm[bases_to_idx[base], ref_coor - exon_start] += 1

                # (b) save fasta to collection of all fasta for that exon
                #fasta.write('>' + gene + '_exon' + exon + '_read' + str(n)\
                #                + '\n' + read.query_sequence)

            # add PFM to Sample data struct
            if not Samples[sample]['PFM'][gene].any():
                Samples[sample]['PFM'][gene] = pfm
            else:
                Samples[sample]['PFM'][gene] = np.concatenate((\
                        Samples[sample]['PFM'][gene], pfm), axis=1)

            # add fasta file to Sample
            #fasta.close()
            #Samples[sample]['FASTA'][exon] = fasta_file
            # (c) call motif caller on that collection
            #subprocess.
            # (d) store motifs in Samples[sample]['motifs'][gene][exon]


##
## Compute position selection parameters
##
for gene in Rh:
    # those that supposedly differentiate genotype
    Rh[gene]['diff_genotype'] = list(pd.read_csv('./variant_data/' \
                                                + gene + '.csv').loc[0,:][2:])
    # whole exome
    Rh[gene]['exomic'] = range(sum([e[1]-e[0] for e in Rh[gene]['exons'].values()]))


##
## Find motifs (1) common and (2) differiential across phenotypes
##
## NOTE: if I do decide to go through with using bioinformatics tools to build
##       features, this will be done in its own stand-alone script
