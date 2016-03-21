'''
Rh_classifier-feature_eng.py

Feature engineering development for Rh classifier project.

Tim Farrell, tmf@bu.edu
20160131
'''
import pysam
import pickle
import pandas as pd
import numpy as np
from multiprocessing import Pool 
# import local utility functions 
exec(open('Rh_classifier-utils.py').read())

## Globals 
datadir = '../../data/Rh_classifier/'
# data structs
Rh = pickle.load(open(datadir + 'Rh.p','rb'))
Samples = pickle.load(open(datadir + 'Samples.p','rb'))
# to fill the Samples struct
def get_sample_alignments(samples): 
    for s in samples: 
        samples[s]['alignment'] = get_alignment(s)
    return samples
Samples = get_sample_alignments(Samples)

## Featurization functions
# abstract function for quantifying an alignment 
# alignment :: pysam.AlignmentFile
# pos :: (int, int)
# quantifier :: 'base' | 'read'
def quantify(alignment, pos, quantifier, reduce_stat=None): 
    if quantifier == 'base': 
        assert reduce_stat, "Must have reductive statistic with base quantifier."
        M = np.asarray(alignment.count_coverage(reference='1', start=pos[0], end=pos[1]),\
                       dtype=np.dtype(np.int32))
        m = np.zeros(M.shape[1])
        for col in range(M.shape[1]): 
            m[col] = reduce_stat(M[:, col])
        return m
    
    if quantifier == 'read': 
        return alignment.count(reference='1', start=pos[0], end=pos[1])
        
def quantify_at_pos_for_samples(pos, index, quantifier, reduce_stat=None, normalize=True):
    global Samples
    df = pd.DataFrame(index=index, columns=list(Samples.keys()))
    for s in list(Samples.keys()): 
        df.loc[index, s] = quantify(Samples[s]['alignment'], pos, quantifier,\
                                          reduce_stat=reduce_stat)
    if normalize: 
        df = df.apply(lambda row: row.apply(lambda x:\
                    float(float(x)-float(row.mean()))/float(row.std())), axis=1)
    return df

# extracts a featureset_type from all samples
# returns a dataframe of features by sample
def extract_featureset(featureset_type, key):
    global Rh, Samples
    
    (selection, metric) = featureset_type
    (genes, withingene) = selection
    (quantifier, reduce_stat) = metric  
    #typeset_str = '-'.join(list(select) + [metric[0], encoding])
    # build positions and index
    positions = [] 
    for g in genes: 
        for n, pos in Rh[g][withingene].items(): 
            positions = positions + [(n, pos)]
    # build feature dataframe 
    features = {} 
    normalize = False if (reduce_stat == np.argmax) else True
    for n, p in positions:
        index = [n]
        if quantifier == 'base': 
            index = list(range(p[0], p[1]))
        features[n] = quantify_at_pos_for_samples(p, index, quantifier,\
                                            reduce_stat=reduce_stat, normalize=normalize)
    df = pd.concat(features, ignore_index=True)
    # add index level for when combining feature typesets 
    df['key'] = key
    df.set_index('key', append=True, inplace=True)
    return df

def extract_pickle(n_ft):
    n, ft = n_ft
    ftstr = ft2str(ft)
    print("Extracting " + ftstr + "...")
    features = extract_featureset(ft, n)
    pickle.dump(features, open('../../data/Rh_classifier/' + ftstr + '_features.p','wb'))
    print("Saved " + ftstr + "_features.p.")
    

## Main
def main(): 
    ## Feature typesets
    Selections = product([['RHD'], ['RHCE']], ['exome', 'snps', 'whole']) 
    Metrics = [('base', np.argmax)] #max, lambda x: np.std(x)**2, np.mean])
    Feature_Typesets = enumerate(product(Selections, Metrics))
    
    # extract and pickle, parallelizing
    #pool = Pool()
    #pool.map(extract_pickle, Feature_Typesets)
    for n_ft in Feature_Typesets: 
        extract_pickle(n_ft)
    # make and pickle Labels dataframe
    #Labels = pd.DataFrame(columns=list(Samples.keys()))
    #for s in Samples: 
    #   for p in Samples[s]['phenotypes']: 
    #      Labels.loc[p[0], s] = p[1]

    # same both to binary
    # pickle.dump(Labels, open('../../data/Labels.p','wb'))

if __name__ == '__main__': 
    main() 