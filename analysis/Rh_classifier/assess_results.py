'''
Rh_classifier-assess_results.py

Assesses the performance of feature selection and classifier techniques for
Rh classifier development project.

Tim Farrell, tmf@bu.edu
20151205
'''
import pickle
import numpy as np
import pandas as pd
from copy import deepcopy
from multiprocessing import Pool
from subprocess import check_output


## Classes
# simple class for result
class Result: 
    def __init__(self, r): 
        self.binary = int(r[0] == r[1])
        self.type_ = ('T' if self.binary else 'F') +\
                     ('P' if r[1] else 'N')
    def __repr__(self): 
        return self.type_
    def getback(self): 
        return {'TP': (1, 1), 'TN': (0, 0),\
                'FP': (0, 1), 'FN': (1, 0)}[self.type_]

    
## Functions 
def safe_divide(a, b): 
    try: 
        return (float(a)/float(b))
    except ZeroDivisionError: 
        return float('Inf')

def indicate(results, binarized=True): 
    l = list(map(lambda (antigen, result):\
                    (antigen, (Result(result) if binarized else int(result[0]==result[1]))),\
                 results))
    l = l + [('DCcEe', int(sum([(r[1].binary if binarized else r[1]) for r in l]) == 5))] 
    return l

def indicate_all(results, binarized=True): 
    return list(map(lambda result: indicate(result, binarized), results))

def indicate_results(results): 
    r_ind = {} 
    for cl in results: 
        r_ind[cl] = indicate_all(results[cl], binarized=('binarize' in cl))
    return r_ind

def dictionize_all(results): 
    antigens = ['D','C','c','E','e']
    antigens = antigens + [''.join(antigens)]
    d = {a: [] for a in antigens}
    for res in results: 
        for a, r in res: 
            d[a].append(r)
    return d

def dictionize_results(results): 
    dictionized = {}
    for cl in results: 
        dictionized[cl] = dictionize_all(results[cl])
    return dictionized

def metricize(cts):
    if type(cts) == dict: 
        recall = safe_divide(cts['TP'], cts['P'])
        precision = safe_divide(cts['TP'], (cts['TP'] + cts['FP']))
        return {'FP-rate': safe_divide(cts['FP'], cts['N']),\
                'precision': precision, 'recall': recall,\
                'accuracy': safe_divide((cts['TP'] + cts['TN']), (cts['P'] + cts['N'])),\
                'F-measure':\
                safe_divide(2, (safe_divide(1, precision) + safe_divide(1, recall)))}
    return {'accuracy': safe_divide(cts, 93)} 
    
def calc_metrics(counts): 
    return {a: metricize(ct) for a, ct in counts.items()}    

def staticise(dictionized, binarized=True):
    if not binarized: 
        return {a: {'accuracy': np.mean(v), 'std': np.std(v)}\
                for a, v in dictionized.items()}
    else: 
        bin_stats = {s: 0 for s in ['TP','TN','FP','FN']}
        cts = {a: (deepcopy(bin_stats) if len(a) < 2 else 0) for a in dictionized} 
        for a, res in dictionized.items(): 
            for r in res: 
                if len(a) < 2:  cts[a][r.type_] += 1
                else:           cts[a] += r
         
        for a, v in cts.items(): 
            if type(v) == dict: 
                cts[a]['P'] = v['TP'] + v['FN']
                cts[a]['N'] = v['FP'] + v['TN']
        
        return calc_metrics(cts)
            
def calc_stats(results):
    dictionized = dictionize_results(indicate_results(results))
    stats = {}
    for cl in dictionized:
        stats[cl] = staticise(dictionized[cl], binarized=('binarize' in cl))
    return stats

# get results, staticize and pickle
def get_stat_pickle(ft):  
    print("Getting results for " + ft + "...")
    results = pickle.load(open('../../results/Rh_classifier/' + ft + '_results.p','rb'))
    print("Calculating stats for " + ft + "...")
    stats = calc_stats(results)
    pickle.dump(stats, open('../../results/Rh_classifier/' + ft + '_stats.p','wb'))
    print("Saved stats for " + ft + ".")

## Main 
def main():  
    result_files = [f for f in check_output(['ls ../../results/Rh_classifier'], shell=True).split('\n')\
                    if '_results.p' in f]
    feature_typesets = [f[:f.rfind('_')] for f in result_files]
    pool = Pool()
    pool.map(get_stat_pickle, feature_typesets)

if __name__ == '__main__': 
    main() 