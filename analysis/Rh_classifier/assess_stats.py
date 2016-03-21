'''
Rh_classifier-assess_stats.py

Collects all stats in dataframe and computes stats on those stats. 

Tim Farrell, tmf@bu.edu
20151205
''' 
import pickle
import pandas as pd
from copy import deepcopy
from subprocess import check_output

# collects all stats in results folder
# puts them in a pandas dataframe for analysis
def collect_stats(): 
    stat_files = [f for f in check_output(['ls ../../results/Rh_classifier'],\
                                            shell=True).split('\n')\
                    if '_stats.p' in f]
    feature_typesets = [f[:f.rfind('_')] for f in stat_files]
    all_stats = pd.DataFrame(columns=['classifier','pheno_enc','gene','select','res','stat',\
                                'antigen','accuracy','F-measure','recall','precision',\
                                'FP-rate','std'])
    #all_stats = [] 
    for ft in feature_typesets: 
        stats = pickle.load(open('../../results/Rh_classifier/' + ft + '_stats.p','rb'))
        gene, select, res, stat = ft.split('-')
        for cl in stats: 
            classifier, pheno_enc = cl.split('-')
            for a in stats[cl]: 
                d = deepcopy(stats[cl][a])
                d.update({'gene':gene,'select':select,'res':res, 'stat':stat,\
                          'classifier':classifier,'pheno_enc':pheno_enc,'antigen':a})
                all_stats = all_stats.append(pd.DataFrame(d, columns=list(d.keys()), index=[1]),\
                                             ignore_index=True) 
                #all_stats = all_stats.append(r)
    return all_stats


def main(): 
    # get all stats and save
    S = collect_stats()
    pickle.dump(S, open('../../results/Rh_classifier/All_Stats.p','wb'))
    
    # rank by composite phenotype call accuracy
    R = rank_by(S, 'accuracy', (antigen,'DCcEe'))

if __name__ == '__main__': 
    main() 