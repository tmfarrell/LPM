'''
Rh_classifier-train_test.py

Tim Farrell, tmf@bu.edu
20151205
'''
import pickle
from sklearn.svm import SVC
from multiprocessing import Pool
from subprocess import check_output
from sklearn.cross_validation import LeavePOut
from sklearn.preprocessing import OneHotEncoder
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier

# import local utils
exec(open('Rh_classifier-utils.py'))

## Classes
# barebones class for encoders 
class Encoder:
    def __init__(self, name, encode_func): 
        self.name = name
        self.encode = encode_func
        
## Globals
Labels = pickle.load(open('../../data/Rh_classifier/Labels.p','rb'))

# init phenotype encoding function types 
Phenotypes = ['D','C','c','E','e']
Phenotype_encoders = [Encoder('identity', lambda x: x),\
                      Encoder('binarize', lambda x: int(x > 0))] 

# init statistical models 
Classifiers = [LogisticRegression(penalty='l1'), SVC(kernel='linear'),\
               DecisionTreeClassifier()]

## Functions 
def fold_cl_results(cl_results): 
    # make dict values into lists of key, item pairs for item in value
    for k, v in cl_results.items(): 
        cl_results[k] = [(k,p) for p in v]
    # get those values 
    clr = [cl_results[k] for k in cl_results.keys()]
    # transpose
    folded = []
    for j in range(len(clr[0])):
        l = []
        for i in range(len(clr)):
            l.append(clr[i][j])
        folded.append(l)
    return folded 

def train_test_pickle(n_ft):
    global Classifiers, Phenotypes, Phenotype_encoders
    n, ft = n_ft
    if 'base' in ft and 'argmax' in ft: 
        encode = True
    else: 
        encode = False
    Results = {} 
    features = pickle.load(open('../../data/Rh_classifier/' + ft + '_features.p','rb'))
    if encode: 
        enc = OneHotEncoder()
        enc.fit(features.as_matrix())
    samples = features.columns
    pout_cv = LeavePOut(len(samples), p=2)
    print("Feature typeset " + str(n) + ": " + ft + "...")
    for cl, pe in product(Classifiers, Phenotype_encoders):
        cl_pe_str = str(cl)[:str(cl).find('(')] + '-' + pe.name
        print("Starting " + cl_pe_str + "...")
        if ('base' in ft) and ('Gradient' in cl_pe_str): 
            continue 
        
        results = {} 
        for p in Phenotypes: 
            train_test_label_pairs = []
            for train_idx, test_idx in pout_cv: 
                # train_model
                train_samples = samples[train_idx]
                train_labels = [pe.encode(Labels.loc[p, s]) for s in train_samples]
                if encode: 
                    train_features = [enc.transform(features.loc[:, s].as_matrix()).toarray()\
                                      for s in train_samples]
                else: 
                    train_features = [features.loc[:, s].as_matrix() for s in train_samples]
                model = cl.fit(train_features, train_labels)

                # test model
                test_sample = samples[test_idx]
                test_label = [pe.encode(Labels.loc[p, s]) for s in test_sample]
                if encode: 
                    test_feature = [enc.transform(features.loc[:, s].as_matrix()).toarray()\
                                    for s in test_sample]
                else: 
                    test_feature = [features.loc[:, s].as_matrix() for s in test_sample]
                train_test_label_pairs = train_test_label_pairs +\
                                        zip(test_label, model.predict(test_feature))
            results[p] = train_test_label_pairs
        print("Finished results for " + cl_pe_str + ".")
        Results[cl_pe_str] = fold_cl_results(results)

    # save results for each feature_typeset in binary 
    pickle.dump(Results, open('../../data/Rh_classifier/' + ft + '_results.p','wb'))
    print("Saved " + ft + "_results.p.")


## Main 
def main(): 
    # load feature_typesets and Labels
    feature_files = [f for f in check_output(['ls ../../data/Rh_classifier'],\
                                            shell=True).split('\n')\
                    if ('_features.p' in f and 'read' in f)]
    feature_typesets = [f[:f.rfind('_')] for f in feature_files]
    
    #for n_ft in enumerate(feature_typesets): 
     #   train_test_pickle(n_ft)
    
    # add combinational features typesets
    #for ft0, ft1 in product(list(Features.keys()), list(Features.keys())): 
     #   Features[ft0 + '-' + ft1] = pd.concat([Features[ft0], Features[ft1]])

    # train test pickle, parallelizing
    pool = Pool()
    pool.map(train_test_pickle, enumerate(feature_typesets))

if __name__ == '__main__': 
    main() 
