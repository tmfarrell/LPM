'''
Rh_classifier-assess_perf.py

Assesses the performance of feature selection and classifier techniques for
Rh classifier development project.

Tim Farrell, tmf@bu.edu
20151205
'''
import matplotlib
from sklearn import svm
from sklearn import preprocessing
from sklearn import cross_validation
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier

# To save plots from command-line
matplotlib.use('Agg')
import matplotlib.pyplot as plt


###############################
##         FUNCTIONS         ##
###############################
# for labeling values on bar plots
def autolabel(rects, stds):
    # attach some text labels
    for rect, std in zip(rects, stds):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.85*height,
                '%.2f+/-%.2f' % (height, std),
                ha='center', va='bottom')

def plot_features(sample_feature_set):
    return 


###############################
##           MAIN            ##
###############################
##
##  Get user input for:
##      (1) whether to classify antigen phenotypes compositionally ('composite')
##          or individually ('individual')
##      (2) numbers of iterations
##
classify_by = sys.argv[1]
assert (classify_by == 'composite' or classify_by == 'individual'),\
        "classify_by arg must be either composite or individual"

iterations = range(int(sys.argv[2]))


##
## Build feature sets from PFMs and assess performance of classifier
##
# get feature sets <- seq of integers representing bases according to bases_to_idx
### can parallelize all of this
antigens = ['D', 'C', 'c', 'E', 'e']
composite_phenotyping_success, antigen_phenotyping = {}, {}
feature_typesets = [(metric, encoding, selection)\
                        for metric in ['max_pos_freq', 'mean_pos_freq']\
                        for encoding in ['encoding']\
                        for selection in ['exomic', 'diff_genotype']]
#'nonencoding', , 'var_pos_freq'

for n, (metric, encoding, selection) in enumerate(feature_typesets):
    ft = metric + '\n' + selection + '\n' + encoding
    #composite_phenotyping_success[ft], antigen_phenotyping[ft] = 10 * [0], 10 * [0]
    antigen_phenotyping[ft] = {}
    print "Starting feature typeset "+ str(n) + ": " + ft.replace('\n', ' ') + "..."

    # feature extraction
    for sample in Samples:
        Samples[sample]['features'] = {}
        for gene in Rh:
            for exon in Rh[gene]['exons']:
                #if metric == 'base':
                #    Samples[sample]['features'][gene] = \
                #            [np.argmax(Samples[sample]['PFM'][gene][:,pos]) \
                #                                for pos in Rh[gene][selection]]
                if metric == 'max_pos_freq':
                    Samples[sample]['features'][gene] = \
                            [max(Samples[sample]['PFM'][gene][:,pos]) \
                                                for pos in Rh[gene][selection]]
                if metric == 'mean_pos_freq':
                    Samples[sample]['features'][gene] = \
                            [np.mean(Samples[sample]['PFM'][gene][:,pos]) \
                                                for pos in Rh[gene][selection]]
                if metric == 'var_pos_freq':
                    Samples[sample]['features'][gene] = \
                            [(np.std(Samples[sample]['PFM'][gene][:,pos]))**2 \
                                                for pos in Rh[gene][selection]]

    # composite antigen phenotype classification via DecisionTreeClassifier
    if classify_by == 'composite':
        if encoding == 'encoding':
            # encode
            enc = preprocessing.OneHotEncoder()
            enc.fit([Samples[sample]['features']['RHCE']+\
                     Samples[sample]['features']['RHD']\
                        for sample in Samples])
            features = [enc.transform(\
                        Samples[sample]['features']['RHCE']+\
                        Samples[sample]['features']['RHD']).toarray().tolist()[0]\
                        for sample in Samples]
        else:
            features = [Samples[sample]['features']['RHCE']+\
                        Samples[sample]['features']['RHD'] for sample in Samples]

        phenotypes = [[ph[1] for ph in Samples[sample]['phenotype']]\
                                          for sample in Samples]
        classifier = DecisionTreeClassifier()

        kfold_cv = cross_validation.KFold(len(features), n_folds=10)
        for iteration in iterations:
            print "\t...iteration: " + str(iteration) + "..."
            test_v_predict_sets  = []
            for traincv, testcv in kfold_cv:
                model = classifier.fit(\
                            [features[t] for t in traincv],\
                            [phenotypes[t] for t in traincv])
                test_features   = [features[t] for t in testcv]
                test_phenotypes = [phenotypes[t] for t in testcv]

                test_v_predict_sets.append(zip(test_phenotypes,\
                                map(lambda l: map(int, l),\
                                model.predict(test_features).tolist())))

    # individual antigen phenotype classification via SVM
    elif classify_by == 'individual':
        features, phenotypes = {}, {}
        for antigen in antigens:
            '''if encoding == 'encoding':
                # encode
                enc = preprocessing.OneHotEncoder()
                enc.fit([Samples[sample]['features']['RHCE']+\
                         Samples[sample]['features']['RHD']\
                            for sample in Samples])
                features = [enc.transform(\
                            Samples[sample]['features']['RHCE']+\
                            Samples[sample]['features']['RHD']).toarray().tolist()[0]\
                            for sample in Samples]
            else:
                features = [Samples[sample]['features']['RHCE']+\
                            Samples[sample]['features']['RHD'] for sample in Samples]
            '''#print "\t...antigen: " + antigen + "..."
            gene2feature = 'RHD' if antigen == 'D' else 'RHCE'
            # encode
            if encoding == 'encoding':
                enc = preprocessing.OneHotEncoder()
                enc.fit([Samples[sample]['features'][gene2feature]\
                        for sample in Samples])
                features[antigen] = [enc.transform(\
                    Samples[sample]['features'][gene2feature]).toarray().tolist()[0]\
                    for sample in Samples]
            else:
                features[antigen] = [Samples[sample]['features'][gene2feature]\
                                    for sample in Samples]

            phenotypes[antigen] = [ph[1] for sample in Samples\
                                         for ph in Samples[sample]['phenotype']\
                                         if ph[0] == antigen]

        classifier = DecisionTreeClassifier()   #svm.SVC()
        test_v_predict_sets  = {a: [] for a in antigens}
        for antigen in antigens:
            kfold_cv = cross_validation.KFold(len(features), n_folds=5)
            for traincv, testcv in kfold_cv:
                model = classifier.fit([features[antigen][t] for t in traincv], \
                                        [phenotypes[antigen][t] for t in traincv])
                test_features   = [features[antigen][t] for t in testcv]
                test_phenotypes = [phenotypes[antigen][t] for t in testcv]

                test_v_predict_sets[antigen].append(zip(test_phenotypes,\
                                            model.predict(test_features).tolist()))

        # transform individual test_v_predict_sets to composite
        t = map(lambda x: zip(x[0],x[1],x[2],x[3],x[4]),\
            zip(test_v_predict_sets['D'], test_v_predict_sets['C'],\
                test_v_predict_sets['c'], test_v_predict_sets['E'],\
                test_v_predict_sets['e']))

        print t

        comp_test_v_predict_sets = [map(lambda x:\
            ([x[0][0], x[1][0], x[2][0], x[3][0], x[4][0]],\
             [x[0][1], x[1][1], x[2][1], x[3][1], x[4][1]]), s) for s in t]

        antigen_phenotyping[ft] = map(np.mean, [[int(test == predict) for test, predict in s]\
                                    for s in comp_test_v_predict_sets])

        '''antigen_phenotyping[ft][iteration] = [[int(test[i] == predict[i]) for i in range(5)]\
                                    for test_v_predict_set in test_v_predict_sets\
                                    for test, predict in test_v_predict_set]
        # can parallelize
        composite_phenotyping_success[ft][iteration] = \
                            (np.mean(map(all, antigen_phenotyping[ft][iteration])),\
                             np.std(map(all, antigen_phenotyping[ft][iteration])))


        antigen_phenotyping[ft][antigen] = [0] * len(iterations)
        kfold_cv = cross_validation.KFold(len(features), n_folds=10)
        for iteration in iterations:
            print "\t...iteration: " + str(iteration) + "..."
            test_v_predict_sets  = []
            for traincv, testcv in kfold_cv:



            antigen_phenotyping[ft][antigen][iteration] = \
                        [int(test == predict) \
                        for test_v_predict_set in test_v_predict_sets\
                        for test, predict in test_v_predict_set]

                ## to do: - compute composite phenotyping data struct
                ##        from individual predictions;
                ##        - tricky, may have to restructure'''

'''
print "done."
print "Calculating antigen-specific accuracy..."
## calc antigen-specific phenotyping success and plot
antigen_phenotyping_success= {}
for typeset in antigen_phenotyping:
    antigen_phenotyping_success[typeset] = {}
    for i, antigen in enumerate(antigens):
        antigen_phenotyping_success[typeset][antigen] = [0] * 10
        for iteration in range(10):
            by_antigen = [l[i] for l in antigen_phenotyping[typeset][iteration]]
            antigen_phenotyping_success[typeset][antigen][iteration] = \
                                    (np.mean(by_antigen), np.std(by_antigen))
print "done."

## plot composite phenotyping result
for typeset in ['base', 'max_pos_freq', 'mean_pos_freq', 'var_pos_freq']:
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    phenotyping_succ = [ft for ft in composite_phenotyping_success if typeset in ft]
    phenotyping_succ.sort()
    idx = np.arange(len(phenotyping_succ)); width = 0.4;
    bars = ax.bar(idx + width,\
                map(np.mean,[composite_phenotyping_success[ft] for ft in phenotyping_succ]),\
                width, color='0.65',\
                yerr=map(np.std,[composite_phenotyping_success[ft] for ft in phenotyping_succ]))
    ax.set_xticks(idx + 1.5*width)
    ax.set_xticklabels(phenotyping_succ)
    autolabel(bars)
    plt.ylim((0.0,1.0))
    plt.ylabel("Mean Composite Phenotyping Success Rate")
    plt.title(typeset + " feature-typesets (10-fold cross-validation)")
    plt.savefig(typeset + '_composite_phenotyping_perf.png')

print "Generating figures..."
for typeset in antigen_phenotyping:
    typeset_str = typeset.replace('\n', ' ')
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    antigen_stats = [antigen_phenotyping_success[typeset][a] for a in antigens]
    antigen_means = map(lambda xs: [mean for mean,std in xs], antigen_stats)
    composite_means = [mean for mean,std in composite_phenotyping_success[typeset]]
    idx = np.arange(len(antigens) + 1); width = 0.4;
    stds = map(np.std, antigen_means) + [np.std(composite_means)]
    bars = ax.bar(idx + width, map(np.mean, antigen_means) +\
                [np.mean(composite_means)], width, color='0.65', yerr=stds)
    ax.set_xticks(idx + 1.5*width)
    ax.set_xticklabels(antigens + ['composite'])
    autolabel(bars, stds)
    plt.ylim((0.0,1.0))
    plt.ylabel("Mean Phenotyping Success Rate")
    plt.title(typeset_str)
    plt.savefig('phenotyping_perf-' + typeset_str +'.png')
print "done."'''
