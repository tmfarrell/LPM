'''
Rh_classifier-feature_eng.py

Feature engineering development for Rh classifier project.

Tim Farrell, tmf@bu.edu
20160131
'''
import pylab
import scipy
import numpy
import matplotlib
matplotlib.use('Agg')               # To save plots from command-line
from sklearn import preprocessing
import scipy.cluster.hierarchy as sch

###############################
##         FUNCTIONS         ##
###############################
def plot_hierarchical_clustering(matrix, typeset, metric, file):
# plots hierarchical clustering of N x N matrix
    if metric == 'distance':
        n = len(matrix)
        D = numpy.zeros((n, n))
        for i, j in [(i, j) for i in range(n) for j in range(n)]:
            u, v = matrix[i], matrix[j]
            D[i][j] = scipy.spatial.distance.hamming(u, v)

    elif metric == 'correlation':
        D = numpy.corrcoef(matrix)

    # first dendrogram
    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = sch.linkage(D, method='centroid')
    T1 = sch.fcluster(Y, 5)
    Z1 = sch.dendrogram(Y, orientation='left')
    ax1.set_xticks([])
    ax1.set_yticks([])

    # second dendrogram
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='single')
    T2 = sch.fcluster(Y, 5)
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # distance matrix
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    file.write(','.join([typeset, 'centroid', str(len(set(idx1))) + '\n']))
    file.write(','.join([typeset, 'single', str(len(set(idx2))) + '\n']))
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # colorbar
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    fig.savefig('./plots/' + typeset + '_' + metric + '.png')


###############################
##           MAIN            ##
###############################
##
## Build feature sets from PFMs and assess performance of classifier
##
f = open('typeset_cluster_stats.csv', 'w')
f.write(','.join(['Typeset,ClusteringMetric,NumClusters\n']))

feature_typesets = [(metric, encoding, selection)\
                        for metric in ['max_pos_freq', 'mean_pos_freq']\
                        for encoding in ['nonencoding', 'encoding']\
                        for selection in ['exomic', 'diff_genotype']]
#'nonencoding', , 'var_pos_freq'

for n, (metric, encoding, selection) in enumerate(feature_typesets):
    ft = metric + '\n' + selection + '\n' + encoding
    typeset = ft.replace('\n', '-')
    print "Starting feature typeset "+ str(n) + ": " + typeset + "..."

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

    if encoding == 'encoding':
        enc = preprocessing.OneHotEncoder()
        enc.fit([Samples[sample]['features']['RHCE']+\
                 Samples[sample]['features']['RHD'] for sample in Samples])
        features = [enc.transform(\
                    Samples[sample]['features']['RHCE']+\
                    Samples[sample]['features']['RHD']).toarray().tolist()[0]\
                    for sample in Samples]
        #plot_hierarchical_clustering(features, typeset, 'distance', f)

    else:
        features = [Samples[sample]['features']['RHCE']+\
                    Samples[sample]['features']['RHD'] for sample in Samples]
        #plot_hierarchical_clustering(features, typeset, 'correlation', f)


phenotypes = [[ph[1] for ph in Samples[sample]['phenotype']]\
                     for sample in Samples]
#plot_hierarchical_clustering(phenotypes, 'phenotypes', 'distance', f)

phenotypes_bin = map(lambda p: map(lambda x: int(x > 0), p), phenotypes)
plot_hierarchical_clustering(phenotypes_bin, 'phenotypes_bin', 'distance', f)
f.close()
