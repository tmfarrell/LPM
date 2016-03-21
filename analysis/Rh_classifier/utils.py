# 
# Rh_classifier_utils.py 
# 
# Tim Farrell, tmf@bu.edu 
# 20160308
# 
## 
## Functions 
## 
# converts coordinate positions to samtools region string
def pos2chr1region(pos): 
    start, end = pos
    return ':'.join(['1', str(start), str(end)]) 

# converts samtools region string to coordinate positions
def chr1region2pos(region): 
    chr_, start, end = region.split(':')
    return (start, end)

# computes the crossproduct of 2 sets|lists
def product(X, Y): 
    return [(x, y) for x in X for y in Y]

# compute the crossproducts of 3 sets|lists
def product3(X, Y, Z): 
    return [(x,y,z) for x in X for y in Y for z in Z]

# converts a sample phenotype of the form 
# [('D',0|1),('C',0|1),('c',0|1),('E',0|1),('e',0|1)]
# to binary identifier str
def phenotypes2bin(phenotypes): 
    return ''.join([str(int(int(p[1]) > 0)) for p in phenotypes])

def ft2str(ft): 
    ((genes, withingene), (quant, stat_func)) = ft
    g = genes[0] if (len(genes) == 1) else (genes[0]+ '-' +genes[1])
    return '-'.join([g, withingene, quant, str(stat_func)])

def get_alignment(s):
    from sys import modules
    assert ('pysam' in modules.keys()), "pysam module must be loaded"
    return pysam.AlignmentFile('/Volumes/Data/work/projects/lpm/data/ngs_rhd_rhce/'\
                               + s + '-extract_for_RHD_RHCE.bam', 'rb')

'''
# plots hierarchical clustering of N x N matrix
def plot_hierarchical_clustering(matrix, typeset, metric, file_):
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
    file_.write(','.join([typeset, 'centroid', str(len(set(idx1))) + '\n']))
    file_.write(','.join([typeset, 'single', str(len(set(idx2))) + '\n']))
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # colorbar
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    fig.savefig('./plots/' + typeset + '_' + metric + '.png')
'''