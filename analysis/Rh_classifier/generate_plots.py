##
## Rh_classifier-generate_plots.py
##
## Tim Farrell, tmf@bu.edu
## 20160120
##
import pickle
import matplotlib
matplotlib.use('Agg')
from numpy import arange
import matplotlib.pyplot as plt
from multiprocessing import Pool
from subprocess import check_output

## Functions 
# labels values on bar plots
def autolabel(ax, rects, stds):
    # attach some text labels
    for rect, std in zip(rects, stds):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.85*height,
                '%.2f+/-%.2f' % (height, std),
                ha='center', va='bottom')

def quad_area(p1, p2, p3, p4): 
    (x1, y1), (x2, y2), (x3, y3), (x4, y4) = p1, p2, p3, p4 
    return 0.5 * abs(x1*y2 + x2*y3 + x3*y4 + x4*y1\
                     - x2*y1 - x3*y2 - x4*y3 - x1*y4)
    
def calc_auc(fp_rate, recall): 
    return quad_area((0,0), (fp_rate, recall), (1,1), (1,0))

def plot_roc(ft, cl, cl_stats): 
    aucs = {} 
    fig = plt.figure() 
    for a in cl_stats: 
        if len(a) < 2: 
            x, y = cl_stats[a]['FP-rate'], cl_stats[a]['recall']
            plt.fill([0, x, 1, 1], [0, y, 1, 0], 'bo-', alpha=0.1)
            plt.annotate(a, xy=(x, y), xycoords='data', xytext=(0, -15),\
                         textcoords='offset points')
            aucs[a] = calc_auc(x, y)
    plt.xlim((0.0, 1.0))
    plt.ylim((0.0, 1.0))
    plt.ylabel('TP rate')
    plt.xlabel('FP rate')
    plt.title("ROC for " + ft + '-' + cl)
    plt.annotate('\n'.join(['AUCs:'] + [k + '  ' + str(v) for k, v in aucs.items()] +\
                           ['', 'Accuracy:', 'DCcEe  ' + str(cl_stats['DCcEe']['accuracy'])]),\
                 xy=(0.6, 0.3), xycoords='data', xytext=(0, -15), textcoords='offset points')
    plt.savefig('../../results/Rh_classifier/plots/' + ft + '-' + cl + '_ROC.png')
    
def plot_accuracy(ft, cl, cl_stats): 
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    antigens = list(cl_stats.keys())
    idx = arange(len(antigens)); width = 0.4;
    accs = [cl_stats[a]['accuracy'] for a in antigens]
    #stds = [cl_stats[a]['std'] for a in antigens]
    bars = ax.bar(idx + width, accs, width, color='0.65') #, yerr=stds)
    ax.set_xticks(idx + 1.5*width)
    ax.set_xticklabels(antigens)
    autolabel(ax, bars, [0]*len(bars)) #stds)
    plt.ylim((0.0, 1.0))
    plt.xlabel("Antigen")
    plt.title("Accuracy for " + ft + '-' + cl)
    plt.savefig('../../results/Rh_classifier/plots/' + ft + '-' + cl + '_ACC.png')
    
def plot(ft):
    ft_stats = pickle.load(open('../../results/Rh_classifier/' + ft + '_stats.p','rb'))
    print("Plotting for " + ft + "...")
    for cl in ft_stats: 
        if 'binarize' in cl: 
            plot_roc(ft, cl, ft_stats[cl])
            plot_accuracy(ft, cl, ft_stats[cl])
        else: 
            plot_accuracy(ft, cl, ft_stats[cl])
    print("Saved all for " + ft + ".")
            

## Main 
def main(): 
    feature_files = [f for f in check_output(['ls ../../results/Rh_classifier'],\
                                            shell=True).split('\n')\
                    if '_stats.p' in f]
    feature_typesets = [f[:f.rfind('_')] for f in feature_files]
    # load and plot, parallelizing
    pool = Pool()
    pool.map(plot, feature_typesets)
    
if __name__ == '__main__': 
    main() 