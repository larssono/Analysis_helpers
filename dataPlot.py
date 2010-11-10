import matplotlib
import numpy as np
import pylab

cdict = {'red': ((0.0, 1.0, 1.0),
                 (0.5, 0.0, 0.0),
                 (1.0, 0.0, 0.0)),
         'green': ((0.0, 1.0, 1.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, 0.0, 0.0),
                  (0.5, 0.0, 0.0),
                  (1.0, 1.0, 1.0))}
blueyellow = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    
def bar_w_error(y, errors=np.array([]), legend=[]):
    colors='bgrcy'
    y=y[...,np.newaxis] if y.ndim==1 else y
    errors=errors[...,np.newaxis] if errors.ndim==1 else errors
    nGroups,nBars = y.shape
    x = np.arange(nGroups)  # the x locations for the groups
    width = 0.8/nBars       # the width of the bars
    
    rects=[]
    for i in range(nBars):
        rects.append(pylab.bar(x+i*width, y[:,i], width, color=colors[i%len(colors)]))
        if len(errors)>0:
            pylab.errorbar(x+width/2+i*width, y[:,i], errors[:,i], fmt=None, ecolor='k')

    pylab.xticks(x+(width*nBars)/2, np.arange(nGroups)+1)

    ax=pylab.gca()
    if len(legend)==nBars:
        ax.legend( [rect[0] for rect in rects], legend, loc=0 )

def qqplot(y, x='norm'):
    """Displays a Quantile-Quantile plot of x versus theoretical
    quantiles (based on y or default normal distribution).

    @param y:   array of pvalues

    @param x: either an array of pvalues or one of the keywords
              {'norm', 'linear'} where 'norm' means y are normally
              distributed and 'linear' means y is linearly distributed
              default: 'norm'
    """
    nLen=len(y)
    if x=='norm':
        x = np.sort(np.random.randn(nLen))
    elif x=='linear':
        x = (np.arange(1,nLen+1)-.5)/nLen
    pylab.loglog(x, np.sort(y), '.', basex=10, basey=10)
    xmin=min(pylab.gca().get_xlim())
    
    pylab.setp(pylab.gca(), 'xlim', [1, xmin])
    pylab.setp(pylab.gca(), 'ylim', [1, min(pylab.gca().get_ylim())])
    #pylab.loglog([1,xmin], [1, xmin], 'k')
    pylab.xlabel('Expected P value')
    pylab.ylabel('Observed P value')
    
def matlabBoxplot(x, groups):
    """Displays box plots of multiple data samples.
    
    BOXPLOT(X,G) specifies one or more grouping variables G, producing a
    separate box for each set of X values sharing the same G value or
    values.  Grouping variables must have one row per element of X, or one
    row per column of X. Specify a single grouping variable in G by using a
    vector of strings; specify multiple grouping variables in G by using a
    array of these vectors, such as [G1 G2 G3], or by using a
    matrix.  If multiple grouping variables are used, they must all be the
    same length.  Groups that contain a NaN or an empty string ('') in a
    grouping variable are omitted, and are not counted in the number of
    groups considered by other parameters.
    
    Arguments:
    - `x`: vector of values
    - `groups`: 
    """
    #split x into the groups
    y=[]; legends=[]
    labels=[]
    groups=map(np.asarray, groups)
    for g in groups:
        labels.append(list(set(g)))
    labelLengths=map(len, labels)
    for i in range(np.prod(labelLengths)):
        idx=np.ones(len(x), np.bool_)
        legend=''
        for j, label in enumerate(labels):
            labelPos=int(i/np.prod(labelLengths[j+1:])%labelLengths[j])
            idx=idx&(groups[j]==label[labelPos])
            legend+='-'+str(label[labelPos])
            #print j, label[int(i/np.prod(labelLengths[j+1:])%labelLengths[j])],
        #print groups[0][idx], groups[1][idx]
        y.append(x[idx])
        legends.append(legend)
    legends=[s.strip('-') for s in  legends]

    pylab.boxplot(y)
    pylab.xticks(range(1,len(y)+1), legends)

