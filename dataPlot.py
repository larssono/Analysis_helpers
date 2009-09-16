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

    pylab.xticks(x+(width*nBars)/2, np.arange(nBars)+1)

    ax=pylab.gca()
    if len(legend)==nBars:
        ax.legend( [rect[0] for rect in rects], legend, loc=0 )

