from numpy.oldnumeric import *
from scipy.linalg import svd, pinv
import scipy.stats as stats

import numpy as np
import pylab
import numbers

import tensor
import dataPlot

#--------------------------------------------------------------            
# Methods for Filtering out missing values
#--------------------------------------------------------------
def filterNaN(data, cutoff):
    """Keeps rows that have a percentage of non missing values
    greater than or equal to cuttoff.

    @param cutoff: percentage of experiments that have to be good
       to keep.
    """
    nGenes, nExps = data.shape
    indexes = nonzero(sum(~np.isnan(data), 1)/float(nExps) >= cutoff)
    data=take(data, indexes)

def filterNaNExp(self, cutoff):
    """Keeps Columns that have a percentage of non missing values
    less than or equal to cuttoff

    @param cutoff: percentage of experiments that have to be good
       to keep.
    """
    nGenes, nExps = data.shape
    indexes = nonzero(sum(~np.isnan(data), 0)/float(nGenes) >= cutoff)
    data=take(data, indexes, 1)



#--------------------------------------------------------------            
#  Interpolate missing values
#--------------------------------------------------------------
def replaceNaNMean(data):
    """Replaces missing values in the matrix with the average
    value across all column for each row.

    Example:
        [1, 2, 1, 3]              [1, 2, 1,    3]
        [2, 2, nan, 3]    ==>     [2, 2, 2.33, 3]
        [1, 4, 0, 2]              [1, 4, 0,    2]
    """
    nGenes, nExps = data.shape
    for i in range(nGenes):
        mdata = np.ma.masked_array(data[i,:],np.isnan(data[i,:]))
        mm = np.mean(mdata)
        data[i,:][mdata.mask]=mm


def replaceNaNInterp(data):
    """Replaces missing values in the matrix by linearly
    interpolating between existing values.

    Example:
        [1, 2, 1,   3]            [1, 2, 1,   3]  
        [2, 2, nan, 3]   ==>      [2, 2, 2.5, 3]
        [1, 4, 0,   2]            [1, 4, 0,   2]  
    """
    nGenes, MAX = data.shape
    for i in range(nGenes):
        pylab.plot(data[i,:], '-o');
        indxs = np.nonzero(np.isnan(data[i, :]))[0]
        x1=0
        print indxs
        if len(indxs) > 0:
            for index in indxs:
                if index < x1: continue
                x0=index-1
                x1=index+1

                while np.isnan(data[i,x0]): x0 -=1     #Find first good point
                while np.isnan(data[i,x1%MAX]): x1 +=1 #Find last good point
                delta = (data[i,x1%MAX] - data[i,x0])/(x1-x0)
                for j in range(x0+1, x1):
                    data[i,j%MAX] = data[i,(j-1)%MAX] + delta


def replaceNaNSVD(data, L):
    """Replaces missing values in the matrix by estimates them as
    a least square supperposition of the L top eigenvectors of the
    row space.

    Ref: Alter, O et. al.  PNAS 100 (6), pp. 3351-3356 (March 2003)

    The data,  when callled as ma.replaceNaNSVD(2), becomes:
        [1, 2, 1,   3]             [1, 2, 1,     3]
        [2, 2, nan, 3]    ==>      [2, 2, 1.087, 3]
        [1, 4, 0,   2]             [1, 4, 0,     2]

    @param L: Number of singular vectors to use in estimating
            missing values.  Has to fullfill L<=min(nGenes,
            nExps).
    """
    nGenes, nExps = data.shape
    indexes = nonzero(sum(~np.isnan(data), 1)/float(nExps) >=1)
    missingRows=set(range(nGenes)).difference(indexes)
    fullmatrix=take(data, indexes)
    U, S, V = svd(fullmatrix, 0)
    V=transpose(V)
    V=V[:,0:L]
    for geneN in missingRows:
        idx2 = nonzero(~np.isnan(data[geneN, :]))  #Places not missing data
        notidx2 = nonzero(np.isnan(data[geneN, :])) #places of missing data
        coeffs = np.dot(pinv(take(V, idx2)),take(data[geneN,:],idx2))
        for idx in notidx2:
            data[geneN, idx] = np.dot(coeffs, V[idx,:])
    


#--------------------------------------------------------------            
#  Normalization methods         
#--------------------------------------------------------------
def scale(data, center=True, scale=True):
    """Normalizes columns of data so that the average is 0 and std is 1 
    That is T_{:j} = T_{:j} - mean(T_{:,j}) and possibly  standard deviation=1:
    That is T_{:j} = T_{:j}/sqrt(T_{:,j} \cdot T_{:,j}).
 
    data: a numeric matrix
    center: logical value indicating weather to mean center the rows of x
    scale: logical value indicating weather to scale the standard deviation of the rows
    """
    if center:
        data = data - data.mean(0)
    if scale:
        data = data/data.std(axis=0, ddof=1)
    return data
    

def normFrobenius(data):
    """Normalizes the whole dataset such that the Frobenius norm is 1.
    That is $T=T/||T||_F$

    B{Example:} ma.normFrobenius()
    """
    from numpy.linalg import norm
    data=data/norm(data)


def QaD_correlation(values, classes, isFactor=None):
    """Given a set of values or class assignments determines correlation/enrichment 
    in list of classificiations.

    Uses, correlation or fisher tests
    
    Arguments:
    - `values`:  Vector of values (either numeric or factors)
    - `classes`: list of lists of classifications (either numeric or factors)
    - `isFactor`: list of same length as classes of T/F indicating weather each 
                  class is a factor or continuous variable (defaults to None) meaning
                  strings are used factors and numbers as continuous

    Returns a list of p-values one for each classification in classes
    """

    isFactor = [None]*len(classes) if isFactor==None else isFactor
    pVals = []
    for classification, factor in zip(classes, isFactor):
        classification = np.asarray(classification)
        #If it is a factor perform ANOVA across groups
        if ((classification.dtype.type in [np.string_,  np.object_, np.bool_, np.unicode_]) or 
            factor):
            groups = [values[classification==l] for l in list(set(classification))]
            f_val, p_val = stats.f_oneway(*groups) 
            pVals.append(p_val)
            #print list(set(classification)), p_val
        #If it is not a factor perform correlation analysis
        else:
            m, b, r_val, p_val, slope_std_error = stats.linregress(values, classification)
            pVals.append(p_val)
    return pVals


def QaD_SVD(d, colorLabels=None, labels=None):
    "d is data matrix and colors is used for color coding the dots."
    u,s,vt = svd(d,0)
    fracs = s**2/np.sum(s**2)
    entropy = -sum(fracs*log(fracs))/np.log(np.min(vt.shape));
    if labels is None:
        labels=range(1, vt.shape[1]+1)
    nGenes, nExps= d.shape
    #Plot Standard SVD plot
    pylab.figure(figsize=(12,4))

    pylab.subplot(1,3,1)
    pylab.imshow(vt, cmap=dataPlot.blueyellow, interpolation='nearest')
    pylab.ylabel('Eigengenes')
    pylab.title('(a) Arrays')
    pylab.xlabel('Arrays')
    pylab.yticks(np.arange(vt.shape[0]), range(1, vt.shape[0]+1))
    pylab.xticks(np.arange(vt.shape[1]), labels)
    pylab.setp(pylab.gca().get_xticklabels(), rotation=45, fontsize=8)

    pylab.subplot(1,3,2)
    pylab.bar(range(1,min(10,nExps)+1), fracs[:min(10,nExps)], width=.8); 
    pylab.ylabel('% Variance')
    pylab.xlabel('Singular Value')
    pylab.xticks(np.arange(1, 11)+.4, np.arange(1, 11))
    pylab.title('(b) Eigenexpression Fraction d=%0.2g' % entropy)

    pylab.subplot(1,3,3)
    pylab.plot(vt[:min(4,vt.shape[0]),:].T, '-o');
    pylab.title('(c) EigenGenes')
    pylab.xlabel('Arrays')
    pylab.ylabel('Expression Level')
    pylab.grid('on')
    pylab.legend(range(1, min(4, vt.shape[0])+1))
    pylab.xticks(np.arange(vt.shape[1]), labels)
    pylab.setp(pylab.gca().get_xticklabels(), rotation=45, fontsize=8)

    pylab.subplots_adjust(left=.07, bottom=None, right=.95, top=None, wspace=.22, hspace=None)

    #Plot Standard PCA plot
    pylab.figure(figsize=(14,14))
    for i in range(4):
        pylab.subplot(2,2,i+1)
        __pcPlot(vt, fracs, colorLabels, i, i+1)

    return u, s, vt
    
     
def __pcPlot(vt, fracs, colorLabels, ax1, ax2):
    """Plots a PC plot base on most enriched colorLabels
    
    Arguments:
    - `vt`:
    - `fracs`:
    - `colorLabels`:
    """
    #Determine the colors of spots
    if  colorLabels is None:
        colors = 'b'
    else: #It is a list or list of lists  #Make sure it works for dataFrames
        if np.asarray(colorLabels, np.object).ndim == 2: #It contains multiple classifications
            pvals = QaD_correlation(vt[ax1,:], colorLabels)
            i = np.nanargmin(pvals)
            #Pick the correct one:
            colorLabels = colorLabels[i]
        #colorLabels is now a list 
        colorTextLabels = sorted(list(set(colorLabels)))
        colors = np.asarray([colorTextLabels.index(val) for val in colorLabels])
    ax = pylab.scatter(vt[ax1,:], vt[ax2,:], c=colors, linewidth=0, s=50, alpha=.7)
    pylab.xlabel('PC%i (%2.1f%%)' %(ax1+1, fracs[ax1]*100))
    pylab.ylabel('PC%i (%2.1f%%)' %(ax2+1, fracs[ax2]*100))
    if colorLabels is not None:
        lines=[]
        for c, label in enumerate(colorTextLabels):
             lines.append(pylab.Rectangle((0, 0), 1, 1, fc=ax.get_cmap()(ax.norm(c))))
        #if len(colorTextLabels)<10:
        pylab.legend(lines, colorTextLabels, loc=4, fontsize=7)
        #else:
        #    pylab.colorbar()




def QaD_HOSVD(t, colors=None):
    """Peforms a quick and dirty hosvd generating several graphics
    
    Arguments:
    - `t`: tensor
    """
    (Z, Un, Sn, Vn) = tensor.hosvd(t)

    nGenes, nExps, nTreats= t.shape
    fractionsn=[]
    entropyn=[]
    for i in range(t.ndim):
        fractionsn.append(Sn[i]**2/np.sum(Sn[i]**2))
        entropyn.append(-np.sum(fractionsn[i]*np.log(fractionsn[i]))/np.log(np.min(size(Vn[i]))))
    
    #Determine ordering
    idx = np.argsort(-np.abs(Z), axis=None)
    zSort=Z.flatten()[idx]
    fractions = zSort**2./sum(zSort**2);
    entropy = -np.sum(fractions*np.log(fractions))/np.log(len(fractions))

    indexes = np.unravel_index(idx, Z.shape)
    labels = []
    for i in range(min(20, len(fractions))):
        labels.append('')
        for j in range(len(indexes)):
            labels[-1]+='%1i,'% indexes[j][i]
        labels[-1]=labels[-1][:-1]
        print '%2i\t(%s) \t%5.3g\t%5.3g%%' %(i, labels[i], Z.flatten()[idx[i]], fractions[i]*100)

    #Plot modes
    for i in range(t.ndim):
        pylab.figure(figsize=(12,8))
        pylab.subplot(2,3,1)
        pylab.imshow(Un[i].T, cmap=dataPlot.blueyellow, interpolation='nearest')
        pylab.axis('tight')
        pylab.ylabel('Eigengenes'); pylab.xticks([]); pylab.yticks([])
        pylab.title('Arrays U_%i' % (i+1))
        pylab.subplot(2,3,4); 
        pylab.plot(Un[i][:,:3], '.-'); pylab.xticks([])

        pylab.subplot(1,3,2); 
        pylab.bar(range(1,min(10,len(fractionsn[i]))+1), fractionsn[i][:min(10,len(fractionsn[i]))], width=.8)
        pylab.ylabel('% Variance'); pylab.xlabel('Singular Value')
        pylab.xticks(np.arange(1, min(10,len(fractionsn[i]))+1)+.4, np.arange(1, min(10,len(fractionsn[i]))+1))
        pylab.title('Fraction d=%0.2g' % entropyn[i])

        pylab.subplot(2,3,3)
        pylab.imshow(Vn[i], cmap=dataPlot.blueyellow, interpolation='nearest')
        pylab.axis('tight'); pylab.xticks([])
        pylab.title('V_%i' %(i+1))
        pylab.subplot(2,3,6); 
        pylab.plot(Vn[i][:,:3], '.-')
        
        pylab.suptitle('Mode %i' % (i+1))
        pylab.subplots_adjust(left=0.1, right=0.95)

    #Plot top patterns
    pylab.figure(figsize=(12,4))
    N=min(10,len(fractions))
    pylab.subplot(1,t.ndim+1,1)
    pylab.barh(range(1,N+1), fractions[:N], height=.8)
    pylab.yticks(np.arange(1, N+1)+.4, labels[:N])
    pylab.gca().set_ylim(pylab.gca().get_ylim()[::-1])
    pylab.grid('on')
    pylab.xlabel('% Variance')

    for i in range(t.ndim):
        pylab.subplot(1,t.ndim+1,i+2)
        pylab.plot(Un[i][:,:3], '.-'); pylab.xticks([])
        pylab.title('U_%i' %(i+1))
    pylab.subplots_adjust(left=0.05, right=0.97)
 
    return (Z, Un, Sn, Vn, indexes)


if __name__ == '__main__':
    data = np.loadtxt('/Users/lom/Dropbox/Sage/log/20121102/yeast_rnaseq_data.csv', delimiter='\t', dtype='|S')
    geneId = data[1:,0]
    geneNames = data[1:,1]
    expLabels= data[0,:][2:11]
    data = np.asarray(data[1:, 2:11], dtype='float')

    print data.shape, expLabels

    data=data[:, [0,1,2,4,5,7,8]]
    expLabels=expLabels[[0,1,2,4,5,7,8]]
    data =  scale(data, center=True, scale=True)
    assert np.all(np.diag(np.dot(data,data.T))-1 < 1e-12) and np.all(np.mean(data, 1)<1e-12)

    #Test output of SVD
    u, s, vt = QaD_SVD(data, [1,1,1,2,2,2,3,3,3], [s.split('.')[0] for s in expLabels])

    #Convert to tensor and perform HOSVD
    t=data.reshape((-1, 3,3))  #genes x repeats x treatment
    (Z, Un, Sn, Vn, indexes) = QaD_HOSVD(t)

