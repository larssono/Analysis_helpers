from numpy.oldnumeric import *
import numpy as np

        

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
    from numpy.linalg import svd, pinv 
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
        coeffs = dot(pinv(take(V, idx2)),take(data[geneN,:],idx2))
        for idx in notidx2:
            data[geneN, idx] = dot(coeffs, V[idx,:])
    


#--------------------------------------------------------------            
#  Normalization methods         
#--------------------------------------------------------------
def scale(data, center=True, scale=True):
    """Normalizes each array so that the average expression is 0.
    That is T_{:j} = T_{:j} - mean(T_{:,j}) and possibly  standard deviation=1:
    That is T_{:j} = T_{:j}/sqrt(T_{:,j} \cdot T_{:,j}).
 
    data: a numeric matrix
    center: logical value indicating weather to mean center the rows of x
    scale: logical value indicating weather to scale the standard deviation of the rows
    """
    if center:
        data=data-np.mean(data, 0);
    if scale:
        data=data/sqrt(diag(dot(transpose(data),data)))
    

def normFrobenius(data):
    """Normalizes the whole dataset such that the Frobenius norm is 1.
    That is $T=T/||T||_F$

    B{Example:} ma.normFrobenius()
    """
    from numpy.linalg import norm
    data=data/norm(data)

