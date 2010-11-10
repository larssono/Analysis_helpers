from scipy.linalg import svd
import numpy as np


def flatten(T, n):
    """D = FLATTEN(T, N) makes a matrix out of a tensor.
    Such that the fibers along dimension N are aligned along the
    columms of D.
    
    SYNOPSIS:  
        T - a tensor
        N - The mode in which to flatten T into.   
    
    AUTHOR Larsson Omberg lom@larssono.com
    DATE   05-June-2009"""

    n=n-1
    if n>T.ndim | n<=0:
        raise ValueError('n has to be between 1 and the number of dimensions of T')
    
    order=[n,]
    order.extend(range(n+1, T.ndim))
    order.extend(range(0,n))

    nrows=T.shape[n]
    ncols=np.prod(T.shape)/nrows

    return np.reshape(np.transpose(T, order),(nrows,ncols), order='FORTRAN')

def hosvd(T, saveSpace=False):
    """HOSVD    N-mode SVD decomposition of N-way tensor 
    (Z, Un, Sn, Vn) = HOSVD(T) decomposes the N-way tensor D into N
    orthogonal matrices stored in Un and a core tensor Z.  Also
    returns the n Tucker1 results in Sn and Vn if saveSpace=False (default)
    
    Author: Larsson Omberg <lom@larssono.com>
    Date:   11-NOV-2007, 28-February-2008, 5-June-2009, 4-Nov-2010"""
    #Find orthogonal matrices
    Un=[]; Vn=[]; Sn=[]
    for n in range(T.ndim):
        Tn = flatten(T, n+1);  #FIX 
        if Tn.shape[1] < Tn.shape[0]:
            [U,S,V] = svd(Tn,0);
            V=V.T
        else:
            [V,S,U] = svd(Tn.T,0);
            U=U.T
        Un.append(U);
        if not saveSpace:
            Vn.append(V);
            Sn.append(S);
    Z=T.copy()
    for n in range(len(Un)):
        Z = nmodmult(Z, Un[n].T, n+1);
    if not saveSpace:
        return [Z, Un, Sn, Vn]
    return [Z, Un]


def unflatten(T, sizes, n):
    """D = UNFLATTEN(T, SIZES, N) remakes a tensor out of a matrix.
    Such that the mode N flattening of the tensor will return the
    matrix T.
    
    SYNOPSIS:  
       T - a matrix of size mxn
       SIZES - size of the output matrix prod(sizes) must equal mxn
       N - The mode in which to unflatten t into.   
       
    AUTHOR Larsson Omberg lom@larssono.com
    DATE   21-January-2005, 5-june-2009"""

    if np.prod(sizes) != np.prod(T.shape) or T.ndim !=2:
        raise ValueError('matrix and output tensor must have same number of elements')
    if n>len(sizes) | n<=0 :
        raise ValueError('n has to be between 1 and the number of dimensions of T')

    order=[n-1,]
    order.extend(range(n, len(sizes)))
    order.extend(range(n-1))

#ndx = [order ones(1,ndims(b)-length(order))];
#ndx(ndx) = 1:max(length(order),ndims(b));  % Inverse permutation order

    ndx = np.arange(len(order))
    ndx[order] = np.arange(len(order))  # Inverse permutation order

    d=np.reshape(T,sizes[order], order='FORTRAN');
    return np.transpose(d, ndx);  #np.transpose(...)


def nmodmult(A, B, n):
    """NMODMULT  peforms the n-mode multiplication of a tensor and matrix
    T = NMODMULT(A, B, N)
      A is a tensor of order < N
      B is a matrix to be multiplied
      N is the order of the multiplication
      return:  tensor product of A x_N B"""
    Asize = np.asarray(A.shape);
    
    Asize[n-1] = B.shape[0];
    
    An=flatten(A,n)
    T = unflatten(np.dot(B,An), Asize, n);
    return T


def test():
    T=np.reshape(np.arange(1,121).T, (3,4,5,2), order='FORTRAN')
    z=np.zeros((3,4,5,2))
    z[:,:,0,0]=[[759.2225,-0.7251,-0.0000, 0.0000],[-0.1844,-0.3383, 0.0000, 0.0000],[-0.0000, 0.0000,-0.0000,-0.0000]]
    z[:,:,1,0]=[[0.1619, 7.3968,-0.0000, 0.0000],[ 1.7977, 0.1504, 0.0000,-0.0000],[-0.0000, 0.0000, 0.0000,-0.0000]]
    z[:,:,2,0]=1.0e-13*np.asarray([[-0.4852, 0.1526,-0.1303, 0.1039],[ 0.0188,-0.0283,-0.0208, 0.0976],[ 0.0109,-0.0099, 0.0548, 0.0638]])
    z[:,:,3,0]=1.0e-13*np.asarray([[-0.5339, 0.0717, 0.1782,-0.0199],[ 0.1186,-0.0234, 0.0284,-0.0507],[ 0.0138,-0.0078,-0.0339,-0.0368]])
    z[:,:,4,0]=1.0e-13*np.asarray([[ 0.3636,-0.0987,-0.1282, 0.0898],[-0.0152,-0.0330, 0.0273, 0.0185],[-0.0194,-0.0093,-0.0181, 0.0594]])
    z[:,:,0,1]=[[ 0.1028,15.1213, 0.0000, 0.0000],[ 3.6749, 0.3076, 0.0000,-0.0000],[-0.0000,-0.0000, 0.0000,-0.0000]]
    z[:,:,1,1]=[[ -80.3200,-6.7990,-0.0000, 0.0000],[-1.6515,-0.1041, 0.0000,-0.0000],[ 0.0000,-0.0000, 0.0000, 0.0000]]
    z[:,:,2,1]=1.0e-13*np.asarray([[-0.3285,-0.0554,-0.0873, 0.0506],[-0.0126,-0.0161,-0.0233, 0.0329],[ 0.0099,-0.0317,-0.0096, 0.0146]])
    z[:,:,3,1]=1.0e-13*np.asarray([[-0.1950,-0.0022, 0.1887, 0.0545],[ 0.0669,-0.0149, 0.0019, 0.0006],[ 0.0032, 0.0422, 0.0003,-0.0039]])
    z[:,:,4,1]=1.0e-13*np.asarray([[ 0.4543,-0.0320,-0.0597, 0.0600],[ 0.0024, 0.0189,-0.0007, 0.0089],[-0.0012, 0.0078,-0.0160, 0.0089]])
    u=[]
    u.append(np.asarray([[-0.5701,-0.7129,0.4082],[-0.5773,-0.0059,-0.8165],[-0.5845,0.7012,0.4082]]))
    u.append(np.asarray([[-0.4715,-0.6912, 0.5477, 0.0054],[-0.4902,-0.2443,-0.7262,-0.4154],[-0.5089, 0.2025,-0.1906, 0.8147],[-0.5276, 0.6493, 0.3691,-0.4046]]))
    u.append(np.asarray([[-0.2959, 0.7158,-0.4263, 0.0037, 0.4672],[-0.3660, 0.4075, 0.0218,-0.0036,-0.8364],[-0.4361, 0.0991, 0.7542, 0.4063, 0.2570],[-0.5062,-0.2093, 0.1315,-0.8165, 0.1265],[-0.5763,-0.5176,-0.4812, 0.4101,-0.0143]]))
    u.append(np.asarray([[-0.3431, 0.9393],[-0.9393,-0.3431]]))

    Z,Un,Sn,Vn = hosvd(T)

    print 'Max Error: z=%0.2g,\tU1=%0.2g,\tU2=%0.2g,\tU3=%0.2g,\tU4=%0.2g' %(abs(Z-z).max(), abs(Un[0]-u[0]).max(), abs(Un[1][:2,:2]-u[1][:2,:2]).max(),abs(Un[2][:2,:2]-u[2][:2,:2]).max(),abs(Un[3]-u[3]).max())
    return Z, Un

