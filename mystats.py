import scipy.io
import numpy as np

def anova1(y, x, returnB=False, returnR2=False):
    from scipy.linalg import lstsq, inv
    import scipy.stats
    
    y = y[...,np.newaxis] if y.ndim==1 else y;
    x = x[...,np.newaxis] if x.ndim==1 else x;

    #Remove missing values i.e. NaNs
    idx=np.logical_not(np.logical_or(np.isnan(x).sum(1),np.isnan(y).sum(1)))
    x=x[idx,:]
    y=y[idx]
    n, m=x.shape

    #this is under the assumption that the variables are continuous and not categorical
    x=np.hstack([np.ones((n, 1)), x]) 

    #Solve for the parameters b
    b, SSE, rank, s = lstsq(x, y)
    if rank<(m+1):
        return np.array([float('NaN')]*m)
    SST = np.dot(np.dot(y.T, np.eye(n)-1.0/n*np.ones((n, n))), y)
    MSE=SSE/(n-m-1)
    s2b=MSE*inv(np.dot(x.T,x))
    t=map(lambda x,y: (x/y)[0], b, np.sqrt(np.diag(s2b)))
    p=scipy.stats.t.sf(np.abs(t),n-m-1)*2  #sf = 1-cdf
    if returnB:
        if returnR2:
            return p[1:], b, 1-SSE/SST
        else:
            return p[1:], b
    return p[1:]



def removeNaNs(data, groups):
    data=np.asarray(data)
    groups=np.asarray(groups)
    idx=~np.isnan(data)
    data=data[idx]
    groups=groups[:,idx]
    return data, groups

def _generateModel(nFactors, level=0):
    if level==0:
        level=nFactors
    q=list(np.ndindex((2,)*nFactors))
    q.sort(lambda x,y: sum(x)-sum(y))
    return filter(lambda x:sum(x)<=level and sum(x)>0, q)


def anovan(data, groups, model=1):
    data, groups = removeNaNs(data, groups)
    if len(groups.shape)==1: 
        groups=[groups]
    if type(model)==type(3):
        model = _generateModel(len(groups), model)
    uniqueGroups=[np.unique(group) for group in groups]
    nClasses = [len(group) for group in uniqueGroups]
    ns=np.zeros(nClasses)
    means=np.zeros(nClasses)
    sse=0
    sst=0
    grandMean=np.mean(data)

    #Find mean and number of measurements for each sub-group 
    for idx in np.ndindex(tuple(nClasses)):
        pos=set(np.nonzero(uniqueGroups[0][idx[0]]==groups[0])[0])
        for dim, i in enumerate(idx[1:]):
            if len(pos)==0: break
            pos=pos.intersection(np.nonzero(uniqueGroups[dim+1][i]==groups[dim+1])[0])
        if len(pos)==0: continue
        pos = list(pos)
        means[idx]=np.mean(data[pos])
        ns[idx]=len(data[pos])
        sse+=np.sum(np.power((data[pos]-means[idx]),2))
        sst+=np.sum(np.power((data[pos]-grandMean),2))
    
    #Start calculating terms and cross terms ss
    #Type 1 SS
    for dim in range(len(nClasses)):
        ssa=np.sum(ns.sum(dim)*(means.mean(dim) - grandMean)**2)
        print 'X%g\t%0.4g\t%g' % (dim, ssa, ns.shape[dim]-1)
    
    print '%s\t%0.4g' % ('Err', sse)
    print '%s\t%0.4g' % ('Tot', sst)
    
    return means

         

def main():

#     vals=scipy.io.loadmat('/home/lom/bin/matlab/stats/carbig.mat')
#     anovan(vals['MPG'], [vals['Origin'],  vals['Model_Year'], vals['Cylinders']])
#     anovan(vals['MPG'], vals['Cylinders'])
    print
    classes = [['f']*20+['m']*20,['n']*10+['c']*10+['n']*10+['c']*10]
    vals=[6,8,5,6,7,8,4,9,6,5,6,3,4,2,4,6,4,3,5,3,7,9,8,7,9,10,6,6,10,8,5,3,2,4,4,3,3,1,5,4]
    anovan(vals, classes, model=2)


    return vals

if __name__ == "__main__":
    vals = main()
