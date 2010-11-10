import numpy as np
import scipy.stats as stats

def encode(data):
    """Converts a matrix of genotypes where the individuals are in rows
    and SNPs are in columns.  Each SNP is represented by two columns for
    the two alleles.  This is accomplished in two steps:

    1) generate coding (break if more that two alleles
    2) Second go through and return -1,0,1 mapping
    """
    nSNPs=data.shape[1]/2
    nIndividuals=data.shape[0]
    code=[[0,0]]*(nSNPs)
    for i in range(0, nSNPs*2, 2):
        code[i/2]=list(np.unique(data[:,i:i+2].flatten()))
        try:
            code[i/2].remove('0')
            code[i/2].remove(0.)
        except:
            pass
        if len(code[i/2]) >2:
            print "Too many alleles {%s} to handle with encoding for SNP %g" %((i/2 +1), `code[i/2]`)
            raise ValueError()
    out=np.zeros((nIndividuals, nSNPs), dtype=float);
    for i in range(nIndividuals):
        for j in range(0,nSNPs*2,2):
            if data[i,j]=='0' or data[i, j+1]=='0': #Missing value case
                out[i,j/2]=float('nan')
            elif data[i,j]==data[i, j+1]: 
                if data[i,j] == code[j/2][0]:
                    out[i,j/2]=-1
                else:
                    out[i,j/2]=1
    return out


#------------------------------ Find top Changing Genes ------------------------------------
def topChangingExpressionProbes_ttest(data1, data2, pVal=0.05, absDiff=0, benjamini=True):
    """Finds probes that are significantly differently expressed between
    data1 and data2.

    Parameters:
       data1:  a matrix n x m1 of expression data for n probes and m1 samples
       data2:   a matrix n x m2 of expression data for n probes and m2 samples
       pVal:   (default 0.05) pvalue cutoff for significance
       absDiff:(default 0)  Minimum fold change in expression between data1 and data1,
            Log_2( |mean(data1[i,:])/mean(data2[i,:])| ) >= absDiff
       benjamini: (T/F) wheter to use Benjamini-Hochtberg correction.

    Return:
         idx:  array of size n with T/F for each gene passing criteria

    """
    (nGenes, nExps)=data1.shape
    (t, p)=stats.ttest_ind(data1, data2, axis=1)

    #Get significant genes by Benjamini-Hodgberg
    if benjamini:
        minP=mystats.benjamini_hochberg(p, pVal)
        idx=p<=minP
    else:
        idx = p<pVal
    if absDiff != 0:
        foldChange=data1.mean(1)/data2.mean(1)
        idx2 = np.logical_or(foldChange >= absDiff, foldChange<=1/absDiff)
        return np.logical_and(idx, idx2)
    return idx
