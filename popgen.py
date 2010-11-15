import  sys,pyximport,  numpy as np
sys.path.append('..')
pyximport.install()
import fileReader

def nchoose2(n):
    return n*(n-1)/2

def fst(vals0, vals1, isNorm=True):
    nij=np.asarray([vals0.shape[1], vals1.shape[1]], np.float)
    nAll=nij.sum()
    xi=((vals0==-1).sum(1)+(vals1==-1).sum(1))/nAll
    xij=[(vals0==-1).sum(1)/nij[0], (vals1==-1).sum(1)/nij[1]]
    if isNorm:
        top=0; bottom=0
        for j in range(2):
            top+=nchoose2(nij[j])*np.sum(2*nij[j]/(nij[j]-1)*xij[j]*(1-xij[j]))
            bottom+=nchoose2(nij[j])
        return 1-top/bottom/np.sum(2*nAll/(nAll-1)*xi*(1-xi))
    return 1-(np.sum((1-xij[0])*xij[0])+np.sum((1-xij[1])*xij[1]))/2/ np.sum((1-xi)*xi)

def fstFile(file1, file2, isNorm=True):
    """Calculates fst based on data stored in two tab delimited files.
    Each file must contain the following columns: snp names; snp position; snp values. """

    files=fileReader.concurrentFileReader(file1, file2, key=1)
    subjects=files.next()
    vals=[]
    for (snpName, snpLocation, snps) in files:
        vals.append(fileReader.nucleotides2Haplotypes(sum(snps, [])))
    vals=np.asarray(vals, np.float)
    nSamples=len(subjects[0])
    return fst(vals[:,:nSamples], vals[:,nSamples:], isNorm)


if __name__ == '__main__':
    print fstFile(*sys.argv[1:])
