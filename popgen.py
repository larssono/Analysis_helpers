import  sys,pyximport,  numpy as np
sys.path.append('..')
pyximport.install()
import fileReader

def nchoose2(n):
    return n*(n-1)/2

def fst(file1, file2):
    files=fileReader.concurrentFileReader(file1, file2, key=1)
    subjects=files.next()
    vals=[]
    for (snpName, snpLocation, snps) in files:
        vals.append(fileReader.nucleotides2Haplotypes(sum(snps, [])))
    vals=np.asarray(vals)
    nAll=float(vals.shape[1])
    nij=np.asarray(map(len, subjects), np.float)
    popIdx=[range(nij[0]), range(nij[0],vals.shape[1])] 

    xi=((vals==-1).sum(1)/nAll)
    xij=[]

    for i, n in enumerate(nij):
        idx=popIdx[i]
        xij.append((vals[:,idx]==-1).sum(1)/n)

    top=0; bottom=0
    for j in range(2):
        top+=nchoose2(nij[j])*np.sum(2*nij[j]/(nij[j]-1)*xij[j]*(1-xij[j]))
        bottom+=nchoose2(nij[j])

    return 1-top/bottom/np.sum(2*nAll/(nAll-1)*xi*(1-xi)),  1-(sum((1-xij[0])*xij[0])+sum((1-xij[1])*xij[1]))/2/ sum((1-xi)*xi)

if __name__ == '__main__':
    print fst(*sys.argv[1:])
