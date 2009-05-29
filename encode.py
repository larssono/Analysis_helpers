import numpy as n

def encode(data):
    """Converts a matrix of genotypes where the individuals are in rows
and SNPs are in columns.  Each SNP is represented by two columns for
the two alleles.  This is accomplished in two steps:

1) generate coding (break if more that two alleles
2)Second go through and return -1,0,1 mapping
    """
    nSNPs=data.shape[1]/2
    nIndividuals=data.shape[0]
    code=[[0,0]]*(nSNPs)
    for i in range(0, nSNPs*2, 2):
        code[i/2]=list(n.unique(data[:,i:i+2].flatten()))
        try:
            code[i/2].remove(0.)
        except:
            pass
        if len(code[i/2]) >2:
            print "Too many alleles to handle with encoding for SNP %g" %(i/2 +1)
            exit()
    out=n.zeros((nIndividuals, nSNPs), dtype=int);
    for i in range(nIndividuals):
        for j in range(0,nSNPs*2,2):
            if data[i,j]==0 or data[i, j+1]==0 or data[i,j]!=data[i, j+1]:
                continue
            elif data[i,j]==data[i, j+1]: 
                if data[i,j] == code[j/2][0]:
                    out[i,j/2]=-1
                else:
                    out[i,j/2]=1
            else:
                print 'WRONG COMBINATION OF GENOTYPES'
                raise ValueError()
    return out
