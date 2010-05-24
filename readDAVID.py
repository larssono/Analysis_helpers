import numpy as np
import gzip 

class enrichmentClass(object):
    """Given one line from a functional clustering output file stores all data
    """

    
    def __init__(self, str):
        """Given a line from a function classification file stores values
        """
        vals=str.strip().split('\t')
        self.category  =vals[0]
        self.term      =vals[1]
        self.count     =vals[2]
        self.percent   =vals[3]
        self.pvalue    =float(vals[4])
        self.genes     =vals[5]
        self.list_total=int(vals[6])
        self.pop_hits  =int(vals[7])
        self.pop_total =int(vals[8])
        self.fold_enrichment=float(vals[9])
        self.bonferroni	=float(vals[10])
        self.benjamini	=float(vals[11])
        self.fdr        =float(vals[12])


def readDavidFile(filename, pvalFilter=1, isGzip=False):
    enrichments={}
    if isGzip:
        fp=gzip.open(filename)
    else:
        fp=open(filename)
    for l in fp:
        try:
            if l.startswith('Category	Term	Count	%'):
                continue
            davidLine=enrichmentClass(l)
            if davidLine.pvalue<=pvalFilter:
                enrichments[davidLine.term]=davidLine
        except IndexError:
            if l.startswith('Functional Group') or l=='\n':
                pass
            else:
                print 'WARNING!  There are unkown structured line in the file'
    return enrichments

if __name__ == '__main__':
    hosvdDown=readDavidFile('david_enrichment_hosvd_down.csv', 1e-3)
    hosvdUp=readDavidFile('david_enrichment_hosvd_up.csv.gz', 1e-3, True)
    anovaDown=readDavidFile('david_enrichment_anova_down.csv.gz', 1e-3, True)
    anovaUp=readDavidFile('david_enrichment_anova_up.csv.gz', 1e-3, True)
    defaultEnrichment=enrichmentClass('none\tnone\t0\t0\t1\t\t0\t0\t0\t0\t1\t1\t1\t1')
    
    downKeys=list(set(hosvdDown.keys()).union(anovaDown.keys()))
    upKeys=list(set(hosvdUp.keys()).union(anovaUp.keys()))

    def printClasses(keys, hosvd, anova):
        for k in keys:
            hosvdP=hosvd.get(k,defaultEnrichment).pvalue
            anovaP=anova.get(k, defaultEnrichment).pvalue
            hosvdCount=hosvd.get(k,defaultEnrichment).count
            anovaCount=anova.get(k,defaultEnrichment).count 
            category=hosvd.get(k, defaultEnrichment).category
            print '%s\t%40s\t%0.2e\t%s\t%0.2e\t%s' %(category,k,
                                                     hosvdP, hosvdCount,
                                                     anovaP, anovaCount)
    print 'Category\tName\tHOSVD pvalue\tHOSVD count\tANOVA pvalue\tANOVA count'
    #printClasses(downKeys, hosvdDown, anovaDown)
    printClasses(upKeys, hosvdUp, anovaUp)

#print enrichments.keys()
