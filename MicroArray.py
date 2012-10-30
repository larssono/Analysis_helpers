from numpy.oldnumeric import *
    
__author__ = "Larsson Omberg <lom@larssono.com>"
__date__ = "7 March 2005"


class MicroArray(object):
    """ Contains data from set of microarray experiments.
    Allows manipulations of data.

    B{Example Usage:}
    
    >>> #Read in the data
    >>> ma = MicroArray('path/to/tab/delimited/file', 2)

    >>> #See size of data
    >>> print ma.data.shape       

    >>> #Extract multiple experiments by specifying
    >>> # the slide names or experiment names 
    >>> ma=ma.getExps(['yB12n099','yB12n100','yB12n138'])

    >>> #if experimental names contain multiple conditions
    >>> # we can split those labels 
    >>> ma.splitExpNames(r'(.*) vs. (.*) (P[0-9]*) (.*?) ([0-9]*\.?[0-9]*).*hr \\
    ...    set *([0-9]*) *(.*)', 7)

    >>> #Filter out all the genes or rows with missing values
    >>> ma.filterNaN(1.0)  

    >>> #Normalize so that each array has the same range of expression values.
    >>> ma.normArrayScale()

    >>> #Save the output as a matlab readable and tab delimited file
    >>> ma.saveMatlab('filename.mat')
    >>> ma.saveTab('filaname.csv')
    
    """
    __NONUMB    = -999.9999

    #-------------------------------------------------------------------
    # Constructors
    #-------------------------------------------------------------------

    def __init__(self, p1='', p2='', p3='', p4=''):
        """Constructor called with two arguments: filename and nGeneAnots;

        @param p1: B{filename} - points to a tab delimeted file in which the
           first row contains columm headers and the first nGeneAnots
           columns contain row headers.

        @param p2: B{nGeneAnots} -  Number of columns of rowHeaders
        """
        import types
        if type(p1)==types.StringType:
            self.__fromFile(p1, p2)
        else:
            self.__fromData(p1,p2,p3,p4)


    def __fromData(self, data, expNames, geneNames, expClass):
        self.data = array(data);
        self.expNames=array(expNames, PyObject)
        self.expClass=array(expClass, PyObject)
        self.geneNames=array(geneNames, PyObject)
        self.nGenes=len(self.geneNames)
        self.nExps=len(self.expNames)
        
        self.data.shape =(self.nGenes, self.nExps)
        self.expNames.shape=(self.nExps,)
        self.geneNames.shape=(self.nGenes,)
        
    def __fromFile(self, filename, nGeneAnots):
        f = open(filename)
        firstLine=f.readline()
        lines=f.readlines()
        firstLine = firstLine.strip().split('\t')

        self.expNames = array(firstLine[nGeneAnots:], PyObject)
        self.nGenes = len(lines)
        self.nExps = len(self.expNames)
        self.data=ones((self.nGenes, self.nExps), 'f')*self.__NONUMB
        self.geneNames=zeros(self.nGenes, PyObject)
        self.expClass=zeros(self.nExps, PyObject)
        i = 0;
        for line in lines:
            line = line.strip()
            cols = line.split('\t')
            self.geneNames[i] = cols[0]#cols[0:nGeneAnots]
            for j in range(nGeneAnots, len(cols)):
                if cols[j]!='':
                    self.data[i][j-nGeneAnots]=float(cols[j])
            i=i+1
        f.close()

    #-------------------------------------------------------------------
    # Manipulation routines
    #-------------------------------------------------------------------

    def splitExpNames(self, regexp, size):
        """splitExpNames - Splits experimental name up into
        classifications

        Each experimental name lists many criterea that was fullfilled
        for that experiments such as time point, cell type, chemical
        environment etc.  In order to be able to sort experiments and
        filter simillar experiments this method takes a regular
        expression and populates MicroArray.expClass with the
        classifications in each column.

        B{Example}
        
        Assume the experiment names look like::
        
              [UHR vs. WI-38 P7 PDGF 4hr set 3 
               UHR vs. WI-38 P8 PDGF 8hr set 1 
               UHR vs. WI-38 P8 PDGF 8hr set 2 
               UHR vs. WI-38 P7 PDGF 8hr set 3 
               UHR vs. WI-38 P6 Serum 0 hr set 1 
               UHR vs. WI-38 P6 Serum 0 hr set 1 hyb2]

        Then the regular expression::
        
        r'(.*) vs. (.*) (P[0-9]*) (.*?) ([0-9]*\.?[0-9]*).*hr set *([0-9]*) *(.*)'

        will split these names into 7 classifications.

        @param regexp: the regular expression with () groupings around
           each type of classification

        @param size: the number of classifciations


        """
        import re
        classifications=[]
        r = re.compile(regexp)
        for heading in self.expNames:
            match = r.match(heading);
            if match == None:
                print "ERROR:   ", heading, " didn't parse"
            else:
                classifications.append(match.groups())
        self.expClass = array(classifications, PyObject)
        self.nClassifications = self.expClass.shape[1]

    def sortExp(self, aOrdering):
        """Sorts Microarrays by ordering experiments according to
        classifications given in sequence nOrdering.  I.e. we sort
        first by classification in nOrdering[0] then without changing
        order of these we sort according to nOrdering[1] etc.

        @param aOrdering : sequence of classifications to sort after.
        """            
        aStrSort = array(('',)*self.nExps, PyObject);
        for i in aOrdering:
            for j in range(len(aStrSort)):
                aStrSort[j]=aStrSort[j]+self.expClass[j, i]+'-';
        sortIndex = argsort(aStrSort)
        self.data = take(self.data, sortIndex,1);
        self.expNames = take(self.expNames, sortIndex)
        self.expClass = take(self.expClass, sortIndex)

    def __isSameExperiment(self, exp1, exp2, requiredSame):
        same = True
        for cl in requiredSame:
            same = same and (self.expClass[exp1, cl]==self.expClass[exp2, cl])
        return same

    def __findSameExperiments(self, requiredSame):
        currExp=0
        uniqueExps =[currExp,]
        while currExp<self.nExps:
            otherExp = currExp
            while self.__isSameExperiment(currExp, otherExp, requiredSame):
                otherExp += 1;
                if otherExp==self.nExps:
                    break
            currExp=otherExp;
            uniqueExps.append(currExp)
        return uniqueExps
        
    def averageRepeats(self, requiredSame):
        """Averages repeated experiments that are same according to
        classifications in requiredSame.

        @param requiredSame: sequence object of integers containing
           classifications that will have to be same.  i.e. [1,3,4]
        """
    
        expCol = 0
        uniqueExps = self.__findSameExperiments(requiredSame)
        
        firstExp = uniqueExps[0];
        for lastExp in uniqueExps[1:]:
            nSameExps = lastExp-firstExp;
            nNumbs    = sum(self.data[:,firstExp:lastExp]!=self.__NONUMB, 1)
            totSum    = sum(self.data[:,firstExp:lastExp], 1) -\
                        (nSameExps-nNumbs)*self.__NONUMB
            for row in xrange(self.nGenes):
                if nNumbs[row] !=0 :
                    self.data[row, expCol] = totSum[row]/nNumbs[row]
                else:
                    self.data[row, expCol] = self.__NONUMB
                    
            self.expNames[expCol] = '%s_'*len(requiredSame) %\
                                    tuple(take(self.expClass, requiredSame, 1)[firstExp])
            self.expClass[expCol,:] = self.expClass[firstExp,:]
            expCol += 1
            firstExp=lastExp
        self.nExps=len(uniqueExps)-1
        self.data=take(self.data, range(self.nExps),1)
        self.expNames=take(self.expNames, range(self.nExps))
        self.expClass=take(self.expClass, range(self.nExps))

    def filterNaN(self, cutoff):
        """Keeps Genes that have a percentage of non missing values
        greater than or equal to cuttoff

        @param cutoff: percentage of experiments that have to be good
           to keep.
        """
        indexes = nonzero(sum(self.data!=self.__NONUMB, 1)/float(self.nExps) >= cutoff)
        self.data=take(self.data, indexes)
        self.geneNames=take(self.geneNames, indexes)
        self.nGenes=len(self.geneNames)
        
    def filterNaNExp(self, cutoff):
        """Keeps Columns that have a percentage of non missing values
        less than or equal to cuttoff

        @param cutoff: percentage of experiments that have to be good
           to keep.
        """
        x = asarray(self.data!=self.__NONUMB, Int)
        indexes = nonzero(sum(x, 0)/float(self.nGenes) >= cutoff)
        self.data=take(self.data, indexes, 1)
        self.expNames=take(self.geneNames, indexes)
        self.nExps=len(self.expNames)


    def saveMatlab(self, filename):
        """Saves microarray class as matlab readable file.

        @param filename: output file that will store tab delimited
           file.
        """
        from scipy.io import savemat
        outputDict={}
 	mask=(self.data==self.__NONUMB)
        putmask(self.data, mask, float('nan'))
        outputDict['data'] = self.data[:, :self.nExps]
        outputDict['geneNames'] = array(list(self.geneNames[:self.nGenes]))
        outputDict['expNames'] = array(list(self.expNames[:self.nExps]))
        savemat(filename, outputDict)
        putmask(self.data, mask, self.__NONUMB)

    def saveTab(self, filename):
        """Saves microarray class as tab delimited file.

        @param filename: output file that will store tab delimited file
        """
	mask = (self.data==self.__NONUMB)
        putmask(self.data, mask, float('nan'))
        fp = open(filename, 'w')
        fp.write(self.__str__())
        putmask(self.data, mask, self.__NONUMB)
        return
    
    def replaceNaNMean(self):
        """Replaces missing values in the matrix with the average
        expression value across all arrays for a specific gene.

        B{Example}:  The microarray::

            [1, 2, 1, 3]
            [2, 2, nan, 3]
            [1, 4, 0, 2]

        becomes::
            [1, 2, 1, 3]
            [2, 2, 2.33, 3]
            [1, 4, 0, 2]
        """
        for i in range(self.nGenes):
            indexes = nonzero(self.data[i, :]==self.__NONUMB)
            mean = (sum(self.data[i, :])-len(indexes)*self.__NONUMB)/(self.nExps-len(indexes))
            for j in indexes:
                self.data[i,j]=mean

    def replaceNaNInterp(self):
        """Replaces missing values in the matrix by linearly
        interpolating between existing values.

        B{Example:}  The microarray::

            [1, 2, 1, 3]
            [2, 2, nan, 3]
            [1, 4, 0, 2]

        becomes::
            [1, 2, 1, 3]
            [2, 2, 2.5, 3]
            [1, 4, 0, 2]
        """
        ##import pylab
        MAX = self.nExps
        for i in range(self.nGenes):
            ##pylab.plot(self.data[i,:], '-o');
            indxs = nonzero(self.data[i, :]==self.__NONUMB)
            x1=0
            if len(indxs) > 0:
                for index in indxs:
                    if index < x1: continue
                    x0=index-1
                    x1=index+1
                
                    while self.data[i,x0] == self.__NONUMB: x0 -=1     #Find first good point
                    while self.data[i,x1%MAX] == self.__NONUMB: x1 +=1 #Find last good point
                    delta = (self.data[i,x1%MAX] - self.data[i,x0])/(x1-x0)
                    for j in range(x0+1, x1):
                        self.data[i,j%MAX] = self.data[i,(j-1)%MAX] + delta

            ##pylab.plot(indxs, take(self.data[i,:], indxs), 'or'); pylab.axis([0, 40,0,5]); pylab.show()

    def replaceNaNSVD(self, L):
        """Replaces missing values in the matrix by estimates them as
        a least square supperposition of the L top eigenvectors of the
        row space.

        Ref: Alter, O et. al.  PNAS 100 (6), pp. 3351-3356 (March 2003)
        
        B{Example:}  The microarray::
            [1, 2, 1,   3]
            [2, 2, nan, 3]
            [1, 4, 0,   2]
                                  
        when callled as ma.replaceNaNSVD(2) becomes::
        
            [1, 2, 1,     3]
            [2, 2, 1.087, 3]
            [1, 4, 0,     2]

        @param L: Number of singular vectors to use in estimating
                missing values.  Has to fullfill L<=min(nGenes,
                nExps).
        """
        from numpy.linalg import svd, pinv 
        indexes = nonzero(sum(self.data!=self.__NONUMB, 1)/float(self.nExps) >=1)
        missingRows=set(range(self.nGenes)).difference(indexes)
        fullmatrix=take(self.data, indexes)
        U, S, V = svd(fullmatrix)
        V=transpose(V)
        V=V[:,0:L]
        for geneN in missingRows:
            idx2 = nonzero(self.data[geneN, :]!=self.__NONUMB)  #Places not missing data
            notidx2 = nonzero(self.data[geneN, :]==self.__NONUMB) #places of missing data
            coeffs = dot(pinv(take(V, idx2)),take(self.data[geneN,:],idx2))
            for idx in notidx2:
                self.data[geneN, idx] = dot(coeffs, V[idx,:])

    #--------------------------------------------------------------            
    #  Normalization methods         
    #--------------------------------------------------------------
    def normArrayCenter(self):
        """Normalizes each array so that the average expression is 0.
        That is T_{:j} = T_{:j} - mean(T_{:,j}).

        B{Example:}  ma.normArrayCener()
        """
        self.data=self.data-mean(self.data, 0);

    def normArrayScale(self):
        """Normalizes each array so that the sum of expression squared is 1.
        That is T_{:j} = T_{:j}/sqrt(T_{:,j} \cdot T_{:,j}).

        B{Example:}  ma.normArrayScale()
        """
        self.datas=self.data/sqrt(diag(dot(transpose(self.data),self.data)))

    def normFrobenius(self):
        """Normalizes the whole dataset such that the Frobenius norm is 1.
        That is $T=T/||T||_F$

        B{Example:} ma.normFrobenius()
        """
        from numpy.linalg import norm
        self.data=self.data/norm(self.data)
        
    #--------------------------------------------------------------            
    #  Sequence methods         
    #--------------------------------------------------------------

    def deleteExps(self, idx):
        """Removes experiments i.e. columns

        @param idx: Indexes of columns (experiments) to be removed.
        """

        
        indexes = [i for i in range(self.nExps) if i not in idx]
        self.data=take(self.data, indexes, 1)
        self.expNames=take(self.expNames, indexes)
        self.expClass=take(self.expClass, indexes)
        self.nExps=len(self.expNames)
        

    def getExpIdx(self, names):
        """Searches through all expNames to find matches to NAMES.
        Where NAMES is either a
        
          1. string, in which case the experiments which have names containing
             this string are found.
          2. Sequence, in which case all experiments with names containing 
             any of the strings in NAMES are found. 

        @return: a sequence if indexes
        """
        if type(names)==type(' '):
            idx=[i for i in range(self.nExps)
                 if self.expNames[i] == names]
        else:
            idx=[i for name in names
                 for i in range(self.nExps)
                 if self.expNames[i] == name]
        return idx

    def getExps(self, names):
        """Searches through all expNames to find matches to NAMES.
        Where NAMES is either a 
          1. string, in which case the experiments which have names containing
             this string are found.
          2. Sequence, in which case all experiments with names containing 
             any of the strings in NAMES are found. 

        @return: a new MicroArray object.
        """
        idx = self.getExpIdx(names)
        return MicroArray(take(self.data,idx, 1), take(self.expNames, idx),
                          self.geneNames, take(self.expClass, idx))

    def getGeneIdx(self, names):
        """Searches through all geneNames to find matches to NAMES.
        Where NAMES is either a 
          1. string, in which case the genes which have names containing
             this string are found.
          2. Sequence, in which case all genes with names containing 
             any of the strings in NAMES are found. 

        @return: a sequence if indexes
        """
        if type(names)==type(' '):
            idx=[i for i in range(self.nGenes)
                 if self.geneNames[i]==names]
        else:
            idx=[i for name in names
                 for i in range(self.nGenes)
                 if self.geneNames[i]==name]
        return idx

    def getGenes(self, names):
        """Searches through all geneNames to find matches to NAMES.
        Where NAMES is either a 
          1. string, in which case the experiments which have names containing
             this string are found.
          2. Sequence, in which case all experiments with names containing 
             any of the strings in NAMES are found. 

        @return: a new MicroArray object.
        """
        idx = self.getGeneIdx(names)
        return MicroArray(take(self.data,idx, 0), self.expNames,
                          take(self.geneNames, idx), self.expClass)
            
    def __str__(self):
        str = 'ORF\t';
        str+= ('%s\t'*(self.nExps-1)+'%s') \
	      % tuple(self.expNames[0:self.nExps])
        str+='\n'
        for i in xrange(self.nGenes):
            str += '%s\t' % self.geneNames[i]
            str += ('%8.4e\t'*(self.nExps-1)+'%8.4e') \
                    % tuple(self.data[i, 0:self.nExps])
            str+='\n'
        return str


    def __getitem__(self, k):
        ###FIX so that ma[1] and ma[1,2] works.  Now self.expNames[1]
        #cannot create array object in constructor
        
        import types
#        print k, self.expClass.shape
        if type(k) == types.TupleType:
            return MicroArray(self.data[k], self.expNames[k[1]],
                              self.geneNames[k[0]], self.expClass[k[1],...])
        elif isinstance(k, types.IntType):
            return MicroArray(self.data[k], self.expNames[:],
                              self.geneNames[k], self.expClass[:])
        else:
            print 'Should handle STRING!!!'
        
    def __setitem__(self, i, value):
        self.data[i] = value

    def __getslice__(self, i,j):
        print 'The array has two Dimensions!!!'
        return None

    #def __setslice__(self, i, j, seq):
    #def __delslice__(self, i,j):
    #def __delitem__(self, i):

