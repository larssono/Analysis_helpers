import gzip, numpy as np

def openfile(file):
    """Returns a file object by opening either a gz or uncompressed file. """
    if file.endswith('.gz'):
        return gzip.open(file)
    else:
        return open(file)


def csv2Dict(file, sep='\t'):
    """Given a tab delimited file returns the first column as list and
    remaining columns as dictionary of lists.
    Parameters:
    `file` input file
    `sep` - optional parameter incase file is not delimited by tab
    """    
    fp=openfile(file)
    head=fp.readline().rstrip().split(sep)
    head=[l.strip('"') for l in head]
    head= dict(zip(head[1:], range(1, len(head)))) #translation between names and columns
    vals=[line.split(sep) for line in fp] 
    phenotypes={}
    subjects=[val[0].strip('"') for val in vals] 
    for key, i in head.items(): 
        phenotypes[key]=np.asarray([val[i].strip('"') for val in vals]) 
    return subjects, phenotypes

def csv2DictDict(file, sep='\t'):
    """Given a tab delimited file returns a dict of dicts where the
      first key references a column and the second key the row.

      Example:"  id       age    gender    
                 NML-001   24     M  
                 NML-002   48     F"
              
          fileDict=csv2DictDict(file)
          fileDict['NML-001']['age']  -> 24
          fileDict['NML-002']['gender'] -> F"""    
    fp=openfile(file)
    head=fp.readline().rstrip().split(sep)
    head=[l.strip('"') for l in head]
    head= dict(zip(head[1:], range(1, len(head)))) #translation between names and columns
    phenotypes={}
    for line in fp:
        line = line.split(sep) 
        subject=line[0].strip('"')
        phenotypes[subject]={}
        for key, i in head.items(): 
            phenotypes[subject][key]=line[i].strip('"')
    return phenotypes
