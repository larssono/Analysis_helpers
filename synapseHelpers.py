import synapseclient
import os
import inspect
import sys
import fileReader
import numpy as np

SYNAPSE_PROPERTIES  = ['benefactorId', 'nodeType', 'concreteType', 'createdByPrincipalId', 
                       'createdOn', 'createdByPrincipalId', 'eTag', 'id', 'modifiedOn', 
                       'modifiedByPrincipalId', 'noteType', 'versionLabel', 'versionComment', 
                       'versionNumber', 'parentId', 'description', 's3Token']

syn=synapseclient.Synapse()
syn.login()

def thisCodeInSynapse(parentId, file=None, description=''):
    """Determines the name of the file that the code is called from
    and uploads that to Synapse returning the synapseId of created codeObject.
    """
    #print inspect.getfile(inspect.currentframe())
    #print os.path.abspath(inspect.getfile(inspect.currentframe()))
    file = inspect.getfile(sys._getframe(1)) if file==None else file
    #Make sure unallowed characters are striped out for the name
    code= synapseclient.File(file, name=os.path.split(file)[-1], parent=parentId, description=description)
    codeEntity = syn.store(code)
    return codeEntity

def query2df(queryContent, filterSynapseFields=True):
    """Converts the returned query object from Synapse into a Pandas DataFrame
    
    Arguments:
    - `queryContent`: content returned from query or chunkedQuery
    - `filterSynapseFields`: Removes Synapse properties of the entity returned with "select * ..."
                             defaults to True
    """
    import pandas as pd
    try:
        queryContent = queryContent['results']
    except TypeError: #If 'results' not present it is a geneartor
        queryContent = list(queryContent)
    for row in queryContent:
        for key in row.keys():
            new_key = '.'.join(key.split('.')[1:])
            item = row.pop(key)
            row[new_key] = item[0] if type(item) is list else item
            if filterSynapseFields and new_key in SYNAPSE_PROPERTIES:
                del row[new_key]
    return pd.DataFrame(queryContent)



def df2markdown(df, wikiEntity=None, subPageId=None, syn=None, prefixText=None, suffixText=None):
    """Todo add documentation"""

    df = df.reset_index()
    if prefixText:
        wikiText = "%s\n\n" % prefixText
    else:
        wikiText = ''
    ncols = df.shape[1]
    nrows = df.shape[0]
    mod_colnames = map(lambda x: x.replace('_', '-'), df.columns.values)
    wikiText += "|%s|\n" %  ('|'.join(mod_colnames))
    wikiText += "|%s|\n" %  ( '|'.join(['--'] * ncols))

    for row in df.iterrows():
        values = row[1].values
        wikiText += "|%s|\n" % ('|'.join(map(str,values)))
    if suffixText:
        wikiText += "%s\n" % suffixText

    #just return the text
    if wikiEntity is None and syn is None:
        return wikiText
    else:
        wiki = syn.getWiki(wikiEntity, subpageId=subPageId)
        wiki['markdown'] = wikiText
        syn.store(wiki)
        return wikiText




def readEntityFile(entity):
    return fileReader.csv2Dict(os.path.join(entity['cacheDir'], entity['files'][0]))


def match(seq1, seq2):
    """Finds the index locations of seq1 in seq2"""
    return [ np.nonzero(seq2==x)[0][0] for x in seq1  if x in seq2 ]


def cleanUpUUIDprojects():
    """
    """
    uuidPattern = re.compile('/[a-f0-9]{8}-[a-f0-9]{4}-4[a-f0-9]{3}-[89aAbB][a-f0-9]{3}-[a-f0-9]{12}/')
    #for each project see if the name matches the 
    myProjects = syn.chunkedQuery('select name from project where createdByPrincipalId=="%s" %syn. )

if __name__ == '__main__':
    thisCodeInSynapse('syn537704')
    print 'done'


