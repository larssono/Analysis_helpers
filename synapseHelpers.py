import synapseclient
import os
import inspect
import sys
import numpy as np

SYNAPSE_PROPERTIES  = ['benefactorId', 'nodeType', 'concreteType', 'createdByPrincipalId', 
                       'createdOn', 'createdByPrincipalId', 'eTag', 'id', 'modifiedOn', 
                       'modifiedByPrincipalId', 'noteType', 'versionLabel', 'versionComment', 
                       'versionNumber', 'parentId', 'description', 's3Token', 'name', 
                       'alias', 'projectId', 'fileNameOverride']

syn=synapseclient.Synapse(skip_checks=True)
syn.login(silent=True)

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


def query2df(queryContent, filterSynapseFields=True, savedSynapseFields = ('id', 'name')):
    """Converts the returned query object from Synapse into a Pandas DataFrame
    
    Arguments:
    - `queryContent`: content returned from query or chunkedQuery
    - `filterSynapseFields`: Removes Synapse properties of the entity returned with "select * ..."
                             defaults to True
    - `savedSynapseFields`: list of synapse fields not to filter defaults to id, name
    """
    import pandas as pd
    removedProperties = set(SYNAPSE_PROPERTIES) - set(savedSynapseFields)
    try:
        queryContent = queryContent['results']
    except TypeError: #If 'results' not present it is a geneartor
        queryContent = list(queryContent)
    for row in queryContent:
        for key in row.keys():
            prefix, new_key = key.split('.', 1)
            item = row.pop(key)
            row[new_key] = item[0] if type(item) is list else item
            if filterSynapseFields and new_key in removedProperties:
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


def match(seq1, seq2):
    """Finds the index locations of seq1 in seq2"""
    return [ np.nonzero(seq2==x)[0][0] for x in seq1  if x in seq2 ]


def cleanUpUUIDprojects():
    """
    """
    uuidPattern = re.compile('/[a-f0-9]{8}-[a-f0-9]{4}-4[a-f0-9]{3}-[89aAbB][a-f0-9]{3}-[a-f0-9]{12}/')
    #for each project see if the name matches the 
    myProjects = syn.chunkedQuery('select name from project where createdByPrincipalId=="%s"' %syn.getUserProfile['userId'] )



###Convenience functions for for tables
def tableDeleteWhere(query):
    """The DELETE statement is used to delete rows in a table.

    :param query      A string with a query to delete rows.
    """    
    syn.delete(syn.tableQuery(query).asRowSet())
        


def tableUpdateWhere(tableSchema, whereClause, setDict):
    """ The UPDATE statement is used to update existing rows in a table.
    """
    from synapseclient.table import Table
    from synapseclient.utils import id_of
    import tempfile
    id = id_of(tableSchema)
    query = 'select %s from %s where %s' % (','.join(setDict.keys()), id, whereClause)
    df = syn.tableQuery(query).asDataFrame()
    for key, value in setDict.items():
        df[key] = value
    print(df)
    # df.to_csv('skit.csv')
    return syn.store(Table(id_of(tableSchema), 'skit.csv'))


if __name__ == '__main__':
    thisCodeInSynapse('syn537704')
    print('done')
