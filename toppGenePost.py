import requests

class ToppGeneEnrichement(object):
    """Given a gene list submits it to toppgene for enrichment analysis.  Useful when coupled with ipython notebooks to display inline 
    """
    
    URL = "http://toppgene.cchmc.org/CheckInput.action"
    LABEL_TYPES = set(['HGNC', "HGNC_SYNONYMS", "ENTREZ", "ENSEMBL", "REFSEQ", "UNIPROT", "GENOME_BROWSER"])


    def __init__(self, genes, labelType='HGNC', correction=None, name='View in ToppFun'):
        """Creates a 
        """
        self.genes = genes
        self.labelType = labelType
        self.correction = correction
        self.name = name


    def _repr_html_(self):
        #if self.correction==None:
        #    response =  requests.post(self.URL, {'type': type, 'training_set': training_set})
        #return response.text
        return """<form action="http://toppgene.cchmc.org/CheckInput.action" method="post" target="_blank">
        <input type="hidden" name="query" value="TOPPFUN">
        <input type="hidden" id="type" name="type" value="%s">
        <input type="hidden" name="training_set" id="training_set" value="%s">
        <input type="submit" value="%s">
        </form>""" %(self.labelType, ' '.join(self.genes), self.name)

