import synapseclient
import os
import inspect
import sys
import fileReader
import numpy as np

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
    code= synapseclient.File(name=file, parentId=parentId, description=description)
    codeEntity = syn.store(code)
    return codeEntity



def createEvaluationBoard(name, parentId, markdown='',status='OPEN', contentSource=''):
        """Creates an evaluation where users can submit entities for evaluation by a monitoring process.
        Currently also creates a folder where the results are displayed in a leaderboard.
        
        Arguments:
        - `name`: Name of the Evaluation
        - `parentId` : Location of leaderboard folder
        - `description`: A string describing the evaluation in detail
        - `status`: A string describing the status: one of { PLANNED, OPEN, CLOSED, COMPLETED}

        Returns:
        - `evaluation`: information about the evaluation 
        - `leaderboard`: folder that contains the output leaderboard
        """
        # # Create an evaluation
        evaluation=syn.store(Evaluation(name=name, status=status, contentSource=contentSource))
        
        # # create a wiki to describe the Challenge
        homeWikiPage = Wiki(title=name, markdown=md, owner=Evaluation)
        homeWikiPage = syn.store(homeWikiPage)

        return evaluation, homeWikiPage
        

def updateLeaderboard(leaderboard, evaluation):
    """Goes through all submissions in an evalution and updates the description markdown in entity.
    
    Arguments:
    - `entity`: A folder/project entity that contains the leaderboard
    - `evaluation`: an evaluation object where submissions are made
    """
    leaderboard.markdown = 'Submission Name|Submitter|Submission|Score|Status|Report|\n'
    userNames={}
    for submission in  syn.getSubmissions(evaluation):
        status =  syn.getSubmissionStatus(submission)
        #Extract the username
        userName = userNames.get(submission['userId'], None)
        if userName==None:
            userName = syn.getUserProfile(submission['userId'])['displayName']
            userNames[submission['userId']] = userName
        print submission.get('name', ''), userName, submission.entityId, status.score, status.status, status.report
        leaderboard.markdown+='%s|%s|%s|%f|%s|%s|\n' %(submission.get('name', ''), 
                                                 userName, submission.entityId, status.score,
                                                 status.status, status.report)
    syn.store(leaderboard)



def readEntityFile(entity):
    return fileReader.csv2Dict(os.path.join(entity['cacheDir'], entity['files'][0]))


def match(seq1, seq2):
    """Finds the index locations of seq1 in seq2"""
    return [ np.nonzero(seq2==x)[0][0] for x in seq1  if x in seq2 ]
        

if __name__ == '__main__':
    thisCodeInSynapse('syn537704')
    print 'done'
