import pickle
import numpy as np
from datetime import datetime
from multiprocessing import Pool
import pandas as pd
import statsmodels.api as sm
import sys
import re
import warnings
from scipy.stats import ks_2samp
import glob
from pathlib import Path
import os
import json


#This script checks if U and average dipole values reach equilibrium and writes summary file

#This script checks pkl files


argErrCode=2
sameErrCode=3
missingErrCode=4

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit(argErrCode)


jsonFromSummaryLast=json.loads(sys.argv[1])
jsonDataFromConf=json.loads(sys.argv[2])


TDirRoot=jsonFromSummaryLast["TDirRoot"]

U_dipole_dataDir=jsonFromSummaryLast["U_dipole_dataDir"]
effective_data_num_required=int(jsonDataFromConf["effective_data_num_required"])
N=int(jsonDataFromConf["N"])

summary_U_dipoleFile=TDirRoot+"/summary_U_dipole.txt"

lastFileNum=10
print(f"U_dipole_dataDir={U_dipole_dataDir}")
theta_data_dir=U_dipole_dataDir+"/theta/"
U_data_dir=U_dipole_dataDir+"/U/"

def sort_data_files_by_flushEnd(oneDir):
    """

    :param oneDir: dir containing pkl data
    :return: sorted file names, by the value of flushEnd
    """
    dataFilesAll=[]
    flushEndAll=[]
    # print("entering sort")
    for oneDataFile in glob.glob(oneDir+"/flushEnd*.pkl"):
        # print(oneDataFile)
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"flushEnd(\d+)",oneDataFile)
        if matchEnd:
            indTmp=int(matchEnd.group(1))
            flushEndAll.append(indTmp)
    endInds=np.argsort(flushEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]
    return sortedDataFiles


def parseSummaryU_dipole(summary_U_dipoleFile):
    """

    :param summary_U_dipoleFile: summary file to be parsed
    :return: startingFileInd
    """
    startingFileInd=-1
    summaryFileExists=os.path.isfile(summary_U_dipoleFile)
    if summaryFileExists==False:
        return startingFileInd
    with open(summary_U_dipoleFile,"r") as fptr:
        lines=fptr.readlines()
    for oneLine in lines:
        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))

    return startingFileInd


def auto_corrForOneColumn(colVec,eps):
    """

    :param colVec: a vector of data
    :param eps: correlation truncation value
    :return: auto-correlation length
    """
    same=False
    # eps=5e-2
    NLags=int(len(colVec)*1/4)
    # print("NLags="+str(NLags))
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(colVec,nlags=NLags)
    except Warning as w:
        same=True
    acfOfVecAbs=np.abs(acfOfVec)
    minAutc=np.min(acfOfVecAbs)

    lagVal=-1
    if minAutc<=eps:
        lagVal=np.where(acfOfVecAbs<=eps)[0][0]
    # np.savetxt("autc.txt",acfOfVecAbs[lagVal:],delimiter=',')
    return same,lagVal


def ksTestOneVec(vec,lag):
    """

    :param vec: a vector of data
    :param lag: auto-correlation length
    :return:
    """
    vecSelected=vec[::lag]

    lengthTmp=len(vecSelected)
    if lengthTmp%2==1:
        lengthTmp-=1
    lenPart=int(lengthTmp/2)

    vecToCompute=vecSelected[-lengthTmp:]

    #ks test
    selectedVecPart0=vecToCompute[:lenPart]
    selectedVecPart1=vecToCompute[lenPart:]
    result=ks_2samp(selectedVecPart0,selectedVecPart1)
    return result.pvalue,result.statistic, lenPart*2