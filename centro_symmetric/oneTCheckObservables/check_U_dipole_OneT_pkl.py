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

# print(f"U_dipole_dataDir={U_dipole_dataDir}")
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


def auto_corrForOneColumn(vec,eps):
    """

    :param vec: a vector of data
    :param eps: correlation truncation value
    :return: auto-correlation length
    """
    same=False
    # eps=5e-2
    NLags=int(len(vec)*1/4)
    # print("NLags="+str(NLags))
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(vec,nlags=NLags)
    except Warning as w:
        same=True
    acfOfVecAbs=np.abs(acfOfVec)
    minAutc=np.min(acfOfVecAbs)
    print(f"minAutc={minAutc}")

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



def check_U_data_files_for_one_T(UData_dir,summary_U_dipoleFile,lastFileNum,eps):

    U_sortedDataFilesToRead=sort_data_files_by_flushEnd(UData_dir)

    if len(U_sortedDataFilesToRead)==0:
        print("no data for U.")
        exit(0)

    startingFileInd=parseSummaryU_dipole(summary_U_dipoleFile)
    if startingFileInd<0:
        startingFileInd=len(U_sortedDataFilesToRead)-lastFileNum

    startingFileName=U_sortedDataFilesToRead[startingFileInd]

    with open(startingFileName,"rb") as fptr:
        in_UArrStart=pickle.load(fptr)

    arr_U=in_UArrStart

    #read the rest of the pkl files
    for pkl_file in U_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            in_UArr=pickle.load(fptr)

        arr_U=np.append(arr_U,in_UArr)

    sameUTmp,lagUTmp=auto_corrForOneColumn(arr_U,eps)
    #if one lag==-1, then the auto-correlation is too large

    if sameUTmp==True or lagUTmp==-1:
        return [sameUTmp,lagUTmp,-1,-1,-1,-1]

    pUTmp,statUTmp,lengthUTmp=ksTestOneVec(arr_U,lagUTmp)

    numDataPoints=lengthUTmp

    return [sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPoints,startingFileInd]


def check_dipole_files_for_one_T(theta_data_dir,summary_U_dipoleFile,lastFileNum,eps):


    dipole_sortedDataFilesToRead=sort_data_files_by_flushEnd(theta_data_dir)

    startingFileInd=parseSummaryU_dipole(summary_U_dipoleFile)

    if startingFileInd<0:
        startingFileInd=len(dipole_sortedDataFilesToRead)-lastFileNum

    dipole_starting_file_name=dipole_sortedDataFilesToRead[startingFileInd]
    with open(dipole_starting_file_name,"rb") as fptr:
        theta_inArrStart=pickle.load(fptr)

    theta_arr=theta_inArrStart
    #read the rest of the theta pkl files
    for pkl_file in dipole_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            theta_inArr=pickle.load(fptr)
            theta_arr=np.append(theta_arr,theta_inArr)

    theta_arr=theta_arr.reshape((-1,N**2))
    cos_theta_arr=np.cos(theta_arr)
    sin_theta_arr=np.sin(theta_arr)

    px_mean_all=np.mean(cos_theta_arr,axis=1)
    py_mean_all=np.mean(sin_theta_arr,axis=1)
    # print(px_mean_all)

    same_px_tmp,lag_px_tmp=auto_corrForOneColumn(px_mean_all,eps)

    #large correlation for px
    if same_px_tmp==True or lag_px_tmp==-1:
        return [same_px_tmp,lag_px_tmp,-1,-1,-1,-1]

    same_py_tmp,lag_py_tmp=auto_corrForOneColumn(py_mean_all,eps)

    #large correlation for py
    if same_py_tmp==True or lag_py_tmp==-1:
        return [same_py_tmp,lag_py_tmp,-1,-1,-1,-1]

    lag_pxpy=np.max((lag_px_tmp,lag_py_tmp))
    p_px_tmp,stat_px_tmp,length_px_tmp=ksTestOneVec(px_mean_all,lag_pxpy)

    p_py_tmp,stat_py_tmp,length_py_tmp=ksTestOneVec(py_mean_all,lag_pxpy)

    numDataPoints=length_py_tmp

    return [same_px_tmp,lag_pxpy,p_px_tmp,p_py_tmp,stat_px_tmp,stat_py_tmp,numDataPoints,startingFileInd]





sameVec=[]
lagVec=[]
pVec=[]
statVec=[]
numDataVec=[]

lastFileNum=1
eps=5e-2
sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPoints_U,startingFileInd=check_U_data_files_for_one_T(U_data_dir,summary_U_dipoleFile,lastFileNum, eps)
print(f"lagUTmp={lagUTmp}")
sameVec.append(sameUTmp)
lagVec.append(lagUTmp)
pVec.append(pUTmp)
statVec.append(statUTmp)
numDataVec.append(numDataPoints_U)
item=check_dipole_files_for_one_T(theta_data_dir,summary_U_dipoleFile,lastFileNum,eps)
print(f"item={item}")
#check if lag==-1
if item[-1]==-1 or lagUTmp==-1:
    msg="high correlation"
    with open(summary_U_dipoleFile,"w+") as fptr:
        fptr.writelines(msg)
        exit(0)

same_exist=any(sameVec)
if same_exist==True:
    with open(summary_U_dipoleFile,"w+") as fptr:
        msg="error: same\n"
        fptr.writelines(msg)
        exit(sameErrCode)


same_px_tmp,lag_pxpy,p_px_tmp,p_py_tmp,stat_px_tmp,stat_py_tmp,numDataPoints_pxpy,startingFileInd=item

lagVec.append(lag_pxpy)
pVec.append(p_px_tmp)
pVec.append(p_py_tmp)

statVec.append(stat_px_tmp)
statVec.append(stat_py_tmp)

numDataVec.append(numDataPoints_pxpy)


if np.min(lagVec)<0:
    msg="high correlation"
    with open(summary_U_dipoleFile,"w+") as fptr:
        fptr.writelines(msg)
    exit(0)

same_exist=any(sameVec)
if same_exist==True:
    with open(summary_U_dipoleFile,"w+") as fptr:
        msg="error: same\n"
        fptr.writelines(msg)
        exit(sameErrCode)


lagMax=np.max(lagVec)
statThreshhold=0.1
print("statVec="+str(statVec))
print("pVec="+str(pVec))
numDataPoints=np.min(numDataVec)

if (np.max(statVec)<=statThreshhold or np.min(pVec)>=0.01) and numDataPoints>=200:
    if numDataPoints>=effective_data_num_required:
        newDataPointNum=0
    else:
        newDataPointNum=effective_data_num_required-numDataPoints
    msg="equilibrium\n" \
        +"lag="+str(lagMax)+"\n" \
        +"numDataPoints="+str(numDataPoints)+"\n" \
        +"startingFileInd="+str(startingFileInd)+"\n" \
        +"newDataPointNum="+str(newDataPointNum)+"\n"
    print(msg)
    with open(summary_U_dipoleFile,"w+") as fptr:
        fptr.writelines(msg)
    exit(0)


#continue
continueMsg="continue\n"
if np.max(statVec)>statThreshhold:
    #not the same distribution
    continueMsg+="stat value: "+str(np.max(statVec))+"\n"


if numDataPoints<200:
    #not enough data number
    continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
    continueMsg+="lag="+str(lagMax)+"\n"

print(continueMsg)
with open(summary_U_dipoleFile,"w+") as fptr:
    fptr.writelines(continueMsg)
exit(0)