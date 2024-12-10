import numpy as np
from datetime import datetime
import sys
import re
import glob
import os
import json
from pathlib import Path
import pandas as pd
import pickle
from decimal import Decimal, getcontext

#this script extracts effective data from pkl files
def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

T=float(sys.argv[1])
TStr=format_using_decimal(T)
print("TStr="+TStr)

N=int(sys.argv[2])
NStr=format_using_decimal(N)


dataRoot=f"../dataAll/N{NStr}/T{TStr}/"

U_dataDir=dataRoot+"/U_dipole_dataFiles/U/"

theta_dataDir=dataRoot+"/U_dipole_dataFiles/theta/"

def parseSummary(dataRoot):
    """

    :param dataRoot:
    :return: startingFileInd,lag
    """
    startingFileInd=-1
    lag=-1
    smrFile=dataRoot+"/summary_U_dipole.txt"


    summaryFileExists=os.path.isfile(smrFile)
    if summaryFileExists==False:
        return startingFileInd,lag
    with open(smrFile,"r") as fptr:
        lines=fptr.readlines()

    for oneLine in lines:
        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))
        #match lag
        matchLag=re.search(r"lag=(\d+)",oneLine)
        if matchLag:
            lag=int(matchLag.group(1))
    return startingFileInd,lag

def sort_pkl_data_files_by_flushEnd(pkl_folder):
    dataFilesAll=[]
    flushEndAll=[]
    for oneDataFile in glob.glob(pkl_folder+"/flushEnd*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"flushEnd(\d+)",oneDataFile)
        if matchEnd:
            flushEndAll.append(int(matchEnd.group(1)))


    endInds=np.argsort(flushEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]

    return sortedDataFiles


def U_pkl_data_2_array(U_sortedDataFiles,startingFileInd,lag):

    U_sorted_pkl_fileds_to_read=U_sortedDataFiles[startingFileInd:]

    with open(U_sorted_pkl_fileds_to_read[0],"rb") as fptr:
        in_U_arr=pickle.load(fptr)

    U_arr=in_U_arr

    for j in range(1,len(U_sorted_pkl_fileds_to_read)):
        with open(U_sorted_pkl_fileds_to_read[j],"rb") as fptr:
            UArrTmp=pickle.load(fptr)
            U_arr=np.append(U_arr,UArrTmp)

    U_selected_vec=U_arr[::lag]

    return U_selected_vec


def theta_pkl_2_array(theta_dataDir,startingFileInd,lag):
    theta_sorted_pkl_files_to_read=theta_dataDir[startingFileInd:]

    with open(theta_sorted_pkl_files_to_read[0],"rb") as fptr:
        in_theta_arr=pickle.load(fptr)

    theta_arr=in_theta_arr

    for j in range(1,len(theta_sorted_pkl_files_to_read)):
        with open(theta_sorted_pkl_files_to_read[j],"rb") as fptr:
            thetaArrTmp=pickle.load(fptr)
            theta_arr=np.append(theta_arr,thetaArrTmp)

    theta_arr=theta_arr.reshape((-1,N**2))
    selected_theta_arr=theta_arr[::lag,:]
    return selected_theta_arr


startingFileInd,lag=parseSummary(dataRoot)

U_sortedDataFiles=sort_pkl_data_files_by_flushEnd(U_dataDir)

theta_sortedDataFiles=sort_pkl_data_files_by_flushEnd(theta_dataDir)
# print(U_sortedDataFiles)
# print(theta_sortedDataFiles)
U_selected_vec=U_pkl_data_2_array(U_sortedDataFiles,startingFileInd,lag)

selected_theta_arr=theta_pkl_2_array(theta_sortedDataFiles,startingFileInd,lag)


outCsv_data=np.c_[U_selected_vec,selected_theta_arr]

outCsvDataRoot="../csvOutAll/"

out_csv_folder=outCsvDataRoot+f"/N{NStr}/T{TStr}/"
Path(out_csv_folder).mkdir(exist_ok=True,parents=True)
out_csv_file=out_csv_folder+"/U_theta.csv"
dfToSave=pd.DataFrame(outCsv_data)

dfToSave.to_csv(out_csv_file,index=False)