import sys
import glob
import re
import json
from decimal import Decimal, getcontext
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import pickle
import random
#this script loads previous data
numArgErr=4
valErr=5
if (len(sys.argv)!=3):
    print("wrong number of arguments.")
    exit(numArgErr)

# print("entering load")
jsonDataFromConf =json.loads(sys.argv[1])
jsonFromSummary=json.loads(sys.argv[2])


U_dipole_dataDir=jsonFromSummary["U_dipole_dataDir"]

startingFileInd=jsonFromSummary["startingFileInd"]

NStr=jsonDataFromConf["N"]

N=int(NStr)
if N<=0:
    print("invalid N="+str(N))
    exit(valErr)


#search and read U_dipole files



#search flushEnd
pklFileList=[]
flushEndAll=[]
#read theta files
for file in glob.glob(U_dipole_dataDir+"/theta/flushEnd*.pkl"):
    pklFileList.append(file)
    matchEnd=re.search(r"flushEnd(\d+)",file)
    if matchEnd:
        flushEndAll.append(int(matchEnd.group(1)))
# print(U_dipole_dataDir)
flushLastFile=-1

def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)
#

def create_init_theta(U_dipole_dataDir):
    """
    create initial values for theta
    :return:
    """
    #

    # print("random_perturbations: "+str(random_perturbations))
    outPath=U_dipole_dataDir+"/coordinate/"

    Path(outPath).mkdir(exist_ok=True,parents=True)

    outFileName=outPath+"/theta_init.pkl"

    thetaMat=np.array([random.uniform(0,2*np.pi) for i in range(0,N*N)])
    with open(outFileName,"wb") as fptr:
        pickle.dump(thetaMat,fptr)





def create_loadedJsonData(flushLastFileVal):

    initDataDict={

        "flushLastFile":str(flushLastFileVal)
    }
    # print(initDataDict)
    return json.dumps(initDataDict)

#if no data found, return flush=-1
if len(pklFileList)==0:
    create_init_theta(U_dipole_dataDir)

    out_U_path=U_dipole_dataDir+"/U/"
    Path(out_U_path).mkdir(exist_ok=True,parents=True)
    loadedJsonDataStr=create_loadedJsonData(flushLastFile)
    loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
    print(loadedJsonData_stdout)
    exit(0)



#if found pkl data with flushEndxxxx
sortedEndInds=np.argsort(flushEndAll)
sortedflushEnd=[flushEndAll[ind] for ind in sortedEndInds]
loadedJsonDataStr=create_loadedJsonData(sortedflushEnd[-1])
loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
print(loadedJsonData_stdout)
exit(0)
