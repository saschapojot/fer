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
from decimal import Decimal, getcontext
import matplotlib.pyplot as plt

argErrCode=2
sameErrCode=3
missingErrCode=4

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit(argErrCode)

def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)


T=float(sys.argv[1])
TStr=format_using_decimal(T)

dataRoot=f"./dataAll/T{TStr}/U_dipole_dataFiles/"

U_data_dir=dataRoot+"/U/"



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


def concatenate_arrs(fileNames_vec):

    len_vec=len(fileNames_vec)

    with open(fileNames_vec[0],"rb") as fptr:
        arr=pickle.load(fptr)

    for j in range(1,len_vec):
        with open(fileNames_vec[j],"rb") as fptr:
            arr_in=pickle.load(fptr)
            arr=np.append(arr,arr_in)
    return arr


sorted_U_pkl_files=sort_data_files_by_flushEnd(U_data_dir)


arr_U=concatenate_arrs(sorted_U_pkl_files)

plt.figure()
plt.plot(arr_U,color="red")
plt.title("$U$")

plt.savefig("UTmp.png")
plt.close()