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
import matplotlib.pyplot as plt

#this script computes effective data in 1 csv file
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

csv_DataRoot="./csvOutAll/"

csv_file=csv_DataRoot+f"/N{NStr}/T{TStr}/U_theta.csv"

inArr=np.array(pd.read_csv(csv_file))


U_vec=inArr[:,0]

theta_all_arr=inArr[1:,:]

#each row of cos_theta_all_arr is N*N px array
cos_theta_all_arr=np.cos(theta_all_arr)

#each row of cos_theta_all_arr is N*N py array
sin_theta_all_arr=np.sin(theta_all_arr)

U_avg=U_vec.mean()
px_avg=cos_theta_all_arr.mean()
py_avg=sin_theta_all_arr.mean()

print(f"px_avg={px_avg}")

print(f"py_avg={py_avg}")
print(f"U_avg={U_avg}")

plt.figure()
plt.scatter(range(0,len(U_vec)),U_vec,color="red")
plt.savefig("UTmp.png")