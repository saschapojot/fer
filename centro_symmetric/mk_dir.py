from pathlib import Path
from decimal import Decimal, getcontext
from traceback import format_stack

import numpy as np
import pandas as pd


#This script creates directories and conf files for mc




def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)




TVals=[1,3,5]#unit is K

dataRoot="./dataAll/"
dataOutDir=dataRoot

TDirsAll=[]
TStrAll=[]
J=2
a=1
N=5 #unit cell number

JStr=format_using_decimal(J)
aStr=format_using_decimal(a)
NStr=format_using_decimal(N)
# print(TDirsAll)
for k in range(0,len(TVals)):
    T=TVals[k]
    # print(T)

    TStr=format_using_decimal(T)
    TStrAll.append(TStr)



# print(TStrAll)
#
# print(PStrAll)



def contents_to_conf(k):
    """

    :param k: index of T

    :return:
    """

    contents=[
        "#This is the configuration file for mc computations\n",
        "#System is centro-symmetric\n",
        "\n" ,
        "#parameters of coefficients\n",
        "\n",
        "#Temperature\n",
        "T="+TStrAll[k]+"\n",
        "\n",
        "a="+aStr+"\n",
        "\n",
        "J="+JStr+"\n",
        "\n",
        "N="+NStr+"\n",
        "\n",
        "erase_data_if_exist=False\n",
        "\n",
        "search_and_read_summary_file=True\n"
        "\n",
        "#For the observable name, only digits 0-9, letters a-zA-Z, underscore _ are allowed\n",
        "\n",
        "observable_name=U_dipole\n",
        "\n",
        "effective_data_num_required=1000\n",
        "\n",
        "#this is the data number in each pkl file, i.e., in each flush\n"
        "sweep_to_write=2\n",
        "\n",
        "#within each flush,  sweep_to_write*sweep_multiple mc computations are executed\n",
        "\n",
        "default_flush_num=10\n",
        "\n",
        "h=0.1\n",
        "\n",
        "#the configurations of the system are saved to file if the sweep number is a multiple of sweep_multiple\n",
        "\n",
        "sweep_multiple=3\n",





    ]
    outDir=dataOutDir+"/T"+TStrAll[k]+"/"
    Path(outDir).mkdir(exist_ok=True,parents=True)
    outConfName=outDir+"/run_T"+TStrAll[k]+".mc.conf"
    with open(outConfName,"w+") as fptr:
        fptr.writelines(contents)



for k in range(0,len(TVals)):
    contents_to_conf(k)