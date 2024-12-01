import re
import sys

import json
import os
#this script parse conf file and return the parameters as json data
fmtErrStr="format error: "
fmtCode=1
valueMissingCode=2
paramErrCode=3
fileNotExistErrCode=4
if (len(sys.argv)!=2):
    print("wrong number of arguments.")
    exit(paramErrCode)
inConfFile=sys.argv[1]

def removeCommentsAndEmptyLines(file):
    """

    :param file: conf file
    :return: contents in file, with empty lines and comments removed
    """
    with open(file,"r") as fptr:
        lines= fptr.readlines()

    linesToReturn=[]
    for oneLine in lines:
        oneLine = re.sub(r'#.*$', '', oneLine).strip()
        if not oneLine:
            continue
        else:
            linesToReturn.append(oneLine)
    return linesToReturn
def parseConfContents(file):
    """

    :param file: conf file
    :return:
    """
    file_exists = os.path.exists(file)
    if not file_exists:
        print(file+" does not exist,")
        exit(fileNotExistErrCode)

    linesWithCommentsRemoved=removeCommentsAndEmptyLines(file)
    # print(linesWithCommentsRemoved)
    TStr=""
    JStr=""
    aStr=""
    NStr=""
    eraseData=""
    searchReadSmrFile=""
    obs_name=""
    confFileName=file
    effective_data_num_required=""
    sweep_to_write=""
    default_flush_num=""
    hStr=""
    swp_multiplyStr=""
    # float_pattern = r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?'
    boolean_pattern = r'(true|false)'

    # coefs_pattern = rf'\[({float_pattern}(?:\s*,\s*{float_pattern})*)\]'
    # print(linesWithCommentsRemoved)
    for oneLine in linesWithCommentsRemoved:
        matchLine=re.match(r'(\w+)\s*=\s*(.+)', oneLine)
        if matchLine:
            key = matchLine.group(1).strip()
            value = matchLine.group(2).strip()

            #match T
            if key=="T":
                match_TValPattern=re.match(r"T\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$",oneLine)
                if match_TValPattern:
                    TStr=match_TValPattern.group(1)
                else:
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)

            #match a
            if key=="a":
                match_a_Pattern=re.match(r"a\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$",oneLine)
                if match_a_Pattern:
                    aStr=match_a_Pattern.group(1)
                else:
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)
            # match J
            if key=="J":
                match_J_Pattern=re.match(r"J\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)$",oneLine)
                if match_J_Pattern:
                    JStr=match_J_Pattern.group(1)
                else:
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)

            # match N
            if key=="N":
                match_N_pattern=re.match("N\s*=\s*([1-9]\d*)",oneLine)
                if match_N_pattern:
                    NStr=match_N_pattern.group(1)
                else:
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)

            #match erase_data_if_exist
            if key=="erase_data_if_exist":
                matchErase=re.match(boolean_pattern,value,re.IGNORECASE)
                if matchErase:
                    eraseData=matchErase.group(1)
                    eraseData=eraseData.capitalize()
                else:
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)

            #match search_and_read_summary_file
            if key=="search_and_read_summary_file":
                matchSmr=re.match(boolean_pattern,value,re.IGNORECASE)
                if matchSmr:
                    searchReadSmrFile=matchSmr.group(1)
                    searchReadSmrFile=searchReadSmrFile.capitalize()
                else:
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)

            #match observable_name
            if key=="observable_name":
                #if matching a non word character
                if re.search(r"[^\w]",value):
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)

                obs_name=value




            #match sweep_to_write
            if key=="sweep_to_write":
                if re.search(r"[^\d]",value):
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)
                sweep_to_write=value

            #match default_flush_num
            if key=="default_flush_num":
                if re.search(r"[^\d]",value):
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)
                default_flush_num=value


            #match effective_data_num_required
            if key=="effective_data_num_required":
                if re.search(r"[^\d]",value):
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)
                effective_data_num_required=value


            #match h
            if key=="h":
                match_h=re.match(r'([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',value)
                # print(value)
                if match_h:
                    hStr=match_h.group(1)
                else:
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)
            #match sweep_multiply
            if key=="sweep_multiple":
                match_swpMultiply=re.match(r"(\d+)",value)
                if match_swpMultiply:
                    swp_multiplyStr=match_swpMultiply.group(1)
                else:
                    print(fmtErrStr+oneLine)
                    exit(fmtCode)



        else:
            print("line: "+oneLine+" is discarded.")
            continue

    if TStr=="":
        print("T not found in "+str(file))
        exit(valueMissingCode)
    if aStr=="":
        print("a not found in "+str(file))
        exit(valueMissingCode)
    if JStr=="":
        print("J not found in "+str(file))
        exit(valueMissingCode)
    if NStr=="":
        print("J not found in "+str(file))
        exit(valueMissingCode)
    if eraseData=="":
        print("erase_data_if_exist not found in "+str(file))
        exit(valueMissingCode)
    if searchReadSmrFile=="":
        print("search_and_read_summary_file not found in "+str(file))
        exit(valueMissingCode)

    #do not check if observable exists

    if effective_data_num_required=="":
        print("effective_data_num_required not found in "+str(file))
        exit(valueMissingCode)

    if sweep_to_write=="":
        print("sweep_to_write not found in "+str(file))
        exit(valueMissingCode)

    if default_flush_num=="":
        print("default_flush_num not found in "+str(file))
        exit(valueMissingCode)

    if hStr=="":
        print("h not found in "+str(file))
        exit(valueMissingCode)
    if obs_name=="":
        print("observable_name not found in "+str(file))
        exit(valueMissingCode)

    if swp_multiplyStr=="":
        swp_multiplyStr="1"




    dictTmp={
            "T":TStr,
            "a":aStr,
            "J":JStr,
            "erase_data_if_exist":eraseData,
            "search_and_read_summary_file":searchReadSmrFile,
            "observable_name":obs_name,
            "effective_data_num_required":effective_data_num_required,
            "sweep_to_write":sweep_to_write,
            "default_flush_num":default_flush_num,
            "confFileName":file,
            "h":hStr,
            "sweep_multiple":swp_multiplyStr,



        }
    return dictTmp






jsonDataFromConf=parseConfContents(inConfFile)

confJsonStr2stdout="jsonDataFromConf="+json.dumps(jsonDataFromConf)

print(confJsonStr2stdout)