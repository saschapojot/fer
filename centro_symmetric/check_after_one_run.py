import re
import subprocess
import sys

import json
argErrCode=2

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    print("example: python launch_one_run.py /path/to/mc.conf")
    exit(argErrCode)


confFileName=str(sys.argv[1])
# print("confFileName is "+confFileName)
invalidValueErrCode=1
summaryErrCode=2
loadErrCode=3
confErrCode=4

#################################################
#parse conf, get jsonDataFromConf
confResult=subprocess.run(["python3", "./init_run_scripts/parseConf.py", confFileName], capture_output=True, text=True)
confJsonStr2stdout=confResult.stdout
# print(confJsonStr2stdout)
if confResult.returncode !=0:
    print("Error running parseConf.py with code "+str(confResult.returncode))
    # print(confResult.stderr)
    exit(confErrCode)


match_confJson=re.match(r"jsonDataFromConf=(.+)$",confJsonStr2stdout)
if match_confJson:
    jsonDataFromConf=json.loads(match_confJson.group(1))
else:
    print("jsonDataFromConf missing.")
    exit(confErrCode)
# print(jsonDataFromConf)
################################################

##################################################
#read summary file, get jsonFromSummary
parseSummaryResult=subprocess.run(["python3","./init_run_scripts/search_and_read_summary.py", json.dumps(jsonDataFromConf)],capture_output=True, text=True)
# print(parseSummaryResult.stdout)
if parseSummaryResult.returncode!=0:
    print("Error in parsing summary with code "+str(parseSummaryResult.returncode))
    # print(parseSummaryResult.stdout)
    # print(parseSummaryResult.stderr)
    exit(summaryErrCode)

match_summaryJson=re.match(r"jsonFromSummary=(.+)$",parseSummaryResult.stdout)
if match_summaryJson:
    jsonFromSummary=json.loads(match_summaryJson.group(1))
# print(jsonFromSummary)
##################################################


##########################################################
#statistics
checkU_dipoleErrCode = 5
# print("entering statistics")
# Start the subprocess
# print("jsonFromSummary="+json.dumps(jsonFromSummary))
# print("jsonDataFromConf="+json.dumps(jsonDataFromConf))
checkU_dipoleProcess = subprocess.Popen(
    ["python3", "-u", "./oneTCheckObservables/check_U_dipole_OneT_pkl.py",
     json.dumps(jsonFromSummary), json.dumps(jsonDataFromConf)],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
)
# return_code = checkU_dipoleProcess.poll()
# if return_code is not None:
#     print(f"Process exited immediately with return code: {return_code}")

# Read output in real-time
while True:
    output = checkU_dipoleProcess.stdout.readline()
    if output == '' and checkU_dipoleProcess.poll() is not None:
        break
    if output:
        print(output.strip())

# Collect remaining output and error messages
stdout, stderr = checkU_dipoleProcess.communicate()
# Check if the process was killed
if checkU_dipoleProcess.returncode is not None:
    if checkU_dipoleProcess.returncode < 0:
        # Process was killed by a signal
        print(f"checkU_dipoleProcess was killed by signal: {-checkU_dipoleProcess.returncode}")
    else:
        # Process exited normally
        print(f"checkU_dipoleProcess exited with return code: {checkU_dipoleProcess.returncode}")
else:
    print("checkU_dipoleProcess is still running")
# Print any remaining standard output
if stdout:
    print(stdout.strip())

# Handle errors and print the return code if there was an error
if stderr:
    print(f"checkU_dipoleProcess return code={checkU_dipoleProcess.returncode}")
    print(stderr.strip())

##########################################################