import os,sys

def reportTime(seconds):
    if seconds < 60:
        return str(round(seconds,2)) + " s"
    elif seconds < 3600:
        return str(round(seconds/60,2)) + " m"
    elif seconds < 86400:
        return str(round(seconds/60/60,2)) + " h"
    else:
        return str(round(seconds/60/60/24,2)) + " days"

logFile = os.path.realpath("testDetails.log")
if len(sys.argv) > 1:
    if os.path.isfile(sys.argv[1]):
        logFile = os.path.realpath(sys.argv[1])
    else:
        print "\nFile ", sys.argv[1], " not Found! \n"
        sys.exit(1)

fileH = open(logFile,"r")
lines = fileH.readlines()
fileH.close()

testTimes = []
testNames = []
os.system("clear")
print "%-30s%-16s%-16s" % ("Test Name","Test time","Total elapsed time")
print "-" * 64
looking4Time = False
for line in lines:
    if "+-+-+-+-+-+-+-+" in line:
        testNames.append((line.replace("+","").replace("-","")).strip())
        if looking4Time:
            testTimes.append(0.0)
        else:
            looking4Time = True
    splitedLine = line.strip().split()
    if len(splitedLine) == 5 and splitedLine[0] == "Ran" and "test" in splitedLine[2] and splitedLine[3] == "in" and splitedLine[4][-1] == "s":
        if looking4Time:
            testTimes.append(float(splitedLine[4].strip()[:-1]))
        else:
            print "%-30s%-16s%-16s" % ("Found unrelated test time", line.strip()[line.index("in"):-1].strip()+" s","-")
        looking4Time = False
        
if looking4Time:
    testTimes.append(0.0)

orderedList = []
for idx in range(len(testNames)):
    orderedList.append((testTimes[idx],testNames[idx]))
orderedList.sort()

runningSum = 0
for iTime in orderedList:
    runningSum += iTime[0]
    print "%-30s%-16s%-16s" % (iTime[1], reportTime(iTime[0]), reportTime(runningSum))
print "-" * 64
print "Total tests running time: " + str(int(sum(testTimes))) + " s"
print "                          (" + str(round(sum(testTimes)/60,3)) + " m)"
print "                          ("+ str(round(sum(testTimes)/60/60,3)) + " h)"
        
print ""
