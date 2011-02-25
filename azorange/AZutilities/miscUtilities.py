from datetime import datetime
from AZOrangeConfig import SCRATCHDIR
from time import strptime
import orange
import random
import types
import os,time
import commands
import fcntl
from xml.dom.minidom import parseString,Document,Element,Node,_write_data
import zipfile
from shutil import copy,copytree
from StringIO import StringIO
from glob import glob

version = 2
verbose = 0
##scPA


def getText(nodelist):
    rc = ""
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc = rc + node.data
    return rc

def balanceJobs():
    c=0
    (code,output)=commands.getstatusoutput("qstat -xml")
    if code == 0: 
        xml=parseString(output)
        jobs=xml.getElementsByTagName("job_list")
        for i in range(0,len(jobs)):
           state=getText(jobs[i].getElementsByTagName("state")[0].childNodes)
           jobname=getText(jobs[i].getElementsByTagName("JB_name")[0].childNodes)
           if "appspack_mpi" in jobname:
               if state == "r":
                   c=c+1
           if "runTest-" in jobname:
               c=c-1
        return c
    else:
        print "Qstat error code" + str(code)
        return -100

class LockError(IOError):
        """LockError"""

def lockFile(filename, mode='r', bufsize=-1, timeout=5, step=0.1):
        file = open(filename, mode, bufsize)
        locked=False
        while locked == False:
            try:
                fcntl.flock(file.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
                locked=True
            except IOError:
                time.sleep(step)
                timeout=timeout-step
                if timeout <= 0:
                        raise LockError
        return file 

def autoValidateRSAKey(machine, user=None):
    """ SSH to the destination machine and add it's fingerprint to the ~/.ssh/known_hosts
        It will exit imediatly from the ssg connection
        returns True is success.
    """
    tmpDir = createScratchDir(desc="SSHAutoValidateRSA")
    tmpScript=os.path.join(tmpDir,"validateRSAkey")
    if user:
        target = user+"@"+machine        
    else:
        target = machine

    tmpfile = open(tmpScript,"w")
    tmpfile.write("spawn ssh "+target+" exit\n")   
    tmpfile.write("expect *yes*no*\n")
    tmpfile.write('send "yes\\r"\n')
    tmpfile.write("expect eof\n")
    tmpfile.close()

    status,out = commands.getstatusoutput("expect "+tmpScript)
    if "Offending key in" in out:
        known_hostsFile = out[out.index("Offending key in")+16:out.find("\n",out.index("Offending key in"))].strip().split(":")[0]
        offendedKey = int(out[out.index("Offending key in"):out.find("\n",out.index("Offending key in"))].strip().split(":")[1])
	print "There are invalid keys in known_hosts.\nDelleting keyHost "+str(offendedKey)+" and trying again..."
        os.system("cp " + known_hostsFile + " " + known_hostsFile+"_save")
        file = open(known_hostsFile+"_save","r")
        lines = file.readlines()
        file.close()

        os.system("head -n "+str(offendedKey-1)+" "+known_hostsFile+"_save" + " > " +known_hostsFile)
        os.system("tail -n "+str(len(lines)-offendedKey) + " "+known_hostsFile+"_save" + " >> " +known_hostsFile)

        #file = open(known_hostsFile,"w")
        #for idx,line in enumerate(lines):
        #    if idx != offendedKey-1:
        #        file.write(line.strip().strip()+"\n")
        #file.close()
        status,out = commands.getstatusoutput("expect "+tmpScript)	

    removeDir(tmpDir)
    retValue = False
    if "Permanently added" in out and "known hosts" in out and status == 0:
        if verbose > 0: print "Added RSA fingerprint to known hosts"
        retValue = True
    elif "validateRSAkey" in out  and "not known" not in out:
        if verbose > 0: print "RSA key already in known hosts"
        retValue = True
    elif "Name or service not known" in out:
        print "Unable to connect"
        return False
    elif "command not found" in out:
        print "The command 'expect' is unavailable. You have to ssh a first time to "+\
                target+" manually or install 'expect'"
        return False

    if "password:" in out:
        print "SSH server is waiting for password. Please configure SSH to use DSA key authentication."
        retValue = False 
    return retValue

def testInterfacelessSSH(machine, user=None, timeout = 5):
    """Tests ssh to the user@machine without user interface.
       Use ONLY for debuging/testing since it adds an overhead of at least 5 seconds
       timeout in seconds
    """
    if user:
        target = user+"@"+machine
    else:
        target = machine

    #Test SSH connection to localhost
    Exec = "ssh"
    args = [Exec, target, "echo \"SSH ok\""]
    if verbose > 0: print "Testing SSH connection to "+machine+"..."
    PID = os.spawnvpe(os.P_NOWAIT, Exec , args, os.environ)
    if PID == 0:
         print "Cannot start SSH"
         return False

    #Pooling checking ih the ssh command finished by itself
    waiting = True
    for n in range(timeout*10):
        psLine = os.popen( "ps -p "+str(PID) ).read()
        if "ssh" in psLine and "<defunct>" not in psLine:
            waiting = True
        else:
            waiting = False
            break
        time.sleep(0.1)
    if waiting:
        print "\nSSH is waiting for user response. Please configure SSH to use DSA key authentication."
        os.kill(PID,9)
        return False
    #Double checking
    status = commands.getstatusoutput("ssh "+target+" echo \"SSH ok\"")
    if status[0] != 0 or status[1].split("\n")[-1] != "SSH ok":
        print "\nSSH did not respond correctly. Response received:",status
        return False

    return True 



def createScratchDir(desc = "", rmFirst = True, baseDir = SCRATCHDIR):
    """Creates a unique scratch dir with an optional description on dir name passed with keyword desc.
       if the scratch dir exists, remove it first (if flag rmFirst is true)
       A baseDir can be specified using he parameter baseDir
       Returns the created scratchDir or None if it was not possible to create it"""
    randNr = random.randint(0,10000)
    scratchdir = os.path.realpath(os.path.join(baseDir, "scratchdir"+str(desc)+str(time.time()).replace(".","")+"_"+str(randNr)))
    try:
        # Make sure we will not never delete the entire base dir
        if os.path.exists(scratchdir) and rmFirst and os.path.realpath(scratchdir) != os.path.realpath(baseDir):
            os.system("rm -rf " + scratchdir)
        os.system("mkdir -p " + scratchdir)
        if os.path.isdir(scratchdir):
            return scratchdir
        else:
            print "ERROR: Could not create ",scratchdir
            return None
    except:     
        if os.path.isdir(scratchdir):
            print "WARNING: The scratch dir specified exists, but something wrong happened while removing/creating it."
            return scratchdir
        else:
            print "ERROR: Could not create the directory: ",scratchdir
            return None

def removeDir(dirToRem):
    """Securely remove a directory (maily for use with scraatch dirs)
       If the dir exists, remove it with -rf options, if not, do not rise any error
       Returns True if the dir existed, otherwise, returns False"""
    if dirToRem and os.path.isdir(dirToRem):
        os.system("rm -rf "+dirToRem)
        return True
    else:
        return False

def isNumber(inStr):
    """Returns True if the string inStr can be converted to float, i.e. it ts a number (int or float)"""
    try:
        num = float(inStr)
        return True
    except:
        return False

def isInt(inStr):
    """Returns True if the string inStr can be converted to int"""
    try:
        num = int(inStr)
        return True
    except:
        return False
##ecPA


def getLatest(dateList):
    """Method which returns the latest date outof a list of strings of the format YYYY-MM-DD"""

    dateTimeList = []
    for date in dateList:
        try: dateTimeList.append(datetime(*strptime(date, "%Y-%m-%d")[0:3]))
        except: pass

    dateTimeList.sort()

    latest = str(dateTimeList[len(dateTimeList)-1])[0:10]
    return latest

def power2Range(first, last, step = 1):
    """
    Returns a list with values separated by "step" starting at "first"
    value until reaches the "last" value.
    unlike the orange range(), the step can be a non-integer value. ex: power2Range(1,10,1)
        result = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    The step can be negative
    """
    list = []
    ##scPA
    if (type(first) != types.IntType and type(first) != types.FloatType):
        if verbose > 0: print "miscUtilities.power2Range: invalid input type! Only alowed int or float"
        return list
    if (type(last) != types.IntType and type(last) != types.FloatType):
        if verbose > 0: print "miscUtilities.power2Range: invalid input type! Only alowed int or float"
        return list
    if (type(step) != types.IntType and type(step) != types.FloatType):
        if verbose > 0: print "miscUtilities.power2Range: invalid input type! Only alowed int or float"
        return list


    if step==0 or (last<first and step>0) or (last>first and step<0) or first==last:
        return [2**first]
    ##ecPA
    list = Range(first, last, step)
 
    for i in range(0,len(list)):
        list[i] = 2**list[i]

    if len(list)==0:
	return [2**first]	
    else:
        return list

##scPA
def Range(first, last, step = 1):
    """
    Returns a list with values separated by "step" starting at "first" 
    value until reaches the "last" value.
    The steps can be negative
    """
    list = []

    if (type(first) != types.IntType and type(first) != types.FloatType):
        if verbose > 0: print "miscUtilitiesp.Range: invalid input type! Only alowed int or float"
        return list
    if (type(last) != types.IntType and type(last) != types.FloatType):
        if verbose > 0: print "miscUtilitiesp.Range: invalid input type! Only alowed int or float"
        return list
    if (type(step) != types.IntType and type(step) != types.FloatType):       
        if verbose > 0: print "miscUtilitiesp.Range: invalid input type! Only alowed int or float"
        return list


    if step==0 or (last<first and step>0) or (last>first and step<0) or first==last:
        return [first]

    list.append(first)
    if last>first:
        while list[-1]+step <= last:
            list.append(list[-1]+step)
    else:
        while list[-1]+step >= last:
            list.append(list[-1]+step)

    if len(list)==0:
        return [first]
    else:
        return list

##ecPA
