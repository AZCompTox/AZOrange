import AZOrangeConfig as AZOC
from AZutilities import miscUtilities
import os

class config():
    def __init__(self):
        self.report = ""

    def ConfigNFSScratchDir(self):
        # Create the AZOrange scratch dir on the nfs system is if does not exist
        if not os.path.exists(AZOC.NFS_SCRATCHDIR):
            os.system("mkdir "+AZOC.NFS_SCRATCHDIR)
            if not os.path.exists(AZOC.NFS_SCRATCHDIR):
                return False
        return True

    def ConfigInterfacelessSSH(self,machine,user=None):
        #Test SSH connection to localhost
        if not miscUtilities.testInterfacelessSSH(machine,user,5):
            print "SSH is not configured for local machine. Trying to fix problem now..."
            if miscUtilities.autoValidateRSAKey(machine,user):
               print "RSA fingerprint added to known_hosts with success"
            else:
               print "Not able to add fingerprint to known_hosts."
            if miscUtilities.testInterfacelessSSH(machine,user,5):
                print "SSH correctly configured for "+machine
            else:
                print "Unable to configure properly interfaceless SSH"
                return False
        return True
                
    def __call__(self):
        self.report = ""
        undefinedError = False
        if not self.ConfigNFSScratchDir():
            self.report += "Unable to create the NFS scratch dir: "+AZOC.NFS_SCRATCHDIR+"\n"
        if not self.ConfigInterfacelessSSH("localhost") or not self.ConfigInterfacelessSSH("127.0.0.1"):
            self.report += "Unable to configure interfaceless ssh.You should take the following actions:\n"
            if os.path.isfile(os.path.join(os.environ["HOME"],".ssh","id_dsa")) and \
               not os.path.isfile(os.path.join(os.environ["HOME"],".ssh","id_dsa.pub")):
                   self.report += "1)In a terminal window create a public SSH key by running the commands:\n"+\
                                  "  On ANY question, just hit ENTER!\n"+\
                                  "      ssh-keygen -y -t dsa > ~/.ssh/id_dsa.pub\n"+\
                                  "      cat ~/.ssh/id_dsa.pub >> ~/.ssh/authorized_keys\n"
            elif not os.path.isfile(os.path.join(os.environ["HOME"],".ssh","id_dsa")):
                   self.report += "2)In a terminal window create private and public SSH keys by running the commands:\n"+\
                                  "  On ANY question, just hit ENTER!\n"+\
                                  "      ssh-keygen -b 1024 -t dsa\n"+\
                                  "      cat ~/.ssh/id_dsa.pub >> ~/.ssh/authorized_keys\n"
            elif not os.path.isfile(os.path.join(os.environ["HOME"],".ssh","authorized_keys")):
                   self.report += "1)In a terminal window run the command:\n"+\
                                  "      cat ~/.ssh/id_dsa.pub >> ~/.ssh/authorized_keys\n"
            else:
                   undefinedError = True
                   self.report += "The problem appears to be on invalid known_hosts or authorized_keys.\n"+\
                                  "1)Try fix it by running in a terminal window run the commands:\n"+\
                                  "  You will probably be asked to take further procedures at next start of AZOrange\n"+\
                                  "      rm ~/.ssh/known_hosts\n"+\
                                  "      rm ~/.ssh/authorized_keys\n"
        if not os.path.isfile(os.path.join(os.environ["HOME"],".ssh","config")):
            if undefinedError:
                self.report = ""
            self.report += "It was detected that there is no '~/.ssh/config' file. \n"+\
                           "1)It is advised to create the SSH config file running for example the commands:\n"+\
                           "      echo \"Host *\" >> ~/.ssh/config\n"+\
                           "      echo \"    StrictHostKeyChecking no\" >> ~/.ssh/config\n"+\
                           "      echo \"    ServerAliveInterval 45\" >> ~/.ssh/config\n"
        else:
            # Read the .ssh/config file
            sshcfg = open(os.path.join(os.environ["HOME"],".ssh","config"),"r") 
            sshcfgText = sshcfg.readlines()
            sshcfg.close()
            #Get the line with the keyword StrictHostKeyChecking
            lineWithKeyCh = None
            for idx,line in enumerate(sshcfgText):
                if "StrictHostKeyChecking" in line:
                    lineWithKeyCh = idx
            #Get the line with the keyword Host *
            lineWithHost = None
            for idx,line in enumerate(sshcfgText):
                if "Host *" in line:
                    lineWithHost = idx
            #Get the line with the keyword KeepAlive
            lineWithKeepAlive = None
            for idx,line in enumerate(sshcfgText):
                if "ServerAliveInterval" in line:
                    lineWithKeepAlive = idx

            #Check ~/.ssh/config is not properly configured
            if (lineWithHost == None or "#") in (sshcfgText[lineWithHost]) or \
               (lineWithKeyCh == None) or ("#" in sshcfgText[lineWithKeyCh]) or ("no" not in sshcfgText[lineWithKeyCh].lower()) or \
               (lineWithKeepAlive == None) or ("#" in sshcfgText[lineWithKeepAlive]) or (" 0" in sshcfgText[lineWithKeepAlive].lower()):
                   if undefinedError:
                       self.report = ""
                   self.report += "It was detected that the file '~/.ssh/config' may not be properly configured. \n"+\
                                  "1)It is advised to edit the file '~/.ssh/config' and make sure that \n"+\
                                  "  you have the the the following configuration it is not commented:\n"+\
                                  "      \n"+\
                                  "      Host *\n"+\
                                  "          StrictHostKeyChecking no\n"+\
                                  "          ServerAliveInterval 45\n"
                   if lineWithHost != None:
                       self.report += "\nTIP: Look in '~/.ssh/config' file arround line number "+str(lineWithHost+1)
                          
        return self.report


if __name__ == "__main__":
    AC = config()
    report = AC()
    if report != "":
        print "AutoConfig report:\n"+report
    else:
        print "AutoConfiguration was done without any errors reported!"


