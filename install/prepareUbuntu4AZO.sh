echo "Installing dependent packages..."
sudo apt-get install -y tcsh
sudo apt-get install -y ssh
sudo apt-get install -y build-essential
sudo apt-get install -y GCC
sudo apt-get install -y g++
sudo apt-get install -y gfortran
sudo apt-get install -y libnspr4-dev
sudo apt-get install -y mpich2  
sudo apt-get install -y libcr-dev
sudo apt-get install -y swig	
sudo apt-get install -y libatlas-base-dev
sudo apt-get install -y python-numpy
sudo apt-get install -y python-sip
sudo apt-get install -y python-networkx 
#Using qt Ver. 4.6.2
sudo apt-get install -y python-qt4
sudo apt-get install -y python-dev
sudo apt-get install -y liblapack-dev
sudo apt-get install -y python-qwt5-qt4
sudo apt-get install -y git-core gitosis
sudo apt-get install -y subversion
sudo apt-get install -y git-core gitosis

#Openbabel
sudo apt-get install -y openbabel
sudo apt-get install -y python-openbabel

#RDKit
sudo apt-get install -y cmake
sudo apt-get install -y bison
sudo apt-get install -y flex
sudo apt-get install -y sqlite3
sudo apt-get install -y libsqlite3-dev
sudo apt-get install -y libboost-all-dev

#CDK
sudo apt-get install -y python-jpype

echo "Updating db..."
sudo updatedb
echo ""
echo "============================================================="
echo "             Finished preparation of Ubuntu"
echo "-------------------------------------------------------------"
echo "Configure Git if you plan to check in changes:"
echo '      git config --global user.name "YouGitUserName"'
echo "      git config --global user.email YourEmail@ServerX"
echo ""

