echo "Installing dependent packages..."
sudo apt-get install -y tcsh
sudo apt-get install -y ssh
#sudo apt-get install -y dpkg-dev
sudo apt-get install -y build-essential
sudo apt-get install -y GCC
sudo apt-get install -y g++
sudo apt-get install -y gfortran
sudo apt-get install -y libnspr4-dev
sudo apt-get install -y mpich2  
sudo apt-get install -y swig	
sudo apt-get install -y libatlas-base-dev
sudo apt-get install -y python-numpy
#sudo apt-get install -y libgtk2.0-dev
sudo apt-get install -y python-sip
#sudo apt-get install -y qt3-dev-tools
#sudo apt-get install -y qt4-dev-tools
#Using qt Ver. 4.6.2
sudo apt-get install -y python-qt4
#sudo apt-get install -y python-qt-dev
sudo apt-get install -y python-dev
#sudo apt-get install -y python2.5-dev
sudo apt-get install -y liblapack-dev
#sudo apt-get install -y python2.5-qt3
sudo apt-get install -y python-qwt5-qt4
sudo apt-get install -y git-core gitosis
sudo apt-get install -y subversion
#sudo apt-get install -y sip4
#Openbabel
sudo apt-get install -y openbabel
sudo apt-get install -y python-openbabel
#WARNING for uninstall: cinfony is installed by the install script in the system default location !

echo "Updating db..."
sudo updatedb

echo "Configure Git if you plan to check in changes:"
echo '      git config --global user.name "YouGitUserName"'
echo "     git config --global user.email YourEmail@ServerX"

