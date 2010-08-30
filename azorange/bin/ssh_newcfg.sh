#!/bin/sh

if [ -d $HOME/.ssh ]; then
	cp -r $HOME/.ssh $HOME/.ssh.backup-`date|tr " " "_"`
	rm -rf $HOME/.ssh
else 
	echo "$HOME/.ssh does not exist. Backup not needed"
fi
ssh-keygen -t dsa -f $HOME/.ssh/id_dsa -N ""
cp $HOME/.ssh/id_dsa.pub $HOME/.ssh/authorized_keys
cat << EOF > $HOME/.ssh/config
Host *
	StrictHostKeyChecking no
        ServerAliveInterval 45
EOF
ssh -o BatchMode=yes localhost true
if [ "$?" != 0 ]; then
	echo "Failed setting up a new working ssh configuration. See your system adm."
        exit 1
fi
ssh -o BatchMode=yes 127.0.0.1 true
echo "Success setting up ssh configuration."

