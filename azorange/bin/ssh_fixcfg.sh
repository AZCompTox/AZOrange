#!/bin/sh

# Backup current .ssh directory if it exists

if [ ! -d $HOME/.ssh ]; then
	echo "$HOME/.ssh does not exist. Backup not needed."
else
	cp -r $HOME/.ssh $HOME/.ssh.backup-`date|tr " " "_"`
fi 

# Remove localhost if it exists, we want to check that ssh will work even if hosts are missing.

if [ -f $HOME/.ssh/known_hosts ]; then
	sed -i '/^localhost/d' $HOME/.ssh/known_hosts
fi

OUTPUT=`ssh -o BatchMode=yes localhost true 2>&1`
echo "$OUTPUT" |grep -q "Host key verification failed"
if [ "$?" == 0 ]; then
        echo "Fix Host"
        if [ -f $HOME/.ssh/config ]; then
                grep -q "^Host" $HOME/.ssh/config
            	if [ "$?" == 0 ]; then
                    grep -q "[sS]trictHostKeyChecking" $HOME/.ssh/config
                     if [ "$?" != 0 ]; then
	             	sed -i 's/\(^Host.*\)/\1\n       StrictHostKeyChecking no/g' $HOME/.ssh/config
                     else
			sed -i 's/^\(.*\)\([sS]trictHostKeyChecking[      ]\)\(.*\)$/\1\2no/g' $HOME/.ssh/config
		     fi
                fi
	else
            cat << EOF > $HOME/.ssh/config
        Host *
        StrictHostKeyChecking no
        ServerAliveInterval 45
EOF
	fi
fi

OUTPUT=`ssh -o BatchMode=yes localhost true 2>&1`
EXITCODE=$?
echo $OUTPUT
echo "$OUTPUT" |grep -q "Permission denied"
if [ "$?" == 0 ]; then
	echo "Fixing keys..."
        if [ -f $HOME/.ssh/id_dsa.pub ] && [ -f $HOME/.ssh/id_dsa ]; then 
		cat $HOME/.ssh/id_dsa.pub >>$HOME/.ssh/authorized_keys
        else 
	        if [ -f $HOME/.ssh/id_rsa.pub ] && [ -f $HOME/.ssh/id_rsa ]; then 
			cat $HOME/.ssh/id_rsa.pub >>$HOME/.ssh/authorized_keys
	        else 
                        mkdir -p $HOME/.ssh
		        ssh-keygen -t dsa -f $HOME/.ssh/id_dsa -N ""
		        cp $HOME/.ssh/id_dsa.pub $HOME/.ssh/authorized_keys
		fi
        fi
fi

# Test the fix above, if still errors exist then print error

OUTPUT=`ssh -o BatchMode=yes localhost true 2>&1`
echo "$OUTPUT" |grep -q "Permission denied"
if [ "$?" == 0 ]; then
	echo "Failed to configure passwordless authentication keys!"
        echo "Error: $OUTPUT"
	exit 1
fi

# Remove localhost if it exists, we want to check that ssh will work even if hosts are missing.

if [ -f $HOME/.ssh/known_hosts ]; then
	sed -i '/^localhost/d' $HOME/.ssh/known_hosts
fi

OUTPUT=`ssh -o BatchMode=yes localhost true 2>&1`
echo "$OUTPUT" |grep -q "Host key verification failed"
if [ "$?" == 0 ]; then
	echo "Failed to configure automatic host acceptance!"
        echo "Error: $OUTPUT"
	exit 1
fi

# Remove localhost if it exists, we want to check that ssh will work even if hosts are missing.

if [ -f $HOME/.ssh/known_hosts ]; then
	sed -i '/^localhost/d' $HOME/.ssh/known_hosts
fi

OUTPUT=`ssh -o BatchMode=yes localhost true 2>&1`
if [ "$?" != 0 ]; then
	echo "Failed to configure ssh. Try manually to fix your ssh config"
        echo "Error: $OUTPUT"
	exit 1
fi
ssh -o BatchMode=yes 127.0.0.1 true 2>&1

echo "Success setting up ssh configuration."
