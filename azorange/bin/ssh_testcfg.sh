#!/bin/sh

errmsg() {
	DIR=$(dirname `readlink /proc/$$/fd/255`)
        echo "--------------------------------------------------------------------------------------------"
	echo "SSH authentication and/or keys not properly configured for AZOrange. Either:"
	echo
	echo "1. Run $DIR/ssh_fixcfg.sh"
	echo "   This creates a backup or your .ssh directory and attempts to fix the current configuration."
	echo "2. Run $DIR/ssh_newcfg.sh"
        echo "   This creates a backup of your .ssh directory and creates a new working ssh configuration."
        echo "--------------------------------------------------------------------------------------------"
}

# Remove localhost if it exists, we want to check that ssh will work even if hosts are missing.

if [ -f $HOME/.ssh/known_hosts ]; then
        sed -i '/^localhost/d' $HOME/.ssh/known_hosts
fi

OUTPUT=`ssh -o BatchMode=yes localhost true 2>&1`
EXITCODE=$?
echo "$OUTPUT" |grep -q "Permission denied"
if [ "$?" -eq 0 ]; then
        errmsg
	exit 1
fi
echo "$OUTPUT" |grep -q "Host key verification failed"
if [ "$?" -eq 0 ]; then
        errmsg
	exit 1
fi

