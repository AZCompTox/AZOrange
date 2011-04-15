#!/bin/sh

# Clean old files and directories from previous executions

NFS_SD="$HOME/AZO_NFS_scratchDir"
if [ ! -d $NFS_SD ]; then
    mkdir $NFS_SD
fi

DIRS="$HOME/qsubDataDebug $NFS_SD"
for DIR in $DIRS; do
        if [ -d $DIR ]; then
            OLDDIRS=`find "$DIR" -mindepth 1 -maxdepth 1 -type d -ctime +14 -exec echo {} \;`
            for OLD in $OLDDIRS; do
                    rm -rf $OLD
            done
        fi
done

OLDFILES=`find "$NFS_SD" -mindepth 1 -maxdepth 1 -type f -ctime +14 -exec echo {} \;`
for OLD in $OLDFILES; do
        rm -f $OLD
done

