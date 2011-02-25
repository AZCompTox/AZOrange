#!/bin/sh

# Archive old directories from previous executions

NFS_SD="$HOME/AZO_NFS_scratchDir"
if [ ! -d $NFS_SD ]; then
    mkdir $NFS_SD
fi

DIRS="$HOME/qsubDataDebug $NFS_SD"
for DIR in $DIRS; do
        if [ -d $DIR ]; then
            mkdir -p $DIR.archive
            OLDDIRS=`find "$DIR" -mindepth 1 -maxdepth 1 -type d -ctime +14 -exec echo {} \;`
            for OLD in $OLDDIRS; do
                    tar czf `dirname $OLD`.archive/`basename $OLD`.tar.gz $OLD 2>/dev/null
                    rm -rf $OLD
            done
        fi
done

#Archive loose files in AZO_NFS_scratchDir
mkdir -p $NFS_SD.archive/looseFiles
mkdir -p $NFS_SD/looseFiles.tmp
FOUNDFILES=0
OLDFILES=`find "$NFS_SD" -mindepth 1 -maxdepth 1 -type f -ctime +14 -exec echo {} \;`
for OLD in $OLDFILES; do
        FOUNDFILES=1
        mv $OLD $NFS_SD/looseFiles.tmp
done
if [ $FOUNDFILES -ne 0 ]; then
    tar czf $NFS_SD.archive/looseFiles/`date +%s`.tar.gz $NFS_SD/looseFiles.tmp 2>/dev/null
fi
rm -rf $NFS_SD/looseFiles.tmp

