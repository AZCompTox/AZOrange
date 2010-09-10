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
