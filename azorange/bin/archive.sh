#!/bin/sh

# Archive old directories from previous executions

DIRS="$HOME/qsubDataDebug $HOME/AZO_NFS_scratchDir"
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
