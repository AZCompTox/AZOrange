#!/bin/sh

patchDir=$1

if [ "$patchDir" = "" ]
  then
    echo "You must specify a patch Directory!"
    exit 1
fi

if [ ! -d $patchDir ]
  then
  echo "$patchDir is not a Directory!"
  exit 1
fi

if [ ! -d $patchDir/install-scripts ]
  then
  echo "The Directory specified does not correspond to an orange2.0 structure:"
  echo "$patchDir/orange is not a valid Directory!"
  exit 1
fi

if [ ! -d $patchDir/source ]
  then
  echo "The Directory specified does not correspond to an orange2.0 structure:"
  echo "$patchDir/source is not a valid Directory!"
  exit 1
fi

if [ ! -d $patchDir/Orange/OrangeCanvas ]
  then
  echo "The Directory specified does not correspond to an orange2.0 structure:"
  echo "$patchDir/OrangeCanvas is not a valid Directory!"
  exit 1
fi

if [ ! -d $patchDir/Orange/OrangeWidgets ]
  then
  echo "The Directory specified does not correspond to an orange2.0 structure:"
  echo "$patchDir/OrangeWidgets is not a valid Directory!"
  exit 1
fi


clear

iniDir=`pwd`

#Copy all the AZO widgets to the right place
#/bin/cp -rf $iniDir/OrangeWidgets/* $patchDir/OrangeWidgets/

cd $patchDir
echo "Patching orange2.0 to AZO in $patchDir..."
echo "==============================================================="

patch --dry-run -fs -p1 -i $iniDir/Orange2_0.patch
STATUS=$?
#patch --dry-run -fs -p1 -i $iniDir/OrangeLayout.patch
#STATUS=`expr $? + $STATUS`
echo "Will Return Status: $STATUS"
if [ $STATUS -ne 0  ]
then
    exit $STATUS
fi

echo "==============================================================="
echo "The outputs before this statement correspond to a simulation"
echo "of the patch to be applied."
echo "Please check if the simulation was successful to continue to"
echo "the real patch."
echo "Example of success patched files:"
echo "    patching file orngTest.py"
echo "    patching file setup.py"
echo ""
echo "Continue? (yes/no)"

#change "yes" to "" in order to ask before apply patch
TEST="yes"

while [ "$TEST" != "yes" ]
do
        read  TEST
        if [ "$TEST" = "no" ]
          then
            echo "Patch canceled by user."
            exit 1
        fi
        if [ "$TEST" != "yes" ]
          then
            echo "Please choose yes/no"
        fi
done

echo "----------------------------------------------------------"
echo "Applying Patch..."
echo "----------------------------------------------------------"
#Apply the AZO patch to original orange files
patch -fs -p1 -i $iniDir/Orange2_0.patch
STATUS=$?
#patch -fs -p1 -i $iniDir/OrangeLayout.patch
#STATUS=`expr $? + $STATUS`

# Patch for removing duplicated files in Prototypes which would conflict with the ones in other tabs:
rm -f OrangeWidgets/Prototypes/OWMergeData.py
rm -f OrangeWidgets/Prototypes/OWConcatenate.py 

cd $iniDir
echo "Return Status: $STATUS"
exit $STATUS


