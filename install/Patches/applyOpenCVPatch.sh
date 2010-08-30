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

if [ ! -d $patchDir/src ]
  then
  echo "The Directory specified does not correspond to an PLearn structure:"
  echo "$patchDir/src is not a valid Directory!"
  exit 1
fi

if [ ! -d $patchDir/interfaces ]
  then
  echo "The Directory specified does not correspond to an PLearn structure:"
  echo "$patchDir/interfaces is not a valid Directory!"
  exit 1
fi


clear

iniDir=`pwd`

cd $patchDir
echo "Patching PLearn in $patchDir..."
echo "==============================================================="

patch --dry-run -fs -p4 -i $iniDir/opencv.patch
STATUS=$?
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
patch -fs -p4 -i $iniDir/opencv.patch
STATUS=$?
cd $iniDir
echo "Return Status: $STATUS"
exit $STATUS


