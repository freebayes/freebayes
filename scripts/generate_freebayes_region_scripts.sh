#!/bin/bash

regionfile=$1
stub=$2
outputdir=$3
scriptdir=$4

if [ $# -ne 4 ]
then
    echo "usage: $0 [region list] [stub script] [output directory] [script directory]"
    echo "Writes single-region scripts using the list of regions in [region list]."
    echo "The scripts will write their output to [output directory]."
    echo "The stub script has REGION in the place of the region and OUTPUT in the"
    echo "place of the output directory."
    exit
fi

mkdir -p $scriptdir
mkdir -p $outputdir

for region in $(cat $regionfile)
do
    echo writing script for $region with output to $outputdir
    cat $stub | sed -e "s/REGION/$region/g" | sed -e "s%OUTPUT%$outputdir%g" >$scriptdir/$region.sh
done
