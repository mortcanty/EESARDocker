#!/bin/bash
#
# Usage: 
# run_sar_seqQ.sh pattern imdir enl significance 
#

alpha="${@: -1}"
enl=("${@: -2}")
imdir=("${@: -3}")

echo $1

fns=$(ls -l $imdir | grep $1 | \
     grep -v 'sub.tif' |  grep -v 'enl.tif' | awk '{print $9}')
     
echo ${fns//$1/$imdir$1}    
     
python scripts/sar_seqQ.py -s $alpha -m \
                     ${fns//$1/$imdir$1} sarseqQ.tif $enl 