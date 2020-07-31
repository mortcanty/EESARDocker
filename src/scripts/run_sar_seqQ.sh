#!/bin/bash
#
# Usage: 
# run_sar_seqQ.sh imdir enl alpha outfile dims 
#

fns=$(ls -l $1 | grep -v _ | awk '{print $9}')

paths=()
b=' '

for fn in $fns
do
   paths+=$1$fn$b
done
    
python scripts/sar_seqQ.py -s $3 -d $5 -m  $paths $4 $2
