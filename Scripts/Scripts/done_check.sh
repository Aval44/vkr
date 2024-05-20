#!/bin/bash
#$1 –> path to directory with .tsv files $2 –> number of proteins

One=$(wc -l < ../Data/MHC/MHCI.txt) 
Two=$(wc -l < ../Data/MHC/MHCII.txt) 
MHC=$((One+Two)) 
Prot=$2 
Files=$((MHC*Prot))

while true; do

sleep 60

Got=$(ls $1 | wc -l) 
echo "Calculated" $((Got*100/Files)) "%"
done
