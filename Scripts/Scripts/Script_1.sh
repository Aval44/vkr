#!/bin/bash

###===Loop through all files in Data/Viral_Proteomes directory===###
Viruses=($(ls ../Data/Viral_proteomes/ | grep fasta))

for Virus in "${Viruses[@]}"
do
Virus=${Virus::-6}
echo "!!!_Analysing Virus :" $Virus
T=$(date)
echo "TIME: " $T
./Script_2.sh $Virus
echo "––––––––"
echo "––––––––"
done
