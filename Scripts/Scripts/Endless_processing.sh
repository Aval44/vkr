#!/bin/bash

while true; do

Files=$(ls ../Data/Viral_proteomes/Ready | wc -l)

if [ $Files != "1" ]; then
echo "New Viruses found"

###===Loop through all files in Data/Viral_Proteomes directory===###
Viruses=($(ls ../Data/Viral_proteomes/Ready | grep fasta))

for Virus in "${Viruses[@]}"
do
Virusx=${Virus::-6}
echo "!!!_Analysing Virus :" $Virusx
T=$(date)
echo "TIME: " $T

python Raw_data_processing.py --name $Virusx

cp ../Data/Viral_proteomes/Ready/${Virus} ../Data/Viral_proteomes/Ready/Ready_2/${Virus}
rm ../Data/Viral_proteomes/Ready/${Virus}
done

rm tmp/test_pep.fasta
rm human.db
echo "Done!"
echo "–––––"
echo "–––––"

fi

sleep 1200
done




