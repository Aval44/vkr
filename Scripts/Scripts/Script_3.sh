#!/bin/bash
#$1 –> file for cleavage
#$2 –> method [easy, netchop]
#$3 –> if netchop is chosen, provide the threshold
echo "running script 3 –> proteins cleavage $2 mode ...."

# output directories
mkdir tmp/Peptides_I
mkdir tmp/Peptides_II

###===Extracting rows with protein names and starting a loop===###
mapfile -t Start_Index < <(awk '/>/{print NR}' $1)

# adding artificial row (for iteration)
Seq_len=$(wc -l < $1)
NL=$((Seq_len+1))
Start_Index+=($NL)

Check=0
Prev=0
for Row in "${Start_Index[@]}"
do

if [ $Check -ne 0 ]; then

Prot_name=$(awk -v i="$Prev" 'NR==i{print $2}' $1)
end_line=$((Row-1))
start_line=$((Prev+1))
Seq=$(sed -n "${start_line},${end_line}p" $1)

# arrays for cleaved peptides
Peptides_I=()
Peptides_II=()

########### E A S Y ###########
if [ $2 == "easy" ]; then

echo "$Seq" > tmp/seq
Seq=$(tr -d '\n' < tmp/seq)
Len=${#Seq}
rm tmp/seq

#------------------- M H C I ------------------- (8–11)
for (( i=8; i < 12; i++ )); do
for (( j=0; j <= $((Len-i-1)); j++ )); do

Add=${Seq:$j:$i}
if [[ $Add != *"X"* ]]; then
Peptides_I+=($Add)
fi
done
done 

#------------------- M H C I I ------------------- (12-18)
for (( i=12; i < 19; i++ )); do
for (( j=0; j <= $((Len-i-1)); j++ )); do

end=$((i+j))
Add=${Seq:$j:$i}
if [[ $Add != *"X"* ]]; then
Peptides_II+=($Add)
fi
done
done 
fi
########### N E T C H O P ##########
if [ $2 == "netchop" ]; then
echo "> LH44Blessed" > tmp/seq
echo "$Seq" >> tmp/seq

#------------------- M H C I ------------------- (8–10)
./../software/netchop-3.1/netchop -v 0 -t $3 tmp/seq > tmp/chopped_I

Start_row=21
Seq_len=$(wc -l < tmp/chopped_I)
End_row=$((Seq_len-5))
Add="" # storing peptide sequence

for (( k=$Start_row; k <= $End_row; k++ )); do

AAcid=$(awk -v i="$k" 'NR==i{print $2}' tmp/chopped_I)
If_Chopped=$(awk -v i="$k" 'NR==i{print $3}' tmp/chopped_I)

if [[ $AAcid != "X" ]]; then
	Add+=$AAcid
	if [[ $If_Chopped == "S" ]] ; then
		if [ "${#Add}" -gt 7 ]; then
			if [ "${#Add}" -lt 12 ]; then
				Peptides_I+=($Add)
			fi
		fi
	Add=""
	fi
else
	Add=""
fi
done

if [ "${#Add}" -gt 7 ]; then
	if [ "${#Add}" -lt 11 ]; then
		Peptides_I+=($Add)
	fi
fi

#------------------- M H C I I ------------------- (12-18)
./../software/netchop-3.1/netchop -v 1 -t $3 tmp/seq > tmp/chopped_II

Add="" # storing peptide sequence
for (( k=$Start_row; k <= $End_row; k++ )); do

AAcid=$(awk -v i="$k" 'NR==i{print $2}' tmp/chopped_II)
If_Chopped=$(awk -v i="$k" 'NR==i{print $3}' tmp/chopped_II)

if [[ $AAcid != "X" ]]; then
	Add+=$AAcid
	if [[ $If_Chopped == "S" ]]; then
		if [ "${#Add}" -gt 11 ]; then
			if [ "${#Add}" -lt 19 ]; then
				Peptides_II+=($Add)
			fi
		fi
	Add=""
	fi
else
	Add=""
fi
done

if [ "${#Add}" -gt 11 ]; then
	if [ "${#Add}" -lt 19 ]; then
		Peptides_II+=($Add)
	fi
fi

# deleting tmp
rm tmp/seq
rm tmp/chopped_I
rm tmp/chopped_II
fi

###===Writing output===###
touch tmp/Peptides_I/${Prot_name}
touch tmp/Peptides_II/${Prot_name}

for peptide in "${Peptides_I[@]}"
do
echo "$peptide" >> tmp/Peptides_I/"${Prot_name}"
done

for peptide in "${Peptides_II[@]}"
do
echo "$peptide" >> tmp/Peptides_II/"${Prot_name}"
done

fi

Prev=$Row
Check=1
done #loop through all proteins
