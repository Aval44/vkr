#!/bin/bash
#$1 –> Viral name
#$2 –> Allele name
#$3 –> Rank Threshold for strong binding peptides
#$4 –> Rank Threshold for weak binding peptides
#$5 –> Allele type (I or II)
#$6 –> Path to the parent directory of Affinity.txt

###===Getting the Protein names list===###

PepArray=()
for Pep in tmp/Peptides_I/*; do
   Pep=${Pep:15}
   PepArray+=("$Pep")
done

# Allele name for Affinity_$1.txt
echo "> " $2 > tmp/Affinity_$1_$2.txt
MHC_name=$2

###===Checking MHC type===###
############################## M H C I ##############################
if [ $5 == "I" ]; then

## Interating for all peptides
for Prot_name in ${PepArray[@]}
do

###===Running affinity software (netMHCpan)===###

./../software/netMHCpan-4.1/netMHCpan -f tmp/Peptides_I/${Prot_name} -a $2  -rth $3 -rlt $4 -inptype 1 -s > tmp/tmp_af_${2}

###===Parsing===###

## for each protein filling Affinity_$1.txt and a .tsv file with complete data

# Number of row with final information
Line=$(wc -l < tmp/tmp_af_${2})
Line=$((Line-2))

Min_aff=$(awk -v i="$Line" 'NR==i{print $9}' tmp/tmp_af_${2})
Min_aff=${Min_aff:0:-1}
Max_aff=$(awk -v i="$Line" 'NR==i{print $14}' tmp/tmp_af_${2})
Max_aff=${Max_aff:0:-1}
Total_pep=$(awk -v i="$Line" 'NR==i{print $18}' tmp/tmp_af_${2})

#.tsv file
echo "   Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL" > $6/$1/${Prot_name}_${MHC_name}.tsv
end_line=$((Line-3))
sed -n "53,${end_line}p" tmp/tmp_af_${2} >> $6/$1/${Prot_name}_${MHC_name}.tsv

sed -i  's/<= WB//g' $6/$1/${Prot_name}_${MHC_name}.tsv
sed -i  's/<= SB//g' $6/$1/${Prot_name}_${MHC_name}.tsv
sed -i  's/^...//g' $6/$1/${Prot_name}_${MHC_name}.tsv
sed -i   -E 's/ +/\t/g' $6/$1/${Prot_name}_${MHC_name}.tsv

#Affinity_$1.txt
echo "Protein ${Prot_name}" >> tmp/Affinity_$1_$2.txt
echo "high: ${Max_aff}; low: ${Min_aff}; overall: ${Total_pep}" >>  tmp/Affinity_$1_$2.txt

done
fi

############################## M H C II ##############################
if [ $5 == "II" ]; then

## Interating for all peptides
for Prot_name in ${PepArray[@]}
do

###===Running affinity software (netMHCIIpan)===###

./../software/netMHCIIpan-4.3/netMHCIIpan -f tmp/Peptides_II/${Prot_name} -a $2 -rankS 1  -rankW 5 -filter -rankF 5 -inptype 1 -s > tmp/tmp_af_${2}

###===Parsing===###

## for each protein filling Affinity_$1.txt and a .tsv file with complete data

# Number of row with final information
Line=$(wc -l < tmp/tmp_af_${2})
Line=$((Line-1))

Min_aff=$(awk -v i="$Line" 'NR==i{print $10}' tmp/tmp_af_${2})
Max_aff=$(awk -v i="$Line" 'NR==i{print $5}' tmp/tmp_af_${2})
Total_pep=$(wc -l < tmp/Peptides_II/${Prot_name})
end_line=$((Line-2))

#.tsv file
echo "   Pos MHC Peptide Of Core Core_Rel Inverted Identity Score_EL %Rank_EL Exp_Bind" > $6/$1/${Prot_name}_${MHC_name}.tsv

sed -n "16,${end_line}p" tmp/tmp_af_${2} >> $6/$1/${Prot_name}_${MHC_name}.tsv

sed -i  's/<= WB//g' $6/$1/${Prot_name}_${MHC_name}.tsv
sed -i  's/<= SB//g' $6/$1/${Prot_name}_${MHC_name}.tsv
sed -i  's/^[ \t]*//g' $6/$1/${Prot_name}_${MHC_name}.tsv
sed -i   -E 's/ +/\t/g' $6/$1/${Prot_name}_${MHC_name}.tsv

#Affinity_$1.txt
echo "Protein ${Prot_name}" >> tmp/Affinity_$1_$2.txt
echo "high: ${Max_aff}; low: ${Min_aff}; overall: ${Total_pep}" >> tmp/Affinity_$1_$2.txt

done
fi

cat tmp/Affinity_$1_$2.txt >> ${6}/Affinity_$1.txt
rm tmp/Affinity_$1_$2.txt

rm tmp/tmp_af_${2}
