#!/bin/bash
#$1 –> Virus name

###===Folders for output data===###
#Affiniy_$1 –> file describing affinity check results for all alleles
touch ../Data/Affinity_results/Affinity_$1.txt
#Directory with complete affinity results (.tsv files)
mkdir ../Data/Affinity_results/$1

PtF="../Data/Viral_proteomes/"$1".fasta" # path to the sequence file

###===Extracting proteins of interest===###
OLDIFS=$IFS
ProtArray=()
mapfile -t LineArray < <(sed -n '/^>/s/^/ /p' $PtF | cut -d' ' -f2-)

for Line in "${LineArray[@]}"
do
FW=$(echo "$Line" | awk '{print $1}')
IFS="|"
read -a PNs <<< $FW
IFS=$OLDIFS
ProtArray+=("${PNs[1]}""|""${PNs[2]}")
done
###===Loop through selected proteins===###
echo ">" >> $PtF

echo "Proteins:"
touch tmp/proteins.fasta
for Protein in "${ProtArray[@]}"
do
echo $Protein
echo "> "$Protein >> tmp/proteins.fasta
sed -n -e "/$Protein/,/>/p" $PtF > tmp/add
sed -i '1d;$d' tmp/add
cat tmp/add >> tmp/proteins.fasta
rm tmp/add
done

###===Running cleavage script===###
./Script_3.sh tmp/proteins.fasta easy 

###===Allelles of interest===###
mapfile -t AllelesI < ../Data/MHC/MHCI.txt
mapfile -t AllelesII < ../Data/MHC/MHCII.txt

###===Running affinity script===###
echo "running script 4 –> calculating affinity...."

# loop through alleles
function mhci {
	./Script_4.sh $1 $2 0.5 2 I ../Data/Affinity_results
}

function mhcii {
	./Script_4.sh $1 $2 1 5 II ../Data/Affinity_results
}

export -f mhci
export -f mhcii

./done_check.sh "../Data/Affinity_results/${1}" ${#ProtArray[@]} &

parallel  --jobs 46 mhci $1 {} ::: "${AllelesI[@]}"
T=$(date)
echo "MHCI DONE AT: " $T

parallel --jobs 49 mhcii $1 {} ::: "${AllelesII[@]}"
T=$(date)
echo "MHCII DONE AT: " $T

killer=$(ps -ef | grep "done_check.sh" | grep -v "grep" |  awk '{print $2}')
kill $killer

rm -rf tmp/Peptides_I
rm -rf tmp/Peptides_II
rm tmp/proteins.fasta

cp ../Data/Viral_proteomes/${1}.fasta ../Data/Viral_proteomes/Ready/${1}.fasta
rm ../Data/Viral_proteomes/${1}.fasta
