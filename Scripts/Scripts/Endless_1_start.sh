#!/bin/bash

while true; do

Files=$(ls ../Data/Viral_proteomes | wc -l)

if [ $Files != "1" ]; then
echo "New Viruses found"
./Script_1.sh
fi

sleep 1200
done
