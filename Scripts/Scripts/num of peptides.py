###===Importing libraries===###
import argparse
import pandas as pd
pd.options.mode.chained_assignment = None
from pepmatch import Preprocessor, Matcher
# for some reason I cant specify the db file path manually(

###===Reading arguments===###
parser= argparse.ArgumentParser(description= 'Script arguments')
parser.add_argument("-n","--name", type= str, help= "viral name", required= True)

args= parser.parse_args()
VirName= args.name

###===Loading MHC alleles===###
MHCI= open("../Data/MHC/MHCI.txt", "r")
MHCI= MHCI.read().split("\n")
MHCI= MHCI[:-1]

MHCII= open("../Data/MHC/MHCII.txt", "r")
MHCII= MHCII.read().split("\n")
MHCII= MHCII[:-1]

###===Extracting protein names===###
Short= open("../Data/Affinity_results/Affinity_"+VirName+".txt", "r")
Short= Short.read().split(">")[1:]
Shortone= Short[0].split(" ")[2:]

Proteins=[]
for I in range (1, len(Shortone), 6):
    add= Shortone[I].split("\n")[0]
    Proteins+= [add]

###===File with viral human-like presented peptides
VHLPP= VirName+"\n"

###===Viral table for MHC_I===###
ColNames= ["Virus", "MHC", "Protein", "# of high binding", "# of low binding", "sum_score high", "sum_score low", "high normed on protein", "low normed on protein"]
VirDF_I= pd.DataFrame(columns= ColNames)

###===Loop through all MHC_I alleles===###
DDD= 0
HHH= []
for x in MHCI:
    for y in Proteins:
        TN= y+"_"+x+".tsv"
        Panda= pd.read_csv('../Data/Affinity_results/'+VirName+"/"+TN, sep= '\t') # Proten vs HLA table
        if type(list(Panda['Pos'].unique())[0]) == type("str"):
            Panda.reset_index(inplace= True)
            Panda.drop("index",axis= 1,inplace= True)
            Keys= list(Panda.keys())
            Keys= Keys[1:]
            Panda.drop("%Rank_EL",axis= 1,inplace= True)
            Panda.columns= Keys
        else:
            Panda.drop('Pos', axis= 1, inplace= True)
        Panda05= Panda.loc[Panda['%Rank_EL'] <= 0.5] # high affinity
        Panda20= Panda.loc[Panda['%Rank_EL'] <= 2] # high and low affinity
        HHH+= list(Panda05['Peptide'].unique())

###===Viral table for MHC_II===###
VirDF_II= pd.DataFrame(columns= ColNames)

###===Loop through all MHC_II alleles===###
for x in MHCII:
    for y in Proteins:
        TN= y+"_"+x+".tsv"
        Panda= pd.read_csv('../Data/Affinity_results/'+VirName+"/"+TN, sep= '\t') # Proten vs HLA table
        if type(list(Panda['Pos'].unique())[0]) == type("str"):
            Panda.reset_index(inplace= True)
            Panda.drop("index",axis= 1,inplace= True)
            Keys= list(Panda.keys())
            Keys= Keys[1:]
            Panda.drop("%Rank_EL",axis= 1,inplace= True)
            Panda.columns= Keys
        else:
            Panda.drop('Pos', axis= 1, inplace= True)
        Panda10= Panda.loc[Panda['%Rank_EL'] <= 1] # high affinity
        Panda50= Panda.loc[Panda['%Rank_EL'] <= 5] # high and low affinity
        HHH+= list(Panda05['Peptide'].unique())

DDD+= len(list(set(HHH)))

print(VirName+" "+str(DDD))
