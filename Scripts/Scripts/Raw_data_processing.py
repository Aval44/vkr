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
        
        ###===PEPMatch===###
        #Preprocessor("../Data/human.fasta").sql_proteome(k = 5)
        
        # only high binding
        add= ""
        Pep= list(Panda05['Peptide'])
        if len(Pep) != 0:
            for z in Pep:
                add+= "> peptide\n"
                add+= z+"\n"
            f= open("tmp/test_pep.fasta", "w")
            f.write(add)
            f.close()
        
            Match= Matcher('tmp/test_pep.fasta', '../Data/human.fasta', 0, 5, output_format= 'dataframe').match()

            Out= Match.loc[Match['Index start'].notna()]
            Ptd= list(Out['Query Sequence'])
            for z in Ptd:
                VHLPP+= "Allele: "+x+" Protein: "+y+" Peptide: "+z+"\n"
        
            Match= Match.loc[Match['Index start'].isnull()]
            Ptd= list(Match['Query Sequence'])
            Panda05= Panda05.loc[Panda05['Peptide'].isin(Ptd)]

        ###===Adding rows to VirDF_I===###
        Row= [None for f in range(9)]
        Row[0]= VirName
        Row[1]= x
        Row[2]= y
        Row[3]= len(Panda05)
        Row[4]= len(Panda20)
        Row[5]= 0.5*Row[3]-Panda05['%Rank_EL'].sum()
        Row[6]= 2*Row[4]-Panda20['%Rank_EL'].sum()
        New= pd.DataFrame(data= [Row], columns= ColNames)
        VirDF_I= pd.concat([VirDF_I, New], axis= 0)

###===Adding "Raw_Proteome" values===###
Proteins+= ["Raw_Proteome"]

for x in MHCI:
    Temp= VirDF_I.loc[VirDF_I['MHC'] == x]
    Row=[None for f in range(9)]
    Row[0]= VirName
    Row[1]= x 
    Row[2]= "Raw_Proteome"
    Row[3]= Temp['# of high binding'].sum()
    Row[4]= Temp['# of low binding'].sum()
    Row[5]= Temp['sum_score high'].sum()
    Row[6]= Temp['sum_score low'].sum()
    New= pd.DataFrame(data= [Row], columns= ColNames)
    VirDF_I= pd.concat([VirDF_I, New], axis= 0)

###===Adding normed values===###
VirDF_I.reset_index(inplace= True)
VirDF_I.drop("index", axis= 1, inplace= True)

Dict= {}
for x in Proteins:
    Dict[x]= [None, None]
    Temp= VirDF_I.loc[VirDF_I['Protein'] == x]
    Dict[x][0]= Temp["sum_score high"].max()
    Dict[x][1]= Temp["sum_score low"].max()

for I in range(len(VirDF_I)):
    VirDF_I['high normed on protein'][I]= VirDF_I['sum_score high'][I]/Dict[VirDF_I['Protein'][I]][0]
    VirDF_I['low normed on protein'][I]= VirDF_I['sum_score low'][I]/Dict[VirDF_I['Protein'][I]][1]

###===Saving MHC_I data===###
VirDF_I.to_csv("../Data/Affinity_results/Processed_3/"+VirName+"_I.csv", index= True)

Proteins= Proteins[:-1]

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
        
        ###===PEPMatch===###
        
        # only high binding
        add= ""
        Pep= list(Panda10['Peptide'])
        if len(Pep) != 0:
            for z in Pep:
                add+= "> peptide\n"
                add+= z+"\n"
            
            f= open("tmp/test_pep.fasta", "w")
            f.write(add)
            f.close()
        
            Match= Matcher('tmp/test_pep.fasta', '../Data/human.fasta', 0, 5, output_format= 'dataframe').match()

            Out= Match.loc[Match['Index start'].notna()]
            Ptd= list(Out['Query Sequence'])
            for z in Ptd:
                VHLPP+= "Allele: "+x+" Protein: "+y+" Peptide: "+z+"\n"
        
            Match= Match.loc[Match['Index start'].isnull()]
            Ptd= list(Match['Query Sequence'])
            Panda10= Panda10.loc[Panda10['Peptide'].isin(Ptd)]

        ###===Adding rows to VirDF_II===###
        Row= [None for f in range(9)]
        Row[0]= VirName
        Row[1]= x
        Row[2]= y
        Row[3]= len(Panda10)
        Row[4]= len(Panda50)
        Row[5]= 1*Row[3]-Panda10['%Rank_EL'].sum()
        Row[6]= 5*Row[4]-Panda50['%Rank_EL'].sum()
        New= pd.DataFrame(data= [Row], columns= ColNames)
        VirDF_II= pd.concat([VirDF_II, New], axis= 0)

###===Adding "Raw_Proteome" values===###
Proteins+= ["Raw_Proteome"]

for x in MHCII:
    Temp= VirDF_II.loc[VirDF_II['MHC'] == x]
    Row= [None for f in range(9)]
    Row[0]= VirName
    Row[1]= x 
    Row[2]= "Raw_Proteome"
    Row[3]= Temp['# of high binding'].sum()
    Row[4]= Temp['# of low binding'].sum()
    Row[5]= Temp['sum_score high'].sum()
    Row[6]= Temp['sum_score low'].sum()
    New= pd.DataFrame(data= [Row], columns= ColNames)
    VirDF_II= pd.concat([VirDF_II, New], axis= 0)

###===Adding normed values===###
VirDF_II.reset_index(inplace= True)
VirDF_II.drop("index", axis= 1, inplace= True)

Dict= {}
for x in Proteins:
    Dict[x]= [None, None]
    Temp= VirDF_II.loc[VirDF_II['Protein'] == x]
    Dict[x][0]= Temp["sum_score high"].max()
    Dict[x][1]= Temp["sum_score low"].max()

for I in range(len(VirDF_II)):
    VirDF_II['high normed on protein'][I]= VirDF_II['sum_score high'][I]/Dict[VirDF_II['Protein'][I]][0]
    VirDF_II['low normed on protein'][I]= VirDF_II['sum_score low'][I]/Dict[VirDF_II['Protein'][I]][1]

###===Saving MHC_II data===###
VirDF_II.to_csv("../Data/Affinity_results/Processed_3/"+VirName+"_II.csv", index= True)

###===Write to file with human-like peptides===###
f= open( "../Data/Affinity_results/Processed_3/human-like_peptides.txt", "a")
f.write(VHLPP)
f.close()

