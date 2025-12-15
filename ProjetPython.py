#!/usr/bin/python3
#-*- coding : utf-8 -*-

############### IMPORT MODULES ###############

import argparse                                                    #Permet de gerer les options en ligne de commande.
import sys                                                         #Lien entre l'interpreteur Python et des commandes specifiques.
from collections import Counter, defaultdict                       #Compteur pour les lignes de commande.
import re                                                          #Travaille avec des expressions regulieres 
import os                                                          #Interactions entre python et le systeme d'exploitation.


############### FUNCTIONS TO :

## 1/ Check, 

def checkSAM(file): #Verification de l'existence du fichier SAM dans le doc. On verifie dans l'ordre si le chemin est correct, puis l'extension, et on voit si le fichier est vide ou non.
    if not os.path.isfile(file):
        print(f"[ERREUR] Le fichier nomme '{file}' n'existe pas.")
        sys.exit(1)

    if not file.endswith(".sam"):
        print(f"[ERREUR] Le fichier nomme '{file}' doit etre un fichier possedant une extension en .sam")
        sys.exit(1)

    if os.path.getsize(file) == 0:
        print(f"[ERREUR] Le fichier nomme '{file}' est vide.") 
        sys.exit(1) 

    f = open(file, 'r')

    for line in f:
        if line.startswith("@"):
            continue  # Ignore les lignes d'en-tête

        #Analyse de la ligne SAM
        cols = line.rstrip().split("\t")

        #Vérification du nombre de colonnes
        if len(cols) < 11:
            print("[ERREUR] Attention, le nombre de colonnes dans ce fichier est incorrect. Il doit etre au minimum de 11.")
            f.close()
            sys.exit(1)

        #Vérification que le FLAG est un entier
        if not cols[1].isdigit():
            print("[ERREUR] Le FLAG (colonne 2) doit etre un entier.")
            f.close()
            sys.exit(1)

    f.close()
    print("Votre fichier SAM est valide et pret a etre utilise !")

#On passe dès à présent à la lecture du fichier SAM.

## 2/ Read,

def readSAM(file): #Lit un fichier SAM et ignore les lignes d'en-tete commencant par @.
   
    sam_lines = [] #Initialise une liste.

    with open(file, "r") as f:  #Ouverture du fichier en mode lecture
        for line in f:
            if line.startswith("@"): #Si la ligne commence par un @, alors on peut supprimer ce premier caractere pour de la visibilite.
                continue
            sam_lines.append(line.rstrip("\n"))

    return sam_lines

## 3/ Store,

def storeSAM(sam_lines):

    flag_counts = defaultdict(int) #On compte desormais les flags, puis les chromosomes, puis les mapq, on s'assure que les valeurs soient entières avec la fonction int.
    chrom_counts = defaultdict(int)
    mapq_counts = defaultdict(int)

    total_reads = 0  #On definit ces variables a la valeur 0, pour les reads et pour les reads mappes. 
    mapped_reads = 0

    for line in sam_lines:
        cols = line.split("\t")   #Le symbole \t correspond a un separateur de tabulation.

        flag = int(cols[1])
        chrom = cols[2]             #Chque caracteristique est assigne a une colonne.
        mapq = int(cols[4])

        total_reads += 1            #On incremente ici chacune des valeurs successivement.
        flag_counts[flag] += 1
        mapq_counts[mapq] += 1

        if not (flag & 4):          #Read mappe si le bit 4 du flag n'est pas actif.
            mapped_reads += 1
            chrom_counts[chrom] += 1
        
    return total_reads, mapped_reads, flag_counts, chrom_counts, mapq_counts    #On retourne toutes les valeurs.

## 4/ Analyse,

#### Conversion du flag en binaire ####
def flagBinary(flag) :

    flagB = bin(int(flag)) # Transform the integer into a binary.
    flagB = flagB[2:] # Remove '0b' Example: '0b1001101' > '1001101'
    flagB = list(flagB) 
    if len(flagB) < 12: # Size adjustement to 12 (maximal flag size)
        add = 12 - len(flagB) # We compute the difference between the maximal flag size (12) and the length of the binary flag.
        for t in range(add):
            flagB.insert(0,'0') # We insert 0 to complete until the maximal flag size.
    return flagB


#### Analyze the unmapped reads (not paired) ####
def unmapped(sam_line):
    
    unmapped_count = 0
    with open ("only_unmapped.fasta", "a+") as unmapped_fasta, open("summary_unmapped.txt", "w") as summary_file: 
        for line in sam_line:
            col_line = line.split("\t")  #On découpe les colonnes
            flag = flagBinary(col_line[1])

            if int(flag[-3]) == 1:     
                unmapped_count += 1
                unmapped_fasta.write(tostring(line))

        summary_file.write("Total unmapped reads: " + str(unmapped_count) + "\n") 
        return unmapped_count

#### Analyze the partially mapped reads ####
def partiallyMapped(sam_line):
    
    partially_mapped_count = 0

    with open ("only_partially_mapped.fasta", "a+") as partillay_mapped_fasta, open("summary_partially_mapped.txt", "w") as summary_file:
        for line in sam_line:
            col_line = line.split("\t")
            flag = flagBinary(col_line[1]) # We compute the same 

            if int(flag[-2]) == 1: 
                if col_line[5] != "100M":
                    partially_mapped_count += 1
                    partillay_mapped_fasta.write(tostring(line))

        summary_file.write("Total partially mapped reads: " + str(partially_mapped_count) + "\n") 
        return partially_mapped_count


def filterSAM(sam_lines, mapped_only=False, mapq_max=None):
   
    filtered = []
    for line in sam_lines:
        cols = line.split("\t")
        flag = int(cols[1])
        mapq = int(cols[4])

        # Ignore unmapped reads si mapped_only
        if mapped_only and (flag & 4):
            continue

        # Ignore reads avec MAPQ >= mapq_max
        if mapq_max is not None and mapq >= mapq_max:
            continue

        filtered.append(line)

    return filtered

### Analyse the CIGAR = regular expression that summarise each read alignment ###
def readCigar(cigar): 
   
    ext = re.findall('\w',cigar) # split cigar 
    key=[] 
    value=[]    
    val=""

    for i in range(0,len(ext)): # For each numeric values or alpha numeric
        if (ext[i] == 'M' or ext[i] == 'I' or ext[i] == 'D' or ext[i] == 'S' or ext[i] == 'H' or ext[i] == "N" or ext[i] == 'P'or ext[i] == 'X'or ext[i] == '=') :
            key.append(ext[i])
            value.append(val)
            val = ""
        else :
            val = "" + val + ext[i]  # Else concatenate in order of arrival
    
    dico = {}
    n = 0
    for k in key:   # Dictionnary contruction in range size lists              
        if k not in dico.keys():    # for each key, insert int value
            dico[k] = int(value[n])   # if key not exist, create and add value
            n += 1
        else:
            dico[k] += int(value[n])  # inf key exist add value
            n += 1
    return dico

### Analyse the CIGAR = regular expression that summarise each read alignment ###
def percentMutation(dico):
        
    totalValue = 0 # Total number of mutations
    for v in dico :
        totalValue += dico[v]

    mutList = ['M','I','D','S','H','N','P','X','=']
    res = ""
    for mut in mutList : # Calculated percent of mutation if mut present in the dictionnary, else, percent of mut = 0
        if mut in dico.keys() :
            res += (str(round((dico[mut] * 100) / totalValue, 2)) + ";")
        else :
            res += ("0.00" + ";")
    return res

def globalPercentCigar():
    """
      Global representation of cigar distribution.
    """
    
    with open ("outpuTable_cigar.txt","r") as outpuTable, open("Final_Cigar_table.txt", "w") as FinalCigar:
        nbReads, M, I, D, S, H, N, P, X, Egal = [0 for n in range(10)]

        for line in outpuTable :
            mutValues = line.split(";")
            nbReads += 2
            M += float(mutValues[2])+float(mutValues[12])    #Aligment Match
            I += float(mutValues[3])+float(mutValues[13])    #Insertion
            D += float(mutValues[4])+float(mutValues[14])    #Deletion
            S += float(mutValues[5])+float(mutValues[15])    #Soft Clipping (non mappe mais conservees)
            H += float(mutValues[6])+float(mutValues[16])    #Hard Clipping (non mappe et non conservees)
            N += float(mutValues[7])+float(mutValues[17])    #Skipped (region non comptee)
            P += float(mutValues[8])+float(mutValues[18])    #Placeholder
            X += float(mutValues[9])+float(mutValues[19])    #Sequence mismatch
            Egal += float(mutValues[10])+float(mutValues[20]) #Sequence Match

        FinalCigar.write("Global cigar mutation observed :"+"\n"                   #Comptage de chacun de différents évènements survenus sur le read.
                        +"Alignlent Match : "+str(round(M/nbReads,2))+"\n"
                        +"Insertion : "+str(round(I/nbReads,2))+"\n"
                        +"Deletion : "+str(round(D/nbReads,2))+"\n"
                        +"Skipped region : "+str(round(S/nbReads,2))+"\n"
                        +"Soft Clipping : "+str(round(H/nbReads,2))+"\n"
                        +"Hard Clipping : "+str(round(N/nbReads,2))+"\n"
                        +"Padding : "+str(round(P/nbReads,2))+"\n"
                        +"Sequence Match : "+str(round(Egal/nbReads,2))+"\n"
                        +"Sequence Mismatch : "+str(round(X/nbReads,2))+"\n")


#### Summarise the results ####

def Summary(total_reads, mapped_reads, flag_counts, chrom_counts, mapq_counts): #On ecrit le resume de notre analyse.
    
    s = open("ResumeGlobalDeLAnalyse.txt", "w")

    s.write("===== RESUME =====\n")
    s.write("Total reads : " + str(total_reads) + "\n")        #Notes : \n correspond à un saut de ligne.
    s.write("Mapped reads : " + str(mapped_reads) + "\n\n") 

    s.write("----- Reads par FLAG -----\n")
    for fl in sorted(flag_counts.keys()):
        s.write("FLAG " + str(fl) + " : " + str(flag_counts[fl]) + "\n") #Nombre de flags comptés

    s.write("\n----- Reads par chromosome -----\n")
    for chrom in sorted(chrom_counts.keys()):
        s.write(str(chrom) + " : " + str(chrom_counts[chrom]) + "\n")    #Nombre de chromosomes comptés

    s.write("\n----- MAPQ distribution -----\n")
    for q in sorted(mapq_counts.keys()):
        s.write("MAPQ " + str(q) + " : " + str(mapq_counts[q]) + "\n")  #MAPQ comptés

    s.close()
    print("ResumeGlobalDeLAnalyse.txt")
   
#### Main function ####

def main():

    parser = argparse.ArgumentParser(description="SAM File Analyzer")

    # Argument obligatoire
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input SAM file"
    )

    # Filtre : reads entièrement mappés
    parser.add_argument(
        "--mapped-only",
        action="store_true",
        help="Ne conserver que les reads entièrement mappés"
    )

    # Filtre : qualité de mapping
    parser.add_argument(
        "--mapq-max",
        type=int,
        help="Ne conserver que les reads avec MAPQ inférieur à cette valeur"
    )

    # Analyse des arguments (TOUJOURS À LA FIN)
    args = parser.parse_args()
    # 1) CHECK
    checkSAM(args.input)

    # 2) READ
    sam_lines = readSAM(args.input)
    sam_lines = filterSAM(sam_lines, mapped_only=args.mapped_only, mapq_max=args.mapq_max)


    # 3) STORE
    total_reads, mapped_reads, flag_counts, chrom_counts, mapq_counts = storeSAM(sam_lines)

    # 4) ANALYSE
    Summary(total_reads, mapped_reads, flag_counts, chrom_counts, mapq_counts)

############### LAUNCH THE SCRIPT ###############

if __name__ == "__main__":
    main()                
