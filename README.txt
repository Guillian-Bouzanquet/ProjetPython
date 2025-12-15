READ ME.
SAM File Analyzer

Ce programme permet d'analyser des fichiers SAM de grande taille.
Avant de commencer, ce programme fonctionne avec la version 3.X de Python. Vous devez également disposer d'un terminal.

A noter que pour utiliser ce programme, les modules présents dans python seront utilisés :
- argparse 
- sys 
- os 
- re 
- collections (defaultdict, Counter)

Comment initialiser le programme ?
python3 sam_analyzer.py -i <input.sam> [--mapped-only] [--mapq-max N]

Les différentes options disponibles sont :
-i, --input : fichier SAM à analyser (obligatoire)
--mapped-only : ne conserver que les reads entièrement mappés (optionnel)
--mapq-max N : ne conserver que les reads avec score MAPQ inférieur à N (optionnel)


Les données qui sont acceptées au sein du programme sont dans un fichier terminant par une extension .sam
Les fichiers qui ne comportent pas cette extension arrêteront le fonctionnement du programme.

Fichier SAM valide (.sam)

Le programme va vérifier les informations suivantes automatiquement :
- Si le fichier est existant.
- Si l'Extension .sam est bien présente.
- Si le fichier contient des données.

Les fichiers que vous obtiendrez dans le cadre de ce programme : ResumeGlobalDeLAnalyse.txt 
Il contiendra les informations suivantes :
- Le total de reads et reads mappés.
- Le nombre de reads par FLAG.
- Le nombre de reads par chromosome.
- La distribution des scores MAPQ.

D'autres informations optionnelles sont également disponibles comme :
- only_unmapped.fasta (optionnel) : reads non mappés
- summary_unmapped.txt (optionnel) : nombre de reads non mappés
- only_partially_mapped.fasta (optionnel) : reads partiellement mappés
- summary_partially_mapped.txt (optionnel) : nombre de reads partiellement mappés
- Final_Cigar_table.txt (optionnel) : pourcentage global des différents types d’événements dans le CIGAR

License Creative Commons (CC BY 4.0)

Partage et adaptation autorisés, mention de l’auteur obligatoire.

BOUZANQUET Guillian (guillain.bouzanquet@etu.umontpellier.fr)