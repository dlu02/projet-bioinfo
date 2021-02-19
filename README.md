# Projet de Bio-Informatique

- Télécharger la base de données RNANet : https://evryrna.ibisc.univ-evry.fr/evryrna/rnanet (SQLLite3 Database)  
-  Créer un dossier pour le projet, appelons le "Projet"  
-  Dans "Projet", créer deux dossiers : "results" et "src"  
-  Placer le script "projet.py" dans "src"
-  Placer le fichier "RNANet.db" dans "results"  
-  Exécuter le script "projet.py"  
-  10 fichiers au format .csv ont été générés dans "Projet/src" 


## Instructions à suivre
- Créer un dossier pour le projet, appelons le "Projet"
- Dans "Projet", créer un dossier : "src"
- Dans "src", créer un dossier : "data"
- Placer le script "projet.py" dans "src"
- Extraite l'archive contenant les 17120 fichiers au format CSV dans "data"
- Installer les modules python suivants :
        - pip install python-igraph
        - pip install pycairo

## TODO
- couleurs à ajouter sur le graphe obtenu en fonction du type des liaisons (12 couleurs pour 12 liaisons -majorité en bleu (cWW)
 --> fonction auxiliaire à appeler avant draw_graph_from_csv

- recherche de sous-graphe : convertir le motif en graphe puis rechercher celui-ci dans le graphe obtenu

- pb de l'image renvoyée deux fois
