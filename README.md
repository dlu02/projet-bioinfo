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
- Rédiger un manuel utilisateur de 2 ou 3 pages, accompagnés de capture d'écran de l'exécution attendue
- Ce manuel doit comprendre les parties suivantes :  
        1 - Brève introduction (= but du projet)  
        2 - Installation du projet (1/2 page)  
        3 - Explication rapide des fonctions principales (développer surtout sur la recherche des motifs dans find_subgraph)  
        4 - Limite du code (proposer des améliorations, en terme de modélisation de l'ARN ou même de performance du code)  
        5 - Conclusion (rappeler ce qui fonctionne / ce qui ne fonctionne pas)  
