# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:51:07 2021

@author: REMI DECOUTY, DAMIEN LU
"""

import pandas
import glob, os
import igraph as ig
import math


# récupération de la liste de tous les fichiers à étudier

os.chdir("data")
allCsvFiles = glob.glob("*")
os.chdir("..")

# coloration d'une arête, en fonction du type de liaison

def color_edge(node, index):
    
    arrayPairType = []
    
    if (index != -1):                    
        arrayPairType = node.split(',')  
        
    cases = {
        "cWW": "blue",
        
        "tWW" : "red",
        
        "cWH" : "orange",
        "cHW" : "orange",
        
        "tWH": "yellow",
        "tHW": "yellow",

        "cWS" : "brown",
        "cSW" : "brown",

        "tWS": "magenta",
        "tSW": "magenta",
        
        "cHH": "cyan",
        
        "tHH": "purple",

        "cHS": "pink",
        "cSH": "pink",
        
        "tHS" : "green",
        "tSH" : "green",

        "cSS" : "maroon",
        
        "tSS": "gray",
    }
            
    if (index != -1):                    
        return (cases.get(arrayPairType[index]))

    else:    
        return (cases.get(node))


# fonction d'affichage du graphe, à partir du nom d'un fichier au format CSV

def draw_graph_from_csv(file):
    
    global filename
    filename = file
    df = pandas.read_csv("data/" + filename)[['old_nt_resnum','paired','pair_type_LW']]
    
    # print("\nDataset du fichier " + file + "\n\n" + df.to_string() + "\n\n")      # affichage dans la console du fichier sous forme de dataset
    
    global g                    # déclaration d'une variable g globale
    g = ig.Graph()              # initialisation du graphe
    
    for node in df.iterrows():
        g.add_vertices(1)       # ajout des noeuds (1 nucléotide = 1 noeud)
        
    for i in range(len(g.vs)):
        g.vs[i]["label"]= str(i+1)  # label du noeud, correspondant à la valeur index_chain 
        
            
    for i in range(len(g.vs)):
        
        currentNode = df['old_nt_resnum'][i]      # noeud parcouru actuellement
        pairedNode = str(df['paired'][i])       # liste des nucléotides possédant une interaction avec le noeud actuel
        
        
        # liaisons phosphodiester entre les nucléotides
        
        if i != len(g.vs)-1:    
            g.add_edges([(i,i+1)]) 

        # interactions canoniques entre certains nucléotides

        if pairedNode.find(",") != -1:                       # la liste "paired" contient plus d'un élément
            arrayPairedNodes = pairedNode.split(',')         # on récupère chacun des nucléotides de la liste
            for node in arrayPairedNodes:
                if (int(node) in df['old_nt_resnum'].tolist()):         # le nucléotide apparié existe bien dans la chaîne
                    if not g.are_connected(int(node)-1,currentNode-1):  # si l'arête n'existe pas déjà dans l'autre sens, alors on la crée
                        g.add_edges([(currentNode-1,int(node)-1)])      # ajout dans le graphe de la liason canonique
                        res_edge = g.es.find(_source=currentNode-1, _target=int(node)-1)
                        res_edge["color"] = color_edge(df['pair_type_LW'][i], arrayPairedNodes.index(node))  
                        
                        
        else:
            pairedNode = pandas.to_numeric(df['paired'][i])  # cast d'un string en int, afin de vérifier si le nucléotide n'est pas apparié (NaN)
            if not math.isnan(pairedNode) and int(pairedNode) in df['old_nt_resnum'].tolist():
                if not g.are_connected(int(pairedNode)-1,currentNode-1):
                    g.add_edges([(currentNode-1,int(pairedNode)-1)])
                    res_edge = g.es.find(_source=currentNode-1, _target=int(pairedNode)-1)
                    res_edge["color"] = color_edge(df['pair_type_LW'][i],-1)              


    visual_style = {}
    visual_style["edge_width"] = 3
    g.vs["color"] = "white"
    
    # draw_g = ig.plot(g, **visual_style) # affichage du graphe, sous forme de fichier .png dans une fenêtre externe
    
    # si l'affichage n'a pas fonctionné (booléen is_dirty == True), alors on appelle la méthode show()
    # Ce cas de figure peut apparaitre si le code source est éxécuté via un IDE
    
    # if (str(draw_g._is_dirty) == "True"):   
    #      draw_g.show()
    
 
    
# find_subgraph : fonction de recherches des sous-graphes, à partir d'un motif issu de carnaval  
# Le principe est le suivant : la fonction get_subisomorphisms_vf2 de la librairie igraph renvoie 
# la liste des sous-graphes isomorphes entre deux graphes. De plus, la fonction de comparaison compare_edges
# va vérifier pour chaque arête 2 à 2 si leur couleur est identique ou non. 
# Si tel est le cas, alors on aura trouvé un sous-graphe correspondant au motif
# Enfin, subgraph_list, qui contient la liste de tous les sous-graphes trouvés, peut contenir des doublons
# Ex: [0,2,4,1] et [4,2,0,1]
# Pour éviter cela, on va filter subgraph_list dans un set, pour supprimer tous ces doublons
    
def find_subgraph(graph, motif, motif_name):
    
    def compare_edges(g1, g2, i1, i2):
        try:
            result = (g1.es[i1]['color'] == g2.es[i2]['color'])
        except:
            return False
        else:
            return result

    subgraph_list = g.get_subisomorphisms_vf2(motif,edge_compat_fn=compare_edges)
    
    if len(subgraph_list) >= 1:
        results = [tuple(x) for x in set(map(frozenset, subgraph_list))]
        print ("Le motif " + motif_name + " est présent " + str(len(results)) + " fois dans la chaîne d'ARN " + filename + "\n")
        for result in sorted(results):
            print (result)
        print ("\n\n")
        return len(results)

    else:
        print ("Le motif " + motif_name + " est absent de la chaîne d'ARN " + filename + "\n\n")
        return 0


# initialisation des motifs RIN issus de Carnaval, sous forme de graphe
# la couleur des arêtes permet d'identiier le type de liaison pair_type_LW
    
def transorm_RIN_to_graph():
    
    global rin_23, rin_129
    
    rin_23 = ig.Graph()
    rin_129 = ig.Graph()

    rin_129.add_vertices(4)
    rin_129.add_edges([(0,1),(1,2),(2,3),(0,3)])
    rin_129.vs["color"] = "white"

    for i in range(len(rin_129.vs)):
        rin_129.vs[i]["label"]= str(i+1) 
        
    res_edge = rin_129.es.find(_source=1, _target=2)
    res_edge["label"] = "tWW"        # liaison tWW
    res_edge["color"] = "red"

    res_edge = rin_129.es.find(_source=0, _target=3)
    res_edge["label"] = "tHS"        # liaison tHS    (d'après la nomenclature officielle)
    res_edge["color"] = "green"
    
    
    rin_23.add_vertices(8)
    rin_23.add_edges([(0,1),(1,2),(2,3),(3,4),(4,5),(5,7),(1,6),(6,7),(1,4),(0,5)])
    rin_23.vs["color"] = "white"

    for i in range(len(rin_23.vs)):
        rin_23.vs[i]["label"]= str(i+1) 
        
    res_edge = rin_23.es.find(_source=0, _target=5)
    res_edge["label"] = "cWW"
    res_edge["color"] = "blue"

    res_edge = rin_23.es.find(_source=1, _target=4)
    res_edge["label"] = "cWW"
    res_edge["color"] = "blue"
    
    res_edge = rin_23.es.find(_source=2, _target=3)
    res_edge["label"] = "cWW"
    res_edge["color"] = "blue"
    
    res_edge = rin_23.es.find(_source=1, _target=6)
    res_edge["label"] = "cSS"
    res_edge["color"] = "maroon"
    
    res_edge = rin_23.es.find(_source=5, _target=7)
    res_edge["label"] = "tSS"
    res_edge["color"] = "gray"
    
    # draw_rin_23 = ig.plot(rin_23)    # affichage du RIN 23
    # draw_rin_129 = ig.plot(rin_129)  # affichage du RIN 129
    
    # if (str(draw_rin_23._is_dirty) == "True"):   
    #      draw_rin_23.show()
    
    # if (str(draw_rin_129._is_dirty) == "True"):   
    #     draw_rin_129.show()

transorm_RIN_to_graph()
draw_graph_from_csv("1mms_1_C")    # test de la fonction principale sur le fichier nommé "1mms_1_C"
find_subgraph(g,rin_23,"rin_23")
find_subgraph(g,rin_129,"rin_129")    


# Si vous souhaitez parcourir tous les fichiers pour rechercher un motif, il suffit de décommenter
# le code ci-dessous. Mais attention, il faudra absolument mettre en commentaire les divers affichages dans 
# draw_graph_from_csv (plot de la chaine d'ARN + affichage du .csv dans la console), au risque de faire planter 
# le programme. Il est déconseillé d'exécuter le programme ci-dessous si votre machine n'est pas très puissante,
# le programme mettra plusieurs dizaines de minutes à s'éxécuter. 
# Il faudra alors privilégier l'analyse fichier par fichier comme ci-dessus.

cpt_23 = 0
cpt_129 = 0

for csvFile in allCsvFiles:
    draw_graph_from_csv(csvFile)
    cpt_23 = cpt_23 + find_subgraph(g,rin_23,"rin_23")
    cpt_129 = cpt_129 + find_subgraph(g,rin_129,"rin_129")

print("Nombre d'occurences de RIN 23 :" + str(cpt_23))
print("Nombre d'occurences de RIN 129 :" + str(cpt_129))