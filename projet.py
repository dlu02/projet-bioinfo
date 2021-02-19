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


# fonction d'affichage du graphe, à partir du nom d'un fichier au format CSV

def draw_graph_from_csv(filename):
    
    df = pandas.read_csv("data/" + filename)[['index_chain','nt_code','paired','pair_type_LW']]
    
    print(df.to_string())       # affichage dans la console du fichier sous forme de dataset
    
    g = ig.Graph()              # initialisation du graph
    
    for node in df.iterrows():
        g.add_vertices(1)       # ajout des noeuds (1 nucléotide = 1 noeud)
        
    for i in range(len(g.vs)):
        g.vs[i]["label"]= str(i+1)  # label du noeud, correspondant à la valeur index_chain 
        
        # coloration des noeuds selon le nucléotide représenté
        
        if (df['nt_code'][i] == "A"):
            g.vs[i]["color"] = "cyan"
        if (df['nt_code'][i] == "C"):
            g.vs[i]["color"] = "yellow"
        if (df['nt_code'][i] == "G"):
            g.vs[i]["color"] = "green"
        if (df['nt_code'][i] == "U"):
            g.vs[i]["color"] = "pink"
            
    for i in range(len(g.vs)):
        
        currentNode = df['index_chain'][i]      # noeud parcouru actuellement
        pairedNode = str(df['paired'][i])       # liste des nucléotides possédant une interaction avec le noeud actuel
        
        # liaisons phosphodiester entre les nucléotides
        
        if i != len(g.vs)-1:    
            g.add_edges([(i,i+1)]) 

        # interactions canoniques entre certains nucléotides

        if pairedNode.find(",") != -1:                       # la liste "paired" contient plus d'un élément
            arrayPairedNodes = pairedNode.split(',')         # on récupère chacun des nucléotides de la liste
            for node in arrayPairedNodes:
                if int(node) in range (1,len(g.vs)):         # le nucléotide apparié existe bien dans la chaîne
                    if not g.are_connected(int(node)-1,currentNode-1):  # si l'arête n'existe pas déjà dans l'autre sens, alors on la crée
                        g.add_edges([(currentNode-1,int(node)-1)])      # ajout dans le graphe de la liason canonique
                        res_edge = g.es.find(_source=currentNode-1, _target=int(node)-1)
                        res_edge["color"] = "blue"            # mise en rouge des interactions canoniques
        
        else:
            pairedNode = pandas.to_numeric(df['paired'][i])  # cast d'un string en int, afin de vérifier si le nucléotide n'est pas apparié (NaN)
            if not math.isnan(pairedNode) and int(pairedNode) in range (1,len(g.vs)):
                if not g.are_connected(int(pairedNode)-1,currentNode-1):
                    g.add_edges([(currentNode-1,int(pairedNode)-1)])
                    res_edge = g.es.find(_source=currentNode-1, _target=int(pairedNode)-1)
                    res_edge["color"] = "blue"                # mise en rouge des interactions canoniques

    visual_style = {}
    visual_style["edge_width"] = 3
    ig.plot(g, **visual_style).show()       # affichage du graphe, sous forme de fichier .png dans une fenêtre externe
    
    
draw_graph_from_csv("1asz_1_S")     # test de la fonction principale sur un fichier donné