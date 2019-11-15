#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 16:10:51 2019

@author: erlend
"""

import gudhi as gd
import networkx as nx
import matplotlib.pyplot as plt


# cplx is the input simplicial complex, filtrated should be true if cplx is filtered in a reasonable way, print dualization prints the dictionaries for the vertices 
def dualize(cplx : gd.SimplexTree(), filtrated : bool, birth : tuple, death : tuple, printDualization : bool, inftyVertex : bool):
    
    # Sets inftyVertex to true if we want the filtrated case and also retrives filtration values for birthSimplex and deathSimplex
    if filtrated:
        inftyVertex = True
        (bSimplex, bFilt) = birth
        (dSimplex, dFilt) = death
        
        
    
    # Initializing the dual complex
    dualCplx = gd.SimplexTree()
    
    # Initializes the dictionary for finding name of a simplex given the name of a dual simplex:
    fDualDict = {}
    
    # Initializes the dictionary for finding name of a dual simplex given the name of a simplex:
    fCplxDict = {}                          
    
    
    # Adds the infinity vertex if inftyVertex is set to true
    if inftyVertex:                        
        dualCplx.insert([0])
        fDualDict[(0)] = [0]
        fCplxDict[(0,)] = [0]
    
    # A list of tuples of the simplices and their filtration value
    simplices= cplx.get_filtration()        #
    
    # Stores the dimension of the input simplicial complex
    dimension = cplx.dimension()            
    
    # Initializes the capacity dictionary giving flow capacities
    capDict={}
    
    # Initializes a counter i to count and enumerate the dual simplices
    i = 1
    
    
    # A loop that adds vertices for every d simplex in the input complex
    for (simplex, filtration) in simplices:
        
        # Picks out the maximal simplices (facets)
        if len(simplex) == dimension +1:
            
            # Inserts a 1 simplex in the dual complex with the filtration value equal to the filtration value of the input simplex. Note that the filtration values are not valid.
            dualCplx.insert([i],filtration)
            
            # Puts into the dictionary that the vertex i is the dual of the current simplex
            fDualDict[(i)]= simplex
            
            # Puts into the dictionary that the dual of the current simplex is the node i
            fCplxDict[tuple(simplex)]=[i]
            
            #
            if filtrated:
                if dSimplex < filtration:
                    dualCplx.insert( [0,i] )
                    capDict[(0,i)] = 100000000000000000
            
            # Increases the counter
            i += 1
    
    
    # Prints the dualization 
    if printDualization:
        for j in fDualDict:
            print("The vertex " + str([j]) + " is the dual of the simplex " + str(fDualDict[j]) + " in the original simplicial complex")
    
        for j in fCplxDict:
            print("The simplex " + str(list(j))+ " is mapped to the vertex " + str(fCplxDict[j]) + " in the dual simplicial complex")
    
    # A for-loop running that adds the edges (later: triangles, tetrahedra etc.) to the dual complex
    for (simplex, filtration) in simplices:
        
        # Picks out the simplices of dimension d-1       
        if len(simplex) == dimension:
            
            # Gets a list of the cofaces of codimension 1
            cofaces = cplx.get_cofaces(simplex,1)
            
            # As per the algorithm we ignore simplices that have 0 cofaces of codimension 1
            
            # The simplices that have 1 coface of codimension 1 are connected to the infinity vertex by an edge, if the infinity vertex has been added
            if len(cofaces) == 1 and inftyVertex:
                (coface, coFilt) = cofaces[0]
                dualEdge = [0] + fCplxDict[tuple(coface)]
                capDict[tuple(dualEdge)] = 100000000000000000
                # Inserts the edge into the dual simplex
                dualCplx.insert(tuple(dualEdge))
            
            
            # The simplices that have 2 coface of codimension 1 are coded in the dual complex as an edge between the vertices corresponding to the cofaces
            if len(cofaces) == 2:
                
                # Unpacking the two cofaces
                (cofaceA, coFiltA) = cofaces[0]
                (cofaceB, coFiltB) = cofaces[1]
                
                # Makes a list reresenting the dual edge (to be inserted)
                dualEdge = fCplxDict[tuple(cofaceA)] + fCplxDict[tuple(cofaceB)]
                
                # Sorts the list to be on the safe side
                dualEdge.sort()
                
                # Sets the edge capacity to the apropiate value
                capDict[tuple(dualEdge)] = 1
                if filtrated:
                    if bFilt < filtration:
                        capDict[tuple(dualEdge)] = 100000000000000000
                
                # Inserts the edge into the dual simplex
                dualCplx.insert(tuple(dualEdge))
        
    return (dualCplx, capDict, fDualDict, fCplxDict)

    
# cplx is the input simplicial complex, capDict is the capacity for the edges, in the graph, vertexName is a dictionary that should give the name for the vertices in the input simplicial complex
def toGraph(cplx : gd.SimplexTree , capDict : dict, vertexName : dict):
    
    # Initializes output graph
    graph = nx.Graph()
    
    # Gets list of simplices
    simplices = cplx.get_filtration()
    
    # Adds 1 simplices with labels from the dictionary
    for (simplex, filtration) in simplices:
        if len(simplex) == 1:
            graph.add_node(str(vertexName[simplex[0]]))
            
                
    
    # Adds 2-simplices to graph
    for (simplex, filtration) in simplices:
        if len(simplex) == 2:
           graph.add_edge(str(vertexName[simplex[0]]), str(vertexName[simplex[1]]), capacity = capDict[tuple(simplex)])
            
    return graph


def tograph2(cplx : gd.SimplexTree):
    
    graph = nx.Graph()
    simplices = cplx.get_filtration()
    
    for (simplex, filtration) in simplices:
        if len(simplex) == 2:
            #print("hei")
            #print(simplex)
            graph.add_nodes_from(simplex)
            graph.add_edge(simplex[0],simplex[1])
            
    return graph

def makeTorus(n : int, m :int):
    simp = gd.SimplexTree()
    
    counter = 0
    for i in range(0, n*m):
        simp.insert([counter],counter)
        counter += 1   
        
        

            
    for j in range(0,m):
        for i in range(0,n):
            a = (i + j*n) % (m*n)
            b = ((i+1 % n) + j*n) % (m*n)
            c = (i + ((j+1) % m)*n) % (m*n)
            d = ((i+1 % n) + ((j+1) % m)*n) % (m*n)
            #print((a,b,c,d))
            simp.insert([a,b], counter)
            counter += 1
            simp.insert([a,c], counter)
            counter += 1
            simp.insert([a,d], counter)
            counter += 1
            
    for j in range(0,m):
        for i in range(0,n):
            a = (i + j*n) % (m*n)
            b = ((i+1 % n) + j*n) % (m*n)
            c = (i + ((j+1) % m)*n) % (m*n)
            d = ((i+1 % n) + ((j+1) % m)*n) % (m*n)
            #print((a,b,c,d))
            simp.insert([a,b,d], counter)
            counter += 1
            simp.insert([a,c,d], counter)
            counter += 1

    return simp


def sampleSphere(n, radius : int):
    
    print ("hi")
    

    
def main():
    
    simp = makeTorus(5,100)
    
#    
#    #Torus:
#    simp.insert([1,2,4])
#    simp.insert([2,3,5])
#    simp.insert([3,1,6])
#        
#    simp.insert([4,5,7])
#    simp.insert([5,6,8])
#    simp.insert([6,4,9])
#    
#    simp.insert([7,8,1])
#    simp.insert([8,9,2])
#    simp.insert([9,7,3])
#    
#    simp.insert([4,5,2])
#    simp.insert([5,6,3])
#    simp.insert([6,4,1])
#    
#    simp.insert([7,8,5])
#    simp.insert([8,9,6])
#    simp.insert([9,7,4])
#    
#    simp.insert([1,2,8])
#    simp.insert([2,3,9])
#    simp.insert([3,1,7])
#    
    
#    # 6-Cycle
#    simp.insert([1,2])
#    simp.insert([2,3])
#    simp.insert([3,4])
#    simp.insert([4,5])
#    simp.insert([5,6])
#    simp.insert([6,1])



    
    #for skeleton in simp.get_skeleton(2):
    #    print(skeleton)
    graph = tograph2(simp)
    plt.figure(num=None, figsize=(30, 20), dpi=40, facecolor='w', edgecolor='k')
    
    plt.subplot(221)
    nx.draw(graph,pos = nx.spring_layout(graph), with_labels=False, font_weight='bold', node_color='r', edge_color='b', node_size= 10)
    plt.subplot(222)
    nx.draw(graph,pos = nx.kamada_kawai_layout(graph), with_labels=False, font_weight='bold', node_color='r', edge_color='b',node_size= 10)
    
    dualCplx, capDict, fDualDict, fCplxDict = dualize(simp, False, (0,0) , (0,0), False, True)
    
    graph = toGraph(dualCplx, capDict, fDualDict)
    #graph.add_nodes_from([1,2,3,4])
    #graph = nx.petersen_graph()
    plt.subplot(223)
    nx.draw(graph, pos = nx.spring_layout(graph), with_labels=False, font_weight='bold', node_color='r', edge_color='b',node_size= 10)
    
    plt.subplot(224)
    nx.draw(graph, pos=nx.kamada_kawai_layout(graph),  with_labels=False, font_weight='bold', node_color='r', edge_color='b',node_size= 10)
    plt.style.use('dark_background')
    plt.show()
    
    cut = nx.minimum_cut(graph, str([1,11,12]), str([41,51,52]), capacity='capacity')


if __name__=="__main__":
    main()

print("COMPLETED")
    
    
    
