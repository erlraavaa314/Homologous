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
                #print(simplex)
                #print(cofaces[0])
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


def makePlane(n : int, m : int):
    simp = gd.SimplexTree()
    
    counter = 0
    for i in range(0, n*m):
        simp.insert([counter],counter)
        counter += 1   
            
            
            
    for j in range(0,m-1):
        for i in range(0,n-1):
            a = i + j*n
            b = i+1 + j*n
            c = i + (j+1)*n
            d = i+1 + (j+1)*n
            #print((a,b,c,d))
            simp.insert([a,b], counter)
            counter += 1
            simp.insert([a,c], counter)
            counter += 1
            simp.insert([a,d], counter)
            counter += 1
            simp.insert([b,d], counter)
            counter += 1
            simp.insert([c,d], counter)
            counter += 1
            simp.insert([a,b,d], counter)
            counter += 1
            simp.insert([a,c,d], counter)
        
        
        
    return simp


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

def makeKlein(n : int, m :int):
    simp = makePlane(n,m)
    
    counter = 5*(n*m)^2

    for i in range(0, n-1):
        a = n*(m-1) + i
        b = n*(m-1) + 1 + i
        c = n-1-i
        d = n-2-i
        simp.insert([a,c],counter)
        counter += 1
        simp.insert([a,d],counter)
        counter += 1
        simp.insert([a,b,d],counter)
        counter += 1
        simp.insert([a,c,d],counter)
        counter += 1
        
        #print(a,b,c,d)

        
    for i in range(1, m):
        a = n*i-1
        b = n*(i-1)
        c = n*(i+1) -1
        d = n*(i)
        #print(a,b,c,d)
        simp.insert([a,c],counter)
        counter += 1
        simp.insert([a,d],counter)
        counter += 1
        simp.insert([a,b,d],counter)
        counter += 1
        simp.insert([a,c,d],counter)
        counter += 1

    a = n-1
    b = 0
    c = n*(m-1)
    d = n*m-1
    
    simp.insert([a,c],counter)
    counter += 1
    simp.insert([a,d],counter)
    counter += 1
    simp.insert([a,b,d],counter)
    counter += 1
    simp.insert([a,c,d],counter)
    counter += 1
    return simp



def translate(dictionary : dict, keys):
    setBack = set()
    for element in keys:
        if len(eval(element)) > 1:
            
            #print(element)
            setBack.add(eval(element)[0])
            setBack.add(eval(element)[1])
            setBack.add(eval(element)[2])
    return setBack
        
    
def myPlotter(simp):    
    #for skeleton in simp.get_skeleton(2):

    graph1 = tograph2(simp)
    plt.figure(num=None, figsize=(100, 50), dpi=40, facecolor='black', edgecolor='black')
    positionA1 = nx.spring_layout(graph1)
    positionA2 = nx.kamada_kawai_layout(graph1)

    plt.subplot(421, facecolor='black')
    nx.draw_networkx(graph1, pos = positionA1, with_labels=False, font_weight='bold', node_color='g', edge_color='b', node_size= 100)
    plt.subplot(422, facecolor='black')
    nx.draw_networkx(graph1, pos = positionA2, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    dualCplx, capDict, fDualDict, fCplxDict = dualize(simp, False, (0,0) , (0,0), False, True)
    
    graph = toGraph(dualCplx, capDict, fDualDict)
    #graph.add_nodes_from([1,2,3,4])
    #graph = nx.petersen_graph()
    position1 = nx.spring_layout(graph)
    position2 = nx.kamada_kawai_layout(graph)
    
    
    #for node in graph1.nodes :
        #print(node, graph1.degree(node))
    
    
    plt.subplot(425, facecolor='black')
    nx.draw_networkx(graph, pos = position1, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    plt.subplot(426, facecolor='black')
    nx.draw_networkx(graph, pos = position2,  with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)

    
    cutSize, cut = nx.minimum_cut(graph, str([1,11,12]), str([21,31,32]), capacity='capacity')
    
    cut1, cut2 = cut
    
    
    plt.subplot(427, facecolor='black')
    nx.draw_networkx_nodes(graph,nodelist = list(cut1), pos = position1, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    nx.draw_networkx_nodes(graph,nodelist = list(cut2), pos = position1, with_labels=False, font_weight='bold', node_color='r', edge_color='b',node_size= 100)
    
    nx.draw_networkx_edges(graph, pos = position1, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    
    plt.subplot(428, facecolor='black')
    nx.draw_networkx_nodes(graph,nodelist = list(cut1), pos = position2, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    nx.draw_networkx_nodes(graph,nodelist = list(cut2), pos = position2, with_labels=False, font_weight='bold', node_color='r', edge_color='b',node_size= 100)
    
    nx.draw_networkx_edges(graph, pos = position2, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
   
    cutA1 = translate(fCplxDict, cut1)
    
    cutA2 = translate(fCplxDict, cut2)
    #print(cut1)
    #print(len(cut1))
    
    #print(graph1.nodes())
    plt.subplot(423, facecolor='black')
    
    nx.draw_networkx_nodes(graph1,nodelist = list(cutA1), pos = positionA1, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    nx.draw_networkx_nodes(graph1,nodelist = list(cutA2), pos = positionA1, with_labels=False, font_weight='bold', node_color='r', edge_color='b',node_size= 100)
    
    nx.draw_networkx_edges(graph1, pos = positionA1, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    
    plt.subplot(424, facecolor='black')
    nx.draw_networkx_nodes(graph1, nodelist = list(cutA1), pos = positionA2, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    nx.draw_networkx_nodes(graph1, nodelist = list(cutA2), pos = positionA2, with_labels=False, font_weight='bold', node_color='r', edge_color='b',node_size= 100)
    
    nx.draw_networkx_edges(graph1, pos = positionA2, with_labels=False, font_weight='bold', node_color='g', edge_color='b',node_size= 100)
    
    #for node in graph.nodes :
    #    print(graph.degree(node))
    plt.show()


def main():
    
    
    simp1 = makePlane(10,20)
    simp2 = makeTorus(10,20)
    simp3 = makeKlein(10,20)
    
    myPlotter(simp1)
    myPlotter(simp2)
    myPlotter(simp3)
    

if __name__=="__main__":
    main()

print("COMPLETED")
    
    
    
