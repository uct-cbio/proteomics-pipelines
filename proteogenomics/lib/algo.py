#!usr/bin/env python

import math
import numpy as np
from scipy import stats

# Statistics

def log2(val):
    try:
        return math.log(val,2)
    except:
        return np.nan

def impute(df, col, width = 0.3, downshift = 1.8):
    values = df[col].dropna()
    missing = df[col].isnull()
    assert len(values)/num_missing >= 2  # require min 50 % sample not missin
    std = np.std(values, ddof = 1, dtype=np.float64)       # need to subtract ddof?
    mean =  np.mean(values, dtype=np.float64)
    shifted = mean - (downshift * std) 
    for row in df.iterrows():
        ind = row[0]
        if df.loc[ind, col] == np.nan:
            imputed = np.random.normal(loc=shifted, scale = width * std, size = 1)
            df.loc[ind, col] = imputed
    return df
    
def welch(arr1, arr2):
    return stats.ttest_ind(arr1, arr2, equal_var = False)

# Trie graph
class Vertex:
    def __init__(self,key):
        self.id = key
        self.connectedTo = {}
        
    def addNeighbor(self,nbr,edge_label=None):
        self.connectedTo[nbr] = edge_label
        
    def __str__(self):
        return str(self.id) + ' connectedTo: '+ str([x.id for x in self.connectedTo])

    def getConnections(self):
        return self.connectedTo.keys()

    def getId(self):
        return self.id

    def getEdge(self,nbr):
        return self.connectedTo[nbr]

class Graph:
    def __init__(self):
        self.vertList = {}
        self.numVertices = 0

    def addVertex(self,key):
        self.numVertices = self.numVertices + 1
        newVertex = Vertex(key)
        self.vertList[key] = newVertex
        return newVertex

    def getVertex(self,n):
        if n in self.vertList:
            return self.vertList[n]
        else:
            return None

    def __contains__(self,n):
        return n in self.vertList

    def addEdge(self,f,t,cost=0):
        if f not in self.vertList:
            nv = self.addVertex(f)
        if t not in self.vertList:
            nv = self.addVertex(t)
        self.vertList[f].addNeighbor(self.vertList[t], cost)

    def getVertices(self):
        return self.vertList.keys()

    def __iter__(self):
        return iter(self.vertList.values())

def trie_graph(lst):
    assert not isinstance(lst, str)    # make sure that the items you are making a trie_graph from are in list format
    trie = Graph()
    nodes = 0
    trie.addVertex(nodes)
    for word in lst:
        node = 0
        word = word.lstrip().rstrip()+'#' # separates substrings of longer strings also in list
        for i in range(0,len(word)):
            symbol = word[i]
            found = False
            if trie.__contains__(node) == False:
                trie.addEdge(node, nodes, symbol)
                nodes = nodes + 1
                node = nodes
            else:
                for neighbour in trie.getVertex(node).getConnections(): 
                    if found != True:
                        if symbol == Vertex.getEdge(trie.getVertex(node), neighbour):
                            node = neighbour.getId()
                            found = True
                if found == False:
                    nodes = nodes + 1
                    trie.addEdge(node, nodes, symbol)
                    node = nodes
    return trie

def prefix_trie_match(trie, string, start): 
    v = 0
    i = start
    pep = []
    last_edge = None 
    while i < len(string) + 1:
        found = False 
        if len(trie.getVertex(v).getConnections()) == 0:
            return start, i
        if i == len(string) + 1:
            return
        elif i <len(string):
            symbol = string[i]
        
            for child in trie.getVertex(v).getConnections():
                edge = Vertex.getEdge(trie.getVertex(v), child)
                if edge == '#':
                    last_edge = start, i
                if symbol == edge:
                    v = child.getId()
                    i += 1
                    pep.append(edge)
                    found = True
                    break
        if found == False:
            for child in trie.getVertex(v).getConnections():
                edge = Vertex.getEdge(trie.getVertex(v), child) 
                if edge == '#':
                    last_edge = start, i
            return last_edge

def trie_upper(Trie, Text):
    positions = []
    start = 0
    str = []
    pos_dct = {}
    while start < len(Text):
        val = prefix_trie_match(Trie, Text, start) 
        if val != None:
            for i in range(val[0], val[1]):
                pos_dct[i] = Text[i].upper()
            start += 1
        else:
            if start not in pos_dct:
                pos_dct[start] = Text[start].lower()
            start += 1
    newText = ''.join(pos_dct.values())
    
    if len(Text) != len(newText):
        assert newText == ''
        print(newText)
        newText = Text.lower()
    return newText

def trie_matching(Trie, Text):
    positions = []
    start = 0
    while start < len(Text):
        val = prefix_trie_match(Trie, Text, start)
        if val != None:
            positions.append(val[0])
        start += 1
    return positions

