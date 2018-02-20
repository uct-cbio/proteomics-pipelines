#!/usr/bin/env python

# Python graph functions and algorithms
import string
import datrie

def fibonacci( n, a=0, b=1):
    for i in range(0, n):
        a, b = b, a + b
    return a

class Vertex:
    def __init__(self,key):
        self.id = key
        self.connectedTo = {}
    def addNeighbour(self,neighbour,edge_label=None):
        self.connectedTo[neighbour] = edge_label
    def __str__(self):
        return str(self.id) + ' connectedTo: '+ str([x.id for x in self.connectedTo])
    def getConnections(self):
        return self.connectedTo.keys()
    def getId(self):
        return self.id
    def getEdge(self,neighbour):
        return self.connectedTo[neighbour]

class Graph:
    def __init__(self):
        self.vertices = {}
        self.numVertices = 0
    def addVertex(self,key):
        assert key not in self.vertices
        self.numVertices = self.numVertices + 1
        newVertex = Vertex(key)
        self.vertices[key] = newVertex
        return newVertex
    def getVertex(self,n):
        if n in self.vertices:
            return self.vertices[n]
        else:
            return None
    def __contains__(self,n):
        return n in self.vertices 
    def addEdge(self,f,t,edge_label=0):
        if f not in self.vertices:
            self.addVertex(f)
        if t not in self.vertices:
            self.addVertex(t)
        self.vertices[f].addNeighbour(t, edge_label)
    def getVertices(self):
        return self.vertices.keys()
    def __iter__(self):
        return iter(self.vertices.values())


class Trie:
    def __init__(self, lst):
        self.trie = self.trie_graph(lst)

    def trie_graph(self, lst):
        trie = datrie.BaseTrie(string.ascii_uppercase)
        for l in lst:
            trie[l] = 0
        return trie

class TrieMatch:
    def __init__(self, Trie, Text):
        trie = Trie.trie
        self.Text=Text
        self.Text_coordinates = self.trie_coordinates(trie)

    def prefix_trie_match(self, string, start, trie):
        last_edges = []
        prefixes  = trie.prefixes(string[start:])
        for pref in prefixes:
            last_edge = start, start + len(pref)
            last_edges.append(last_edge)
        return last_edges

    def trie_matching(self):
        coordinates = self.Text_coordinates
        positions = [i[0] for i in coordinates]
        return positions

    def trie_upper(self):
        positions = set()
        coords =  self.Text_coordinates
        for coord in coords:
            for indx in range(coord[0], coord[1]):
                positions.add(indx)
        new_Text_list = []
        for indx in range(len(self.Text)):
            if indx not in positions:
                new_Text_list.append(self.Text[indx].lower())
            else:
                new_Text_list.append(self.Text[indx].upper())
        new_Text = ''.join(new_Text_list)
        return new_Text

    def trie_coordinates(self, trie): 
        coordinates = []
        start = 0
        while start < len(self.Text):
            vals = self.prefix_trie_match( self.Text, start, trie)
            for val in vals:
                coordinates.append(val)
            start += 1 
        return coordinates

    def trie_export(self):
        coordinates = self.Text_coordinates
        words = set()
        for coord in coordinates:
            word = self.Text[coord[0]:coord[1]]
            words.add(word)
        return words

    def trie_coverage(self): 
        positions = set()
        coords =  self.Text_coordinates
        for coord in coords:
            for indx in range(coord[0], coord[1]):
                positions.add(indx)
        Text_indxs=set([ pos for pos in range(len(self.Text))]) 
        Text_coverage = len(positions)/float(len(Text_indxs)) * 100
        return Text_coverage
        



