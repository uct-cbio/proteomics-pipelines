#!/usr/bin/env python

# Python graph functions and algorithms

def fibonacci( n, a=0, b=1):
    for i in xrange(0, n):
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
    def __init__(self, lst, Text=None):
        self.trie = self.trie_graph(lst)
        self.Text_coordinates=[]
        self.Text = None
        if Text != None:
            self.Text_coordinates = self.trie_coordinates(Text)
    def trie_graph(self, lst):
        assert not isinstance(lst, basestring)    # make sure that the items you are making a trie_gr    aph from are in list format
        trie = Graph()
        root=0
        mark = 0
        trie.addVertex(root)
        for word in lst:
            node = root
            word = word.lstrip().rstrip() + '#' # separates substrings of longer strings also in list
            for i in range(0,len(word)):
                symbol = word[i]
                found = False
                for neighbour in trie.getVertex(node).getConnections():
                    if symbol == trie.getVertex(node).getEdge(neighbour):
                        node = neighbour
                        found = True
                        break
                if found == False:
                    mark = mark + 1
                    trie.addEdge(node, mark, symbol)
                    node = mark
        return trie
    def prefix_trie_match(self, string, start):
        v = 0
        i = start
        pep = []
        last_edge = None
        last_edges = []
        while i < len(string) + 1:
            found = False
            if len(self.trie.getVertex(v).getConnections()) == 0:
                return start, i
            if i == len(string) + 1:
                return
            if i < len(string):
                symbol = string[i]
                new_child = None
                for child in self.trie.getVertex(v).getConnections():
                    edge = Vertex.getEdge(self.trie.getVertex(v), child)
                    if edge == '#':
                        last_edge = start, i
                        last_edges.append(last_edge)
                    elif symbol == edge:
                        new_child = child
                        pep.append(edge)        
                        found = True
                        i += 1
                if found == True:
                    v = new_child
            elif i == len(string):
                new_child = None
                for child in self.trie.getVertex(v).getConnections():
                    edge = Vertex.getEdge(self.trie.getVertex(v), child)
                    if edge == '#':
                        last_edge = start, i
                        last_edges.append(last_edge)
            if found == False:
                return last_edges
    def trie_matching(self, Text):
        coordinates = self.trie_coordinates(Text)
        positions = [i[0] for i in coordinates]
        return positions
    def trie_upper(self, Text):
        positions = set()
        coords =  self.trie_coordinates(Text)
        for coord in coords:
            for indx in xrange(coord[0], coord[1]):
                positions.add(indx)
        new_Text_list = []
        for indx in xrange(len(Text)):
            if indx not in positions:
                new_Text_list.append(Text[indx].lower())
            else:
                new_Text_list.append(Text[indx].upper())
        new_Text = ''.join(new_Text_list)
        return new_Text
    def trie_coordinates(self, Text): 
        if Text == self.Text:
            coordinates = self.Text_coordinates
        else:
            coordinates = []
            start = 0
            while start < len(Text):
                vals = self.prefix_trie_match( Text, start)
                for val in vals:
                    coordinates.append(val)
                start += 1 
            self.Text=Text
            self.Text_coordinates=coordinates
        return coordinates
    def trie_export(self, Text):
        coordinates = self.trie_coordinates(Text)
        words = set()
        for coord in coordinates:
            word = Text[coord[0]:coord[1]]
            words.add(word)
        return words
    def trie_coverage(self, Text): 
        positions = set()
        coords =  self.trie_coordinates(Text)
        for coord in coords:
            for indx in xrange(coord[0], coord[1]):
                positions.add(indx)
        Text_indxs=set([ pos for pos in xrange(len(Text))]) 
        Text_coverage = len(positions)/float(len(Text_indxs)) * 100
        return Text_coverage
        



