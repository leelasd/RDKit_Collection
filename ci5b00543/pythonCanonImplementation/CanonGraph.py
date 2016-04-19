from __future__ import absolute_import, division, print_function
from collections import defaultdict
from types import ListType
import operator

class CanonGraph():
    def __init__(self,graph,setNodes,setProperties,setNeighbors,useEdgeProperties=False,setNeighborsStructuredByEdgeProperties=None):
        self.nodes = []
        self.nodeNeighbors = defaultdict(list)
        self.nodeProperties = defaultdict(dict)
        self.nodeIndex = defaultdict(int)
        self.useIndexKey = False
        self.useIndexDependentPropKey = False
        self.useEdgeProperties = useEdgeProperties
        self.graph = graph
        self.setNodes = setNodes
        self.setProperties = setProperties
        self.setNeighbors = setNeighbors
        self.setNeighborsStructuredByEdgeProperties = setNeighborsStructuredByEdgeProperties
        self.getKeyNode = self._getKeyNode
        self.calcIndexDependentPropertyKey = None
        self._graphToCanonGraph()
            
    def _graphToCanonGraph(self):
        self.nodes = self.setNodes(self.graph)
        if self.useEdgeProperties and (self.setNeighborsStructuredByEdgeProperties is not None):
            self.nodeNeighbors = self.setNeighborsStructuredByEdgeProperties(self.graph)      
        else:
            self.nodeNeighbors = self.setNeighbors(self.graph)
        self.nodeProperties = self.setProperties(self.graph)
        self.nodeIndex = defaultdict(int, [(n,0) for n in self.nodes])
        
    def _getPropertyKey(self,i):
        key = []
        for k in sorted(self.nodeProperties[i]):
            if type(self.nodeProperties[i][k]) is ListType:
                key.extend(self.nodeProperties[i][k])
            else:
                key.append(self.nodeProperties[i][k])
        return tuple(key)

    def _getIndexKey(self,i):       
        key=0
        if self.useEdgeProperties:
            key = [self.nodeIndex[i]]
        for k in self.nodeNeighbors[i]:
            if not self.useEdgeProperties:
                key+=self.nodeIndex[k]
            else:
                nkeys = []
                for e in k:
                    nkeys.append(self.nodeIndex[e])
                nkeys.sort(reverse=True)
                key.extend(nkeys)
        if self.useEdgeProperties:
            return tuple(key)
        else:
            return tuple([key])
        
    def _getIndexDependentPropertyKey(self,i):
        key=[]
        if self.calcIndexDependentPropertyKey is not None:
            key=self.calcIndexDependentPropertyKey(i, self.nodeProperties, self.nodeIndex, self.nodeNeighbors)
        return tuple(key)
    
    def _getKeyNode(self,i):
        if self.useIndexDependentPropKey:
            return self._getIndexDependentPropertyKey(i)
        elif self.useIndexKey:
            return self._getIndexKey(i)
        else:
            return self._getPropertyKey(i)


def _updateSymClassDict(sortedElements,symClassDict,graph,currentSymClass):
    last_i = sortedElements[0]
    sameSymClass = [last_i]
    newSymClasses = []
    count = currentSymClass
    for i in sortedElements[1:]:
        count+=1
        if graph.getKeyNode(last_i) != graph.getKeyNode(i):
            symClassDict[currentSymClass] = sameSymClass
            currentSymClass=count
            newSymClasses.append(currentSymClass)
            sameSymClass = []
            sameSymClass.append(i)
            last_i = i
        else:
            sameSymClass.append(i)
            last_i = i
        symClassDict[currentSymClass] = sameSymClass
    return newSymClasses


def _updateSymClassesCanonGraph(graph,symClassDict):
    for k in symClassDict:
        for e in symClassDict[k]:
            graph.nodeIndex[e] = k


def createSymClassDict(numElements):
    """ Create a symmetry class node map for the graph """
    symClassElementDict=defaultdict(list)
    symClassElementDict[0]=[i for i in range(numElements)]
    return symClassElementDict


def refinePartitions(graph,symClassDict):
    """ Iteratively refine the partitions of the graph to obtain a canonical ranking """
    partitionsToSort = sorted([i for i in symClassDict.items() if len(i[1])>1],reverse=True)
    while len(partitionsToSort) :
        # get the partition with highest symClass
        index,partitionMembers = partitionsToSort.pop(0)
        # canonical ranking within the current partition
        sortedElements = sorted(partitionMembers,key=graph.getKeyNode)
        #update the symmetry class
        newSymClasses = _updateSymClassDict(sortedElements,symClassDict,graph,index)
        if len(newSymClasses) :
            # update the symmetry classes of the elements
            _updateSymClassesCanonGraph(graph,symClassDict)
            graph.useIndexKey=True
            if graph.useIndexDependentPropKey:
                graph.useIndexKey=False
            addPartitions=[]
            # get the partitions which are affected by the new symmetry classes
            for cls in newSymClasses:
                for e in symClassDict[cls]:
                    for nbr in graph.nodeNeighbors[e]:
                        if graph.useEdgeProperties:
                            for n in nbr:
                                addPartitions.append(graph.nodeIndex[n])
                        else:
                            addPartitions.append(graph.nodeIndex[nbr])
            # add the partitions at the beginning which should be re-analyzed
            for p in sorted(set(addPartitions)):
                if len(symClassDict[p]) > 1 and tuple((p,symClassDict[p])) not in partitionsToSort:
                    partitionsToSort.insert(0,(p,symClassDict[p]))


def breakTies(graph,symClassDict):
    """ Break ties in canon graph and run again the refinement of the graph to obtain a canonical ranking """
    partitionsToSort = sorted([i for i in symClassDict.items() if len(i[1])>1])
    while len(partitionsToSort) :
        # get the partition with highest symClass
        index,partitionMembers = partitionsToSort.pop(0)
        # get the length of the current partition
        length = len(partitionMembers)
        # determine the symclass of the last element in the current partition
        offset = index+length-1
        # set the new symclass of the element
        graph.nodeIndex[partitionMembers[-1]] = offset
        # update the symmetry class dict
        symClassDict[offset] = [partitionMembers.pop()]
        symClassDict[index] = partitionMembers
        # refine partitions again
        refinePartitions(graph,symClassDict)
        # update the partitions list
        partitionsToSort = sorted([i for i in symClassDict.items() if len(i[1])>1])

