from __future__ import absolute_import, division, print_function
from collections import defaultdict

# this invariant is independent of the toolkit used, it just needs the connectivity of the current canon graph
def setSpecialSymmetryPropertiesCanonGraph(canonGraph, isNodeInCycle=None):
    res=defaultdict(dict)
    for i in canonGraph.nodes:
        res_i={}
        res_i["1_neighborNum"]=""
        res_i["2_revistedNeighbors"]=""
        if isNodeInCycle is not None:
            if not isNodeInCycle(canonGraph.graph, i):
                continue
        neighbors=[i]
        visited=set([i])
        numRevisitedNbrs=dict([(j,0) for j in canonGraph.nodes])
        while len(neighbors):
            numLevelNbrs=0
            lastLevelNeighbors=set([])
            currentLevelNeighbors=set([])
            while len(neighbors):
                nbr = neighbors.pop(0)
                if isNodeInCycle is not None:
                    if not isNodeInCycle(canonGraph.graph, nbr):
                        continue
                lastLevelNeighbors.add(nbr)
                visited.add(nbr)
                nodeNeighbors = canonGraph.nodeNeighbors[nbr]
                if canonGraph.useEdgeProperties:
                    nodeNeighbors=[]
                    for n in canonGraph.nodeNeighbors[nbr]: nodeNeighbors.extend(n)
                for nn in nodeNeighbors:
                    if nn not in visited:
                        currentLevelNeighbors.add(nn)
                        numLevelNbrs+=1
                        visited.add(nn)
            for n in currentLevelNeighbors:
                nodeNeighbors = canonGraph.nodeNeighbors[n]
                if canonGraph.useEdgeProperties:
                    nodeNeighbors=[]
                    for j in canonGraph.nodeNeighbors[n]: nodeNeighbors.extend(j)
                for nn in nodeNeighbors:
                    if nn in currentLevelNeighbors or nn in lastLevelNeighbors:
                        numRevisitedNbrs[nn]+=1
            lastLevelNeighbors = currentLevelNeighbors
            neighbors.extend(currentLevelNeighbors)
            currentLevelNeighbors.clear()
            nr = sorted([k for k in numRevisitedNbrs.values() if k > 0])
            revisitedNbrs = ""
            for n in nr:
              revisitedNbrs+=str(n)
            res_i["1_neighborNum"] += str(numLevelNbrs)
            res_i["1_neighborNum"] += '|'
            res_i["2_revistedNeighbors"] += revisitedNbrs
            res_i["2_revistedNeighbors"] += '|'
            numRevisitedNbrs = {k:0 for k in numRevisitedNbrs.keys()}
        res_i["1_neighborNum"] = res_i["1_neighborNum"]         
        res_i["2_revistedNeighbors"] = res_i["2_revistedNeighbors"]
        res[i]=res_i
    return (True,res)

