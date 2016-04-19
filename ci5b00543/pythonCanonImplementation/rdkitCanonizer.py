import time

import CanonGraph as cg
import RDKit_Graph_Invariants as rdk
import SpecialSymmetryGraph_Invariant as sg

def testUseOfSpecialSymmetryInvariant(graph,symClassDict):
    numNodesInCycles=0
    numSymNodesInCycles=0
    nodesInMultipleCycles=True
    nodeSymClassDict={}
    # invert the symClassDict
    for i in symClassDict.items():
        for j in i[1]:
            nodeSymClassDict[j]=i[0]
    # search for nodes within cycles
    for i in graph.nodes:
        if rdk.isRDKitRingAtom(graph.graph, i):
            numNodesInCycles+=1
            # count nodes within cycles which have more than 2 symmetry equivalent nodes
            numNodesInSymClass = len(symClassDict[nodeSymClassDict[i]]) 
            if numNodesInSymClass > 2:
                numSymNodesInCycles+=numNodesInSymClass
            # find nodes which are in more than one cycle which have symmetry equivalent nodes
            if rdk.isRDKitAtomInMultipleRings(graph.graph, i) > 1 and numNodesInSymClass > 1:
                nodesInMultipleCycles=True
    if nodesInMultipleCycles and numNodesInCycles > 0 and float(numSymNodesInCycles)/numNodesInCycles > 0.5:
        return True
    return False

def canonizeGraph(graph, output=None):
    """ Calculate the canonical labeling of a graph/molecule """
    t=time.time()
    symClassDict=cg.createSymClassDict(len(graph.nodes))
    cg.refinePartitions(graph,symClassDict)
    if output is not None:
        output.append(time.time()-t)
    t=time.time()
    if len(symClassDict.keys()) < len(graph.nodes):
        # try more specialized properties before tie breaking
        # first use chirality to get canonical order
        graph.nodeProperties=rdk.setRDKitSpecialChiralityProperties(graph.graph)
        # this special invariant depends on the current index/order therefor we need to provide a special function for comparision
        graph.calcIndexDependentPropertyKey=rdk.calcIndexDependentRDKitChiralityPropertyKey
        # chose the right key for comparision
        graph.useIndexDependentPropKey=True
        cg.refinePartitions(graph,symClassDict)
        if output is not None:
            output.append(time.time()-t)
        t=time.time()
        graph.useIndexDependentPropKey=False
        keysbefore = len(symClassDict.keys())
        if testUseOfSpecialSymmetryInvariant(graph,symClassDict):
            # try a very specialized graph invariant which helps for symmetrical and cyclic graphs
            useInvariant, properties=sg.setSpecialSymmetryPropertiesCanonGraph(graph,isNodeInCycle=rdk.isRDKitRingAtom)
            if useInvariant:
                graph.nodeProperties=properties
                # do not forget to set back the choice of the sorting key for the first round of refinement
                graph.useIndexKey=False   
                cg.refinePartitions(graph,symClassDict)
                if output is not None:
                    output.append(time.time()-t)
        else:
            if output is not None:
                output.append(time.time()-t)
        # finally if there are still ties break them
        graph.useIndexKey=True
        cg.breakTies(graph,symClassDict)
        if output is not None:
            output.append(time.time()-t)

