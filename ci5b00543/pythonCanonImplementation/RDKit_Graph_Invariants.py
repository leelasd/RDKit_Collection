from __future__ import absolute_import, division, print_function
from rdkit import Chem
from collections import defaultdict
from types import ListType
import operator

# a collection of RDKit specific functions

def _getNumRingStereoNeighbors(a):
    count=0
    for nbr in a.GetNeighbors():
        if (nbr.GetChiralTag()==Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW or nbr.GetChiralTag()==Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW) and nbr.HasProp("_ringStereoAtoms"):
            count+=1
    return count

def setRDKitAtomProperties(mol):
    res =defaultdict(dict)
    for a in mol.GetAtoms():
        props = {}
        props["1_degree"] = a.GetDegree()
        props["2_atomicNum"] = a.GetAtomicNum()
        props["3_isotope"] = a.GetIsotope()
        props["4_totalNumHs"] = a.GetTotalNumHs()
        props["5_charge"] = a.GetFormalCharge()
        props["6_cipCode"] = 0
        if a.HasProp("_CIPCode"):
            cip = a.GetProp("_CIPCode")
            if cip == "R":
                props["6_cipCode"] = 2
            else:
                props["6_cipCode"] = 1
        props["7_chiralTag"] = 0
        if a.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            props["7_chiralTag"] = 1
        props["8_numRingStereoNeighbors"] = _getNumRingStereoNeighbors(a)
        props["9_bonds"] = []
        bonds=[]
        for b in a.GetBonds():
            bond = []
            bondType = int(b.GetBondType())
            if b.GetIsAromatic():
                bondType = int(Chem.rdchem.BondType.AROMATIC)
            bond.append(bondType)
            bondStereo = 0
            if b.GetStereo() != Chem.rdchem.BondStereo.STEREOANY and b.GetStereo() != Chem.rdchem.BondStereo.STEREONONE:
                bondStereo = int(b.GetStereo())
            bond.append(bondStereo)
            bonds.append(bond)
        bonds.sort(reverse=True)
        for bo in bonds:
            props["9_bonds"].extend(bo)
        res[a.GetIdx()]=props
    return res

def setRDKitAtomNeighborsStructedByBondProperties(mol):
    neighborDict=defaultdict(list)
    for a in mol.GetAtoms():
        # first get the edge properties
        if len(a.GetBonds()) < 1:
            neighborDict[a.GetIdx()]=[]
            continue
        bonds=[]
        for b in a.GetBonds():
            bond = []
            bondType = int(b.GetBondType())
            if b.GetIsAromatic():
                bondType = int(Chem.rdchem.BondType.AROMATIC)
            bondStereo = 0
            if b.GetStereo() != Chem.rdchem.BondStereo.STEREOANY and b.GetStereo() != Chem.rdchem.BondStereo.STEREONONE:
                bondStereo = int(b.GetStereo())
            bond.append((bondType,bondStereo))
            bond.append(b.GetOtherAtomIdx(a.GetIdx()))
            bonds.append(bond)
        bonds = sorted(bonds, key=operator.itemgetter(0), reverse=True)
        # combine the neighbors by same edge properties
        last_b=bonds[0]
        neighbors=[]
        sameEdgeProps = []
        for b in bonds:
            if b[0] == last_b[0]:
                sameEdgeProps.append(b[1])
            else:
                neighbors.append(sameEdgeProps)
                sameEdgeProps=[b[1]]
            last_b = b
        neighbors.append(sameEdgeProps)
        neighborDict[a.GetIdx()]=neighbors
    return neighborDict

def setRDKitAtomNeighbors(mol):
    return defaultdict(list, [(a.GetIdx(),[x.GetOtherAtomIdx(a.GetIdx()) for x in a.GetBonds()]) for a in mol.GetAtoms()])

def setRDKitAtomNodes(mol):
    return [a.GetIdx() for a in mol.GetAtoms()]

def setRDKitSpecialChiralityProperties(mol):
    res =defaultdict(dict)
    for a in mol.GetAtoms():
        props = {}
        props["1_chiralTag"] = 0
        if a.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
            props["1_chiralTag"] = 1
        elif a.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            props["1_chiralTag"] = 2
        props["2_refNbrs"] = [x.GetOtherAtomIdx(a.GetIdx()) for x in a.GetBonds()]
        res[a.GetIdx()]=props
    return res

def _getNumSwapsToInterconvert(ref,probe):
    numSwaps=0
    if len(ref) != len(probe):
        return numSwaps
    for a,i in enumerate(ref):
        for b,j in enumerate(probe):
            if i != j:
                continue
            else:
                if a == b:
                    continue
                probe[a], probe[b] = probe[b], probe[a]
                numSwaps+=1
    return numSwaps 

def calcIndexDependentRDKitChiralityPropertyKey(i, propDict, indexDict, nbrDict):
    res=[]
    nbrIds = propDict[i]["2_refNbrs"]
    for nbr in nbrIds:
        if propDict[nbr]["1_chiralTag"] == 0:
            res.append((indexDict[nbr],0))
        else:
            ref = propDict[nbr]["2_refNbrs"]
            probe = [i]
            # get probe neighbors which order depending on the current index/sym class
            idSortedNbr = []
            if type(nbrDict[nbr][0]) == ListType:
                for j in nbrDict[nbr]:
                    idSortedStructNbr=[]
                    for e in j:
                        idSortedStructNbr.append((indexDict[e],e))
                    idSortedStructNbr = sorted(idSortedStructNbr, key=operator.itemgetter(0), reverse=True)
                    idSortedNbr.extend(idSortedStructNbr)
            else:
                for j in nbrDict[nbr]:
                    idSortedNbr.append((indexDict[j],j))
                idSortedNbr = sorted(idSortedNbr, key=operator.itemgetter(0), reverse=True)
            probe.extend([n[1] for n in idSortedNbr if n[1] != i])
            # calculate the number of swaps need to interconvert between original order and index-dependent order
            numSwaps = _getNumSwapsToInterconvert(ref, probe)
            if propDict[nbr]["1_chiralTag"] == 1:
                if numSwaps%2:
                    res.append((indexDict[nbr],2))
                else:
                    res.append((indexDict[nbr],1))
            elif propDict[nbr]["1_chiralTag"] == 2:
                if numSwaps%2:
                    res.append((indexDict[nbr],1))
                else:
                    res.append((indexDict[nbr],2))
    res.sort()
    return res

def isRDKitRingAtom(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    return atom.IsInRing()

def isRDKitAtomInMultipleRings(mol, idx):
    ri = mol.GetRingInfo()
    return (ri.NumAtomRings(idx)>1)

