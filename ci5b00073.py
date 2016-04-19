##MOARF, an Integrated Workflow for Multiobjective Optimization: Implementation, Synthesis, and Biological Evaluation
from rdkit import Chem
from rdkit.Chem import AllChem
from sys import argv


BIARYL = Chem.MolFromSmarts('[a]-&!@[a]')
ALKENE = Chem.MolFromSmarts('[#6](=&!@[#6][#6])[#6]')
IODINE = Chem.MolFromSmarts('[R]-[I&D1]')
HETATS = Chem.MolFromSmarts('[!#6&!R]-[!#6&!R]')
GLYCOS = Chem.MolFromSmarts('[!#6]-[A]-&@[#7,#8,#16]')
UNSATH = Chem.MolFromSmarts('[!#6]-&!@[A,a]!-[!#6]')
ALKYLH = Chem.MolFromSmarts('[!a&!#1]-[!#6]-[a&R]')
BENZYL = Chem.MolFromSmarts('[!#6]-[A&!R]-[a&R]')
EXONIT = Chem.MolFromSmarts('[#7&R]-&!@[A,a]')
ENOLIC = Chem.MolFromSmarts('[A,a!R]-[!#6!R]-[#6!R]!-[#6!R]')

orderedCutList = [BIARYL, ALKENE, IODINE, HETATS, GLYCOS, UNSATH, ALKYLH, BENZYL, EXONIT, ENOLIC]

def molFragmenter(mol, cutList = orderedCutList, addDummies = False, numberedDummies = True):
    if(not mol):
        return None
    numRings = mol.GetRingInfo().NumRings()
    atomPairs =[]
    #iterate of each type of cut and make them if they satisify conditions
    for i in cutList:
        matches = mol.GetSubstructMatches(i)
        for j in matches:
            try:
                Chem.Kekulize(mol, clearAromaticFlags = True)
            except:
                return None
            tempMol = Chem.EditableMol(mol)
            wasBond = mol.GetBondBetweenAtoms(j[0], j[1])
            tempMol.RemoveBond(j[0], j[1])
            makeCut = checkCut(tempMol, numRings)
            if(makeCut and wasBond is not None):
                if(addDummies):
                    atomPairs.append([j[0], j[1]])
                mol = tempMol.GetMol()
            try:
                Chem.SanitizeMol(mol)
            except:
                return None
    try:
        Chem.Kekulize(mol, clearAromaticFlags = True)
    except:
        mol = AdjustAromaticNs(mol)
    try:
        Chem.Kekulize(mol, clearAromaticFlags = True)
    except:
        return None

    #This hopefulyl prevents any dodgey explicit H's lingering, works so far...
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in (7,15,16) and atom.GetNumExplicitHs()==1:
            atom.SetNumExplicitHs(0)

    tempMol = Chem.EditableMol(mol)

    dumAtom = Chem.Atom(0)
    for i in range(len(atomPairs)):
        if(numberedDummies):
            dumAtom.SetIsotope(i+1)
        firstIdx = tempMol.AddAtom(dumAtom)
        secondIdx = tempMol.AddAtom(dumAtom)
        tempMol.AddBond(atomPairs[i][0], firstIdx, Chem.BondType.SINGLE)
        tempMol.AddBond(atomPairs[i][1], secondIdx, Chem.BondType.SINGLE)

    mol = tempMol.GetMol()
    try:
        Chem.SanitizeMol(mol)
    #************************************************************************************
    # Use James Davidson and Greg Landrum's function to try and fix any aromatic problems
    # See http://www.rdkit.org/docs/Cookbook.html#cleaning-up-heterocycles
    #************************************************************************************
    except:
        mol = AdjustAromaticNs(mol)
    try:
        Chem.SanitizeMol(mol)
        mol = Chem.GetMolFrags(mol, asMols = True)
    #Add in any raising of errors in here I just return None for my usage
    except:
        return None
    return mol


def checkCut(tempMol, numRings):
    # This function checks that there has been no ring breakage (depreciated now I've swapped to SMARTS)
    # Also check that I've left no isolated heteroatoms 
    mol = tempMol.GetMol()
    mol = AllChem.DeleteSubstructs(mol, Chem.MolFromSmarts('[#0]'))
    Chem.Kekulize(mol, clearAromaticFlags = True)
    try:
        Chem.SanitizeMol(mol)
    except:
        mol = AdjustAromaticNs(mol)
    try:
        Chem.SanitizeMol(mol)
    except:
        pass
    newNumRings = mol.GetRingInfo()
    newNumRings = newNumRings.NumRings()
    if(numRings != newNumRings):
        return False
    if(mol.HasSubstructMatch(Chem.MolFromSmarts('[!#6&!#53&D0]'))):
        return False
    return True


def removeRGroups(inMol):
    # Strips off the dummies and deals with heterocycles
    patt = Chem.MolFromSmiles('[*]')
    outMol = AllChem.DeleteSubstructs(inMol, patt)
    outMol = AdjustAromaticNs(outMol)
    return outMol

#Everything below this comment is James/Greg's code


def _FragIndicesToMol(oMol,indices):
    em = Chem.EditableMol(Chem.Mol())

    newIndices={}
    for i,idx in enumerate(indices):
        em.AddAtom(oMol.GetAtomWithIdx(idx))
        newIndices[idx]=i

    for i,idx in enumerate(indices):
        at = oMol.GetAtomWithIdx(idx)
        for bond in at.GetBonds():
            if bond.GetBeginAtomIdx()==idx:
                oidx = bond.GetEndAtomIdx()
            else:
                oidx = bond.GetBeginAtomIdx()
            # make sure every bond only gets added once:
            if oidx<idx:
                continue
            em.AddBond(newIndices[idx],newIndices[oidx],bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    Chem.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap=newIndices
    return res

def _recursivelyModifyNs(mol,matches,indices=None):
    if indices is None:
        indices=[]
    res=None
    while len(matches) and res is None:
        tIndices=indices[:]
        nextIdx = matches.pop(0)
        tIndices.append(nextIdx)
        nm = Chem.Mol(mol)
        nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
        nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(1)
        cp = Chem.Mol(nm)
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            res,indices = _recursivelyModifyNs(nm,matches,indices=tIndices)
        else:
            indices=tIndices
            res=cp
    return res,indices

def AdjustAromaticNs(m,nitrogenPattern='[n&D2&H0;r5,r6]'):
    """
        default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
        to fix: O=c1ccncc1
        """
    if(not m):
        return None
    Chem.GetSymmSSSR(m)
    m.UpdatePropertyCache(False)

    # break non-ring bonds linking rings:
    em = Chem.EditableMol(m)
    linkers = m.GetSubstructMatches(Chem.MolFromSmarts('[r]!@[r]'))
    plsFix=set()
    for a,b in linkers:
        em.RemoveBond(a,b)
        plsFix.add(a)
        plsFix.add(b)
    nm = em.GetMol()
    for at in plsFix:
        at=nm.GetAtomWithIdx(at)
        if at.GetIsAromatic() and at.GetAtomicNum()==7:
            at.SetNumExplicitHs(1)
            at.SetNoImplicit(True)

    # build molecules from the fragments:
    fragLists = Chem.GetMolFrags(nm)
    frags = [_FragIndicesToMol(nm,x) for x in fragLists]

    # loop through the fragments in turn and try to aromatize them:
    ok=True
    for i,frag in enumerate(frags):
        cp = Chem.Mol(frag)
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
            lres,indices=_recursivelyModifyNs(frag,matches)
            if not lres:
                #print 'frag %d failed (%s)'%(i,str(fragLists[i]))
                ok=False
                break
            else:
                revMap={}
                for k,v in frag._idxMap.iteritems():
                    revMap[v]=k
                for idx in indices:
                    oatom = m.GetAtomWithIdx(revMap[idx])
                    oatom.SetNoImplicit(True)
                    oatom.SetNumExplicitHs(1)
    if not ok:
        return None
    return m

