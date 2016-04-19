from rdkit import Chem
from rdkit.Chem import AllChem
from math import sqrt,pow
import sys


def eucl_dist(mol_x_com, mol_y_com, mol_z_com, x_com, y_com, z_com):
    """Function to find the Euclidean distance between two points
    Takes 6 floats (coords)
    Returns a Euclidean distance"""
    return sqrt(pow((mol_x_com - x_com), 2) + pow((mol_y_com - y_com), 2) + pow((mol_z_com - z_com), 2))


def find_conformations(mol, core, useTethers=True, coreConfId=-1, randomseed=2342, max_iters=200, opt=None):
    """Function to generate conformations. Heavily based on ConstrainedEmbed in the RDKit
    Uses a forcefield (default MMFF) to generate conformations constrained to a
    core smiles. Does energy minimisation. Calculates the RMSD
    Takes an RDKit molecule and a core. Options are to useTethers,
    coreConfId - the conformer ID to use, randomseed - the randomseed to use,
    maxIts - the maximum number of iterations for the minimisation,
    opt -  the forcefield to use.
    Returns an RDKit molecule
    """
    # Re-read molecule prevents failures
    mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), removeHs=False)
    # Find the shared core
    match = mol.GetSubstructMatch(core)
    # If there is no match - print this error
    if not match:
        print "CORE ", Chem.MolToSmiles(core, isomericSmiles=True)
        print "MOL ", Chem.MolToSmiles(mol, isomericSmiles=True)
        raise ValueError, "molecule doesn't match the core"
    # Derive a coordinate map
    coordMap = {}
    coreConf = core.GetConformer(coreConfId)
    # Go through the atoms in the match
    for i, idxI in enumerate(match):
        # Find the atomic position
        corePtI = coreConf.GetAtomPosition(i)
        # Add this to the coordinate map
        coordMap[idxI] = corePtI
    # Embed this using random coords
    ci = AllChem.EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, useRandomCoords=True)
    # If it was unsuccesful - > raise this error
    if ci < 0:
        mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
        ci = AllChem.EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, useRandomCoords=True)
        if ci < 0:
            print Chem.MolToMolBlock(mol)
            print Chem.MolToMolBlock(core)
            raise ValueError, 'Could not embed molecule.'

    # Now make a map of the points to tether
    algMap = [(j, i) for i, j in enumerate(match)]
    if not useTethers:
        # clean up the conformation
        if opt is "MMFF":
            try:
                # Make the unsanitized molecule
                mmff_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), sanitize=False, removeHs=False)
                # Get the forcefield
                myff = Chem.rdForceFieldHelpers.SetupMMFFForceField(mmff_mol, mmffVerbosity=0)
                # Get the forcefield for this molecule
                ff = AllChem.MMFFGetMoleculeForceField(mol, myff, confId=0)
            # Because the newer version of RDKit has this difference
            except AttributeError:
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
        else:
            # Simply get the forcefield 
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
        # Now loop over the atoms in the match
        for i, idxI in enumerate(match):
            for j in range(i + 1, len(match)):
                idxJ = match[j]
                #  Add this constraint
                d = coordMap[idxI].Distance(coordMap[idxJ])
                ff.AddDistanceConstraint(idxI, idxJ, d, d, 100.)
        # Intitialise the forcefeild
        ff.Initialize()
        n = 4
        # Minimise
        more = ff.Minimize()
        while more and n:
            more = ff.Minimize()
            n -= 1
        # rotate the embedded conformation onto the core:
        rms = AllChem.AlignMol(mol, core, atomMap=algMap)
    else:
        # rotate the embedded conformation onto the core:
        rms = AllChem.AlignMol(mol, core, atomMap=algMap)
        if opt is "MMFF":
            try:
                mmff_mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), sanitize=False, removeHs=False)
                myff = Chem.rdForceFieldHelpers.SetupMMFFForceField(mmff_mol, mmffVerbosity=0)
                ff = AllChem.MMFFGetMoleculeForceField(mol, myff, confId=0)
            # Because the newer version of RDKit has this difference
            except AttributeError:
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=0)
        conf = core.GetConformer()
        if ff is None:
            sys.stderr.write("FORCEFIELD IS NONE\n" + Chem.MolToSmiles(mol))
            return None
        # Now go thtough these atoms
        for i in range(core.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
            # Add a costraint to the FF
            ff.AddDistanceConstraint(pIdx, match[i], 0, 0.0, 100.)
        ff.Initialize()
        n = 0
        #Do an energy minimisation
        # Gradient set to deliver multiple conformations
        more = ff.Minimize(maxIts=max_iters, energyTol=1e-4, forceTol=1e-3)
        if opt is "MMFF":
            n = 4
        while more and n:
            more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
            n -= 1
        # realign
        rms = AllChem.AlignMol(mol, core, atomMap=algMap)
    mol.SetProp('EmbedRMS', str(rms))
    return mol


def check_for_clashes(rdmol, protein_coords=None, threshold=0.5):
    """Function to look for clashes between the protein and the RDKit molecule
    rdmol - the RDKit molecule
    protein - a list of geometry points for the protein
    threshold - the min distance between points"""
    # If the protein restrictions do not exist- 
    if not protein_coords:
        return True
    # Find the molecule coords
    conf = rdmol.GetConformer()
    mol_pos = [conf.GetAtomPosition(atm.GetIdx()) for atm in rdmol.GetAtoms()]
    if type(protein_coords) is list:
        # We need the molecule to break ALL of these
        for pos in mol_pos:
            for prot_counter, protein in enumerate(protein_coords):
                len_c = len(protein) - 1
                print "ON NEXT PROTEIN"
                for i, prot_pos in enumerate(protein):
                    if eucl_dist(pos.x, pos.y, pos.z, prot_pos.x, prot_pos.y, prot_pos.z) < threshold:
                        # break from this protein - look for next
                        break
                    if i == len_c:
                        return True
        # If we get to the end here then we return fasle
        return False
    else:
        for pos in mol_pos:
            for prot_pos in protein:
                if eucl_dist(pos.x, pos.y, pos.z, prot_pos.x, prot_pos.y, prot_pos.z) < threshold:
                    return None
    return True


def generate_conformations(rdmol, core_smi, num_confs=100, num_fails=150, max_iters=200, ff="MMFF", prot_pos=None, min_diff=0.35):
    """Function to generate conformations for the activity point molecule.
    Removes overly similar conformations.
    Input RDKit molecule of the ActivityPoint molecule (already given 3D coords),
    core_smi - a smiles string of the shared core or a molecule of it ,
    num_confs - the number of conformations required,
    num_fails - the number of potential fails where the new conformation is too similar,
    max_iters - the maximum number of iterations,
    ff - the forcefield used,
    prot_pos - the coords of steric hindrance regions,
    min_diff - the min RMSD between conformers.
    Output is a list of RDKit molecule in different conformations"""
    import random
    Chem.SanitizeMol(core_smi)
    # Now replace the * with H (if it exists)
    core_smi = add_hydrogens(core_smi)
    if core_smi is None:
        sys.stderr.write("NONE NEW MOLECULE")
        return None
    # Now add Hydrogens
    rdmol = AllChem.AddHs(rdmol, addCoords=True)
    fail_counter = 0
    conf_list = []
    conf_counter = 0
    while True:
        # If the maximum number of fails has been got - move on
        if fail_counter == num_fails:
            break
        # If the maximum number of confs has beend found - move on
        if conf_counter == num_confs:
            break
        # find a conformation in this loop
        retval = find_conformations(rdmol, core_smi, useTethers=True, randomseed=random.randint(1, 10000000), max_iters=max_iters, opt=ff)
        # If there is an error in the contrained embed
        if retval is None:
            return None
        # Check for overly-similar conformations
        if conf_counter > 0:
          #  Check for overly similar conformations
            if min([AllChem.AlignMol(Chem.MolFromMolBlock(item), retval) for item in conf_list]) < min_diff:
                # Add to the fail conunter
                fail_counter += 1
            else:
                # Add to the conf_counter and add the conf
                conf_counter += 1
                # Check for steric clash with the protein (currently None)
                if check_for_clashes(retval, prot_pos):
                  # Add this molecule
                    conf_list.append(Chem.MolToMolBlock(retval))
        else:
            # Or if this is a new one check for clashes only
            if check_for_clashes(retval, prot_pos):
                conf_counter += 1
                # Add a new molecule
                conf_list.append(Chem.MolToMolBlock(retval))
            else:
                # Otherwise add to the fail counter
                fail_counter += 1
    # Return the confs generated
    return conf_list


##### SANIFIX 4 #####################
def fragment_inds_to_mol(oMol, indices):
    em = Chem.EditableMol(Chem.Mol())
    newIndices = {}
    for i, idx in enumerate(indices):
        em.AddAtom(oMol.GetAtomWithIdx(idx))
        newIndices[idx] = i

    for i, idx in enumerate(indices):
        at = oMol.GetAtomWithIdx(idx)
        for bond in at.GetBonds():
            if bond.GetBeginAtomIdx() == idx:
                oidx = bond.GetEndAtomIdx()
            else:
                oidx = bond.GetBeginAtomIdx()
            # make sure every bond only gets added once:
            if oidx < idx:
                continue
            em.AddBond(newIndices[idx], newIndices[oidx], bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    Chem.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap = newIndices
    return res


def _recursivelyModifyNs(mol, matches, indices=None):
    if indices is None:
        indices = []
    res = None
    while len(matches) and res is None:
        tIndices = indices[:]
        nextIdx = matches.pop(0)
        tIndices.append(nextIdx)
        nm = Chem.Mol(mol.ToBinary())
        nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
        nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(1)
        cp = Chem.Mol(nm.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            res, indices = _recursivelyModifyNs(nm, matches, indices=tIndices)
        else:
            indices = tIndices
            res = cp
    return res, indices


def adjust_arom_Ns(m, nitrogenPattern='[n&D2&H0;r5,r6]'):
    """
       default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
       to fix: O=c1ccncc1
    """
    Chem.GetSymmSSSR(m)
    m.UpdatePropertyCache(False)

    # break non-ring bonds linking rings:
    em = Chem.EditableMol(m)
    linkers = m.GetSubstructMatches(Chem.MolFromSmarts('[r]!@[r]'))
    plsFix = set()
    for a, b in linkers:
        em.RemoveBond(a, b)
        plsFix.add(a)
        plsFix.add(b)
    nm = em.GetMol()
    for at in plsFix:
        at = nm.GetAtomWithIdx(at)
        if at.GetIsAromatic() and at.GetAtomicNum() == 7:
            at.SetNumExplicitHs(1)
            at.SetNoImplicit(True)

    # build molecules from the fragments:
    fragLists = Chem.GetMolFrags(nm)
    frags = [fragment_inds_to_mol(nm, x) for x in fragLists]

    # loop through the fragments in turn and try to aromatize them:
    ok = True
    for i, frag in enumerate(frags):
        cp = Chem.Mol(frag.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
            lres, indices = _recursivelyModifyNs(frag, matches)
            if not lres:
                #print 'frag %d failed (%s)'%(i,str(fragLists[i]))
                ok = False
                break
            else:
                revMap = {}
                for k, v in frag._idxMap.iteritems():
                    revMap[v] = k
                for idx in indices:
                    oatom = m.GetAtomWithIdx(revMap[idx])
                    oatom.SetNoImplicit(True)
                    oatom.SetNumExplicitHs(1)
    if not ok:
        return None
    return m

######################### OOMMPPAA FUNCTIONS ###################################

# The following  were developed to deal with a bug in the RDKit that would not 
# allow constrained embed using aromatic sulphurs


def clear_aromatic_S(my_mol):
    """Function to replace an aromatic sulphur with an aromatic oxygen with an unusual isotope
    Takes an RDKit molecule
    Returns the updated molecule"""
    for atom in my_mol.GetAtoms():
        if atom.GetIsAromatic()is True and atom.GetAtomicNum() == 16:
            atom.SetAtomicNum(8)
            atom.SetIsotope(17)
    return my_mol


def clear_isotopic_O(my_mol):
    """Function to replace an aromatic oxygen with an unusual isotope with an aromatic sulphur
    Takes an RDKit molecule
    Returns the updated molecule"""
    for atom in my_mol.GetAtoms():
        if atom.GetIsotope() == 17 and atom.GetAtomicNum() == 8:
            atom.SetAtomicNum(16)
            atom.SetIsotope(0)
    return my_mol


def remove_breaks(my_mol):
    """Function to remove the break point
    Takes an RDKit molecule
    Returns the upadte molecule"""
    for atom in my_mol.GetAtoms():
        if "*" in atom.GetSmarts():
            atom.SetAtomicNum(1)
            atom.SetIsotope(0)
    return my_mol


def canonicalise_context(my_mol):
    """Function to canonicalise the context
    Takes an RDKit molecule
    Returns the upadte molecule"""
    return Chem.MolFromSmiles(Chem.MolToSmiles(my_mol,isomericSmiles=True))


def clean_up_frags(frag_out1, frag_out2=None):
    """Function to clean up fragments that are broken becasue of nH aromaticity issues
    Takes two fragments
    Runs sanfix4 programs
    Returns the fixed fragments"""
    try:
        Chem.SanitizeMol(frag_out1)
    except ValueError:
        nm = adjust_arom_Ns(frag_out1)
        if nm is not None:
            Chem.SanitizeMol(nm)
            sys.stderr.write('Fixed aromaticity:' + Chem.MolToSmiles(nm))
            frag_out1 = nm
        else:
            sys.stderr.write('Aromaticity still broken')
            return None,None
    if frag_out2:
        try:
            Chem.SanitizeMol(frag_out2)
        except ValueError:
            nm = adjust_arom_Ns(frag_out2)
            if nm is not None:
                Chem.SanitizeMol(nm)
                sys.stderr.write('Fixed aromaticity:' + Chem.MolToSmiles(nm))
                frag_out2 = nm
            else:
                sys.stderr.write('Aromaticity still broken')
                return None,None
    return frag_out1, frag_out2


def find_attachment_point(match1, match2, mol1, mol2):
    """Function to find the attachment point for two molecules
    Takes two lists of matches and the two molecules
    Returns the updated matche lists."""
# Now align the fragments as best you can
    listmatch1 = list(match1)
    listmatch2 = list(match2)
    # Find the linking point -> an atom in match but bonded to an atom not in match
    for atm in  mol1.GetAtoms():
        if atm.GetIdx() in match1:
            if len([at for at in atm.GetNeighbors() if at not in match1]) != 0:
                listmatch1.append(atm.GetIdx())
                #Now find the substituent closest to this guy
                if len([at for at in atm.GetNeighbors() if at not in match1]) == 1:
                    listmatch1.append(atm.GetNeighbors()[0].GetIdx())
                else:
                    listmatch2.append(atm.GetNeighbors()[0].GetIdx())
    for atm in  mol2.GetAtoms():
        if atm.GetIdx() in match2:
            if len([at for at in atm.GetNeighbors() if at not in match2]) != 0:
                listmatch2.append(atm.GetIdx())
                if len([at for at in atm.GetNeighbors()if at not in match2]) == 1:
                    listmatch2.append([at for at in atm.GetNeighbors()if at not in match2][0].GetIdx())
                else:
                # Pick one at random
                    listmatch2.append([at for at in atm.GetNeighbors()if at not in match2][0].GetIdx())
    return tuple(listmatch1), tuple(listmatch2)
