#
# $Id$
#
# module to calculate a fingerprint from SMILES

from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Avalon import pyAvalonTools as fpAvalon
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdMolDescriptors

nbits = 1024
longbits = 16384

# helper function for hybrid fingerprint
def getHybridFP(fp1, fp2):
    newfp = DataStructs.ExplicitBitVect(fp1.GetNumBits()+fp2.GetNumBits())
    newfp.SetBitsFromList(fp1.GetOnBits())
    bits2 = [i+nbits for i in fp2.GetOnBits()]
    newfp.SetBitsFromList(bits2)
    return newfp

# dictionary
fpdict = {}
fpdict['ecfp0'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 0, nBits=nbits)
fpdict['ecfp2'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 1, nBits=nbits)
fpdict['ecfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=nbits)
fpdict['ecfp6'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, nBits=nbits)
fpdict['ecfc0'] = lambda m: AllChem.GetMorganFingerprint(m, 0)
fpdict['ecfc2'] = lambda m: AllChem.GetMorganFingerprint(m, 1)
fpdict['ecfc4'] = lambda m: AllChem.GetMorganFingerprint(m, 2)
fpdict['ecfc6'] = lambda m: AllChem.GetMorganFingerprint(m, 3)
fpdict['fcfp2'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 1, useFeatures=True, nBits=nbits)
fpdict['fcfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, useFeatures=True, nBits=nbits)
fpdict['fcfp6'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, useFeatures=True, nBits=nbits)
fpdict['fcfc2'] = lambda m: AllChem.GetMorganFingerprint(m, 1, useFeatures=True)
fpdict['fcfc4'] = lambda m: AllChem.GetMorganFingerprint(m, 2, useFeatures=True)
fpdict['fcfc6'] = lambda m: AllChem.GetMorganFingerprint(m, 3, useFeatures=True)
fpdict['lecfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=longbits)
fpdict['lecfp6'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, nBits=longbits)
fpdict['lfcfp4'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, useFeatures=True, nBits=longbits)
fpdict['lfcfp6'] = lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 3, useFeatures=True, nBits=longbits)
fpdict['maccs'] = lambda m: MACCSkeys.GenMACCSKeys(m)
fpdict['ap'] = lambda m: Pairs.GetAtomPairFingerprint(m)
fpdict['tt'] = lambda m: Torsions.GetTopologicalTorsionFingerprintAsIntVect(m)
fpdict['hashap'] = lambda m: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(m, nBits=nbits)
fpdict['hashap4'] = lambda m: rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(m, nBits=nbits, maxLength=4)
fpdict['hashtt'] = lambda m: rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(m, nBits=nbits)
fpdict['avalon'] = lambda m: fpAvalon.GetAvalonFP(m, nbits)
fpdict['laval'] = lambda m: fpAvalon.GetAvalonFP(m, longbits)
fpdict['rdk5'] = lambda m: Chem.RDKFingerprint(m, maxPath=5, fpSize=nbits, nBitsPerHash=2)
fpdict['rdk6'] = lambda m: Chem.RDKFingerprint(m, maxPath=6, fpSize=nbits, nBitsPerHash=2)
fpdict['rdk7'] = lambda m: Chem.RDKFingerprint(m, maxPath=7, fpSize=nbits, nBitsPerHash=2)
fpdict['m2ap'] = lambda m: getHybridFP(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=nbits), rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(m, nBits=nbits))


def CalculateFP(fp_name, smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError('SMILES cannot be converted to a RDKit molecules:', smiles)

    return fpdict[fp_name](m)
