from __future__ import absolute_import, division, print_function

from rdkit import Chem
from collections import defaultdict
import operator
import sys

from pythonCanonImplementation import rdkitCanonizer
from pythonCanonImplementation import CanonGraph as cg
from pythonCanonImplementation import RDKit_Graph_Invariants as rdk

if len(sys.argv) < 2:
    print("\nusage: python",sys.argv[0],"[molecule as SMILES or mol.sdf]\n\nFor example use: python moleculeCanonizer.py 'c1ccccc1C(=O)O'\n")
    sys.exit(0)
mol=None
if len(sys.argv[1].split('.')) > 1 and sys.argv[1].split('.')[1].strip() == 'sdf':
    mol = Chem.MolFromMolFile(sys.argv[1])
else:
    try:
        mol = Chem.MolFromSmiles(sys.argv[1].strip())
    except:
        print("Molecule could not be read. Please check your input.")
        sys.exit(0)
if mol is None:
    print("Molecule could not be read. Please check your input.")
    sys.exit(0)
    
res_smis = []
frags = Chem.GetMolFrags(mol, asMols=True)

for frag in frags:
    cG2 = cg.CanonGraph(frag,rdk.setRDKitAtomNodes,rdk.setRDKitAtomProperties,rdk.setRDKitAtomNeighbors,useEdgeProperties=True,setNeighborsStructuredByEdgeProperties=rdk.setRDKitAtomNeighborsStructedByBondProperties)
    
    rdkitCanonizer.canonizeGraph(cG2)
        
    frag.SetProp("_canonicalRankingNumbers","True")
    for a in frag.GetAtoms():
        a.SetProp("_canonicalRankingNumber",str(cG2.nodeIndex[a.GetIdx()]))
    res_smis.append(Chem.MolToSmiles(frag,True))
res_smis.sort()
smi = '.'.join(res_smis)
print("\nCanonical SMILES: ",smi,"\n")
sys.exit(0)


