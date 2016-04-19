from __future__ import print_function
from rdkit import Chem
from rdkit.RDLogger import logger
logger=logger()
import gzip,random,sys,os

from pythonCanonImplementation import rdkitCanonizer
from pythonCanonImplementation import CanonGraph as cg
from pythonCanonImplementation import RDKit_Graph_Invariants as rdk

random.seed(42)
maxToDo=2000000
nprop = '_Name'
nDone=0
npf=0
ncf=0
nef=0


def generateCanonSmiles(mol):
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
    return smi
                  

if len(sys.argv) < 2:
    print("Usage: python",sys.argv[0],"[filename.sdf,filename.smi]")
for fname in sys.argv[1:]:
    logger.info('Doing file: %s'%fname)
    basename,ext=os.path.splitext(os.path.basename(fname))
    if ext=='.smi':
        suppl = Chem.SmilesMolSupplier(fname)
    else:
        if ext=='.gz':
            basename=basename.split('.sdf')[0]
            inf=gzip.open(fname,'r')
        else:
            inf=open(fname,'r')
        suppl = Chem.ForwardSDMolSupplier(inf)

    canonf=file('%s.canonfail.smi'%(basename),'w+')
    extraf=file('%s.extrafail.smi'%(basename),'w+')
    errorf=file('%s.errors.smi'%(basename),'w+')

    for i,m in enumerate(suppl):
        if nDone>=maxToDo:
            break
        try:
            if not m:
                npf+=1
                continue
            nm=m.GetProp(nprop)
            csmi = generateCanonSmiles(m)
            m2=Chem.MolFromSmiles(csmi)
            if not m2:
                print(basename,i,csmi,nm,file=errorf)
                continue
            # try renumbering the atoms and testing that we get the same SMILES:
            nReps = 500
            for rep in range(nReps):
                aids = list(range(m2.GetNumAtoms()))
                random.shuffle(aids)
                rm2 = Chem.RenumberAtoms(m2,aids)
                ncsmi = generateCanonSmiles(rm2)
                if ncsmi!=csmi:
                    ncf+=1
                    print(basename,i,ncsmi,csmi,nm,rep,aids,'\n',file=canonf)
                    break
            nDone+=1
            if not nDone%10:
                logger.info('done %d, %d parse, %d canon, %d extra'%(nDone,npf,ncf,nef))
                canonf.flush()
                extraf.flush()
                errorf.flush()
        except KeyboardInterrupt:
            break
        except:
            import traceback
            traceback.print_exc()

    canonf.flush()
    extraf.flush()
    errorf.flush()

    canonf.close()
    extraf.close()
    errorf.close()




