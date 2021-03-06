from __future__ import absolute_import, division, print_function
import sys,os
import gzip
import time

from rdkit import Chem
from rdkit.Chem import SDMolSupplier
from rdkit.RDLogger import logger
logger=logger()

from pythonCanonImplementation import rdkitCanonizer
from pythonCanonImplementation import CanonGraph as cg
from pythonCanonImplementation import RDKit_Graph_Invariants as rdk

count=0
tt=0

if len(sys.argv) < 2:
    print("Usage: python",sys.argv[0],"[filename.sdf/filename.sdf.gz]")
for fname in sys.argv[1:]:
    logger.info('Doing file: %s'%fname)
    basename,ext=os.path.splitext(os.path.basename(fname))
    if ext=='.gz':
        basename=basename.split('.sdf')[0]
        inf=gzip.open(fname,'r')
    else:
        inf=open(fname,'r')
    suppl = Chem.ForwardSDMolSupplier(inf)

    outfile=file('%s.canonTimings_python.smi'%(basename),'w+')

    for i,mol in enumerate(suppl):
        try: 
            name = mol.GetProp('_Name')
            smi=''
            times=[]
            elapsed=0
            if len(Chem.GetMolFrags(mol)) > 1:
                res_smis = []
                frags = Chem.GetMolFrags(mol, asMols=True)
                for frag in frags:
                    t=time.time()
                    cG2 = cg.CanonGraph(frag,rdk.setRDKitAtomNodes,rdk.setRDKitAtomProperties,rdk.setRDKitAtomNeighbors,\
                                useEdgeProperties=True,setNeighborsStructuredByEdgeProperties=rdk.setRDKitAtomNeighborsStructedByBondProperties)

                    rdkitCanonizer.canonizeGraph(cG2,times)
                    elapsed += (time.time() - t)

                    frag.SetProp("_canonicalRankingNumbers","True")
                    for a in frag.GetAtoms():
                        a.SetProp("_canonicalRankingNumber",str(cG2.nodeIndex[a.GetIdx()]))
                    res_smis.append(Chem.MolToSmiles(frag,True))
                res_smis.sort()
                smi = '.'.join(res_smis)
            else:
                t=time.time()
                cG2 = cg.CanonGraph(mol,rdk.setRDKitAtomNodes,rdk.setRDKitAtomProperties,rdk.setRDKitAtomNeighbors,\
                useEdgeProperties=True,setNeighborsStructuredByEdgeProperties=rdk.setRDKitAtomNeighborsStructedByBondProperties)

                rdkitCanonizer.canonizeGraph(cG2,times)
                elapsed = time.time() - t

                mol.SetProp("_canonicalRankingNumbers","True")
                for a in mol.GetAtoms():
                    a.SetProp("_canonicalRankingNumber",str(cG2.nodeIndex[a.GetIdx()]))
                smi = Chem.MolToSmiles(mol,True)
            tt+=elapsed
            print(elapsed,mol.GetNumAtoms(),name,smi,times,'\n',file=outfile)
            count+=1
            if not count%1000:
                print(count, tt)
                logger.info('done %d, time already elapsed %f'%(count,tt))
                outfile.flush()
        except KeyboardInterrupt:
            break
        except:
            import traceback
            traceback.print_exc()

    outfile.flush()
    outfile.close()


