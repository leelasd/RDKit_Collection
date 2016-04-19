from __future__ import print_function
from rdkit import Chem
from rdkit.RDLogger import logger
logger=logger()
import gzip,random,sys,os

random.seed(42)
maxToDo=2000000
nprop = '_Name'
nDone=0
npf=0
ncf=0
nef=0
if len(sys.argv) < 2:
    print("Usage: python",sys.argv[0],"[filename.sdf]")
for fname in sys.argv[1:]:
    logger.info('Doing file: %s'%fname)
    basename,ext=os.path.splitext(os.path.basename(fname))
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
            csmi = Chem.MolToSmiles(m,True)
            m2=Chem.MolFromSmiles(csmi)
            if not m2:
                print(basename,i,csmi,nm,file=errorf)
                continue
            # try renumbering the atoms and testing that we get the same SMILES:
            nReps = 50
            for rep in range(nReps):
                aids = list(range(m2.GetNumAtoms()))
                random.shuffle(aids)
                rm2 = Chem.RenumberAtoms(m2,aids)
                ncsmi = Chem.MolToSmiles(rm2,True)
                if ncsmi!=csmi:
                    ncf+=1
                    print(basename,i,ncsmi,csmi,nm,rep,aids,'\n',file=canonf)
                    break
            nDone+=1
            if not nDone%1000:
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




