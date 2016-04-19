## 00_get_chembl_data.py

# adapted from https://www.ebi.ac.uk/chembldb/index.php/ws#pythonClient
# 
# with great help from Daniel Kuhn
# 

import urllib2
import json
import re

from sys import argv


codes = {'herg_human'  : 'Q12809' }



######

def looks_like_number(x):
    try:
        float(x)
        return True
    except ValueError:
        return False



def QueryChembl(accession=None):
    '''
    Query chembl
    '''
    ########################################################################

    # 1. Use UniProt accession to get target details

    print """
    # =========================================================
    # 1. Use UniProt accession to get target details
    # =========================================================
    """
    if not accession:
        accession = argv[1]
        #
    # test if we have a CHEMBL_TARGET_ID
    if accession.find("CHEMBL") == -1:     
      target_data = urllib2.urlopen("https://www.ebi.ac.uk/chemblws/targets/uniprot/%s.json" % accession).read()
      target_data =  json.loads(target_data)
      
      #print "Target Description: %s with %s" % (target_data['target']['description'],type(target_data['target']['description']))
      #print "Target CHEMBLID:    %s with %s" % (target_data['target']['chemblId'],type(target_data['target']['chemblId']))
      #
    #
    else:
      target_data = {}
      target_data['target'] = {}
      target_data['target']['chemblId'] = accession
    # 2. Get all bioactivties for target CHEMBL_ID

    print """

    # =========================================================
    # 2. Get all bioactivties for target CHEMBL_ID
    # =========================================================
    """

    bioactivity_data = json.loads(urllib2.urlopen("https://www.ebi.ac.uk/chemblws/targets/%s/bioactivities.json" % target_data['target']['chemblId']).read())

    print "Bioactivity Count:           %d" % len(bioactivity_data['bioactivities'])
    print "Bioactivity Count (IC50):    %d" % len([record for record in bioactivity_data['bioactivities'] if record['bioactivity_type'] == 'IC50'] )
    print "Bioactivity Count (Ki):      %d" % len([record for record in bioactivity_data['bioactivities'] if record['bioactivity_type'] == 'Ki' ])

   
    # 3. Get compounds with high binding affinity (IC50 < 100)
    
    print """

    # =========================================================
    # 3. Get compounds with binding affinity (IC50)
    # =========================================================
    """
    ic50_skip=0
    ki_skip=0
    inhb_skip=0

    count=0
    non_homo=0
    dr={}
    #for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('IC50', record['bioactivity_type']) and looks_like_number(record['value']) and float(record['value']) < 100]:
    #for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('IC50', record['bioactivity_type']) and looks_like_number(record['value']) ] :
    for bioactivity in [record for record in bioactivity_data['bioactivities'] if looks_like_number(record['value']) ] :
         
      #for k,v in bioactivity.items():    print "%30s - %s" % (k,v)
      #print
      if bioactivity['organism'] != 'Homo sapiens':
        non_homo+=1
        continue
      if re.search('IC50', bioactivity['bioactivity_type']):
        if bioactivity['units'] != 'nM':
          ic50_skip+=1
          continue
      elif re.search('Ki', bioactivity['bioactivity_type']):
        #for k,v in bioactivity.items():    print "%30s - %s" % (k,v)
        ki_skip+=1
        continue
      elif re.search('Inhibition', bioactivity['bioactivity_type']):
        inhb_skip+=1
      else:
        #print "Got activity type ", bioactivity['bioactivity_type']
        #for k,v in bioactivity.items():    print "%30s - %s" % (k,v)
        continue
        #
      #
      #print bioactivity['ingredient_cmpd_chemblid']
      cmpd_data = json.loads(urllib2.urlopen("https://www.ebi.ac.uk/chemblws/compounds/%s.json" % bioactivity['ingredient_cmpd_chemblid']).read())
      #print cmpd_data.items()
      my_smiles = cmpd_data['compound']['smiles']
      bioactivity['Smiles']=my_smiles
      dr[count] = bioactivity
      count+=1
      #
    #
    print "Skipped %i IC50 values" % ic50_skip
    print "Skipped %i Ki values" % ki_skip
    print "Skipped %i Inhibition values" % inhb_skip
    print "Skipped %i non-homo values" % non_homo
    #print dr
    import csv
    header = dr[0].keys()
    #print header
    
    ofile = "out_%s.csv" % (accession)
    dw = csv.DictWriter(open(ofile,'w'), delimiter='\t', fieldnames=header)
    dw.writerow(dict((fn,fn) for fn in header))
    for cpd in dr.keys():
        dw.writerow(dr[cpd])


    return
    # 4. Get assay details foe Ki actvity types


    print """

    # =========================================================
    # 4. Get assay details foe Ki actvity types
    # =========================================================
    """

    for bioactivity in [record for record in bioactivity_data['bioactivities'] if re.search('Ki', record['bioactivity_type'], re.IGNORECASE)]:

      #print "Assay CHEMBLID: %s" % bioactivity['assay_chemblid']

      assay_data = json.loads(urllib2.urlopen("https://www.ebi.ac.uk/chemblws/assays/%s.json" % bioactivity['assay_chemblid']).read())

      #print "  %s" % assay_data['assay']['assayDescription']


if __name__ == '__main__':
    
    for name,accession in codes.items()[:]:
        print "Checking ", name
        QueryChembl(accession)


## keeplargestfrag.py

import gzip,sys,os
from rdkit import Chem

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]


file_name    = sys.argv[1]+".onlylargestfrag.sdf"
test_cpd_out = Chem.SDWriter(file_name)

for cpd in cpds:
    fragmented = Chem.GetMolFrags(cpd,asMols=True)
    list_cpds_fragsize = []
    for x in fragmented:
        list_cpds_fragsize.append(x.GetNumAtoms())
    largest_frag_index = list_cpds_fragsize.index(max(list_cpds_fragsize))
    largest_frag = fragmented[largest_frag_index]
    test_cpd_out.write(largest_frag)
test_cpd_out.close()

f_in = open(file_name, 'rb')
file_name_gz = file_name+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
os.remove(file_name)
#f_out.flush()


## remove_dupl.py

from rdkit import Chem
import sys,gzip,os

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

uniqF = inF+".uniq.sdf"
duplF = inF+".dupl.sdf"
uniqF_SD = Chem.SDWriter(uniqF)
duplF_SD = Chem.SDWriter(duplF)

all_struct_dict = {}
all_struct_list = []
for cpd in cpds:
    Chem.RemoveHs(cpd)
    cansmi = Chem.MolToSmiles(cpd,canonical=True)
    try:
        all_struct_dict[cansmi].append(cpd)
    except:
        all_struct_dict[cansmi] = [cpd]
    all_struct_list.append(cansmi)
    
for ent in all_struct_dict:
    if len(all_struct_dict[ent]) == 1:
        can_smi = ent
        all_struct_dict[ent][0].SetProp('cansmirdkit',can_smi)
        uniqF_SD.write(all_struct_dict[ent][0])
    else:
        for singleentries in all_struct_dict[ent]:
            can_smi = ent
            singleentries.SetProp('cansmirdkit',can_smi)
            duplF_SD.write(singleentries)
uniqF_SD.close()
duplF_SD.close()

f_in = open(uniqF, 'rb')
file_name_gz = uniqF+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.remove(uniqF)

f_in = open(duplF, 'rb')
file_name_gz = duplF+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.remove(duplF)

## merge_IC50_calcstddev_removeoutliers.py

from rdkit import Chem
from rdkit.Chem import AllChem
import sys,numpy,gzip,os

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

file_out = inF+".avg.sdf"
cpds_out = Chem.SDWriter(file_out)

IC50_dict = {}
for cpd in cpds:
    cansmi = str(cpd.GetProp("cansmirdkit"))
    IC50_dict[cansmi]={}
    
def get_mean_IC50(mol_list):
    IC50 = 0
    IC50_avg = 0
    for bla in mol_list:
        try:
            IC50 +=  float(bla.GetProp("value"))
        except:
            print "no IC50 reported",bla.GetProp("_Name")
    IC50_avg = IC50 / len(mol_list)
    return IC50_avg

def get_stddev_IC50(mol_list):
    IC50_list = []
    for bla in mol_list:
        try:
            IC50_list.append(round(float(bla.GetProp("value")),2))
        except:
            print "no IC50 reported",bla.GetProp("_Name")
    IC50_stddev = numpy.std(IC50_list,ddof=1)
    # http://stackoverflow.com/questions/4575645/different-standard-deviation-for-same-input-from-wolfram-and-numpy
    return IC50_stddev,IC50_list
    
for cpd in cpds:
    cansmi = str(cpd.GetProp("cansmirdkit"))
    try:
        IC50_dict[cansmi].append(cpd)
    except:
        IC50_dict[cansmi] = [cpd]
        
for ent in IC50_dict:
    #print ent,IC50_dict[ent]
    #for x in IC50_dict[ent]:
    #    print x.GetProp("_Name")
    IC50_avg = str(get_mean_IC50(IC50_dict[ent]))
    #IC50_dict[ent][0].SetProp("value_avg",IC50_avg)
    IC50_stddev = get_stddev_IC50(IC50_dict[ent])[0]
    IC50_dict[ent][0].SetProp("value_stddev",str(IC50_stddev))
    IC50_list   = get_stddev_IC50(IC50_dict[ent])[1]
    IC50_dict[ent][0].SetProp("value_avg",str(IC50_list))
    minimumvalue = float(IC50_avg)-3*float(IC50_stddev)
    maximumvalue = float(IC50_avg)+3*float(IC50_stddev)
    IC50_dict[ent][0].SetProp("value",IC50_avg)
    if round(IC50_stddev,1) == 0.0:
        cpds_out.write(IC50_dict[ent][0])
        #print "perfect!", IC50_list, IC50_avg,IC50_stddev
    elif IC50_stddev > float(IC50_avg):
        runawaylist =[]
        for x in IC50_dict[ent]:
            runawaylist.append(x.GetProp("_Name"))
        print "stddev larger than mean", runawaylist, IC50_list, IC50_avg,IC50_stddev
        #print "stddev larger than mean", IC50_list, IC50_avg,IC50_stddev    
    elif numpy.min(IC50_list) < minimumvalue or numpy.max(IC50_list) > maximumvalue:
        pass
        #print "LARGER than 3fold_sigma", IC50_list, IC50_avg,IC50_stddev,minimumvalue,maximumvalue
    else:
        pass
        #print "deviation from mean in reasonable range. " ##,IC50_list, IC50_avg,IC50_stddev,minimumvalue,maximumvalue
        cpds_out.write(IC50_dict[ent][0])

cpds_out.close()

f_in = open(file_out, 'rb')
file_name_gz = file_out+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.remove(file_out)


## combine_2sdgz.py

import gzip,sys,os
from rdkit import Chem

inF1  = sys.argv[1]
inF2  = sys.argv[2]

if inF1.endswith(".sdf.gz") or inF1.endswith(".sd.gz"):
    cpds1 = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds1 = [x for x in  Chem.SDMolSupplier(inF1) if x is not None]
if inF2.endswith(".sdf.gz") or inF2.endswith(".sd.gz"):
    cpds2 = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[2])) if x is not None]
else:
    cpds2 = [x for x in  Chem.SDMolSupplier(inF2) if x is not None]


file_out = inF1.split(".")[0]+".ready2model"
test_cpd_out = Chem.SDWriter(file_out)

for cpd1 in cpds1:
    test_cpd_out.write(cpd1)
for cpd2 in cpds2:
    test_cpd_out.write(cpd2)

test_cpd_out.close()

f_in = open(file_out, 'rb')
file_name_gz = file_out+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
f_out.close()
f_in.close()
os.remove(file_out)


## set_hERG_TL.py

import sys,gzip,os
from rdkit import Chem

inF       = sys.argv[1]
SD_tag    = sys.argv[2]
cutoff    = sys.argv[3]
file01_wr = inF+".2class_hERGTL.CUT_%snM.sdf" %(cutoff)
ms_wr01   = Chem.SDWriter(file01_wr)

def main(argv=[__name__]):
    if len(argv) != 4:
        print (" \n\nUSAGE: %s <SD file> <IC50 flag> <cutoff value in nM>\n" % argv[0])
        sys.exit()
    if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
        cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
    else:
        cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]
    for mol in cpds:
        # highly active cpds -> class1 (smaller than 1 mikroMol)
        if float(mol.GetProp(SD_tag))> float(float(cutoff)):
            mol.SetProp('hERG_TL','0')
        else:
            mol.SetProp('hERG_TL','1')
        ms_wr01.write(mol)

    ms_wr01.close()
    f_in = open(file01_wr, 'rb')
    file_name_gz = file01_wr+".gz"
    f_out = gzip.open(file_name_gz, 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    os.remove(file01_wr)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

## remove_desciptors.py

import sys,gzip,os
from rdkit import Chem

inF            = sys.argv[1]
SDFilename_out = inF+'.removedSDtags.sdf'
file_SDtags =  sys.argv[2]

def main(argv=[__name__]):
    if len(argv) != 3:
        print (" \n\nUSAGE: %s <SD file> <TXT file with SDtags2remove, line-wise!>\n" % argv[0])
    if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
        SDFile = Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1]))
    else:
        SDFile = Chem.SDMolSupplier(inF)    
    SDTags2remove = open(sys.argv[2],'r').readlines()
    SDFileout = Chem.SDWriter(SDFilename_out)
    SDTags2remove_rstripped = []
    for ent in SDTags2remove:
        SDTags2remove_rstripped.append(ent.rstrip())
    for mol in SDFile:
        for tags2remove in SDTags2remove_rstripped:
            if mol.HasProp(tags2remove):
                mol.ClearProp(tags2remove)
        SDFileout.write(mol)
    SDFileout.close()
    f_in = open(SDFilename_out, 'rb')
    file_name_gz = SDFilename_out+".sdf.gz"
    f_out = gzip.open(file_name_gz, 'wb')
    f_out.writelines(f_in)
    os.remove(SDFilename_out)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

# calc_descr.py

import sys,gzip,os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

outF     = inF+".descr.sdf"
cpds_out = Chem.SDWriter(outF)
nms=[x[0] for x in Descriptors._descList]
nms.remove('MolecularFormula')
print len(Descriptors._descList)

calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
for i in range(len(cpds)):
    descrs = calc.CalcDescriptors(cpds[i])
    for x in range(len(descrs)):
        cpds[i].SetProp(str(nms[x]),str(descrs[x]))
    cpds_out.write(cpds[i])
cpds_out.close()
f_in = open(outF, 'rb')
file_name_gz = outF+".sdf.gz"
f_out = gzip.open(file_name_gz, 'wb')
f_out.writelines(f_in)
os.remove(outF)

## 09_train_models.py

from rdkit import Chem
from scipy import interp
import numpy as np
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.cross_validation import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_score,recall_score
from sklearn import preprocessing
import os,sys,csv,cPickle,gzip
from time import time
t0 = time()
inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]
csv_file        = open(sys.argv[1]+".CV.100est.csv","wb")
csv_file_writer = csv.writer(csv_file,delimiter=";",quotechar=' ')
title_line = ["accuracy","MCC","precision","recall","f1","auc","assaytype","threshold","randomseed"]
csv_file_writer.writerow(title_line)
directory = os.getcwd().split("/")[-2:]
dir_string  = ';'.join(directory)
#
hERG_TL_list = []
property_list_list = []
for cpd in cpds:
    property_list = []
    property_name_list = []
    prop_name = cpd.GetPropNames()
    for property in prop_name:
        if property not in ['hERG_TL','value']:
            property_list.append(float(cpd.GetProp(property)))
            property_name_list.append(property)
        elif property == 'hERG_TL':
            hERG_TL_list.append(int(cpd.GetProp(property)))
        else:
            pass
    property_list_list.append(property_list)
dataDescrs_array = np.asarray(property_list_list)
dataActs_array   = np.array(hERG_TL_list)

for randomseedcounter in range(1,11):
    print "seed", randomseedcounter
    X_train,X_test,y_train,y_test = cross_validation.train_test_split(dataDescrs_array,dataActs_array,test_size=.4,random_state=randomseedcounter)
    ##
    # RF
    clf_RF     = RandomForestClassifier(compute_importances=True,n_estimators=100,random_state=randomseedcounter)
    clf_RF     = clf_RF.fit(X_train,y_train)
    cv_counter = 5
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.zero_one_score)
    accuracy_CV = round(scores.mean(),3)
    accuracy_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.matthews_corrcoef)
    MCC_CV = round(scores.mean(),3)
    MCC_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.f1_score)
    scores_rounded = [round(x,3) for x in scores]
    f1_CV = round(scores.mean(),3)
    f1_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.precision_score)
    scores_rounded = [round(x,3) for x in scores]
    precision_CV = round(scores.mean(),3)
    precision_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.recall_score)
    scores_rounded = [round(x,3) for x in scores]
    recall_CV = round(scores.mean(),3)
    recall_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.auc_score)
    scores_rounded = [round(x,3) for x in scores]
    auc_CV = round(scores.mean(),3)
    auc_std_CV = round(scores.std(),3)
    result_string_cut = [str(accuracy_CV)+"_"+str(accuracy_std_CV),
                         str(MCC_CV)+"_"+str(MCC_std_CV),
                         str(precision_CV)+"_"+str(precision_std_CV),
                         str(recall_CV)+"_"+str(recall_std_CV),
                         str(f1_CV)+"_"+str(f1_std_CV),
                         str(auc_CV)+"_"+str(auc_std_CV),
                         dir_string,randomseedcounter]
    csv_file_writer.writerow(result_string_cut)
t1 = time()-t0
print "done in %1.2f seconds" %(t1)


## 010_train_models.Yscrambled.py

from rdkit import Chem
from scipy import interp
import numpy as np
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.cross_validation import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_score,recall_score
from sklearn import preprocessing
import os,sys,csv,cPickle,gzip
from time import time
t0 = time()
inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(sys.argv[1])) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]
csv_file        = open(sys.argv[1]+".CV.100est.csv","wb")
csv_file_writer = csv.writer(csv_file,delimiter=";",quotechar=' ')
title_line = ["accuracy","MCC","precision","recall","f1","auc","assaytype","threshold","randomseed"]
csv_file_writer.writerow(title_line)
directory = os.getcwd().split("/")[-2:]
dir_string  = ';'.join(directory)
#
hERG_TL_list = []
property_list_list = []
for cpd in cpds:
    property_list = []
    property_name_list = []
    prop_name = cpd.GetPropNames()
    for property in prop_name:
        if property not in ['hERG_TL','value']:
            property_list.append(float(cpd.GetProp(property)))
            property_name_list.append(property)
        elif property == 'hERG_TL':
            hERG_TL_list.append(int(cpd.GetProp(property)))
        else:
            pass
    property_list_list.append(property_list)
dataDescrs_array = np.asarray(property_list_list)
dataActs_array   = np.array(hERG_TL_list)

for randomseedcounter in range(1,11):
    print "seed", randomseedcounter
    X_train,X_test,y_train,y_test = cross_validation.train_test_split(dataDescrs_array,dataActs_array,test_size=.4,random_state=randomseedcounter)
    ##
    # RF
    clf_RF     = RandomForestClassifier(compute_importances=True,n_estimators=100,random_state=randomseedcounter)
    clf_RF     = clf_RF.fit(X_train,y_train)
    cv_counter = 5
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.zero_one_score)
    accuracy_CV = round(scores.mean(),3)
    accuracy_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.matthews_corrcoef)
    MCC_CV = round(scores.mean(),3)
    MCC_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.f1_score)
    scores_rounded = [round(x,3) for x in scores]
    f1_CV = round(scores.mean(),3)
    f1_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.precision_score)
    scores_rounded = [round(x,3) for x in scores]
    precision_CV = round(scores.mean(),3)
    precision_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.recall_score)
    scores_rounded = [round(x,3) for x in scores]
    recall_CV = round(scores.mean(),3)
    recall_std_CV = round(scores.std(),3)
    scores = cross_validation.cross_val_score( clf_RF, X_train,y_train, cv=cv_counter,score_func=metrics.auc_score)
    scores_rounded = [round(x,3) for x in scores]
    auc_CV = round(scores.mean(),3)
    auc_std_CV = round(scores.std(),3)
    result_string_cut = [str(accuracy_CV)+"_"+str(accuracy_std_CV),
                         str(MCC_CV)+"_"+str(MCC_std_CV),
                         str(precision_CV)+"_"+str(precision_std_CV),
                         str(recall_CV)+"_"+str(recall_std_CV),
                         str(f1_CV)+"_"+str(f1_std_CV),
                         str(auc_CV)+"_"+str(auc_std_CV),
                         dir_string,randomseedcounter]
    csv_file_writer.writerow(result_string_cut)
t1 = time()-t0
print "done in %1.2f seconds" %(t1)

## 011_predict_oob.py

from rdkit import Chem
from scipy import interp
import numpy as np
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.cross_validation import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_score,recall_score
from sklearn import preprocessing
import os,sys,csv,cPickle,gzip
import matplotlib
import matplotlib.pyplot as plt

from pylab import *

inF = "./bindingdata.ready2model.sdf.gz.2class_hERGTL.CUT_10000nM.sdf.gz.removedSDtags.sdf.sdf.gz.descr.sdf.sdf.gz"
#inF = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

hERG_TL_list = []
property_list_list = []
for cpd in cpds:
  property_list = []
  property_name_list = []
  prop_name = cpd.GetPropNames()
  for property in prop_name:
    if property not in ['hERG_TL','value']:
      property_list.append(float(cpd.GetProp(property)))
      property_name_list.append(property)
    elif property == 'hERG_TL':
      hERG_TL_list.append(int(cpd.GetProp(property)))
    else:
      pass
  property_list_list.append(property_list)

dataDescrs_array = np.asarray(property_list_list)
dataActs_array   = np.array(hERG_TL_list)

# Den Zaehler benoetige ich, wenn ich mehrere train/test splits durchfuehre 
randomseedcounter = 1
X_train,X_test,y_train,y_test = train_test_split(\
dataDescrs_array,dataActs_array,test_size=.4,random_state=randomseedcounter)
##
# RF
clf_RF     = RandomForestClassifier(compute_importances=True,random_state=0,oob_score=True,n_estimators=100)
clf_RF     = clf_RF.fit(X_train,y_train)
#    print "oob_decision_function_", len(clf_RF.oob_decision_function_)
<
nbins = 10
binhits = [0 for i in range(nbins)]
bin_actives = [0 for i in range(nbins)]
for score,cl in zip(clf_RF.oob_decision_function_,y_train):
  bin_idx = int(nbins * score[1]-0.000000001)
  binhits[bin_idx] += 1
  if cl == 1: bin_actives[bin_idx] += 1

y_frac_active = [float(i)/j for i,j in zip(bin_actives,binhits)]

y_frac_inactive = [float(j-i)/j for i,j in zip(bin_actives,binhits)]

x_active   = [float(i)/10 for i in range(10)]
x_inactive = [float(i)/10 for i in range(9,-1,-1)] #range(9,-1,-1)

print "train set: %i inactives || %i actives \n" %(np.bincount(y_train)[0],np.bincount(y_train)[1])

print "bin_actives %s  %s" %(bin_actives,len(bin_actives))
print "bin_hits    %s  %s" %(binhits,len(binhits))
print "x_inactive ", x_inactive
print "x_active   ", x_active
#
#
print "\n\n"
y_frac_active_fake   = [float(i) for i,j in zip(bin_actives,binhits)]
y_frac_inactive_fake = [float(j-i) for i,j in zip(bin_actives,binhits)]
print "y_frac_active_fake   ", y_frac_active_fake
print "y_frac_inactive_fake ", y_frac_inactive_fake

fig = plt.figure(figsize=(10,5))
#http://matplotlib.1069221.n5.nabble.com/title-when-using-subplot-td22151.html
fig.text(.5, .95, 'train set - out of bag probabilties', \
         horizontalalignment='center', fontsize=14)
ax2 = fig.add_subplot(121)
ax2.set_xlim([0,1])
ax2.set_ylim([0,1])
plt.scatter(x_inactive,y_frac_inactive,marker="o",s=40)
ax2.set_xlabel("relative inactives count   [x_inactive]",fontsize=10)
ax2.set_ylabel("fraction of true_inactives [y_frac_inactive]",fontsize=10)
ax1 = fig.add_subplot(122)
ax1.set_xlim([0,1])
ax1.set_ylim([0,1])
plt.scatter(x_active,y_frac_active,marker="o",s=40)
ax1.set_xlabel("relative actives count   [x_active]",fontsize=10)
ax1.set_ylabel("fraction of true_actives [y_frac_active]",fontsize=10)
plt.show()


## 013A_plot_histo_chembl_assayid.py

from rdkit import Chem
import cPickle,gzip
from pylab import *

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

all_struct_dict = {}
all_struct_list = []
for cpd in cpds:
    chemblassayid = cpd.GetProp("assay__chemblid")
    try:
        all_struct_dict[chemblassayid].append(cpd)
    except:
        all_struct_dict[chemblassayid] = [cpd]
    all_struct_list.append(chemblassayid)

new_dict = {}
new_list = []
for ent in all_struct_dict:
    new_dict[ent] = len(all_struct_dict[ent])
    new_list.append(len(all_struct_dict[ent]))

'''
# binding assay
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(111)
plt.hist(new_list,bins=20)
ax1.set_title("measurements \n per CHEMBL assay ID",fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel("number of measurements",fontsize=24)
ax1.set_ylabel("number of compounds",fontsize=24)
plt.show()

'''
    
'''
# functional assay
fig = plt.figure(figsize=(10,10))

ax1 = fig.add_subplot(121)
ax1.set_title("measurements",fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel("number of",fontsize=24)
ax1.set_ylabel("number of compounds",fontsize=24)
bins1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
#lt.yticks([30,40])
#plt.yticks([5,10,15,20,25,30,35,40])
plt.hist(new_list,bins1)
ax2 = fig.add_subplot(122)
ax2.set_title("per CHEMBL assay ID",fontsize=32)
ax2.tick_params(axis='both', which='major', labelsize=20)
ax2.set_xlabel("measurements",fontsize=24)
bins2 = [600,610,620,630,640,650,660,670,680,690,700]
plt.hist(new_list,bins2)
'''
# functional assay
fig = plt.figure(figsize=(10,10))

###ax1 = fig.add_subplot(111)
###bins1=[1,2,3,4,5,6,7,8,9,10]
###ax1.hist(new_list,bins1,rwidth=0.05,align='left')


ax2 = fig.add_subplot(111)
bins2=[673,675]
ax2.hist(new_list,bins2,rwidth=0.05)

plt.show()

print new_list 


## 013B_plot_histo_measurements_per_cpd.py

from rdkit import Chem
import cPickle,gzip
from pylab import *

# CAVE
# run this on uniq.func.combined/func.uniq.dupl.sdf.gz

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]

all_struct_dict = {}
for cpd in cpds:
    cansmi = cpd.GetProp('cansmirdkit')
    try:
        all_struct_dict[cansmi].append(cpd)
    except:
        all_struct_dict[cansmi] = [cpd]

all_struct_list = []
for ent in all_struct_dict:
    for blub in all_struct_dict[ent]:
        all_struct_list.append(len(all_struct_dict[ent]))
        

fig = plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(111)
plt.xticks(range(1,20))
plt.hist(all_struct_list,bins=20)

ax1.set_title("measurements per compound",fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel("compound count",fontsize=24)
ax1.set_ylabel("measurement count",fontsize=24)
show()


## 013C_pIC50_distro.py

from rdkit import Chem
import cPickle,gzip
from pylab import *

# CAVE
# run this on uniq.func.combined/func.uniq.dupl.sdf.gz

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]


hERG_TL_list = []
property_list_list = []
cpd_name_list = []
for cpd in cpds:
    property_list = []
    property_name_list = []
    prop_name = cpd.GetPropNames()
    cpd_name_list.append(cpd.GetProp("_Name"))
    for property in prop_name:
        if property in ['value']:
            pIC50 = -log10(float(cpd.GetProp(property))*1e-9)
        else:
            pass
    property_list_list.append(pIC50)

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(111)
plt.hist(property_list_list)
#ax1.set_title("pIC50 distribution",fontsize=32)
ax1.tick_params(axis='both', which='major', labelsize=20)
ax1.set_xlabel("pIC50",fontsize=24)
ax1.set_ylabel("compound count",fontsize=24)
show()

## 014_compare_bind_func.py

from rdkit import Chem
import sys,gzip,os



bind_inF  = sys.argv[1]
func_inF  = sys.argv[2]
outF      = Chem.SDWriter("merged.sdf")

if bind_inF.endswith(".sdf.gz") or bind_inF.endswith(".sd.gz"):
    bind_cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(bind_inF)) if x is not None]
else:
    bind_cpds = [x for x in  Chem.SDMolSupplier(bind_inF) if x is not None]

if func_inF.endswith(".sdf.gz") or func_inF.endswith(".sd.gz"):
    func_cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(func_inF)) if x is not None]
else:
    func_cpds = [x for x in  Chem.SDMolSupplier(func_inF) if x is not None]


value_bind_dict = {}
for bind_cpd in bind_cpds:
    cansmi = bind_cpd.GetProp("cansmirdkit")
    value_bind = bind_cpd.GetProp("value")
    bind_cpd.SetProp("value_bind",value_bind)
    bind_cpd.ClearProp("value")
    hERG_TL_bind = bind_cpd.GetProp("hERG_TL")
    bind_cpd.SetProp("hERG_TL_bind",hERG_TL_bind)
    bind_cpd.ClearProp("hERG_TL")
    try:
        value_avg_bind = bind_cpd.GetProp("value_avg")
        bind_cpd.SetProp("value_avg_bind",value_avg_bind)
        bind_cpd.ClearProp("value_avg")
    except:
        pass
    try:
        value_stddev_bind = bind_cpd.GetProp("value_stddev")
        bind_cpd.SetProp("value_stddev_bind",value_stddev_bind)
        bind_cpd.ClearProp("value_stddev")
    except:
        pass

    value_bind_dict[cansmi] = bind_cpd
    
value_func_dict = {}
for func_cpd in func_cpds:
    cansmi = func_cpd.GetProp("cansmirdkit")
    if value_bind_dict.has_key(cansmi):
        value_func_dict[cansmi] = value_bind_dict[cansmi]
        value_func = func_cpd.GetProp("value")

        hERG_TL_func = func_cpd.GetProp("hERG_TL")
        value_bind_dict[cansmi].SetProp("hERG_TL_func",hERG_TL_func)
        value_bind_dict[cansmi].SetProp("value_func",value_func)
        
        #func_cpd.SetProp("value_func",value_func)
        #func_cpd.ClearProp("value")
        try:
            value_avg_func = bind_cpd.GetProp("value_avg")
            func_cpd.SetProp("value_avg_func",value_avg_func)
            func_cpd.ClearProp("value_avg")
        except:
            pass
        try:
            value_stddev_func = func_cpd.GetProp("value_stddev")
            func_cpd.SetProp("value_stddev_func",value_stddev_func)
            func_cpd.ClearProp("value_stddev")
        except:
            pass

print len(value_func_dict)

for x in value_func_dict:
    outF.write(value_func_dict[x])
outF.close()


## 015_exemplary_boxplot.py

from rdkit import Chem
import cPickle,gzip
import numpy as np
from matplotlib import pyplot
from pylab import *
from scipy import stats

filename_inactives = sys.argv[1]
filename_actives = sys.argv[2]
cpds_actives = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(filename_actives)) if x is not None]
cpds_inactives = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(filename_inactives)) if x is not None]

hERG_TL_list_actives = []
property_list_list_actives = []
cpd_name_list_actives = []
prop_dict_actives_all = {}
for cpd in cpds_actives:
    property_list_actives = []
    property_name_list_actives = []
    prop_name_actives = cpd.GetPropNames()
    cpd_name_list_actives.append(cpd.GetProp("_Name"))
    for property in prop_name_actives:
            if property not in prop_dict_actives_all:
                prop_dict_actives_all[property] = []
                prop_dict_actives_all[property].append(float(cpd.GetProp(property)))
            else:
                prop_dict_actives_all[property].append(float(cpd.GetProp(property)))
                property_name_list_actives.append(property)
#
hERG_TL_list_inactives = []
property_list_list_inactives = []
cpd_name_list_inactives = []
prop_dict_inactives_all = {}
for cpd in cpds_inactives:
    property_list_inactives = []
    property_name_list_inactives = []
    prop_name_inactives = cpd.GetPropNames()
    cpd_name_list_inactives.append(cpd.GetProp("_Name"))
    for property in prop_name_inactives:
            if property not in prop_dict_inactives_all:
                list_element = cpd.GetProp(property)
                prop_dict_inactives_all[property] = []
                prop_dict_inactives_all[property].append(float(cpd.GetProp(property)))
            else:
                prop_dict_inactives_all[property].append(float(cpd.GetProp(property)))
                property_name_list_inactives.append(property)
intersectionlist = ['NumHAcceptors']
counter = 0
data = []
label_list = []
for relevant_descr in intersectionlist:
    data.append([prop_dict_actives_all[relevant_descr],prop_dict_inactives_all[relevant_descr]])
    label_list.append(["actives","inactives"])
    counter += 1
    
for ii,x in enumerate(data):
    print intersectionlist[ii]
    medianlist = []
    for i,value in enumerate(x):
        print label_list[ii][i]
        print "%1.1f // %1.1f // %1.1f" %(round(stats.scoreatpercentile(x[i],25),1),round(median(x[i]),1),round(stats.scoreatpercentile(x[i],75),1))


fig = plt.figure(figsize=(10,10))

for x,i in enumerate(intersectionlist):
    plotname = "ax"+str(x)
    subplotnumbering = int("11"+str(x+1))
    plotname = fig.add_subplot(subplotnumbering)
    xtickNames = plt.setp(plotname, xticklabels=label_list[x])
    plt.setp(xtickNames, rotation=45, fontsize=24)
    bp = plt.boxplot(data[x], notch=0, sym='+', vert=1, whis=1.5)
    plotname.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)
    plotname.set_title(intersectionlist[x],fontsize=24)
    plotname.tick_params(axis='both', which='major', labelsize=20)

show()


## 016_mcss.py

import sys
import re
from rdkit import Chem
import hacked_fmcs
from time import time

import argparse
parser = argparse.ArgumentParser("Find fragments")
parser.add_argument("filename",
                    help="SMILES filename")
parser.add_argument("--threshold", "-t", type=float, default=0.5,
                    help="minimum match threshold (default: 0.5)")
parser.add_argument("--min-num-atoms", type=int, default=2,
                    help="minimum number of atom for the match (default: 2)")
parser.add_argument("--min-num-bonds", type=int, default=1,
                    help="minimum number of bonds for the match (default: 1)")
parser.add_argument("--complete-rings-only", action="store_true",
                    help="don't match part ring in a fragment")

args = parser.parse_args()
if args.min_num_atoms < 2:
    parser.error("--min-num-atoms must be at least 2")
if args.min_num_bonds < 1:
    parser.error("--min-num-bonds must be at least 2")
if args.threshold <= 0.0:
    parser.error("--threshold must be positive")

table = Chem.GetPeriodicTable()
_atomic_num_to_symbol = dict((str(i), table.GetElementSymbol(i)) for i in range(105))
atom_pat = re.compile("\[#(\d+)\]")
def atomic_num_to_symbol(m):
    return _atomic_num_to_symbol[m.group(1)]
    
def smarts_to_unique_smiles_fragment(smarts):
    fragment_smiles = atom_pat.sub(atomic_num_to_symbol, smarts)
    fragment_smiles = fragment_smiles.replace("!@", "").replace("@", "")
    mol = Chem.MolFromSmiles(fragment_smiles, sanitize=False)
    # Correct the aromaticity information
    for bond in mol.GetBonds():
        if bond.GetIsAromatic():
            bond.GetBeginAtom().SetIsAromatic(1)
            bond.GetEndAtom().SetIsAromatic(1)
    return Chem.MolToSmiles(mol)

## Hack in code to check the results of the complete_rings_only check

only_has_complete_rings = {}
old_check_complete_rings_only = hacked_fmcs.check_complete_rings_only

def new_check_complete_rings_only(smarts, subgraph, enumeration_mol):
    result = old_check_complete_rings_only(smarts, subgraph, enumeration_mol)
    only_has_complete_rings[smarts] = result
    return result

hacked_fmcs.check_complete_rings_only = new_check_complete_rings_only

hacked_fmcs._maximize_options[("all-matches", False)] = (
    hacked_fmcs.no_pruning, hacked_fmcs.SingleBestAtoms)
hacked_fmcs._maximize_options[("all-matches", True)] = (
    hacked_fmcs.no_pruning, hacked_fmcs.SingleBestAtomsCompleteRingsOnly)

##
# PC
##
t0 = time()

##### Read in the input file

mols = []
for lineno, line in enumerate(open(args.filename)):
    fields = line.split()
    if not fields:
        sys.stderr.write("Missing SMILES on line %d" % (lineno+1,))
        continue
    smiles = fields[0]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.stderr.write("Could not process SMILES %s on line %d\n" % (smiles, lineno+1))
        continue
    mols.append(mol)

mcs_result = hacked_fmcs.fmcs(mols,
                              threshold=args.threshold,
                              min_num_atoms=args.min_num_atoms,
                              complete_rings_only=args.complete_rings_only,
                              # fmcs doesn't implement min_num_bonds; use post-processing
                              maximize="all-matches")

# Convert to unique fragment SMILES
matches = []
unique_fragment_smiles = set()
for smarts, is_match in mcs_result.matches.items():
    if not is_match:  # remember, this is hacked code, and a bit ugly!
        continue
    if args.complete_rings_only:
        if not only_has_complete_rings.get(smarts, False):
            continue

    fragment_smiles = smarts_to_unique_smiles_fragment(smarts)
    if fragment_smiles in unique_fragment_smiles:
        continue
    unique_fragment_smiles.add(fragment_smiles)
    pattern = Chem.MolFromSmarts(smarts)
    if pattern.GetNumBonds() < args.min_num_bonds:
        continue
    matches.append( (pattern, smarts, fragment_smiles) )

def order_by_number_of_atoms(term):
    return term[0].GetNumAtoms()
    
matches.sort(key=order_by_number_of_atoms)


print "frequency #atoms #bonds SMILES SMARTS"
for pattern, smarts, fragment_smiles in matches:
    count = sum(1 for mol in mols if mol.HasSubstructMatch(pattern))
    print count, pattern.GetNumAtoms(), pattern.GetNumBonds(), fragment_smiles, smarts

print "done in %0.3fs" % (time() - t0)


## hacked_fmcs.py

#!/usr/bin/env python

# This work was funded by Roche and generously donated to the free
# and open source cheminformatics community.

## Copyright (c) 2012 Andrew Dalke Scientific AB
## Andrew Dalke <dalke@dalkescientific.com>
##
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""FMCS - Find Maximum Common Substructure

This software finds the maximum common substructure of a set of
structures and reports it as a SMARTS strings.

This implements what I think is a new algorithm for the MCS problem.
The core description is:

  best_substructure = None
  pick one structure as the query, and other as the targets
  for each substructure in the query graph:
    convert it to a SMARTS string based on the desired match properties
    if the SMARTS pattern exists in all of the targets:
       then this is a common substructure
       keep track of the maximum such common structure,

The SMARTS string depends on the desired match properties. For
example, if ring atoms are only allowed to match ring atoms then an
aliphatic ring carbon in the query is converted to the SMARTS "[C;R]",
and the double-bond ring bond converted to "=;@" while the respectice
chain-only version are "[C;!R]" and "=;!@".

The algorithm I outlined earlier will usually take a long time. There
are several ways to speed it up.

== Bond elimination ==

As the first step, remove bonds which obviously cannot be part of the
MCS.

This requires atom and bond type information, which I store as SMARTS
patterns. A bond can only be in the MCS if its canonical bond type is
present in all of the structures. A bond type is string made of the
SMARTS for one atom, the SMARTS for the bond, and the SMARTS for the
other atom. The canonical bond type is the lexographically smaller of
the two possible bond types for a bond.

The atom and bond SMARTS depend on the type comparison used.

The "ring-matches-ring-only" option adds an "@" or "!@" to the bond
SMARTS, so that the canonical bondtype for "C-C" becomes [#6]-@[#6] or
[#6]-!@[#6] if the bond is in a ring or not in a ring, and if atoms
are compared by element and bonds are compared by bondtype. (This
option does not add "R" or "!R" to the atom SMARTS because there
should be a single bond in the MCS of c1ccccc1O and CO.)

The result of all of this atom and bond typing is a "TypedMolecule"
for each input structure.

I then find which canonical bondtypes are present in all of the
structures. I convert each TypedMolecule into a
FragmentedTypedMolecule which has the same atom information but only
those bonds whose bondtypes are in all of the structures. This can
break a structure into multiple, disconnected fragments, hence the
name.

(BTW, I would like to use the fragmented molecules as the targets
because I think the SMARTS match would go faster, but the RDKit SMARTS
matcher doesn't like them. I think it's because the new molecule
hasn't been sanitized and the underlying data structure the ring
information doesn't exist. Instead, I use the input structures for the
SMARTS match.)

== Use the structure with the smallest largest fragment as the query ==
== and sort the targets by the smallest largest fragment             ==

I pick one of the FragmentedTypedMolecule instances as the source of
substructure enumeration. Which one?

My heuristic is to use the one with the smallest largest fragment.
Hopefully it produces the least number of subgraphs, but that's also
related to the number of rings, so a large linear graph will product
fewer subgraphs than a small fused ring system. I don't know how to
quantify that.

For each of the fragmented structures, I find the number of atoms in
the fragment with the most atoms, and I find the number of bonds in
the fragment with the most bonds. These might not be the same
fragment.

I sort the input structures by the number of bonds in the largest
fragment, with ties broken first on the number of atoms, and then on
the input order. The smallest such structure is the query structure,
and the remaining are the targets.

== Use a breadth-first search and a priority queue to    ==
== enumerate the fragment subgraphs                      ==

I extract each of the fragments from the FragmentedTypedMolecule into
a TypedFragment, which I use to make an EnumerationMolecule. An
enumeration molecule contains a pair of directed edges for each atom,
which simplifies the enumeration algorithm.

The enumeration algorithm is based around growing a seed. A seed
contains the current subgraph atoms and bonds as well as an exclusion
set of bonds which cannot be used for future grown. The initial seed
is the first bond in the fragment, which may potentially grow to use
the entire fragment. The second seed is the second bond in the
fragment, which is excluded from using the first bond in future
growth. The third seed starts from the third bond, which may not use
the first or second bonds during growth, and so on.


A seed can grow along bonds connected to an atom in the seed but which
aren't already in the seed and aren't in the set of excluded bonds for
the seed. If there are no such bonds then subgraph enumeration ends
for this fragment. Given N bonds there are 2**N-1 possible ways to
grow, which is just the powerset of the available bonds, excluding the
no-growth case.

This breadth-first growth takes into account all possibilties of using
the available N bonds so all of those bonds are added to the exclusion
set of the newly expanded subgraphs.

For performance reasons, the bonds used for growth are separated into
'internal' bonds, which connect two atoms already in the subgraph, and
'external' bonds, which lead outwards to an atom not already in the
subgraph.

Each seed growth can add from 0 to N new atoms and bonds. The goal is
to maximize the subgraph size so the seeds are stored in a priority
queue, ranked so the seed with the most bonds is processed first. This
turns the enumeration into something more like a depth-first search.


== Prune seeds which aren't found in all of the structures ==

At each stage of seed growth I check that the new seed exists in all
of the original structures. (Well, all except the one which I
enumerate over in the first place; by definition that one will match.)
If it doesn't match then there's no reason to include this seed or any
larger seeds made from it.

The check is easy; I turn the subgraph into its corresponding SMARTS
string and use RDKit's normal SMARTS matcher to test for a match.

There are three ways to generate a SMARTS string: 1) arbitrary, 2)
canonical, 3) hybrid.

I have not tested #1. During most of the development I assumed that
SMARTS matches across a few hundred structures would be slow, so that
the best solution is to generate a *canonical* SMARTS and cache the
match information.

Well, it turns out that my canonical SMARTS match code takes up most
of the FMCS run-time. If I drop the canonicalization step then the
code averages about 5-10% faster. This isn't the same as #1 - I still
do the initial atom assignment based on its neighborhood, which is
like a circular fingerprint of size 2 and *usually* gives a consistent
SMARTS pattern, which I can then cache.

However, there are times when the non-canonical SMARTS code is slower.
Obviously one is if there are a lot of structures, and another if is
there is a lot of symmetry. I'm still working on characterizing this.


== Maximize atoms? or bonds? ==

The above algorithm enumerates all subgraphs of the query and
identifies those subgraphs which are common to all input structures.

It's trivial then to keep track of the current "best" subgraph, which
can defined as having the subgraph with the most atoms, or the most
bonds. Both of those options are implemented.

It would not be hard to keep track of all other subgraphs which are
the same size.

== --complete-ring-only implementation ==

The "complete ring only" option is implemented by first enabling the
"ring-matches-ring-only" option, as otherwise it doesn't make sense.

Second, in order to be a "best" subgraph, all bonds in the subgraph
which are ring bonds in the original molecule must also be in a ring
in the subgraph. This is handled as a post-processing step.

(Note: some possible optimizations, like removing ring bonds from
structure fragments which are not in a ring, are not yet implemented.)


== Prune seeds which have no potential for growing large enough  ==

Given a seed, its set of edges available for growth, and the set of
excluded bonds, figure out the maximum possible growth for the seed.
If this maximum possible is less than the current best subgraph then
prune.

This requires a graph search, currently done in Python, which is a bit
expensive. To speed things up, I precompute some edge information.
That is, if I know that a given bond is a chain bond (not in a ring)
then I can calculate the maximum number of atoms and bonds for seed
growth along that bond, in either direction. However, precomputation
doesn't take into account the excluded bonds, so after a while the
predicted value is too high.

Again, I'm still working on characterizing this, and an implementation
in C++ would have different tradeoffs.
"""

__version__ = "1.1"
__version_info = (1, 1, 0)

import sys

try:
    from rdkit import Chem
except ImportError:
    sys.stderr.write("Please install RDKit from http://www.rdkit.org/\n")
    raise


import argparse
import copy
import itertools
import re
import weakref
from heapq import heappush, heappop, heapify
from itertools import chain, combinations, product
import collections
from collections import defaultdict
import time


### A place to set global options
# (Is this really useful?)

class Default(object):
    timeout = None
    timeout_string = "none"
    maximize = "bonds"
    atom_compare = "elements"
    bond_compare = "bondtypes"
    match_valences = False
    ring_matches_ring_only = False
    complete_rings_only = False


####### Atom type and bond type information #####

# Lookup up the atomic symbol given its atomic number
_get_symbol = Chem.GetPeriodicTable().GetElementSymbol

# Lookup table to get the SMARTS for an atom given its element
# This uses the '#<n>' notation for atoms which may be aromatic.
# Eg, '#6' for carbon, instead of 'C,c'.
# Use the standard element symbol for atoms which can't be aromatic.
class AtomSmartsNoAromaticity(dict):
    def __missing__(self, eleno):
        value = _get_symbol(eleno)
        self[eleno] = value
        return value
_atom_smarts_no_aromaticity = AtomSmartsNoAromaticity()
# Initialize to the ones which need special treatment
# RDKit supports b, c, n, o, p, s, se, and te.
# Daylight and OpenSMILES don't 'te' but do support 'as'
# I don't want 'H'-is-hydrogen to get confused with 'H'-as-has-hydrogens.
# For better portability, I use the '#' notation for all of them.
for eleno in (1, 5, 6, 7, 8, 15, 16, 33, 34, 52):
    _atom_smarts_no_aromaticity[eleno] = "#" + str(eleno)
assert _atom_smarts_no_aromaticity[6] == "#6"
assert _atom_smarts_no_aromaticity[2] == "He"

# Match any atom
def atom_typer_any(atoms):
    return ["*"] * len(atoms)

# Match atom by atomic element; usually by symbol
def atom_typer_elements(atoms):
    return [_atom_smarts_no_aromaticity[atom.GetAtomicNum()] for atom in atoms]

# Match atom by isotope number. This depends on the RDKit version
if hasattr(Chem.Atom, "GetIsotope"):
    def atom_typer_isotopes(atoms):
        return ["%d*" % atom.GetIsotope() for atom in atoms]
else:
    # Before mid-2012, RDKit only supported atomic mass, not isotope.
    # [12*] matches atoms whose mass is 12.000 +/-  0.5/1000
    # This generally works, excepting elements which have no
    #    Tc, Pm, Po, At, Rn, Fr, Ra, Ac, Np, Pu, Am, Cm,
    #    Bk, Cf, Es, Fm, Md, No, Lr
    # natural abundance; [98Tc] is the same as [Tc], etc.
    # This leads to problems because I don't have a way to
    # define the SMARTS for "no defined isotope." In SMILES/SMARTS
    # that's supposed to be through isotope 0.
    # The best I can do is force the non-integer masses to 0 and
    # use isotope 0 to match them. That's clumsy, but it gives
    # the expected result.
    def atom_typer_isotopes(atoms):
        atom_smarts_types = []
        for atom in atoms:
            mass = atom.GetMass()
            int_mass = int(round(mass * 1000))
            if int_mass % 1000 == 0:
                # This is close enough that RDKit's match will work
                atom_smarts = "%d*" % (int_mass//1000)
            else:
                # Probably in natural abundance. In any case,
                # there's no SMARTS for this pattern, so force
                # everything to 0.
                atom.SetMass(0.0)  # XX warning; in-place modification of the input!
                atom_smarts = "0*"
            atom_smarts_types.append(atom_smarts)
        return atom_smarts_types
    

# Match any bond
def bond_typer_any(bonds):
    return ["~"] * len(bonds)


# Match bonds based on bond type, including aromaticity

def bond_typer_bondtypes(bonds):
    # Aromaticity matches are important
    bond_smarts_types = []
    for bond in bonds:
        bond_term = bond.GetSmarts()
        if not bond_term:
            # The SMILES "", means "single or aromatic" as SMARTS.
            # Figure out which one.
            if bond.GetIsAromatic():
                bond_term = ':'
            else:
                bond_term = '-'
        bond_smarts_types.append(bond_term)

    return bond_smarts_types


atom_typers = {
    "any": atom_typer_any,
    "elements": atom_typer_elements,
    "isotopes": atom_typer_isotopes,
    }

bond_typers = {
    "any": bond_typer_any,
    "bondtypes": bond_typer_bondtypes,
    }
default_atom_typer = atom_typers[Default.atom_compare]
default_bond_typer = bond_typers[Default.bond_compare]


####### Support code for handling user-defined atom classes

# User-defined atom classes are handled in a round-about fashion.  The
# fmcs code doesn't know atom classes, but it can handle isotopes.
# It's easy to label the atom isotopes and do an "isotopes" atom
# comparison. The hard part is if you want to get the match
# information back using the original structure data, without the
# tweaked isotopes.

# My solution uses "save_isotopes" and "save_atom_classes" to store
# the old isotope information and the atom class assignments (both
# ordered by atom position), associated with the molecule.

# Use "restore_isotopes()" to restore the molecule's isotope values
# from the saved values. Ise "get_selected_atom_classes" to get the
# atom classes used by specified atom indices.


if hasattr(Chem.Atom, "GetIsotope"):
    def get_isotopes(mol):
        return [atom.GetIsotope() for atom in mol.GetAtoms()]
    def set_isotopes(mol, isotopes):
        if mol.GetNumAtoms() != len(isotopes):
            raise ValueError("Mismatch between the number of atoms and the number of isotopes")
        for atom, isotope in zip(mol.GetAtoms(), isotopes):
            atom.SetIsotope(isotope)

else:
    # Backards compatibility. Before mid-2012, RDKit only supported atomic mass, not isotope.
    def get_isotopes(mol):
        return [atom.GetMass() for atom in mol.GetAtoms()]
    def set_isotopes(mol, isotopes):
        if mol.GetNumAtoms() != len(isotopes):
            raise ValueError("Mismatch between the number of atoms and the number of isotopes")
        for atom, isotope in zip(mol.GetAtoms(), isotopes):
            atom.SetMass(isotope)
    
_isotope_dict = weakref.WeakKeyDictionary()
_atom_class_dict = weakref.WeakKeyDictionary()

def save_isotopes(mol, isotopes):
    _isotope_dict[mol] = isotopes

def save_atom_classes(mol, atom_classes):
    _atom_class_dict[mol] = atom_classes

def get_selected_atom_classes(mol, atom_indices):
    atom_classes = _atom_class_dict.get(mol, None)
    if atom_classes is None:
        return None
    return [atom_classes[index] for index in atom_indices]

def restore_isotopes(mol):
    try:
        isotopes = _isotope_dict[mol]
    except KeyError:
        raise ValueError("no isotopes to restore")
    set_isotopes(mol, isotopes)


def assign_isotopes_from_class_tag(mol, atom_class_tag):
    try:
        atom_classes = mol.GetProp(atom_class_tag)
    except KeyError:
        raise ValueError("Missing atom class tag %r" % (atom_class_tag,))
    fields = atom_classes.split()
    if len(fields) != mol.GetNumAtoms():
        raise ValueError("Mismatch between the number of atoms (#%d) and the number of atom classes (%d)" % (
            mol.GetNumAtoms(), len(fields)))
    new_isotopes = []
    for field in fields:
        if not field.isdigit():
            raise ValueError("Atom class %r from tag %r must be a number" % (field, atom_class_tag))
        isotope = int(field)
        if not (1 <= isotope <= 10000):
            raise ValueError("Atom class %r from tag %r must be in the range 1 to 10000" % (field, atom_class_tag))
        new_isotopes.append(isotope)

    save_isotopes(mol, get_isotopes(mol))
    save_atom_classes(mol, new_isotopes)
    set_isotopes(mol, new_isotopes)



### Different ways of storing atom/bond information about the input structures ###

# A TypedMolecule contains the input molecule, unmodified, along with
# atom type, and bond type information; both as SMARTS fragments. The
# "canonical_bondtypes" uniquely charactizes a bond; two bonds will
# match if and only if their canonical bondtypes match. (Meaning:
# bonds must be of equivalent type, and must go between atoms of
# equivalent types.)


class TypedMolecule(object):
    def __init__(self, rdmol, rdmol_atoms, rdmol_bonds, atom_smarts_types,
                 bond_smarts_types, canonical_bondtypes):
        self.rdmol = rdmol

        # These exist as a performance hack. It's faster to store the
        # atoms and bond as a Python list than to do GetAtoms() and
        # GetBonds() again. The stage 2 TypedMolecule does not use
        # these.
        
        self.rdmol_atoms = rdmol_atoms
        self.rdmol_bonds = rdmol_bonds

        # List of SMARTS to use for each atom and bond
        self.atom_smarts_types = atom_smarts_types
        self.bond_smarts_types = bond_smarts_types

        # List of canonical bondtype strings
        self.canonical_bondtypes = canonical_bondtypes

        # Question: Do I also want the original_rdmol_indices?  With
        # the normal SMARTS I can always do the substructure match
        # again to find the indices, but perhaps this will be needed
        # when atom class patterns are fully implemented.


# Start with a set of TypedMolecules. Find the canonical_bondtypes
# which only exist in all them, then fragment each TypedMolecule to
# produce a FragmentedTypedMolecule containing the same atom
# information but containing only bonds with those
# canonical_bondtypes.
        
class FragmentedTypedMolecule(object):
    def __init__(self, rdmol, rdmol_atoms, orig_atoms, orig_bonds,
                 atom_smarts_types, bond_smarts_types, canonical_bondtypes):
        self.rdmol = rdmol
        self.rdmol_atoms = rdmol_atoms
        self.orig_atoms = orig_atoms
        self.orig_bonds = orig_bonds
        # List of SMARTS to use for each atom and bond
        self.atom_smarts_types = atom_smarts_types
        self.bond_smarts_types = bond_smarts_types

        # List of canonical bondtype strings
        self.canonical_bondtypes = canonical_bondtypes

# A FragmentedTypedMolecule can contain multiple fragments. Once I've
# picked the FragmentedTypedMolecule to use for enumeration, I extract
# each of the fragments as the basis for an EnumerationMolecule.
        
class TypedFragment(object):
    def __init__(self, rdmol,
                 orig_atoms, orig_bonds,
                 atom_smarts_types, bond_smarts_types, canonical_bondtypes):
        self.rdmol = rdmol
        self.orig_atoms = orig_atoms
        self.orig_bonds = orig_bonds
        self.atom_smarts_types = atom_smarts_types
        self.bond_smarts_types = bond_smarts_types
        self.canonical_bondtypes = canonical_bondtypes



# The two possible bond types are
#    atom1_smarts + bond smarts + atom2_smarts
#    atom2_smarts + bond smarts + atom1_smarts
# The canonical bond type is the lexically smaller of these two.

def get_canonical_bondtypes(rdmol, bonds, atom_smarts_types, bond_smarts_types):
    canonical_bondtypes = []
    for bond, bond_smarts in zip(bonds, bond_smarts_types):
        atom1_smarts = atom_smarts_types[bond.GetBeginAtomIdx()]
        atom2_smarts = atom_smarts_types[bond.GetEndAtomIdx()]
        if atom1_smarts > atom2_smarts:
            atom1_smarts, atom2_smarts = atom2_smarts, atom1_smarts
        canonical_bondtypes.append("[%s]%s[%s]" % (atom1_smarts, bond_smarts, atom2_smarts))
    return canonical_bondtypes
    

# Create a TypedMolecule using the element-based typing scheme

# TODO: refactor this. It doesn't seem right to pass boolean flags.

def get_typed_molecule(rdmol, atom_typer, bond_typer, match_valences = Default.match_valences,
                       ring_matches_ring_only = Default.ring_matches_ring_only):
    atoms = list(rdmol.GetAtoms())
    atom_smarts_types = atom_typer(atoms)

    # Get the valence information, if requested
    if match_valences:
        new_atom_smarts_types = []
        for (atom, atom_smarts_type) in zip(atoms, atom_smarts_types):
            valence = atom.GetImplicitValence() + atom.GetExplicitValence()
            valence_str = "v%d" % valence
            if "," in atom_smarts_type:
                atom_smarts_type += ";" + valence_str
            else:
                atom_smarts_type += valence_str
            new_atom_smarts_types.append(atom_smarts_type)
        atom_smarts_types = new_atom_smarts_types
        

    # Store and reuse the bond information because I use it twice.
    # In a performance test, the times went from 2.0 to 1.4 seconds by doing this.
    bonds = list(rdmol.GetBonds())
    bond_smarts_types = bond_typer(bonds)
    if ring_matches_ring_only:
        new_bond_smarts_types = []
        for bond, bond_smarts in zip(bonds, bond_smarts_types):
            if bond.IsInRing():
                if bond_smarts == ":":
                    # No need to do anything; it has to be in a ring
                    pass
                else:
                    if "," in bond_smarts:
                        bond_smarts += ";@"
                    else:
                        bond_smarts += "@"
            else:
                if "," in bond_smarts:
                    bond_smarts += ";!@"
                else:
                    bond_smarts += "!@"
                    
            new_bond_smarts_types.append(bond_smarts)
        bond_smarts_types = new_bond_smarts_types

    canonical_bondtypes = get_canonical_bondtypes(rdmol, bonds, atom_smarts_types, bond_smarts_types)
    return TypedMolecule(rdmol, atoms, bonds, atom_smarts_types, bond_smarts_types, canonical_bondtypes)


# Create a TypedMolecule using the user-defined atom classes (Not implemented!)

def get_specified_types(rdmol, atom_types, ring_matches_ring_only):
    raise NotImplementedError("not tested!")
    # Make a copy because I will do some destructive edits
    rdmol = copy.copy(rdmol)
    
    atom_smarts_types = []
    atoms = list(mol.GetAtoms())
    for atom, atom_type in zip(atoms, atom_types):
        atom.SetAtomicNum(0)
        atom.SetMass(atom_type)
        atom_term = "%d*" % (atom_type,)
        if ring_matches_ring_only:
            if atom.IsInRing():
                atom_term += "R"
            else:
                atom_term += "!R"
        atom_smarts_types.append('[' + atom_term + ']')

    bonds = list(rdmol.GetBonds())
    bond_smarts_types = get_bond_smarts_types(mol, bonds, ring_matches_ring_only)
    canonical_bondtypes = get_canonical_bondtypes(mol, bonds, atom_smarts_types, bond_smarts_types)

    return TypedMolecule(mol, atoms, bonds, atom_smarts_types, bond_smarts_types, canonical_bondtypes)


def convert_input_to_typed_molecules(mols, atom_typer, bond_typer, match_valences, ring_matches_ring_only):
    typed_mols = []
    for molno, rdmol in enumerate(mols):
        typed_mol = get_typed_molecule(rdmol, atom_typer, bond_typer,
                                       match_valences=match_valences, ring_matches_ring_only=ring_matches_ring_only)
        typed_mols.append(typed_mol)

    return typed_mols

def _check_atom_classes(molno, num_atoms, atom_classes):
    if num_atoms != len(atom_classes):
        raise ValueError("mols[%d]: len(atom_classes) must be the same as the number of atoms" % (molno,))
    for atom_class in atom_classes:
        if not isinstance(atom_class, int):
            raise ValueError("mols[%d]: atom_class elements must be integers" % (molno,))
        if not (1 <= atom_class < 1000):
            raise ValueError("mols[%d]: atom_class elements must be in the range 1 <= value < 1000" %
                             (molno,))

#############################################

# This section deals with finding the canonical bondtype counts and
# making new TypedMolecule instances where the atoms contain only the
# bond types which are in all of the structures.

# In the future I would like to keep track of the bond types which are
# in the current subgraph. If any subgraph bond type count is ever
# larger than the maximum counts computed across the whole set, then
# prune. But so far I don't have a test set which drives the need for
# that.

# Return a dictionary mapping iterator item to occurence count
def get_counts(it):
    d = defaultdict(int)
    for item in it:
        d[item] += 1
    return dict(d)

# Merge two count dictionaries, returning the smallest count for any
# entry which is in both.
def intersect_counts(counts1, counts2):
    d = {}
    for k, v1 in counts1.iteritems():
        if k in counts2:
            v = min(v1, counts2[k])
            d[k] = v
    return d


# Figure out which canonical bonds SMARTS occur in every molecule
def get_canonical_bondtype_counts(typed_mols):
    overall_counts = defaultdict(list)
    for typed_mol in typed_mols:
        bondtype_counts = get_counts(typed_mol.canonical_bondtypes)
        for k,v in bondtype_counts.items():
            overall_counts[k].append(v)
    return overall_counts

# If I know which bondtypes exist in all of the structures, I can
# remove all bonds which aren't in all structures. RDKit's Molecule
# class doesn't let me edit in-place, so I end up making a new one
# which doesn't have unsupported bond types.

def remove_unknown_bondtypes(typed_mol, supported_canonical_bondtypes):
    emol = Chem.EditableMol(Chem.Mol())

    # Copy all of the atoms, even those which don't have any bonds. 
    for atom in typed_mol.rdmol_atoms:
        emol.AddAtom(atom)

    # Copy over all the bonds with a supported bond type.
    # Make sure to update the bond SMARTS and canonical bondtype lists.
    orig_bonds = []
    new_bond_smarts_types = []
    new_canonical_bondtypes = []
    for bond, bond_smarts, canonical_bondtype in zip(typed_mol.rdmol_bonds, typed_mol.bond_smarts_types,
                                                     typed_mol.canonical_bondtypes):
        if canonical_bondtype in supported_canonical_bondtypes:
            orig_bonds.append(bond)
            new_bond_smarts_types.append(bond_smarts)
            new_canonical_bondtypes.append(canonical_bondtype)
            emol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())

    new_mol = emol.GetMol()
    return FragmentedTypedMolecule(new_mol, list(new_mol.GetAtoms()),
                                   typed_mol.rdmol_atoms, orig_bonds,
                                   typed_mol.atom_smarts_types, new_bond_smarts_types,
                                   new_canonical_bondtypes)

# The molecule at this point has been (potentially) fragmented by
# removing bonds with unsupported bond types. The MCS cannot contain
# more atoms than the fragment of a given molecule with the most
# atoms, and the same for bonds. Find those upper limits. Note that
# the fragment with the most atoms is not necessarily the one with the
# most bonds.

def find_upper_fragment_size_limits(rdmol, atoms):
    max_num_atoms = max_twice_num_bonds = 0
    for atom_indices in Chem.GetMolFrags(rdmol):
        num_atoms = len(atom_indices)
        if num_atoms > max_num_atoms:
            max_num_atoms = num_atoms

        # Every bond is connected to two atoms, so this is the
        # simplest way to count the number of bonds in the fragment.
        twice_num_bonds = 0
        for atom_index in atom_indices:
            # XXX Why is there no 'atom.GetNumBonds()'?
            twice_num_bonds += sum(1 for bond in atoms[atom_index].GetBonds())
        if twice_num_bonds > max_twice_num_bonds:
            max_twice_num_bonds = twice_num_bonds

    return max_num_atoms, max_twice_num_bonds // 2


####### Convert the selected TypedMolecule into an EnumerationMolecule

# I convert one of the typed fragment molecules (specifically, the one
# with the smallest largest fragment score) into a list of
# EnumerationMolecule instances. Each fragment from the typed molecule
# gets turned into an EnumerationMolecule.

# An EnumerationMolecule contains the data I need to enumerate all of
# its subgraphs.

# An EnumerationMolecule contains a list of 'Atom's and list of 'Bond's.
# Atom and Bond indices are offsets into those respective lists.
# An Atom has a list of "bond_indices", which are offsets into the bonds.
# A Bond has a 2-element list of "atom_indices", which are offsets into the atoms.

EnumerationMolecule = collections.namedtuple("Molecule", "rdmol atoms bonds directed_edges")
Atom = collections.namedtuple("Atom", "real_atom atom_smarts bond_indices is_in_ring")
Bond = collections.namedtuple("Bond", "real_bond bond_smarts canonical_bondtype atom_indices is_in_ring")

# A Bond is linked to by two 'DirectedEdge's; one for each direction.
# The DirectedEdge.bond_index references the actual RDKit bond instance.
# 'end_atom_index' is the index of the destination atom of the directed edge
# This is used in a 'directed_edges' dictionary so that
#     [edge.end_atom_index for edge in directed_edges[atom_index]]
# is the list of all atom indices connected to 'atom_index'
DirectedEdge = collections.namedtuple("DirectedEdge",
                                      "bond_index end_atom_index")

# A Subgraph is a list of atom and bond indices in an EnumerationMolecule
Subgraph = collections.namedtuple("Subgraph", "atom_indices bond_indices")

def get_typed_fragment(typed_mol, atom_indices):
    rdmol = typed_mol.rdmol
    rdmol_atoms = typed_mol.rdmol_atoms

    # I need to make a new RDKit Molecule containing only the fragment.
    # XXX Why is that? Do I use the molecule for more than the number of atoms and bonds?

    # Copy over the atoms
    emol = Chem.EditableMol(Chem.Mol())
    atom_smarts_types = []
    atom_map = {}
    for i, atom_index in enumerate(atom_indices):
        atom = rdmol_atoms[atom_index]
        emol.AddAtom(atom)
        atom_smarts_types.append(typed_mol.atom_smarts_types[atom_index])
        atom_map[atom_index] = i

    # Copy over the bonds.
    orig_bonds = []
    bond_smarts_types = []
    new_canonical_bondtypes = []
    for bond, orig_bond, bond_smarts, canonical_bondtype in zip(
                rdmol.GetBonds(), typed_mol.orig_bonds,
                typed_mol.bond_smarts_types, typed_mol.canonical_bondtypes):
        begin_atom_idx = bond.GetBeginAtomIdx()
        end_atom_idx = bond.GetEndAtomIdx()
        count = (begin_atom_idx in atom_map) + (end_atom_idx in atom_map)
        # Double check that I have a proper fragment
        if count == 2:
            bond_smarts_types.append(bond_smarts)
            new_canonical_bondtypes.append(canonical_bondtype)
            emol.AddBond(atom_map[begin_atom_idx], atom_map[end_atom_idx], bond.GetBondType())
            orig_bonds.append(orig_bond)
        elif count == 1:
            raise AssertionError("connected/disconnected atoms?")
    return TypedFragment(emol.GetMol(),
                         [typed_mol.orig_atoms[atom_index] for atom_index in atom_indices],
                         orig_bonds,
                         atom_smarts_types, bond_smarts_types, new_canonical_bondtypes)


def fragmented_mol_to_enumeration_mols(typed_mol, min_num_atoms=2):
    if min_num_atoms < 2:
        raise ValueError("min_num_atoms must be at least 2")

    fragments = []
    for atom_indices in Chem.GetMolFrags(typed_mol.rdmol):
        # No need to even look at fragments which are too small.
        if len(atom_indices) < min_num_atoms:
            continue

        # Convert a fragment from the TypedMolecule into a new
        # TypedMolecule containing only that fragment.

        # You might think I could merge 'get_typed_fragment()' with
        # the code to generate the EnumerationMolecule. You're
        # probably right. This code reflects history. My original code
        # didn't break the typed molecule down to its fragments.
        typed_fragment = get_typed_fragment(typed_mol, atom_indices)
        rdmol = typed_fragment.rdmol
        atoms = []
        for atom, orig_atom, atom_smarts_type in zip(rdmol.GetAtoms(), typed_fragment.orig_atoms,
                                                typed_fragment.atom_smarts_types):
            bond_indices = [bond.GetIdx() for bond in atom.GetBonds()]
            #assert atom.GetSymbol() == orig_atom.GetSymbol()
            atom_smarts = '[' + atom_smarts_type + ']'
            atoms.append(Atom(atom, atom_smarts, bond_indices, orig_atom.IsInRing()))

        directed_edges = collections.defaultdict(list)
        bonds = []
        for bond_index, (bond, orig_bond, bond_smarts, canonical_bondtype) in enumerate(
                zip(rdmol.GetBonds(), typed_fragment.orig_bonds,
                    typed_fragment.bond_smarts_types, typed_fragment.canonical_bondtypes)):
            atom_indices = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
            bonds.append(Bond(bond, bond_smarts, canonical_bondtype, atom_indices, orig_bond.IsInRing()))

            directed_edges[atom_indices[0]].append(DirectedEdge(bond_index, atom_indices[1]))
            directed_edges[atom_indices[1]].append(DirectedEdge(bond_index, atom_indices[0]))

        fragment = EnumerationMolecule(rdmol, atoms, bonds, dict(directed_edges))
        fragments.append(fragment)

    # Optimistically try the largest fragments first
    fragments.sort(key = lambda fragment: len(fragment.atoms), reverse=True)
    return fragments


####### Canonical SMARTS generation using Weininger, Weininger, and Weininger's CANGEN

# CANGEN "combines two separate algorithms, CANON and GENES.  The
# first stage, CANON, labels a molecualr structure with canonical
# labels. ... Each atom is given a numerical label on the basis of its
# topology. In the second stage, GENES generates the unique SMILES
# ... . [It] selects the starting atom and makes branching decisions
# by referring to the canonical labels as needed."


# CANON is based on the fundamental theorem of arithmetic, that is,
# the unique prime factorization theorem. Which means I need about as
# many primes as I have atoms.

# I could have a fixed list of a few thousand primes but I don't like
# having a fixed upper limit to my molecule size. I modified the code
# Georg Schoelly posted at http://stackoverflow.com/a/568618/64618 .
# This is one of many ways to generate an infinite sequence of primes.
def gen_primes():
    d = defaultdict(list)
    q = 2
    while 1:
        if q not in d:
            yield q
            d[q*q].append(q)
        else:
            for p in d[q]:
                d[p+q].append(p)
            del d[q]
        q += 1

_prime_stream = gen_primes()

# Code later on uses _primes[n] and if that fails, calls _get_nth_prime(n)
_primes = []

def _get_nth_prime(n):
    # Keep appending new primes from the stream until I have enough.
    current_size = len(_primes)
    while current_size <= n:
        _primes.append(next(_prime_stream))
        current_size += 1
    return _primes[n]

# Prime it with more values then will likely occur
_get_nth_prime(1000)

###

# The CANON algorithm is documented as:
#  (1) Set atomic vector to initial invariants. Go to step 3.
#  (2) Set vector to product of primes corresponding to neighbors' ranks.
#  (3) Sort vector, maintaining stability over previous ranks.
#  (4) Rank atomic vector.
#  (5) If not invariants partitioning, go to step 2.
#  (6) On first pass, save partitioning as symmetry classes [fmcs doesn't need this]
#  (7) If highest rank is smaller than number of nodes, break ties, go to step 2
#  (8) ... else done.


# I track the atom information as a list of CangenNode instances.

class CangenNode(object):
    # Using __slots__ improves get_initial_cangen_nodes performance by over 10%
    # and dropped my overall time (in one benchmark) from 0.75 to 0.73 seconds
    __slots__ = ["index", "atom_smarts", "value", "neighbors", "rank", "outgoing_edges"]
    def __init__(self, index, atom_smarts):
        self.index = index
        self.atom_smarts = atom_smarts  # Used to generate the SMARTS output
        self.value = 0
        self.neighbors = []
        self.rank = 0
        self.outgoing_edges = []

# The outgoing edge information is used to generate the SMARTS output
# The index numbers are offsets in the subgraph, not in the original molecule
OutgoingEdge = collections.namedtuple("OutgoingEdge",
                                      "from_atom_index bond_index bond_smarts other_node_idx other_node")

# Convert a Subgraph of a given EnumerationMolecule into a list of
# CangenNodes. This contains the more specialized information I need
# for canonicalization and for SMARTS generation.
def get_initial_cangen_nodes(subgraph, enumeration_mol, atom_assignment, do_initial_assignment=True):
    # The subgraph contains a set of atom and bond indices in the enumeration_mol.
    # The CangenNode corresponds to an atom in the subgraph, plus relations
    # to other atoms in the subgraph.
    # I need to convert from offsets in molecule space to offset in subgraph space.

    # Map from enumeration mol atom indices to subgraph/CangenNode list indices
    atom_map = {}

    cangen_nodes = []
    atoms = enumeration_mol.atoms
    canonical_labels = []
    for i, atom_index in enumerate(subgraph.atom_indices):
        atom_map[atom_index] = i
        cangen_nodes.append(CangenNode(i, atoms[atom_index].atom_smarts))
        canonical_labels.append([])

    # Build the neighbor and directed edge lists
     
    for bond_index in subgraph.bond_indices:
        bond = enumeration_mol.bonds[bond_index]
        from_atom_index, to_atom_index = bond.atom_indices
        from_subgraph_atom_index = atom_map[from_atom_index]
        to_subgraph_atom_index = atom_map[to_atom_index]

        from_node = cangen_nodes[from_subgraph_atom_index]
        to_node = cangen_nodes[to_subgraph_atom_index]
        from_node.neighbors.append(to_node)
        to_node.neighbors.append(from_node)

        canonical_bondtype = bond.canonical_bondtype
        canonical_labels[from_subgraph_atom_index].append(canonical_bondtype)
        canonical_labels[to_subgraph_atom_index].append(canonical_bondtype)

        from_node.outgoing_edges.append(
            OutgoingEdge(from_subgraph_atom_index, bond_index, bond.bond_smarts,
                         to_subgraph_atom_index, to_node))
        to_node.outgoing_edges.append(
            OutgoingEdge(to_subgraph_atom_index, bond_index, bond.bond_smarts,
                         from_subgraph_atom_index, from_node))

    if do_initial_assignment:
        # Do the initial graph invariant assignment. (Step 1 of the CANON algorithm)
        # These are consistent only inside of the given 'atom_assignment' lookup.
        for atom_index, node, canonical_label in zip(subgraph.atom_indices, cangen_nodes, canonical_labels):
            # The initial invariant is the sorted canonical bond labels
            # plus the atom smarts, separated by newline characters.
            #
            # This is equivalent to a circular fingerprint of width 2, and
            # gives more unique information than the Weininger method.
            canonical_label.sort()
            canonical_label.append(atoms[atom_index].atom_smarts)
            label = "\n".join(canonical_label)

            # The downside of using a string is that I need to turn it
            # into a number which is consistent across all of the SMARTS I
            # generate as part of the MCS search. Use a lookup table for
            # that which creates a new number of the label wasn't seen
            # before, or uses the old one if it was.
            node.value = atom_assignment[label]

    return cangen_nodes


# Rank a sorted list (by value) of CangenNodes
def rerank(cangen_nodes):
    rank = 0     # Note: Initial rank is 1, in line with the Weininger paper
    prev_value = -1
    for node in cangen_nodes:
        if node.value != prev_value:
            rank += 1
            prev_value = node.value
        node.rank = rank

# Given a start/end range in the CangenNodes, sorted by value,
# find the start/end for subranges with identical values
def find_duplicates(cangen_nodes, start, end):
    result = []
    prev_value = -1
    count = 0
    for index in xrange(start, end):
        node = cangen_nodes[index]
        if node.value == prev_value:
            count += 1
        else:
            if count > 1:
                # New subrange containing duplicates
                result.append( (start, index) )
            count = 1
            prev_value = node.value
            start = index
    if count > 1:
        # Last elements were duplicates
        result.append( (start, end) )
    return result

#@profile 
def canon(cangen_nodes):
    # Precondition: node.value is set to the initial invariant
    # (1) Set atomic vector to initial invariants (assumed on input)
    
    # Do the initial ranking
    cangen_nodes.sort(key = lambda node: node.value)
    rerank(cangen_nodes)

    # Keep refining the sort order until it's unambiguous
    master_sort_order = cangen_nodes[:]

    # Find the start/end range for each stretch of duplicates
    duplicates = find_duplicates(cangen_nodes, 0, len(cangen_nodes))

    PRIMES = _primes # micro-optimization; make this a local name lookup
    
    while duplicates:
        # (2) Set vector to product of primes corresponding to neighbor's ranks
        for node in cangen_nodes:
            try:
                node.value = PRIMES[node.rank]
            except IndexError:
                node.value = _get_nth_prime(node.rank)
        for node in cangen_nodes:
            # Apply the fundamental theorem of arithmetic; compute the
            # product of the neighbors' primes
            p = 1
            for neighbor in node.neighbors:
                p *= neighbor.value
            node.value = p
            

        # (3) Sort vector, maintaining stability over previous ranks
        # (I maintain stability by refining ranges in the
        # master_sort_order based on the new ranking)
        cangen_nodes.sort(key = lambda node: node.value)

        # (4) rank atomic vector
        rerank(cangen_nodes)

        # See if any of the duplicates have been resolved.
        new_duplicates = []
        unchanged = True  # This is buggy? Need to check the entire state XXX
        for (start, end) in duplicates:
            # Special case when there's only two elements to store.
            # This optimization sped up cangen by about 8% because I
            # don't go through the sort machinery
            if start+2 == end:
                node1, node2 = master_sort_order[start], master_sort_order[end-1]
                if node1.value > node2.value:
                    master_sort_order[start] = node2
                    master_sort_order[end-1] = node1
            else:
                subset = master_sort_order[start:end]
                subset.sort(key = lambda node: node.value)
                master_sort_order[start:end] = subset

            subset_duplicates = find_duplicates(master_sort_order, start, end)
            new_duplicates.extend(subset_duplicates)
            if unchanged:
                # Have we distinguished any of the duplicates?
                if not (len(subset_duplicates) == 1 and subset_duplicates[0] == (start, end)):
                    unchanged = False

        # (8) ... else done
        # Yippee! No duplicates left. Everything has a unique value.
        if not new_duplicates:
            break
            
        # (5) If not invariant partitioning, go to step 2
        if not unchanged:
            duplicates = new_duplicates
            continue
        
        duplicates = new_duplicates
        
        # (6) On first pass, save partitioning as symmetry classes
        pass # I don't need this information
        
        # (7) If highest rank is smaller than number of nodes, break ties, go to step 2
        # I follow the Weininger algorithm and use 2*rank or 2*rank-1.
        # This requires that the first rank is 1, not 0.
        for node in cangen_nodes:
            node.value = node.rank * 2

        # The choice of tie is arbitrary. Weininger breaks the first tie.
        # I break the last tie because it's faster in Python to delete
        # from the end than the beginning.
        start, end = duplicates[-1]
        cangen_nodes[start].value -= 1
        if end == start+2:
            # There were only two nodes with the same value. Now there
            # are none. Remove information about that duplicate.
            del duplicates[-1]
        else:
            # The first N-1 values are still duplicates.
            duplicates[-1] = (start+1, end)
        rerank(cangen_nodes)

    # Restore to the original order (ordered by subgraph atom index)
    # because the bond information used during SMARTS generation
    # references atoms by that order.
    cangen_nodes.sort(key=lambda node: node.index)


def get_closure_label(bond_smarts, closure):
    if closure < 10:
        return bond_smarts + str(closure)
    else:
        return bond_smarts + "%%%02d" % closure

# Precompute the initial closure heap. *Overall* performance went from 0.73 to 0.64 seconds!
_available_closures = range(1, 101)
heapify(_available_closures)

# The Weininger paper calls this 'GENES'; I call it "generate_smiles."

# I use a different algorithm than GENES. It's still use two
# passes. The first pass identifies the closure bonds using a
# depth-first search. The second pass builds the SMILES string.

def generate_smarts(cangen_nodes):
    start_index = 0
    best_rank = cangen_nodes[0].rank
    for i, node in enumerate(cangen_nodes):
        if node.rank < best_rank:
            best_rank = node.rank
            start_index = i
        node.outgoing_edges.sort(key=lambda edge: edge.other_node.rank)

    visited_atoms = [0] * len(cangen_nodes)
    closure_bonds = set()

    ## First, find the closure bonds using a DFS
    stack = []
    atom_idx = start_index
    stack.extend(reversed(cangen_nodes[atom_idx].outgoing_edges))
    visited_atoms[atom_idx] = True

    while stack:
        edge = stack.pop()
        if visited_atoms[edge.other_node_idx]:
            closure_bonds.add(edge.bond_index)
        else:
            visited_atoms[edge.other_node_idx] = 1
            for next_edge in reversed(cangen_nodes[edge.other_node_idx].outgoing_edges):
                if next_edge.other_node_idx == edge.from_atom_index:
                    # Don't worry about going back along the same route
                    continue
                stack.append(next_edge)


    available_closures = _available_closures[:]
    unclosed_closures = {}

    # I've identified the closure bonds.
    # Use a stack machine to traverse the graph and build the SMARTS.
    # The instruction contains one of 4 instructions, with associated data
    #   0: add the atom's SMARTS and put its connections on the machine
    #   1: add the bond's SMARTS and put the other atom on the machine
    #   3: add a ')' to the SMARTS
    #   4: add a '(' and the bond SMARTS
    
    smiles_terms = []
    stack = [(0, (start_index, -1))]
    while stack:
        action, data = stack.pop()
        if action == 0:
            # Add an atom.

            # The 'while 1:' emulates a goto for the special case
            # where the atom is connected to only one other atom.  I
            # don't need to use the stack machinery for that case, and
            # can speed up this function by about 10%.
            while 1:
                # Look at the bonds starting from this atom
                num_neighbors = 0
                atom_idx, prev_bond_idx = data
                smiles_terms.append(cangen_nodes[atom_idx].atom_smarts)
                outgoing_edges = cangen_nodes[atom_idx].outgoing_edges
                for outgoing_edge in outgoing_edges:
                    bond_idx = outgoing_edge.bond_index

                    # Is this a ring closure bond?
                    if bond_idx in closure_bonds:
                        # Have we already seen it before?
                        if bond_idx not in unclosed_closures:
                            # This is new. Add as a ring closure.
                            closure = heappop(available_closures)
                            smiles_terms.append(get_closure_label(outgoing_edge.bond_smarts, closure))
                            unclosed_closures[bond_idx] = closure
                        else:
                            closure = unclosed_closures[bond_idx]
                            smiles_terms.append(get_closure_label(outgoing_edge.bond_smarts, closure))
                            heappush(available_closures, closure)
                            del unclosed_closures[bond_idx]
                    else:
                        # This is a new outgoing bond.
                        if bond_idx == prev_bond_idx:
                            # Don't go backwards along the bond I just came in on
                            continue
                        if num_neighbors == 0:
                            # This is the first bond. There's a good chance that
                            # it's the only bond. 
                            data = (outgoing_edge.other_node_idx, bond_idx)
                            bond_smarts = outgoing_edge.bond_smarts
                        else:
                            # There are multiple bonds. Can't shortcut.
                            if num_neighbors == 1:
                                # Capture the information for the first bond
                                # This direction doesn't need the (branch) characters.
                                stack.append((0, data))
                                stack.append((1, bond_smarts))
                            
                            # Add information for this bond
                            stack.append((3, None))
                            stack.append((0, (outgoing_edge.other_node_idx, bond_idx)))
                            stack.append((4, outgoing_edge.bond_smarts))

                        num_neighbors += 1
                if num_neighbors != 1:
                    # If there's only one item then goto action==0 again.
                    break
                smiles_terms.append(bond_smarts)
        elif action == 1:
            # Process a bond which does not need '()'s
            smiles_terms.append(data) # 'data' is bond_smarts
            continue
            
        elif action == 3:
            smiles_terms.append(')')
        elif action == 4:
            smiles_terms.append('(' + data)  # 'data' is bond_smarts
        else:
            raise AssertionError

    return "".join(smiles_terms)


# Full canonicalization is about 5% slower unless there are well over 100 structures
# in the data set, which is not expected to be common.
# Commented out the canon() step until there's a better solution (eg, adapt based
# in the input size.)
def make_canonical_smarts(subgraph, enumeration_mol, atom_assignment):
    cangen_nodes = get_initial_cangen_nodes(subgraph, enumeration_mol, atom_assignment, True)
    #canon(cangen_nodes)
    smarts = generate_smarts(cangen_nodes)
    return smarts

## def make_semicanonical_smarts(subgraph, enumeration_mol, atom_assignment):
##     cangen_nodes = get_initial_cangen_nodes(subgraph, enumeration_mol, atom_assignment, True)
##     # There's still some order because of the canonical bond typing, but it isn't perfect
##     #canon(cangen_nodes)
##     smarts = generate_smarts(cangen_nodes)
##     return smarts

def make_arbitrary_smarts(subgraph, enumeration_mol, atom_assignment):
    cangen_nodes = get_initial_cangen_nodes(subgraph, enumeration_mol, atom_assignment, False)
    # Use an arbitrary order
    for i, node in enumerate(cangen_nodes):
        node.value = i
    smarts = generate_smarts(cangen_nodes)
    return smarts

        
############## Subgraph enumeration ##################

# A 'seed' is a subgraph containing a subset of the atoms and bonds in
# the graph. The idea is to try all of the ways in which to grow the
# seed to make a new seed which contains the original seed.

# There are two ways to grow a seed:
#   - add a bond which is not in the seed but where both of its
#            atoms are in the seed
#   - add a bond which is not in the seed but where one of its
#            atoms is in the seed (and the other is not)

# The algorithm takes the seed, and finds all of both categories of
# bonds. If there are N total such bonds then there are 2**N-1
# possible new seeds which contain the original seed. This is simply
# the powerset of the possible bonds, excepting the case with no
# bonds.

# Generate all 2**N-1 new seeds. Place the new seeds back in the
# priority queue to check for additional growth.

# I place the seeds in priority queue, sorted by score (typically the
# number of atoms) to preferentially search larger structures first. A
# simple stack or deque wouldn't work because the new seeds have
# between 1 to N-1 new atoms and bonds.

    
# Some useful preamble code
    
# Taken from the Python documentation
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# Same as the above except the empty term is not returned
def nonempty_powerset(iterable):
    "nonempty_powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    it = chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    it.next()
    return it


# Call this to get a new unique function. Used to break ties in the
# priority queue.
tiebreaker = itertools.count().next

### The enumeration code


# Given a set of atoms, find all of the ways to leave those atoms.
# There are two possibilities:
#   1) bonds; which connect two atoms which are already in 'atom_indices'
#   2) directed edges; which go to atoms that aren't in 'atom_indices'
#     and which aren't already in visited_bond_indices. These are external
#     to the subgraph.
# The return is a 2-element tuple containing:
#  (the list of bonds from (1), the list of directed edges from (2))
def find_extensions(atom_indices, visited_bond_indices, directed_edges):
    internal_bonds = set()
    external_edges = []
    for atom_index in atom_indices:
        for directed_edge in directed_edges[atom_index]:
            # Skip outgoing edges which have already been evaluated
            if directed_edge.bond_index in visited_bond_indices:
                continue

            if directed_edge.end_atom_index in atom_indices:
                # case 1: This bond goes to another atom which is already in the subgraph.
                internal_bonds.add(directed_edge.bond_index)
            else:
                # case 2: This goes to a new (external) atom
                external_edges.append(directed_edge)

    # I don't think I need the list()
    return list(internal_bonds), external_edges


# Given the 2-element tuple (internal_bonds, external_edges),
# construct all of the ways to combine them to generate a new subgraph
# from the old one. This is done via a powerset.
# This generates a two-element tuple containing:
#   - the set of newly added atom indices (or None)
#   - the new subgraph

def all_subgraph_extensions(enumeration_mol, subgraph, visited_bond_indices, internal_bonds, external_edges):
    #print "Subgraph", len(subgraph.atom_indices), len(subgraph.bond_indices), "X", enumeration_mol.rdmol.GetNumAtoms()
    #print "subgraph atoms", subgraph.atom_indices
    #print "subgraph bonds", subgraph.bond_indices
    #print "internal", internal_bonds, "external", external_edges
    # only internal bonds
    if not external_edges:
        #assert internal_bonds, "Must have at least one internal bond"
        it = nonempty_powerset(internal_bonds)
        for internal_bond in it:
            # Make the new subgraphs
            bond_indices = set(subgraph.bond_indices)
            bond_indices.update(internal_bond)
            yield None, Subgraph(subgraph.atom_indices, frozenset(bond_indices)), 0, 0
        return

    # only external edges
    if not internal_bonds:
        it = nonempty_powerset(external_edges)
        exclude_bonds = set(chain(visited_bond_indices, (edge.bond_index for edge in external_edges)))
        for external_ext in it:
            new_atoms = frozenset(ext.end_atom_index for ext in external_ext)
            atom_indices = frozenset(chain(subgraph.atom_indices, new_atoms))
            bond_indices = frozenset(chain(subgraph.bond_indices,
                                            (ext.bond_index for ext in external_ext)))
            num_possible_atoms, num_possible_bonds = find_extension_size(
                enumeration_mol, new_atoms, exclude_bonds, external_ext)

            #num_possible_atoms = len(enumeration_mol.atoms) - len(atom_indices)
            #num_possible_bonds = len(enumeration_mol.bonds) - len(bond_indices)
            yield new_atoms, Subgraph(atom_indices, bond_indices), num_possible_atoms, num_possible_bonds
        return

    # Both internal bonds and external edges
    internal_powerset = list(powerset(internal_bonds))
    external_powerset = powerset(external_edges)

    exclude_bonds = set(chain(visited_bond_indices, (edge.bond_index for edge in external_edges)))

    for external_ext in external_powerset:
        if not external_ext:
            # No external extensions. Must have at least one internal bond.
            for internal_bond in internal_powerset[1:]:
                bond_indices = set(subgraph.bond_indices)
                bond_indices.update(internal_bond)
                yield None, Subgraph(subgraph.atom_indices, bond_indices), 0, 0
        else:
            new_atoms = frozenset(ext.end_atom_index for ext in external_ext)
            atom_indices = frozenset(chain(subgraph.atom_indices, new_atoms))
            #            no_go_bond_indices = set(chain(visited_bond_indices, extern

            bond_indices = frozenset(chain(subgraph.bond_indices,
                                            (ext.bond_index for ext in external_ext)))
            num_possible_atoms, num_possible_bonds = find_extension_size(
                enumeration_mol, atom_indices, exclude_bonds, external_ext)
            #num_possible_atoms = len(enumeration_mol.atoms) - len(atom_indices)
            for internal_bond in internal_powerset:
                bond_indices2 = frozenset(chain(bond_indices, internal_bond))
                #num_possible_bonds = len(enumeration_mol.bonds) - len(bond_indices2)
                yield new_atoms, Subgraph(atom_indices, bond_indices2), num_possible_atoms, num_possible_bonds


def find_extension_size(enumeration_mol, known_atoms, exclude_bonds, directed_edges):
    num_remaining_atoms = num_remaining_bonds = 0
    visited_atoms = set(known_atoms)
    visited_bonds = set(exclude_bonds)
    #print "start atoms", visited_atoms
    #print "start bonds", visited_bonds
    #print "Along", [directed_edge.bond_index for directed_edge in directed_edges]
    for directed_edge in directed_edges:
        #print "Take", directed_edge
        stack = [directed_edge.end_atom_index]

        # simple depth-first search search
        while stack:
            atom_index = stack.pop()
            for next_edge in enumeration_mol.directed_edges[atom_index]:
                #print "Visit", next_edge.bond_index, next_edge.end_atom_index
                bond_index = next_edge.bond_index
                if bond_index in visited_bonds:
                    #print "Seen bond", bond_index
                    continue
                num_remaining_bonds += 1
                visited_bonds.add(bond_index)
                #print "New BOND!", bond_index, "count", num_remaining_bonds

                next_atom_index = next_edge.end_atom_index
                if next_atom_index in visited_atoms:
                    #print "Seen atom"
                    continue
                num_remaining_atoms += 1
                #print "New atom!", next_atom_index, "count", num_remaining_atoms
                visited_atoms.add(next_atom_index)

                stack.append(next_atom_index)
                
    #print "==>", num_remaining_atoms, num_remaining_bonds
    return num_remaining_atoms, num_remaining_bonds

# Check if a SMARTS is in all targets.
# Uses a dictionary-style API, but please only use matcher[smarts]
# Caches all previous results.

class CachingTargetsMatcher(dict):
    def __init__(self, targets, required_match_count=None):
        self.targets = targets
        if required_match_count is None:
            required_match_count = len(targets)
        self.required_match_count = required_match_count
        self._num_allowed_errors = len(targets) - required_match_count
        super(dict, self).__init__()

    def shift_targets(self):
        assert self._num_allowed_errors >= 0, (self.required_match_count, self._num_allowed_errors)
        self.targets = self.targets[1:]
        self._num_allowed_errors = len(self.targets) - self.required_match_count
        
    def __missing__(self, smarts):
        num_allowed_errors = self._num_allowed_errors
        if num_allowed_errors < 0:
            raise AssertionError("I should never be called")
            self[smarts] = False
            return False
        
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            raise AssertionError("Bad SMARTS: %r" % (smarts,))

        num_allowed_errors = self._num_allowed_errors
        for target in self.targets:
            if not MATCH(target, pat):
                if num_allowed_errors == 0:
                    # Does not match. No need to continue processing
                    self[smarts] = False
                    return False
                num_allowed_errors -= 1
        # Matches enough structures, which means it will always
        # match enough structures. (Even after shifting.)
        self[smarts] = True
        return True

class VerboseCachingTargetsMatcher(object):
    def __init__(self, targets, required_match_count=None):
        self.targets = targets
        if required_match_count is None:
            required_match_count = len(targets)
        self.cache = {}
        self.required_match_count = required_match_count
        self._num_allowed_errors = len(targets) - required_match_count
        self.num_lookups = self.num_cached_true = self.num_cached_false = 0
        self.num_search_true = self.num_search_false = self.num_matches = 0

    def shift_targets(self):
        assert self._num_allowed_errors >= 0, (self.required_match_count, self._num_allowed_errors)
        if self._num_allowed_errors > 1:
            self.targets = self.targets[1:]
            self._num_allowed_errors = len(self.targets) - self.required_match_count
            
    def __getitem__(self, smarts, missing=object()):
        self.num_lookups += 1
        x = self.cache.get(smarts, missing)
        if x is not missing:
            if x:
                self.num_cached_true += 1
            else:
                self.num_cached_false += 1
            return x
        
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            raise AssertionError("Bad SMARTS: %r" % (smarts,))

        for i, target in enumerate(self.targets):
            if not MATCH(target, pat):
                # Does not match. No need to continue processing
                self.num_search_false += 1
                self.num_matches += i+1
                self.cache[smarts] = False
                N = len(self.targets)
                return False
                # TODO: should I move the mismatch structure forward
                # so that it's tested earlier next time?
        # Matches everything
        self.num_matches += i+1
        self.num_search_true += 1
        self.cache[smarts] = True
        return True
        
    def report(self):
        print >>sys.stderr, "%d tests of %d unique SMARTS, cache: %d True %d False, search: %d True %d False (%d substructure tests)" % (self.num_lookups, len(self.cache), self.num_cached_true, self.num_cached_false, self.num_search_true, self.num_search_false, self.num_matches)
        

##### Different maximization algorithms ######
def prune_maximize_bonds(subgraph, mol, num_remaining_atoms, num_remaining_bonds, best_sizes):
    # Quick check if this is a viable search direction
    num_atoms = len(subgraph.atom_indices)
    num_bonds = len(subgraph.bond_indices)
    best_num_atoms, best_num_bonds = best_sizes

    # Prune subgraphs which are too small can never become big enough
    diff_bonds = (num_bonds + num_remaining_bonds) - best_num_bonds
    if diff_bonds < 0:
        return True
    elif diff_bonds == 0:
        # Then we also maximize the number of atoms
        diff_atoms = (num_atoms + num_remaining_atoms) - best_num_atoms
        if diff_atoms <= 0:
            return True

    return False
    

def prune_maximize_atoms(subgraph, mol, num_remaining_atoms, num_remaining_bonds, best_sizes):
    # Quick check if this is a viable search direction
    num_atoms = len(subgraph.atom_indices)
    num_bonds = len(subgraph.bond_indices)
    best_num_atoms, best_num_bonds = best_sizes

    # Prune subgraphs which are too small can never become big enough
    diff_atoms = (num_atoms + num_remaining_atoms) - best_num_atoms
    if diff_atoms < 0:
        return True
    elif diff_atoms == 0:
        diff_bonds = (num_bonds + num_remaining_bonds) - best_num_bonds
        if diff_bonds <= 0:
            return True
    else:
        #print "Could still have", diff_atoms
        #print num_atoms, num_remaining_atoms, best_num_atoms
        pass

    return False

##### Callback handlers for storing the "best" information #####x

class _SingleBest(object):
    def __init__(self, timer, verbose):
        self.best_num_atoms = self.best_num_bonds = -1
        self.best_smarts = None
        self.sizes = (-1, -1)
        self.timer = timer
        self.verbose = verbose

    def _new_best(self, num_atoms, num_bonds, smarts):
        self.best_num_atoms = num_atoms
        self.best_num_bonds = num_bonds
        self.best_smarts = smarts
        self.sizes = sizes = (num_atoms, num_bonds)
        self.timer.mark("new best")
        if self.verbose:
            dt = self.timer.mark_times["new best"] - self.timer.mark_times["start fmcs"]
            sys.stderr.write("Best after %.1fs: %d atoms %d bonds %s\n" % (dt, num_atoms, num_bonds, smarts))
        return sizes

    def get_result(self, completed, matcher):
        return MCSResult(self.best_num_atoms, self.best_num_bonds, self.best_smarts, completed, matcher)

class MCSResult(object):
    def __init__(self, num_atoms, num_bonds, smarts, completed, matches):
        self.num_atoms = num_atoms
        self.num_bonds = num_bonds
        self.smarts = smarts
        self.completed = completed
        self.matches = matches
    def __nonzero__(self):
        return self.smarts is not None
        

class SingleBestAtoms(_SingleBest):
    def add_new_match(self, subgraph, mol, smarts):
        sizes = self.sizes
        
        # See if the subgraph match is better than the previous best
        num_subgraph_atoms = len(subgraph.atom_indices)
        if num_subgraph_atoms < sizes[0]:
            return sizes

        num_subgraph_bonds = len(subgraph.bond_indices)
        if num_subgraph_atoms == sizes[0]:
            if num_subgraph_bonds <= sizes[1]:
                return sizes

        return self._new_best(num_subgraph_atoms, num_subgraph_bonds, smarts)

class SingleBestBonds(_SingleBest):
    def add_new_match(self, subgraph, mol, smarts):
        sizes = self.sizes
        
        # See if the subgraph match is better than the previous best
        num_subgraph_bonds = len(subgraph.bond_indices)
        if num_subgraph_bonds < sizes[1]:
            return sizes

        num_subgraph_atoms = len(subgraph.atom_indices)
        if num_subgraph_bonds == sizes[1] and num_subgraph_atoms <= sizes[0]:
            return sizes
        return self._new_best(num_subgraph_atoms, num_subgraph_bonds, smarts)



### Check if there are any ring atoms; used in --complete-rings-only

# This is (yet) another depth-first graph search algorithm

def check_complete_rings_only(smarts, subgraph, enumeration_mol):
    #print "check", smarts, len(subgraph.atom_indices), len(subgraph.bond_indices)

    atoms = enumeration_mol.atoms
    bonds = enumeration_mol.bonds

    # First, are any of bonds in the subgraph ring bonds in the original structure?
    ring_bonds = []
    for bond_index in subgraph.bond_indices:
        bond = bonds[bond_index]
        if bond.is_in_ring:
            ring_bonds.append(bond_index)

    #print len(ring_bonds), "ring bonds"
    if not ring_bonds:
        # No need to check .. this is an acceptable structure
        return True

    if len(ring_bonds) <= 2:
        # No need to check .. there are no rings of size 2
        return False

    # Otherwise there's more work. Need to ensure that
    # all ring atoms are still in a ring in the subgraph.

    confirmed_ring_bonds = set()
    subgraph_ring_bond_indices = set(ring_bonds)
    for bond_index in ring_bonds:
        #print "start with", bond_index, "in?", bond_index in confirmed_ring_bonds
        if bond_index in confirmed_ring_bonds:
            continue
        # Start a new search, starting from this bond
        from_atom_index, to_atom_index = bonds[bond_index].atom_indices

        # Map from atom index to depth in the bond stack
        atom_depth = {from_atom_index: 0,
                      to_atom_index: 1}
        bond_stack = [bond_index]
        backtrack_stack = []
        prev_bond_index = bond_index
        current_atom_index = to_atom_index

        while 1:
            # Dive downwards, ever downwards
            next_bond_index = next_atom_index = None
            this_is_a_ring = False
            for outgoing_edge in enumeration_mol.directed_edges[current_atom_index]:
                if outgoing_edge.bond_index == prev_bond_index:
                    # Don't loop back
                    continue
                if outgoing_edge.bond_index not in subgraph_ring_bond_indices:
                    # Only advance along ring edges which are in the subgraph
                    continue

                if outgoing_edge.end_atom_index in atom_depth:
                    #print "We have a ring"
                    # It's a ring! Mark everything as being in a ring
                    confirmed_ring_bonds.update(bond_stack[atom_depth[outgoing_edge.end_atom_index]:])
                    confirmed_ring_bonds.add(outgoing_edge.bond_index)
                    if len(confirmed_ring_bonds) == len(ring_bonds):
                        #print "Success!"
                        return True
                    this_is_a_ring = True
                    continue

                # New atom. Need to explore it.
                #print "we have a new bond", outgoing_edge.bond_index, "to atom", outgoing_edge.end_atom_index
                if next_bond_index is None:
                    # This will be the immediate next bond to search in the DFS
                    next_bond_index = outgoing_edge.bond_index
                    next_atom_index = outgoing_edge.end_atom_index
                else:
                    # Otherwise, backtrack and examine the other bonds
                    backtrack_stack.append(
                        (len(bond_stack), outgoing_edge.bond_index, outgoing_edge.end_atom_index) )

            if next_bond_index is None:
                # Could not find a path to take. Might be because we looped back.
                if this_is_a_ring:
                    #assert prev_bond_index in confirmed_ring_bonds, (prev_bond_index, confirmed_ring_bonds)
                    # We did! That means we can backtrack
                    while backtrack_stack:
                        old_size, prev_bond_index, current_atom_index = backtrack_stack.pop()
                        if bond_index not in confirmed_ring_bonds:
                            # Need to explore this path.
                            # Back up and start the search from here
                            del bond_stack[old_size:]
                            break
                    else:
                        # No more backtracking. We fail. Try next bond?
                        # (If it had been sucessful then the
                        #    len(confirmed_ring_bonds) == len(ring_bonds)
                        # would have return True)
                        break
                else:
                    # Didn't find a ring, nowhere to advance
                    return False
            else:
                # Continue deeper
                bond_stack.append(next_bond_index)
                atom_depth[next_atom_index] = len(bond_stack)
                prev_bond_index = next_bond_index
                current_atom_index = next_atom_index

        # If we reached here then try the next bond
        #print "Try again"


class SingleBestAtomsCompleteRingsOnly(_SingleBest):
    def add_new_match(self, subgraph, mol, smarts):
        sizes = self.sizes
        
        # See if the subgraph match is better than the previous best
        num_subgraph_atoms = len(subgraph.atom_indices)
        if num_subgraph_atoms < sizes[0]:
            return sizes

        num_subgraph_bonds = len(subgraph.bond_indices)
        if num_subgraph_atoms == sizes[0] and num_subgraph_bonds <= sizes[1]:
            return sizes

        if check_complete_rings_only(smarts, subgraph, mol):
            return self._new_best(num_subgraph_atoms, num_subgraph_bonds, smarts)
        return sizes
        

class SingleBestBondsCompleteRingsOnly(_SingleBest):
    def add_new_match(self, subgraph, mol, smarts):
        sizes = self.sizes
        
        # See if the subgraph match is better than the previous best
        num_subgraph_bonds = len(subgraph.bond_indices)
        if num_subgraph_bonds < sizes[1]:
            return sizes
        
        num_subgraph_atoms = len(subgraph.atom_indices)
        if num_subgraph_bonds == sizes[1] and num_subgraph_atoms <= sizes[0]:
            return sizes

        if check_complete_rings_only(smarts, subgraph, mol):
            return self._new_best(num_subgraph_atoms, num_subgraph_bonds, smarts)
        return sizes

def no_pruning(subgraph, mol, num_remaining_atoms, num_remaining_bonds, best_sizes):
    return False

_maximize_options = {
    ("atoms", False): (prune_maximize_atoms, SingleBestAtoms),
    ("atoms", True): (prune_maximize_atoms, SingleBestAtomsCompleteRingsOnly),
    ("bonds", False): (prune_maximize_bonds, SingleBestBonds),
    ("bonds", True): (prune_maximize_bonds, SingleBestBondsCompleteRingsOnly),
    ("all-matches", False): (no_pruning, SingleBestAtoms),
    # XXX need to implement more here
    }


###### The engine of the entire system. Enumerate subgraphs and see if they match. #####

def enumerate_subgraphs(enumeration_mols, prune, atom_assignment, matches_all_targets, hits, timeout,
                        heappush, heappop):
    if timeout is None:
        end_time = None
    else:
        end_time = time.time() + timeout

    seeds = []
        
    best_sizes = (0, 0)
    # Do a quick check for the not uncommon case where one of the input fragments
    # is the largest substructure or one off from the largest.
    for mol in enumeration_mols:
        atom_range = range(len(mol.atoms))
        bond_set = set(range(len(mol.bonds)))
        subgraph = Subgraph(atom_range, bond_set)
        if not prune(subgraph, mol, 0, 0, best_sizes):
            # Micro-optimization: the largest fragment SMARTS doesn't
            # need to be canonicalized because there will only ever be
            # one match. It's also unlikely that the other largest
            # fragments need canonicalization.
            smarts = make_arbitrary_smarts(subgraph, mol, atom_assignment)
            if matches_all_targets[smarts]:
                best_sizes = hits.add_new_match(subgraph, mol, smarts)
            
    
    for mol in enumeration_mols:
        directed_edges = mol.directed_edges
        # Using 20001 random ChEMBL pairs, timeout=15.0 seconds
        #  1202.6s with original order
        #  1051.9s sorting by (bond.is_in_ring, bond_index)
        #  1009.7s sorting by (bond.is_in_ring + atom1.is_in_ring + atom2.is_in_ring)
        #  1055.2s sorting by (if bond.is_in_ring: 2; else: -(atom1.is_in_ring + atom2.is_in_ring))
        #  1037.4s sorting by (atom1.is_in_ring + atom2.is_in_ring)
        sorted_bonds = list(enumerate(mol.bonds))
        def get_bond_ring_score((bond_index, bond), atoms=mol.atoms):
            a1, a2 = bond.atom_indices
            return bond.is_in_ring + atoms[a1].is_in_ring + atoms[a2].is_in_ring
        sorted_bonds.sort(key = get_bond_ring_score)

        visited_bond_indices = set()
        num_remaining_atoms = len(mol.atoms)-2
        num_remaining_bonds = len(mol.bonds)
        for bond_index, bond in sorted_bonds: #enumerate(mol.bonds): #
            #print "bond_index", bond_index, len(mol.bonds)
            visited_bond_indices.add(bond_index)
            num_remaining_bonds -= 1
            subgraph = Subgraph(bond.atom_indices, frozenset([bond_index]))

            # I lie about the remaining atom/bond sizes here.
            if prune(subgraph, mol, num_remaining_atoms, num_remaining_bonds, best_sizes):
                continue
            # bond.canonical_bondtype doesn't necessarily give the same
            # SMARTS as make_canonical_smarts, but that doesn't matter.
            # 1) I know it's canonical, 2) it's faster, and 3) there is
            # no place else which generates single-bond canonical SMARTS.
            #smarts = make_canonical_smarts(subgraph, mol, atom_assignment)
            smarts = bond.canonical_bondtype
            if matches_all_targets[smarts]:
                best_sizes = hits.add_new_match(subgraph, mol, smarts)
            else:
                # This can happen if there's a threshold
                #raise AssertionError("This should never happen: %r" % (smarts,))
                continue

            a1, a2 = bond.atom_indices
            outgoing_edges = [e for e in (directed_edges[a1] + directed_edges[a2])
                  if e.end_atom_index not in bond.atom_indices and e.bond_index not in visited_bond_indices]

            empty_internal = frozenset()
            if not outgoing_edges:
                pass
            else:
                # The priority is the number of bonds in the subgraph, ordered so
                # that the subgraph with the most bonds comes first. Since heapq
                # puts the smallest value first, I reverse the number. The initial
                # subgraphs have 1 bond, so the initial score is -1.
                heappush(seeds, (-1, tiebreaker(), subgraph,
                                 visited_bond_indices.copy(), empty_internal, outgoing_edges,
                                 mol, directed_edges))
    
    # I made so many subtle mistakes where I used 'subgraph' instead
    # of 'new_subgraph' in the following section that I finally
    # decided to get rid of 'subgraph' and use 'old_subgraph' instead.
    del subgraph

    while seeds:
        if end_time:
            if time.time() >= end_time:
                return False
            
        #print "There are", len(seeds), "seeds", seeds[0][:2]
        score, _, old_subgraph, visited_bond_indices, internal_bonds, external_edges, mol, directed_edges = heappop(seeds)

        new_visited_bond_indices = visited_bond_indices.copy()
        new_visited_bond_indices.update(internal_bonds)
        ## for edge in external_edges:
        ##     assert edge.bond_index not in new_visited_bond_indices
        new_visited_bond_indices.update(edge.bond_index for edge in external_edges)

        for new_atoms, new_subgraph, num_remaining_atoms, num_remaining_bonds in \
               all_subgraph_extensions(mol, old_subgraph, visited_bond_indices, internal_bonds, external_edges):
            if prune(new_subgraph, mol, num_remaining_atoms, num_remaining_bonds, best_sizes):
                #print "PRUNE", make_canonical_smarts(new_subgraph, mol, atom_assignment)
                continue
            smarts = make_canonical_smarts(new_subgraph, mol, atom_assignment)
            if matches_all_targets[smarts]:
                #print "YES", smarts
                best_sizes = hits.add_new_match(new_subgraph, mol, smarts)
            else:
                #print "NO", smarts
                continue

            if not new_atoms:
                continue

            new_internal_bonds, new_external_edges = find_extensions(
                new_atoms, new_visited_bond_indices, directed_edges)

            if new_internal_bonds or new_external_edges:
                # Rank so the subgraph with the highest number of bonds comes first
                heappush(seeds, (-len(new_subgraph.bond_indices), tiebreaker(), new_subgraph,
                                 new_visited_bond_indices, new_internal_bonds, new_external_edges,
                                 mol, directed_edges))

    return True


# Assign a unique identifier to every unique key
class Uniquer(dict):
    def __init__(self):
        self.counter = itertools.count().next
    def __missing__(self, key):
        self[key] = count = self.counter()
        return count


# This is here only so I can see it in the profile statistics
def MATCH(mol, pat):
    return mol.HasSubstructMatch(pat)

class VerboseHeapOps(object):
    def __init__(self, trigger, verbose_delay):
        self.num_seeds_added = 0
        self.num_seeds_processed = 0
        self.verbose_delay = verbose_delay
        self._time_for_next_report = time.time() + verbose_delay
        self.trigger = trigger
        
    def heappush(self, seeds, item):
        self.num_seeds_added += 1
        return heappush(seeds, item)
    
    def heappop(self, seeds):
        if time.time() >= self._time_for_next_report:
            self.trigger()
            self.report()
            self._time_for_next_report = time.time() + self.verbose_delay
        self.num_seeds_processed += 1
        return heappop(seeds)

    def trigger_report(self):
        self.trigger()
        self.report()

    def report(self):
        print >>sys.stderr, "  %d subgraphs enumerated, %d processed" % (
            self.num_seeds_added, self.num_seeds_processed)

def compute_mcs(fragmented_mols, typed_mols, min_num_atoms, threshold_count=None, maximize = Default.maximize,
                complete_rings_only = Default.complete_rings_only,
                timeout = Default.timeout,
                timer = None, verbose=False, verbose_delay=1.0):
    assert timer is not None
    assert 0 < threshold_count <= len(fragmented_mols), threshold_count
    assert len(fragmented_mols) == len(typed_mols)
    assert len(fragmented_mols) >= 2
    if threshold_count is None:
        threshold_count = len(fragmented_mols)
    else:
        assert threshold_count >= 2, threshold_count
    
    atom_assignment = Uniquer()
    if verbose:
        if verbose_delay < 0.0:
            raise ValueError("verbose_delay may not be negative")
        matches_all_targets = VerboseCachingTargetsMatcher(typed_mols[1:], threshold_count-1)
        heapops = VerboseHeapOps(matches_all_targets.report, verbose_delay)
        push = heapops.heappush
        pop = heapops.heappop
        end_verbose = heapops.trigger_report
    else:
        matches_all_targets = CachingTargetsMatcher(typed_mols[1:], threshold_count-1)
        push = heappush
        pop = heappop
        end_verbose = lambda: 1

    try:
        prune, hits_class = _maximize_options[(maximize, bool(complete_rings_only))]
    except KeyError:
        raise ValueError("Unknown 'maximize' option %r" % (maximize,))

    hits = hits_class(timer, verbose)

    remaining_time = None
    if timeout is not None:
        stop_time = time.time() + timeout
    
    for query_index, fragmented_query_mol in enumerate(fragmented_mols):
        enumerated_query_fragments = fragmented_mol_to_enumeration_mols(fragmented_query_mol, min_num_atoms)
        
        targets = typed_mols
        if timeout is not None:
            remaining_time = stop_time - time.time()
        success = enumerate_subgraphs(enumerated_query_fragments, prune, atom_assignment, matches_all_targets, hits,
                                      remaining_time, push, pop)
        if query_index + threshold_count >= len(fragmented_mols):
            break
        if not success:
            break
        matches_all_targets.shift_targets()
        
    end_verbose()

    result = hits.get_result(success, matches_all_targets)
    if result.num_atoms < min_num_atoms:
        return MCSResult(-1, -1, None, result.completed, dict(matches_all_targets))
    return result
        
########## Main driver for the MCS code

class Timer(object):
    def __init__(self):
        self.mark_times = {}
    def mark(self, name):
        self.mark_times[name] = time.time()

def _update_times(timer, times):
    if times is None:
        return
    for (dest, start, end) in ( ("fragment", "start fmcs", "end fragment"),
                                ("select", "end fragment", "end select"),
                                ("enumerate", "end select", "end fmcs"),
                                ("best_found", "start fmcs", "new best"),
                                ("mcs", "start fmcs", "end fmcs") ):
        try:
            diff = timer.mark_times[end] - timer.mark_times[start]
        except KeyError:
            diff = None
        times[dest] = diff

def _get_threshold_count(num_mols, threshold):
    if threshold is None:
        return num_mols

    x = num_mols * threshold
    threshold_count = int(x)
    if threshold_count < x:
        threshold_count += 1
    
    if threshold_count < 2:
        # You can specify 0.00001 or -2.3 but you'll still get
        # at least one *common* substructure.
        threshold_count = 2

    return threshold_count


def fmcs(mols, min_num_atoms=2,
         maximize = Default.maximize,
         atom_compare = Default.atom_compare,
         bond_compare = Default.bond_compare,
         threshold = 1.0,
         match_valences = Default.match_valences,
         ring_matches_ring_only = False,
         complete_rings_only = False,
         timeout=Default.timeout,
         times=None,
         verbose=False,
         verbose_delay=1.0,
         ):

    timer = Timer()
    timer.mark("start fmcs")

    if min_num_atoms < 2:
        raise ValueError("min_num_atoms must be at least 2")
    if timeout is not None:
        if timeout <= 0.0:
            raise ValueError("timeout must be None or a positive value")

    threshold_count = _get_threshold_count(len(mols), threshold)
    if threshold_count > len(mols):
        # Threshold is too high. No possible matches.
        return MCSResult(-1, -1, None, 1, {})
        
    if complete_rings_only:
        ring_matches_ring_only = True

    try:
        atom_typer = atom_typers[atom_compare]
    except KeyError:
        raise ValueError("Unknown atom_compare option %r" % (atom_compare,))
    try:
        bond_typer = bond_typers[bond_compare]
    except KeyError:
        raise ValueError("Unknown bond_compare option %r" % (bond_compare,))


    # Make copies of all of the molecules so I can edit without worrying about the original
    typed_mols = convert_input_to_typed_molecules(mols, atom_typer, bond_typer,
                                                  match_valences = match_valences,
                                                  ring_matches_ring_only = ring_matches_ring_only)
    bondtype_counts = get_canonical_bondtype_counts(typed_mols)
    supported_bondtypes = set()
    for bondtype, count_list in bondtype_counts.items():
        if len(count_list) >= threshold_count:
            supported_bondtypes.add(bondtype)
            # For better filtering, find the largest count which is in threshold
            # Keep track of the counts while building the subgraph.
            # The subgraph can never have more types of a given count.

    
    fragmented_mols = [remove_unknown_bondtypes(typed_mol, bondtype_counts) for typed_mol in typed_mols]
    timer.mark("end fragment")

    sizes = []
    max_num_atoms = fragmented_mols[0].rdmol.GetNumAtoms()
    max_num_bonds = fragmented_mols[0].rdmol.GetNumBonds()
    ignored_count = 0
    for tiebreaker, (typed_mol, fragmented_mol) in enumerate(zip(typed_mols, fragmented_mols)):
        num_atoms, num_bonds = find_upper_fragment_size_limits(fragmented_mol.rdmol,
                                                               fragmented_mol.rdmol_atoms)
        if num_atoms < min_num_atoms:
            # This isn't big enough to be in the MCS
            ignored_count += 1
            if ignored_count + threshold_count > len(mols):
                # I might be able to exit because enough of the molecules don't have
                # a large enough fragment to be part of the MCS
                timer.mark("end select")
                timer.mark("end fmcs")
                _update_times(timer, times)
                return MCSResult(-1, -1, None, True, {})
        else:
            if num_atoms < max_num_atoms:
                max_num_atoms = num_atoms
            if num_bonds < max_num_bonds:
                max_num_bonds = num_bonds
            sizes.append( (num_bonds, num_atoms, tiebreaker, typed_mol, fragmented_mol) )

    if len(sizes) < threshold_count:
        timer.mark("end select")
        timer.mark("end fmcs")
        _update_times(timer, times)
        return MCSResult(-1, -1, None, True, {})

    assert min(size[1] for size in sizes) >= min_num_atoms

    # Sort so the molecule with the smallest largest fragment (by bonds) comes first.
    # Break ties with the smallest number of atoms.
    # Break secondary ties by position.
    sizes.sort()
    #print "Using", Chem.MolToSmiles(sizes[0][4].rdmol)

    timer.mark("end select")

    # Extract the (typed mol, fragmented mol) pairs.
    fragmented_mols = [size_info[4] for size_info in sizes]  # used as queries
    typed_mols = [size_info[3].rdmol for size_info in sizes]    # used as targets

    timer.mark("start enumeration")
    mcs_result = compute_mcs(fragmented_mols, typed_mols, min_num_atoms,
                             threshold_count=threshold_count, maximize=maximize,
                             complete_rings_only=complete_rings_only, timeout=timeout,
                             timer=timer, verbose=verbose, verbose_delay=verbose_delay)
    timer.mark("end fmcs")
    _update_times(timer, times)
    return mcs_result


######### Helper functions to generate structure/fragment output given an MCS match

# Given a Subgraph (with atom and bond indices) describing a
# fragment, make a new molecule object with only that fragment

def subgraph_to_fragment(mol, subgraph):
    emol = Chem.EditableMol(Chem.Mol())
    atom_map = {}
    for atom_index in subgraph.atom_indices:
        emol.AddAtom(mol.GetAtomWithIdx(atom_index))
        atom_map[atom_index] = len(atom_map)

    for bond_index in subgraph.bond_indices:
        bond = mol.GetBondWithIdx(bond_index)
        emol.AddBond(atom_map[bond.GetBeginAtomIdx()],
                     atom_map[bond.GetEndAtomIdx()],
                     bond.GetBondType())

    return emol.GetMol()
    

# Convert a subgraph into a SMILES
def make_fragment_smiles(mcs, mol, subgraph, args=None):
    fragment = subgraph_to_fragment(mol, subgraph)
    new_smiles = Chem.MolToSmiles(fragment)
    return "%s %s\n" % (new_smiles, mol.GetProp("_Name"))


def _copy_sd_tags(mol, fragment):
    fragment.SetProp("_Name", mol.GetProp("_Name"))
    # Copy the existing names over
    for name in mol.GetPropNames():
        if name.startswith("_"):
            continue
        fragment.SetProp(name, mol.GetProp(name))

def _MolToSDBlock(mol):
    # Huh?! There's no way to get the entire SD record?
    mol_block = Chem.MolToMolBlock(mol, kekulize=False)
    tag_data = []
    for name in mol.GetPropNames():
        if name.startswith("_"):
            continue
        value = mol.GetProp(name)
        tag_data.append("> <" + name + ">\n")
        tag_data.append(value + "\n")
        tag_data.append("\n")
    tag_data.append("$$$$\n")
    return mol_block + "".join(tag_data)

def _save_other_tags(mol, fragment, mcs, orig_mol, subgraph, args):
    if args.save_counts_tag is not None:
        if not mcs:
            line = "-1 -1 -1"
        elif mcs.num_atoms == 0:
            line = "0 0 0"
        else:
            line = "1 %d %d" % (mcs.num_atoms, mcs.num_bonds)
        mol.SetProp(args.save_counts_tag, line)

    if args.save_smiles_tag is not None:
        if mcs and mcs.num_atoms > 0:
            smiles = Chem.MolToSmiles(fragment)
        else:
            smiles = "-"
        mol.SetProp(args.save_smiles_tag, smiles)

    if args.save_smarts_tag is not None:
        if mcs and mcs.num_atoms > 0:
            smarts = mcs.smarts
        else:
            smarts = "-"
        mol.SetProp(args.save_smarts_tag, smarts)
    

# Convert a subgraph into an SD file
def make_fragment_sdf(mcs, mol, subgraph, args):
    fragment = subgraph_to_fragment(mol, subgraph)
    Chem.FastFindRings(fragment)
    _copy_sd_tags(mol, fragment)

    if args.save_atom_class_tag is not None:
        output_tag = args.save_atom_class_tag
        atom_classes = get_selected_atom_classes(mol, subgraph.atom_indices)
        if atom_classes is not None:
            fragment.SetProp(output_tag, " ".join(map(str, atom_classes)))

    _save_other_tags(fragment, fragment, mcs, mol, subgraph, args)

    return _MolToSDBlock(fragment)

# 
def make_complete_sdf(mcs, mol, subgraph, args):
    fragment = copy.copy(mol)
    _copy_sd_tags(mol, fragment)

    if args.save_atom_indices_tag is not None:
        output_tag = args.save_atom_indices_tag
        s = " ".join(str(index) for index in subgraph.atom_indices)
        fragment.SetProp(output_tag, s)

    _save_other_tags(fragment, subgraph_to_fragment(mol, subgraph), mcs, mol, subgraph, args)

    return _MolToSDBlock(fragment)


structure_format_functions = {
    "fragment-smiles": make_fragment_smiles,
    "fragment-sdf": make_fragment_sdf,
    "complete-sdf": make_complete_sdf,
    }
def make_structure_format(format_name, mcs, mol, subgraph, args):
    try:
        func = structure_format_functions[format_name]
    except KeyError:
        raise ValueError("Unknown format %r" % (format_name,))
    return func(mcs, mol, subgraph, args)



def parse_num_atoms(s):
    num_atoms = int(s)
    if num_atoms < 2:
        raise argparse.ArgumentTypeError("must be at least 2, not %s" % s)
    return num_atoms

def parse_threshold(s):
    try:
        import fractions
    except ImportError:
        threshold = float(s)
        one = 1.0
    else:
        threshold = fractions.Fraction(s)
        one = fractions.Fraction(1)
    if not (0 <= threshold <= one):
        raise argparse.ArgumentTypeError("must be a value between 0.0 and 1.0, not %s" % s)
    return threshold

def parse_timeout(s):
    if s == "none":
        return None
    timeout = float(s)
    if timeout < 0.0:
        raise argparse.ArgumentTypeError("Must be a non-negative value, not %r" % (s,))
    return timeout

class starting_from(object):
    def __init__(self, left):
        self.left = left
    def __contains__(self, value):
        return self.left <= value 

range_pat = re.compile(r"(\d+)-(\d*)")
value_pat = re.compile("(\d+)")
def parse_select(s):
    ranges = []
    start = 0
    while 1:
        m = range_pat.match(s, start)
        if m is not None:
            # Selected from 'left' to (and including) 'right'
            # Convert into xrange fields, starting from 0
            left = int(m.group(1))
            right = m.group(2)
            if not right:
                ranges.append(starting_from(left-1))
            else:
                ranges.append( xrange(left-1, int(right)) )
        else:
            # Selected a single value
            m = value_pat.match(s, start)
            if m is not None:
                val = int(m.group(1))
                ranges.append( xrange(val-1, val) )
            else:
                raise argparse.ArgumentTypeError("Unknown character at position %d of %r" %(
                    start+1, s))
        start = m.end()
        # Check if this is the end of string or a ','
        t = s[start:start+1]
        if not t:
            break
        if t == ",":
            start += 1
            continue
        raise argparse.ArgumentTypeError("Unknown character at position %d of %r" % (
            start+1, s))
    return ranges



parser = argparse.ArgumentParser(description="Find the maximum common substructure of a set of structures",
     epilog = "For more details on these options, see https://bitbucket.org/dalke/fmcs/")
parser.add_argument("filename", nargs=1,
                    help="SDF or SMILES file")

parser.add_argument("--maximize", choices=["atoms", "bonds", "all-matches"],
                    default=Default.maximize,
                    help="Maximize the number of 'atoms' or 'bonds' in the MCS. (Default: %s)" % (Default.maximize,))
parser.add_argument("--min-num-atoms", type=parse_num_atoms, default=2,
                    metavar="INT",
                    help="Minimimum number of atoms in the MCS (Default: 2)")


compare_shortcuts = {
    "topology": ("any", "any"),
    "elements": ("elements", "any"),
    "types": ("elements", "bondtypes"),
    }
class CompareAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        atom_compare_name, bond_compare_name = compare_shortcuts[value]
        namespace.atom_compare = atom_compare_name
        namespace.bond_compare = bond_compare_name


parser.add_argument("--compare", choices = ["topology", "elements", "types"],
                    default=None, action=CompareAction, help=
                    "Use 'topology' as a shorthand for '--atom-compare any --bond-compare any', "
                    "'elements' is '--atom-compare elements --bond-compare any', "
                    "and 'types' is '--atom-compare elements --bond-compare bondtypes' "
                    "(Default: types)")

parser.add_argument("--atom-compare", choices=["any", "elements", "isotopes"],
                    default=None, help=(
                        "Specify the atom comparison method. With 'any', every atom matches every "
                        "other atom. With 'elements', atoms match only if they contain the same element. "
                        "With 'isotopes', atoms match only if they have the same isotope number; element "
                        "information is ignored so [5C] and [5P] are identical. This can be used to "
                        "implement user-defined atom typing. "
                        "(Default: elements)"))

parser.add_argument("--bond-compare", choices=["any", "bondtypes"],
                    default="bondtypes", help=(
                        "Specify the bond comparison method. With 'any', every bond matches every "
                        "other bond. With 'bondtypes', bonds are the same only if their bond types "
                        "are the same. (Default: bondtypes)"))

parser.add_argument("--threshold", default="1.0", type=parse_threshold, help=
                    "Minimum structure match threshold. A value of 1.0 means that the common "
                    "substructure must be in all of the input structures. A value of 0.8 finds "
                    "the largest substructure which is common to at least 80%% of the input "
                    "structures. (Default: 1.0)")

parser.add_argument("--atom-class-tag", metavar="TAG", help=
                    "Use atom class assignments from the field 'TAG'. The tag data must contain a space "
                    "separated list of integers in the range 1-10000, one for each atom. Atoms are "
                    "identical if and only if their corresponding atom classes are the same. Note "
                    "that '003' and '3' are treated as identical values. (Not used by default)")

## parser.add_argument("--match-valences", action="store_true",
##                     help=
##                     "Modify the atom comparison so that two atoms must also have the same total "
##                     "bond order in order to match.")

                        
parser.add_argument("--ring-matches-ring-only", action="store_true",
                    help=
                    "Modify the bond comparison so that ring bonds only match ring bonds and chain "
                    "bonds only match chain bonds. (Ring atoms can still match non-ring atoms.) ")

parser.add_argument("--complete-rings-only", action="store_true",
                    help=
                    "If a bond is a ring bond in the input structures and a bond is in the MCS "
                    "then the bond must also be in a ring in the MCS. Selecting this option also "
                    "enables --ring-matches-ring-only.")

parser.add_argument("--select", type=parse_select, action="store",
                    default="1-",
                    help=
                    "Select a subset of the input records to process. Example: 1-10,13,20,50- "
                    "(Default: '1-', which selects all structures)")

parser.add_argument("--timeout", type=parse_timeout, metavar="SECONDS",
                    default=Default.timeout,
                    help=
                    "Report the best solution after running for at most 'timeout' seconds. "
                    "Use 'none' for no timeout. (Default: %s)" % (Default.timeout_string,))

parser.add_argument("--output", "-o", metavar="FILENAME",
                    help="Write the results to FILENAME (Default: use stdout)")

parser.add_argument("--output-format", choices = ["smarts", "fragment-smiles", "fragment-sdf", "complete-sdf"],
                    default="smarts", help=
                    "'smarts' writes the SMARTS pattern including the atom and bond criteria. "
                    "'fragment-smiles' writes a matching fragment as a SMILES string. "
                    "'fragment-sdf' writes a matching fragment as a SD file; see --save-atom-class for "
                    "details on how atom class information is saved. "
                    "'complete-sdf' writes the entire SD file with the fragment information stored in "
                    "the tag specified by --save-fragment-indices-tag. (Default: smarts)")

parser.add_argument("--output-all", action="store_true",
                    help=
                    "By default the structure output formats only show an MCS for the first input structure. "
                    "If this option is enabled then an MCS for all of the structures are shown.")

parser.add_argument("--save-atom-class-tag", metavar="TAG", help=
                    "If atom classes are specified (via --class-tag) and the output format is 'fragment-sdf' "
                    "then save the substructure atom classes to the tag TAG, in fragment atom order. By "
                    "default this is the value of --atom-class-tag.")

parser.add_argument("--save-counts-tag", metavar="TAG", help=
                    "Save the fragment count, atom count, and bond count to the specified SD tag as "
                    "space separated integers, like '1 9 8'. (The fragment count will not be larger than "
                    "1 until fmcs supports disconnected MCSes.)")

parser.add_argument("--save-atom-indices-tag", metavar="TAG", help=
                    "If atom classes are specified and the output format is 'complete-sdf' "
                    "then save the MCS fragment atom indices to the tag TAG, in MCS order. "
                    "(Default: mcs-atom-indices)")

parser.add_argument("--save-smarts-tag", metavar="TAG", help=
                    "Save the MCS SMARTS to the specified SD tag. Uses '-' if there is no MCS")

parser.add_argument("--save-smiles-tag", metavar="TAG", help=
                    "Save the fragment SMILES to the specified SD tag. Uses '-' if there is no MCS")


parser.add_argument("--times", action="store_true",
                    help="Print timing information to stderr")
parser.add_argument("-v", "--verbose", action="count", dest="verbosity",
                    help="Print progress statistics to stderr. Use twice for higher verbosity.")
parser.add_argument("--version", action="version", version="%(prog)s " + __version__)


# RDKit's match function only returns the atom indices of the match.
# To get the bond indices, I need to go through the pattern molecule.
def _get_match_bond_indices(pat, mol, match_atom_indices):
    bond_indices = []
    for bond in pat.GetBonds():
        mol_atom1 = match_atom_indices[bond.GetBeginAtomIdx()]
        mol_atom2 = match_atom_indices[bond.GetEndAtomIdx()]
        bond = mol.GetBondBetweenAtoms(mol_atom1, mol_atom2)
        assert bond is not None
        bond_indices.append(bond.GetIdx())
    return bond_indices

def main(args=None):
    args = parser.parse_args(args)

    filename = args.filename[0]
    fname = filename.lower()
    if fname.endswith(".smi"):
        try:
            reader = Chem.SmilesMolSupplier(filename, titleLine=False)
        except IOError:
            raise SystemExit("Unable to open SMILES file %r" % (filename,))
    elif fname.endswith(".sdf"):
        try:
            reader = Chem.SDMolSupplier(filename)
        except IOError:
            raise SystemExit("Unable to open SD file %r" % (filename,))
    elif fname.endswith(".gz"):
        raise SystemExit("gzip compressed files not yet supported")
    else:
        raise SystemExit("Only SMILES (.smi) and SDF (.sdf) files are supported")

    if args.min_num_atoms < 2:
        parser.error("--min-num-atoms must be at least 2")

    if args.atom_compare is None:
        if args.atom_class_tag is None:
            args.atom_compare = "elements" # Default atom comparison
        else:
            args.atom_compare = "isotopes" # Assing the atom classes to the isotope fields
    else:
        if args.atom_class_tag is not None:
            parser.error("Cannot specify both --atom-compare and --atom-class-tag fields")

    # RDKit uses special property names starting with "_"
    # It's dangerous to use some of them directly
    for name in ("atom_class_tag", "save_atom_class_tag", "save_counts_tag", "save_atom_indices_tag",
                 "save_smarts_tag", "save_smiles_tag"):
        value = getattr(args, name)
        if value is not None:
            if value.startswith("_"):
                parser.error("--%s value may not start with a '_': %r" % (name.replace("_", "-"), value))

    # Set up some defaults depending on the output format
    atom_class_tag = args.atom_class_tag
    if args.output_format == "fragment-sdf":
        if atom_class_tag is not None:
            if args.save_atom_class_tag is None:
                args.save_atom_class_tag = atom_class_tag

    if args.output_format == "complete-sdf":
        if (args.save_atom_indices_tag is None and
            args.save_counts_tag is None and
            args.save_smiles_tag is None and
            args.save_smarts_tag is None):
            parser.error("Using --output-format complete-sdf is useless without at least one "
                         "of --save-atom-indices-tag, --save-smarts-tag, --save-smiles-tag, "
                         "or --save-counts-tag")

    t1 = time.time()
    structures = []
    if args.verbosity > 1:
        sys.stderr.write("Loading structures from %s ..." % (filename,))
        
    for molno, mol in enumerate(reader):
        if not any(molno in range_ for range_ in args.select):
            continue
        if mol is None:
            print >>sys.stderr, "Skipping unreadable structure #%d" % (molno+1,)
            continue
        if atom_class_tag is not None:
            try:
                assign_isotopes_from_class_tag(mol, atom_class_tag)
            except ValueError, err:
                raise SystemExit("Structure #%d: %s" % (molno+1, err))
        structures.append(mol)
        if args.verbosity > 1:
            if len(structures) % 100 == 0:
                sys.stderr.write("\rLoaded %d structures from %s ..." % (len(structures), filename))
                sys.stderr.flush() # not needed; it's stderr. But I'm cautious.

    if args.verbosity > 1:
        sys.stderr.write("\r")

    times = {"load": time.time()-t1}
    
    if args.verbosity:
        print >>sys.stderr, "Loaded", len(structures), "structures from", filename, "    "

    if len(structures) < 2:
        raise SystemExit("Input file %r must contain at least two structures" % (filename,))

    mcs = fmcs(structures,
               min_num_atoms = args.min_num_atoms,
               maximize = args.maximize,
               atom_compare = args.atom_compare,
               bond_compare = args.bond_compare,
               threshold = args.threshold,
               #match_valences = args.match_valences,
               match_valences = False, # Do I really want to support this?
               ring_matches_ring_only = args.ring_matches_ring_only,
               complete_rings_only = args.complete_rings_only,
               timeout = args.timeout,
               times = times,
               verbose = args.verbosity > 1,
               verbose_delay = 1.0,
        )

    msg_format = "Total time %(total).2f seconds: load %(load).2f fragment %(fragment).2f select %(select).2f enumerate %(enumerate).2f"
    times["total"] = times["mcs"] + times["load"]

    if mcs and mcs.completed:
        msg_format += " (MCS found after %(best_found).2f)"

    del mol

    if args.output:
        outfile = open(args.output, "w")
    else:
        outfile = sys.stdout

    if args.output_format == "smarts":
        if not mcs:
            outfile.write("No MCS found\n")
        else:
            if mcs.completed:
                status = "(complete search)"
            else:
                status = "(timed out)"
            outfile.write("%s %d atoms %d bonds %s\n" % (
                mcs.smarts, mcs.num_atoms, mcs.num_bonds, status))
            outfile.write("== Successful matches ==\n")
            for smarts, in_all in sorted(mcs.matches.iteritems(), key=lambda x: len(x[0])):
                if in_all:
                    outfile.write(smarts+"\n")

    else:
        if mcs.smarts is None:
            # There is no MCS. Use something which can't match.
            pat = Chem.MolFromSmarts("[CN]")
        else:
            # Need to make a structure output
            pat = Chem.MolFromSmarts(mcs.smarts)
        for structure in structures:
            atom_indices = structure.GetSubstructMatch(pat)
            if not atom_indices:
                # The only time that a SMARTS shouldn't match an input
                # structure is if there's a threshold cutoff and this
                # structure didn't make it.
                assert args.threshold < 1, "No indices but should have matched everything!"
                continue
            bond_indices = _get_match_bond_indices(pat, structure, atom_indices)
            subgraph = Subgraph(atom_indices, bond_indices)
            if atom_class_tag:
                restore_isotopes(structure)

            outfile.write(make_structure_format(args.output_format, mcs, structure, subgraph, args))

            if not args.output_all:
                break

    if args.output:
        outfile.close()

    if args.times or args.verbosity:
        print >>sys.stderr, msg_format % times

if __name__ == "__main__":
    main(sys.argv[1:])


## 017_COOcount.py

from rdkit import Chem
import cPickle,gzip
from pylab import *

# CAVE
# run this on uniq.func.combined/func.uniq.dupl.sdf.gz

inF  = sys.argv[1]
if inF.endswith(".sdf.gz") or inF.endswith(".sd.gz"):
    cpds = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(inF)) if x is not None]
else:
    cpds = [x for x in  Chem.SDMolSupplier(inF) if x is not None]


filename_inactives = sys.argv[1]
filename_actives = sys.argv[2]
cpds_actives = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(filename_actives)) if x is not None]
cpds_inactives = [x for x in  Chem.ForwardSDMolSupplier(gzip.open(filename_inactives)) if x is not None]

hERG_TL_list_actives = []
property_list_list_actives = []
cpd_name_list_actives = []
prop_dict_actives_all = {}
for cpd in cpds_actives:
    property_list_actives = []
    property_name_list_actives = []
    prop_name_actives = cpd.GetPropNames()
    cpd_name_list_actives.append(cpd.GetProp("_Name"))
    for property in prop_name_actives:
            if property not in prop_dict_actives_all:
                prop_dict_actives_all[property] = []
                prop_dict_actives_all[property].append(float(cpd.GetProp(property)))
            else:
                prop_dict_actives_all[property].append(float(cpd.GetProp(property)))
                property_name_list_actives.append(property)
#
hERG_TL_list_inactives = []
property_list_list_inactives = []
cpd_name_list_inactives = []
prop_dict_inactives_all = {}
for cpd in cpds_inactives:
    property_list_inactives = []
    property_name_list_inactives = []
    prop_name_inactives = cpd.GetPropNames()
    cpd_name_list_inactives.append(cpd.GetProp("_Name"))
    for property in prop_name_inactives:
            if property not in prop_dict_inactives_all:
                list_element = cpd.GetProp(property)
                prop_dict_inactives_all[property] = []
                prop_dict_inactives_all[property].append(float(cpd.GetProp(property)))
            else:
                prop_dict_inactives_all[property].append(float(cpd.GetProp(property)))
                property_name_list_inactives.append(property)

intersectionlist = ['fr_COO']


counter = 0
data = []
label_list = []
for relevant_descr in intersectionlist:
    data.append([prop_dict_actives_all[relevant_descr],prop_dict_inactives_all[relevant_descr]])
    label_list.append(["actives","inactives"])
    counter += 1
fig = plt.figure(figsize=(10,10))

master_dict_list = [prop_dict_inactives_all,prop_dict_actives_all]

title_list = ["inactives","actives"]

for x,i in enumerate(master_dict_list):
    plotname = "ax"+str(x)
    subplotnumbering = int("12"+str(x+1))
    plotname = fig.add_subplot(subplotnumbering)
    plotname.hist(master_dict_list[x]['fr_COO'])
    plotname.tick_params(axis='both', which='major', labelsize=14)
    plotname.set_title(title_list[x],fontsize=30)
    plotname.set_xlabel("COO count",fontsize=20)
    plotname.set_ylabel("compound count",fontsize=20)
show()
