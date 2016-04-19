from IOhandle.models import Compound,Molecule,Protein,ActivityPoint,Project,Target
from MMPMaker.functions import index_hydrogen_change, make_mol_mmp, make_mmp_database, make_list_mmps, act_mmp_3d, make_ph4_difference_points, find_pharma_changes
from MMPMaker.models import MMPDiffMap, MMPComparison, ActPharmaPoint, MMPFrag
from Group.models import Group
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from django.core.exceptions import ValidationError
import sys, os
from IOhandle.functions import add_new_comp
import csv
import uuid
# Global variable to define what things can be used possibly
__POSS_LIST__ = ["IC50", "Ki"]


def check_smiles_or_SD(file_path):
    """Function to check if a is a valid smiles or SD file
    Takes a file path
    Returns a lits of Molecules or None"""
    if os.path.isfile(file_path):
        pass
    else:
        print "NO SUCH FILE: ", file_path
        return None
    try:
        mols = Chem.SDMolSupplier(file_path)
    except IOError:
        mols = None
    try:
        mols = Chem.SmilesMolSupplier(file_path)
    except IOError:
        mols = None

    if mols == None:
        print "NOT VALID SMILE OR SD FILE: ", file_path
    else:
        return mols


def read_CSV(file_in):
    """Function to read a csv file.
    Takes a file path
    Returns a csv dict object"""
    return csv.DictReader(file_in)


def add_new_mol(rdmol, target):
    """Function to add a new bound Molecule object
    Takes an RDKit molecule and a Target
    Returns None"""
    new_mol = Molecule()
    rdProps = rdmol.GetProp("_Name").split("_")
    comp_ref = add_new_comp(rdmol)
    if comp_ref is None:
        return None
    # To get rid of the .pdb suffix
    pdb_id = rdProps[0].split(".")[0]
    # Check that the name is unique and more than 3 characters long
    # and doesn't contain the target.title
    if len(pdb_id) < 4 or len(Protein.objects.filter(code=pdb_id)) > 0 or target.title in pdb_id:
        # Make a new uniqid
        # First up check that this molecule has not been added before
        mols = [[Chem.MolFromMolBlock(str(x.sdf_info)), x.pk] for x in Molecule.objects.filter(prot_id__target_id=target, cmpd_id=comp_ref)]
        [x[0].SetProp("_Name", "N") for x in mols]
        rdmol.SetProp("_Name", "N")
        # If this is actually a duplicate then continue
        sd_block = Chem.MolToMolBlock(Chem.MolFromMolBlock(Chem.MolToMolBlock(rdmol)))
        matches = [x for x in mols if sd_block == Chem.MolToMolBlock(x[0])]
        if len(matches) > 0:
            # Just return
            return
        else:
            # We have not put this EXACT mol into the database
            # Now lets make an ID
            molid = uuid.uuid4().hex + "_" + pdb_id
            rdmol.SetProp("_Name", molid)
    else:
        molid = pdb_id
    # Make a protein object by which it is related in the DB
    new_mol.prot_id = Protein.objects.get_or_create(code=molid, target_id=target)[0]
    new_mol.sdf_info = Chem.MolToMolBlock(rdmol)
    new_mol.smiles = Chem.MolToSmiles(rdmol, isomericSmiles=True)
    try:
        new_mol.lig_id = rdProps[1]
        new_mol.chain_id = rdProps[2]
        new_mol.occupancy = float(rdProps[3])
    except IndexError:
        new_mol.lig_id = "UNL"
        new_mol.chain_id = "Z"
        new_mol.occupancy = 0.0
    # Add this to the compound list -> make sure this passes in for the
    # correct molecule. I.e. if it fails where does it go???
    # Now link that compound back
    new_mol.cmpd_id = comp_ref
    try:
        new_mol.validate_unique()
        new_mol.save()
    except ValidationError:
        pass


def initialise_dummys():
    """Function to initialise all the dummy database objects
    Takes no args
    Returns None"""
    # Make the dummy objects required for other objects
    d_targ = Target.objects.get_or_create(title="DUMMY", uniprot_id="DUMMY")[0]
    d_prot = Protein.objects.get_or_create(code="DUMMY", target_id=d_targ)[0]
    d_cmpd = Compound.objects.get_or_create(smiles="DUMMY", inchi="DUMMY",
                                            num_added=0, mol_log_p=0.0,
                                            mol_wt=0.0, heavy_atom_count=0,
                                            heavy_atom_mol_wt=0, nhoh_count=0,
                                            no_count=0, num_h_acceptors=0,
                                            num_h_donors=0, num_het_atoms=0,
                                            num_rot_bonds=0, num_val_electrons=0,
                                            ring_count=0, tpsa=0.0)[0]
    Molecule.objects.get_or_create(prot_id=d_prot, cmpd_id=d_cmpd, rscc=0.0,
                                   lig_id="DUMMY", chain_id="X", occupancy=0.0,
                                   x_com=0.0, y_com=0.0, z_com=0.0, sdf_info="DUMMY")
    Group.objects.get_or_create(group_number=1, group_text="DUMMY", method_used="DUMMY")
    ActivityPoint.objects.get_or_create(confidence=0, cmpd_id=d_cmpd, target_id=d_targ,
                                        activity=0.0, units="DUMMY", internal_id="DUMMY",
                                        source="DUMMY")


def do_oommppaa_proc(target_id=1):
    """Function to generate the full MMP database and find differences.
    Make the mpps, index h change, find all matched pairs
    make 3D coords and make and find differences
    Takes a target_id
    Returns None"""
    initialise_dummys()
    # Make them for activitv points
    make_mmp_database(target_id, option="TWOD")
    # Make them for 3D molecules
    make_mmp_database(target_id, option=None)
    # Do index H change
    index_hydrogen_change(target_id)
    # Find the mmps
    out_mmps = make_list_mmps(target_id, option="ACT", max_size=10,
                              ratio=3, use_ratio=False)
    # Now generate the  three-d mmps from this
    act_mmp_3d(out_mmps, target_id)
    # Find differences between molecules
    make_ph4_difference_points(target_id, opt_type=__POSS_LIST__)
    # Record them to the database
    find_pharma_changes(target_id, my_type=str(__POSS_LIST__))


def find_3d_mmps(target_id=1):
    """Function to fragment all Molecules for a target
    Takes a target id
    Returns None"""
    make_mmp_database(target_id, option=None)


def loading_index_h_change(target_id=1):
    """Function to index hydrogen change for a target
    Takes a target id
    Returns None"""
    index_hydrogen_change(target_id)


def delete_target(target_id):
    """Function to delete all the molecules, 
    proteins and MMP data for a given target
    Takes a target id
    Returns None
    """
    # We need to delete them on by one for SQLite
    # Get the mols, prots and mmpfragss
    mols = Molecule.objects.filter(prot_id__target_id=target_id)
    prots = Protein.objects.filter(target_id=target_id)
    mmpfrags = MMPFrag.objects.filter(mol_id__prot_id__target_id=target_id)
    MMPDiffMap.objects.filter(target_id=target_id).delete()
    acts = ActivityPoint.objects.filter(target_id=target_id)
    ap = ActPharmaPoint.objects.filter(target_id=target_id)
    mmpcomps = MMPComparison.objects.filter(target_id=target_id)
    # Now delete all of them
    for a in acts:
        a.delete()
    for a in ap:
        a.delete()
    for m in mols:
        m.delete()
    for p in prots:
        p.delete()
    for m in mmpfrags:
        m.delete()
    for m in mmpcomps:
        m.delete()
    ActPharmaPoint.objects.filter(target_id=target_id).delete()
    # Now clean up by deleting the target
    Target.objects.get(pk=target_id).delete()


def list_targets():
    """Function to list all the targets
    Takes no args
    Returns None"""
    targs = Target.objects.all()
    for t in targs:
        print t.title


def make_3d_confs(target_id=1):
    """Function to make 3D coords for a target
    Takes a target_id
    Returns None"""
    out_mmps = make_list_mmps(target_id, option="ACT", max_size=10, ratio=3, use_ratio=False)
    # Now generate the  three-d mmps from this
    act_mmp_3d(out_mmps, target_id)


def add_new_act(comp_ref, target, activity, units, chembl_id, source, operator):
    """Function to register a new activity point
    Takes a Compound object a Target object, a measured activity, measured units,
    an id for the compound and a source, e.g. IC50
    Returns None"""
    new_act = ActivityPoint()
    try:
        new_act.activity = float(activity)
    except ValueError:
        print "ACTIVITY DATA NOT A NUMBER -> SKIPPING"
        return
    new_act.units = units
    new_act.source = source
    new_act.cmpd_id = comp_ref
    new_act.target_id = target
    new_act.internal_id = chembl_id
    new_act.operator = operator
    try:
        new_act.validate_unique()
        new_act.save()
    except ValidationError:
        pass


def load_mols(target, file_path):
    """Function to load in the 3D molecules
    Takes a Target object and a file path
    Returns None"""
    mols = Chem.SDMolSupplier(file_path)
    if not mols:
        return
    tot = len(mols)
    if tot == 0:
        print "No molecules given"
        return
    old = -1
    print "Adding molecules..."
    for i, m in enumerate(mols):
        # Catch none molecules
        if m is None:
            print "None molecule", sys.exit()
        # Print the progress
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
        Chem.AssignAtomChiralTagsFromStructure(m)
        add_new_mol(m, target)
    old = 100
    sys.stdout.write("\r%d%%" % old)
    sys.stdout.flush()
    print "\nAdding molecules complete"


def load_activity_data(target, file_path):
    """Function to load in a CSV file of activity data
    Takes a Target object and a file path
    Returns None"""
    # Read the file into a CSV dict
    in_d = read_CSV(open(file_path))
    # Fields looking for
    all_fields = ["smiles", "Activity", "ID", "operator"]
    mandatory_fields = ["smiles", "Activity"]
    # Check to see if fields are missing
    missing_fields = [x for x in mandatory_fields if x not in in_d.fieldnames]
    if len(missing_fields) != 0:
        print " ".join(missing_fields), " fields required"
        sys.exit()
    if len([x for x in all_fields if x not in in_d.fieldnames]) != 0:
        print " ".join([x for x in all_fields if x not in in_d.fieldnames]), " fields missing"
    tot = len(open(file_path).readlines()) - 1
    if tot == 0:
        print "No activity data"
        return
    old = -1
    print "Loading activity data"
    for i, l in enumerate(in_d):
        # Do the percent clock
        if i * 100 / tot != old:
            old = i * 100 / tot
            sys.stdout.write("\r%d%% complete..." % old)
            sys.stdout.flush()
        m = Chem.MolFromSmiles(str(l["smiles"]).decode('string-escape'))
        if m is None:
            print "Error None molecule", l["smiles"]
            continue
        comp_ref = add_new_comp(m)
        if comp_ref is None:
            continue
        # Now add the required information if no column is entered
        units = l.get("units")
        if units is None:
            units = "pnM"
        chid = l.get("ID")
        if chid is None:
            chid = "NONE"
        source = l.get("Source")
        if source is None:
            source = "IC50"
        operator = l.get("operator")
        if operator is None:
            operator = "NA"
        add_new_act(comp_ref, target, l["Activity"], units, chid, source, operator)
    old = 100
    sys.stdout.write("\r%d%%" % old)
    sys.stdout.flush()
    print "\nAdding activity data complete"
    return None


def load_protein(target, file_path):
    """Function to load in a pdb file
    Takes a Target object and a file path
    Returns None"""
    new_prot = Protein()
    new_prot.code = target.title + "TEMP"
    new_prot.pdb_info = open(file_path).read()
    new_prot.target_id = target
    try:
        new_prot.validate_unique()
        new_prot.save()
    except ValidationError:
        my_prot = Protein.objects.get(code=target.title + "TEMP")
        print "TEMPLATE PROTEIN ALREADY ENTERED"
        print "OVERWRITING PROTEIN"
        my_prot.pdb_info = open(file_path).read()
        my_prot.save()


def load_compounds(file_path):
    """Function to load compounds and make the MMPs
    Takes a file path
    Returns None"""
    mols = Chem.SDMolSupplier(file_path)
    counter = 0
    for m in mols:
        if m is None:
            print "NONE MOL"
            continue
        counter +=1
        print counter
        # add the new compound to the database
        comp_ref = add_new_comp(m)
        if comp_ref is None:
            continue
        new_m = Chem.MolFromSmiles(str(comp_ref.smiles))
        # Filter too big molecules
        if Descriptors.ExactMolWt(new_m) > 560:
            continue
        make_mol_mmp(new_m, id="cmp" + str(comp_ref.pk), target_id=None)


def refresh_mmp_maps(target_id):
    """Function to refresh maps for a target
    Takes a target id
    Returns None"""
    MMPDiffMap.objects.filter(target_id=target_id).delete()
    try:
        MMPComparison.objects.filter(target_id=target_id).delete()
    except:
        # Fix for SQLite
        mols = MMPComparison.objects.filter(target_id=target_id)
        for m in mols:
            m.delete()
    ActPharmaPoint.objects.filter(target_id=target_id).delete()
    print "RECHARGING"
    make_ph4_difference_points(target_id, opt_type=__POSS_LIST__)
    find_pharma_changes(target_id, my_type=str(__POSS_LIST__))
