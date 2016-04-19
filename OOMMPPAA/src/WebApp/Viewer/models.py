import matplotlib
# Required for the Django app
matplotlib.use('Agg')
from django.db import models
from IOhandle.models import Protein,Molecule,Target,ActivityPoint,Compound
from Pharmacophore.models import PharmaPoint
from MMPMaker.models import MMPDiffMap,MMPComparison,ActMapPoint,MMPFrag,MMP,ActPharmaPoint
from MMPMaker.functions import find_points
#from MakeMap import main as makemap
from math import exp,sqrt,pow
import ast,string,random,StringIO,gzip,sys,numpy,math,os,tempfile

from base64 import encodestring
from rdkit.Chem import RDConfig,ChemicalFeatures,AllChem,Draw,MCS
from rdkit.Chem.Draw import SimilarityMaps
from rdkit import DataStructs,Chem
from pylab import cm,colorbar,subplot,get_cmap,ones,imshow,arange
try:
    import Image,ImageDraw,ImageFont
except ImportError:
    from PIL import Image,ImageDraw,ImageFont


def eucl_dist(mol_1, mol_2):
    """Function to find the euclidean distance between two points.
    Takes two objects with attribute x_com, y_com and z_com
    Returns a float"""
    x_squared = pow((mol_1.x_com - mol_2.x_com), 2)
    y_squared = pow((mol_1.y_com - mol_2.y_com), 2)
    z_squared = pow((mol_1.z_com - mol_2.z_com), 2)
    dist = sqrt(x_squared + y_squared + z_squared)
    return dist


class ViewHandler():
    """Class to deal with functions for viewing molecules/proteins in 
    different ways.
    Consists of a dict relating calls from the web to functions here.
    Arguments to functions are handled in views.py"""
    # These are the functions to be used here for IOhandling
    def __init__(self):
        # Routines for loading files in as PDBs
        # All these routines take an option of PDB_code(either 1, a list or all) and return a PDB file
        self.pdbroutine_list = {"CHECKPOINTS": check_points,
                                "VIEWPROTEIN": view_protein,
                                "VIEWALLMOLS": view_all_mols,
                                "VIEWMOL": view_mol, "CHECKMOL": check_mol,
                                "GETMMP": get_mmp, "GETMAP": get_map,
                                "MOLS": get_mol, "2DMOL": view_2dmol,
                                "MAKESIM": make_similarity_map, "REFRESH": refresh_up,
                                "VIEWLOG": view_log}


def view_log(option=None, maps=None, out_put=None, target_id=None, extra=None):
    """Function to view a log file"""
    # Get the target we're referring to
    t = Target.objects.get(pk=target_id)
    print option
    gen_path = os.path.split(sys.argv[0])[0]
    # Get the path to the log file
    if int(option) == 1:
      log_path = os.path.join(os.path.split(sys.argv[0])[0], t.title + ".std")
      print log_path
      if os.path.isfile(log_path):
          return "<strong> log files available:</strong> \n " + os.path.join(gen_path,t.title+".std") + "\n" + os.path.join(gen_path, "out.stderr")
      else:
          return "No log file available"
#    elif int(option) == 2:
#      # Get the path for this logfikle
#      log_path = os.path.join(os.path.split(sys.argv[0])[0], "out.std")
#      print log_path
#      return "<strong> Other log files available: " + log_path + "</strong> \n\n" + open(log_path).read()
#    elif int(option) == 3:
#      log_path = os.path.join(os.path.split(sys.argv[0])[0], "out.stderr")
#      print log_path
#      return "<strong> Other log files available: " + log_path + "</strong> \n\n" + open(log_path).read()


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


def refresh_up(option=None, maps=None, out_put=None, target_id=None, extra=None):
    """Function to refresh the upload page AND to delete any targets that need deleting"""
    if option == "DELETE":
        # I
        delete_target(target_id)
    # First delete what needs deleting 
    ts = Target.objects.exclude(title="DUMMY")
    out_s = """<li ><a onclick="document.getElementById('buttons-one').style.display='none';document.getElementById('input-id').style.display='';$this.addClass('active');";>Add new target</a></li>"""
    for target in ts:
        out_s += '''<li class="dropdown">
    <a class="dropdown-toggle" data-toggle="dropdown" href="#">
      ''' + target.title + '''<span class="caret"></span>
    </a>
    <ul class="dropdown-menu">
      <li><a onclick="document.getElementById('inputDefault').value='''+"'"+target.title+"'"+''';document.getElementById('buttons-one').style.display='';document.getElementById('input-id').style.display='none';" >Edit</a></li>
      <li><a onclick="refresh_info('''+"'"+str(target.pk)+"'"+''','delete','''+"'"+target.title+"'"+''')">Delete</a></li>
      <li class="divider"></li>
      <li><a href="#">Log file</a></li>
    </ul>
  </li>'''
    return out_s


def make_similarity_maps(mol, weights, colorMap=cm.PiYG, scale=-1, size=(250, 250), sigma=None,coordScale=1.5, step=0.01, colors='k', contourLines=10, alpha=0.5, **kwargs):
    """Function to calculate similarity maps
    Heavily based on the similarity map function in the RDKit. A few changes
    to deal with exceptions and change rendering,
    Takes an RDKit molecule and a list of atom-based weights.
    Returns an image."""
    if mol.GetNumAtoms() < 2:
        raise ValueError("too few atoms")
    fig = Draw.MolToMPL(mol, coordScale=coordScale, size=size, **kwargs)
    if sigma is None:
        if mol.GetNumBonds() > 0:
            bond = mol.GetBondWithIdx(0)
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            sigma = 0.3 * math.sqrt(sum([(mol._atomPs[idx1][i] - mol._atomPs[idx2][i]) ** 2 for i in range(2)]))
        else:
            sigma = 0.3 * math.sqrt(sum([(mol._atomPs[0][i] - mol._atomPs[1][i]) ** 2 for i in range(2)]))
        sigma = round(sigma, 2)
    x, y, z = Draw.calcAtomGaussians(mol, sigma, weights=weights, step=step)
    # scaling
    if scale <= 0.0:
        maxScale = max(math.fabs(numpy.min(z)), math.fabs(numpy.max(z)))
    else:
        maxScale = scale
    # coloring
    cax = fig.axes[0].imshow(z, cmap=colorMap, interpolation='bilinear', origin='lower', extent=(0,1,0,1), vmin=-maxScale, vmax=maxScale)
    cbar = fig.colorbar(cax, shrink=.75, pad=.02,ticks=[-maxScale, 0, maxScale], orientation='vertical')
    cbar.ax.set_yticklabels(['', '', ''])  # contour lines
    fig.axes[0].contour(x, y, z, contourLines, colors=colors, alpha=alpha, **kwargs)
    return fig


def view_2dmol(option, maps=None, out_put=None, target_id=None, legends=None):
    """Function to render a mol image from a smiles.
    The input (option) could be 1) a list of smiles 2) a smiles 3) pdb_code
    Returns a molecule image as data"""
    option = str(option)
    print option
    try:
        option = ast.literal_eval(option)
    except:
        pass
    if type(option) is list:
        mols = [Chem.MolFromSmiles(str(x)) for x in option]
        p = Chem.MolFromSmarts(MCS.FindMCS(mols).smarts)
        AllChem.Compute2DCoords(p)
        [AllChem.GenerateDepictionMatching2DStructure(x, p) for x in mols]
        image = Draw.MolsToGridImage(mols, 2, legends=legends)
    elif Chem.MolFromSmiles(str(option)) is None and type(option) is str:
        mol = Chem.MolFromSmiles(str(Molecule.objects.filter(prot_id__code=option)[0].cmpd_id.smiles))
        image = Draw.MolToImage(mol)
    elif type(option) is str:
        mol = Chem.MolFromSmiles(str(option))
        image = Draw.MolToImage(mol)
    else:
        print "NOT VALID TYPE"
        return "NOT VALID TYPE"
    output = StringIO.StringIO()
    image.save(output, format="PNG")
    contents = output.getvalue()
    return contents


def view_protein(option, maps=None, out_put=None, target_id=None):
    """Function to return the PDB file for a protein given it's code
    Takes a code for a Protein object
    Returns an open file handle"""
    my_protein = Protein.objects.get(code=option)
    try:
        my_val = open(str(my_protein.pdb_info.file))
    except IOError:
        my_val = str(my_protein.pdb_info)
    return my_val


def view_all_mols(option, maps=None, out_put=None, target_id=None):
    """Function to show the SD block of all the molecules for a target
    Takes a Target as input
    Returns an SD block of all the molecules for that target.
    """
    my_mols = Molecule.objects.filter(prot_id__target_id=int(option))
    new_mol = ""
    for mol in my_mols:
        new_mol += (str(mol.sdf_info))+"\n\n$$$$\n"
    return new_mol


def get_mmp(option, maps=None, out_put=None, target_id=None):
    """Function to find the MMP comparison
    Takes a pk of the MMP comparison
    Returns the SD block of that data"""
    mmp = MMPComparison.objects.get(pk=option)
    return mmp.sdf_info


def get_map(option, maps=None, out_put=None, target_id=None):
    """Function to show the information for a given map
    Takes the pk for the map
    Returns the PDB information for that map as a GZIP"""
    mymap = MMPDiffMap.objects.get(pk=option)
    out = StringIO.StringIO()
    f = gzip.GzipFile(fileobj=out, mode='w')
    f.write(str(mymap.pdb_info))
    f.close()
    contents = out.getvalue()
    return contents


def check_mol(option, maps=None, out_put=None, target_id=None):
    """Function to check whether an input is either a valid smiles or a valid 
    protein code
    Takes a string and a Target
    Returns an answer to be used by jquery"""
    my_mols = Molecule.objects.filter(prot_id__code__icontains=option)
    target = Target.objects.get(pk=target_id)
    if len(my_mols) == 0:
        tmpmol = Chem.MolFromSmiles(str(option))
        if tmpmol is None:
            return "None molecule"
        # Now do a similarity search on this against all the molecules
        cmps = [Chem.MolFromSmiles(str(x)) for x in Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title).exclude(cmpd_id__smiles="DUMMY").values_list("cmpd_id__smiles", flat=True)]
        fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024) for x in cmps]
        sims = DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprintAsBitVect(tmpmol, 2, nBits=1024), fps)
        ind = max(enumerate(sims), key=lambda x: x[1])[0]
        mycmp = cmps[ind]
        my_mols = Molecule.objects.filter(cmpd_id__smiles=Chem.MolToSmiles(mycmp, isomericSmiles=True))
    # Now return the appropriate PDBcode
    return my_mols[0].prot_id.code


def check_points(option, maps=None, out_put=None, target_id=None):
    """Function to check how many points correspond to a selection
    Takes the min and max activity change and the min and max number of pharmacophore point differences in option
    Returns the number of points"""
    spl_opt = option.split(",")
    points = ActPharmaPoint.objects.filter(target_id=target_id, in_diff_15=True, num_diff_15__gte=spl_opt[0], num_diff_15__lte=spl_opt[1], act_change__gte=spl_opt[2], act_change__lte=spl_opt[3])
    return len(points)


def view_mol(option, maps=None, out_put=None, target_id=None):
    """Function to render the 3D coordinates of a molecule
    Takes a PDB code as input
    Returns an SD block"""
    my_mols = Molecule.objects.filter(prot_id__code=option, prot_id__target_id=target_id)
    new_mol = ""
    for mol in my_mols:
        new_mol += (str(mol.sdf_info)) + "\n\n$$$$\n"
    return new_mol


def find_act_map_points(option, maps=None, out_put=None, target_id=None):
    """Function to return a series of activity map points from a map and a JSON dict
    Option is a json dict of points and maps is the indices of maps to use
    Returns the points that apply."""
    maps = MMPDiffMap.objects.filter(pk__in=[int(x) for x in maps.split(",")])
    act = maps[0]
    inact = maps[1]
    shape = maps[2]
    my_d = ast.literal_eval(option)
    if "my_res" in my_d:
        mols = my_d["my_res"].split("|")
        mol_d = {}
        # Now make a dictionary of the inds for the maps
        for mol in mols:
            if mol.split("_")[1].split(".//")[0] not in ["act", "inact", "shape"]:
                continue
            mol_d[mol.split("_")[1].split(".//")[0]] = get_list_inds(mol.split("_")[1].split(".//")[1])
        my_points = []
        if "act" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=act, point_id__in=mol_d["act"]))
        if "inact" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=inact, point_id__in=mol_d["inact"]))
        if "shape" in mol_d:
            my_points.extend(ActMapPoint.objects.filter(map_id=shape, point_id__in=mol_d["shape"]))
        # Find things with a CofM within these ranges
        return my_points
    elif "my_var" in my_d:
        points = [float(x) for x in my_d["my_var"]]
        return points[3], points[0], points[4], points[1], points[5], points[2]


def render_act(act):
    """Function to render activity data
    Takes a float of activity data and the operator for that data 
    Returns the data rendered to 2DP with activity associated"""
    # Deal with inactive data
    if act.operator == "<":
        return act.source + " Inactive: under {0:.2f}".format(act.activity)
    else:
        return act.source + " Activity: {0:.2f}".format(act.activity)


def draw_acts(mols, legends, subs, molsPerRow=2, subImgSize=(200, 200), **kwargs):
    """Function to draw a grid of activity data. 
    Based very strongly on the RDKit Draw.MolsToGridImage()
    Takes a list of RDKit molecules, assoicated legends, the associated 
    substituions, the numnber of mols per row and the size
    Returns a PIL image"""
    try:
        import Image, ImageDraw
    except ImportError:
        from PIL import Image, ImageDraw
    nRows = len(mols) // molsPerRow
    # Add the extra one if they are divisble
    if len(mols) % molsPerRow:
        nRows += 1
    res = Image.new("RGB", (molsPerRow * subImgSize[1], nRows * subImgSize[1]), (255, 255, 255))
    for i, mol in enumerate(mols):
        highlightAtoms = mol.GetSubstructMatch(subs[i // 2])
        row = i // molsPerRow
        col = i % molsPerRow
        highlightMap = {}
        for idx in highlightAtoms:
            highlightMap[idx] = (0.7, 0.7, 0.7)
        molimage = Draw.MolToImage(mol, subImgSize, legend=legends[i], highlightMap=highlightMap)
        res.paste(molimage, (col * subImgSize[0], row * subImgSize[1]))
    draw = ImageDraw.Draw(res)
# draw.line((res.size[0]/2, 0,res.size[0]/2,res.size[1]), width=2,fill=(0,0,0))
    for i in range(len(mols) / 2):
        if i == 0:
            draw.line((0, 0, res.size[0], 0), width=2, fill=(0, 0, 0))
        else:
            draw.line((0, res.size[1] * i / nRows, res.size[0], res.size[1] * i / nRows ), width=2, fill=(0, 0, 0))
    draw.line((0, res.size[1], res.size[0], res.size[1]), width=2, fill=(0, 0, 0))
    return res


def get_mol(option, maps=None, out_put=None, target_id=None):
    """Function to return an image or SD block of a molecule
    Option is a dict => values to find, Maps is a list of map pks to cross reference back,
    out_put indicates whether to return SD block or image
    Returns an image or an SD block of the molecules that have been selected."""
    # First get my maps
    my_points = find_act_map_points(option, maps)
    # Now get all the mols that are associated to this
    mmpcomps = MMPComparison.objects.filter(actmappoint__in=my_points).distinct()
    #Now make the molecule SD information
    new_mol = ""
    i = 0
    # Now make the image to look at
    if out_put == "images":
        rdmols = []
        acts = []
        subs = []
        smsubs = []
        donemols = []
        for m in mmpcomps:
            # Ensure that this is not just the same comparison
            mol1 = Chem.MolFromMolBlock((str(m.xtal_mol.sdf_info)))
            mol2 = Chem.MolFromMolBlock((str(m.chembl_mol.sdf_info)))
            if [m.xtal_mol.cmpd_id.pk, m.chembl_mol.pk] in donemols or [m.chembl_mol.cmpd_id.pk, m.xtal_mol.pk] in donemols:
                # Don't do the same comparison twice
                continue
            else:
                donemols.append([m.xtal_mol.cmpd_id.pk, m.chembl_mol.cmpd_id.pk])

            # Set the molecule name for the 3D display
            mol1.SetProp("_Name", "inactMOL" + ''.join(random.sample(string.letters * 5, 5)))
            acts.append(render_act(m.xtal_act))
            acts.append(render_act(m.chembl_act))
            mol2.SetProp("_Name", "actMOL" + ''.join(random.sample(string.letters * 5, 5)))
            # Generate the two-d depictions after canonicalising the smiles
            mol1 = Chem.MolFromSmiles(Chem.MolToSmiles(mol1, isomericSmiles=True))
            mol2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol2, isomericSmiles=True))
            smp = MCS.FindMCS([mol1, mol2], completeRingsOnly=True, ringMatchesRingOnly=True, timeout=0.5).smarts
            p = Chem.MolFromSmarts(smp)
            subs.append(p)
            smsubs.append(smp)
            AllChem.Compute2DCoords(p)
            AllChem.GenerateDepictionMatching2DStructure(mol1, p, acceptFailure=True)
            AllChem.GenerateDepictionMatching2DStructure(mol2, p, acceptFailure=True)
            rdmols.extend([mol1, mol2])
        # So now we have the mols in a list with actvity information in a list
        # Order this list of molecules based on scaffold (p)
        # Get a list of the indices of rdmols to rearrange
        myinds = sorted(range(len(smsubs)), key=lambda x:smsubs[x])
        nmols = []
        nacts = []
        nsubs = []
        # Now rearrange everthing to suit
        for ind_m in myinds:
            nmols.extend([rdmols[ind_m * 2], rdmols[ind_m * 2 + 1]])
            nacts.extend([acts[ind_m * 2], acts[ind_m * 2 + 1]])
            nsubs.append(subs[ind_m])
        image = draw_acts(nmols, nacts, nsubs)
        output = StringIO.StringIO()
        image.save(output, format="PNG")
        contents = output.getvalue()
        return contents
    elif out_put == "sds":
        for mol in mmpcomps:
            i += 1
            rd_mol = Chem.MolFromMolBlock((str(mol.xtal_mol.sdf_info)))
            rd_mol.SetProp("_Name", "MOL" + str(i))
            new_mol += Chem.MolToMolBlock(rd_mol) + "\n\n$$$$\n"
            i += 1
            rd_mol = Chem.MolFromMolBlock((str(mol.chembl_mol.sdf_info)))
            rd_mol.SetProp("_Name", "MOL" + str(i))
            new_mol += Chem.MolToMolBlock(rd_mol) + "\n\n$$$$\n"
        return new_mol


def get_empty_image(error_message="Invalid Molecule"):
    """Function to return an empyty image if molecule rendering fails
    Takes an error message
    Returns an empty image with that message."""
    try:
        import Image, ImageDraw, ImageFont
    except ImportError:
        from PIL import Image, ImageDraw, ImageFont
    out_map = Image.new("RGB", (240, 240), (255, 255, 255))
    draw = ImageDraw.Draw(out_map)
    draw.text((out_map.size[0] / 2, out_map.size[1] / 2), error_message, (0, 0, 0))
    output = StringIO.StringIO()
    out_map.save(output, format="PNG")
    contents = output.getvalue()
    return contents


def make_similarity_map(option, maps=1, out_put=None, target_id=None):
    """Function to make a sim map for a given PDB code relating to a molecule
    Takes option as the code of the molecules protein, choice indicates the 
    option for the map, target_id is the Target
    Returns a png image  as data"""
    # Find the molecule from the option i
    dj_mols = Molecule.objects.filter(prot_id__code=option)
    if len(dj_mols) != 0:
        dj_mol = dj_mols[0]
    else:
        # Assume it is a smiles string
        tmpmol = Chem.MolFromSmiles(str(option))
        if tmpmol is None:
            # If it is not then return an empty image
            #Except if the values are all nought!
            return get_empty_image()
        # Now do a similarity search on this against all the molecules
        cmps = [Chem.MolFromSmiles(str(x)) for x in Molecule.objects.exclude(cmpd_id__smiles="DUMMY").filter(prot_id__target_id=target_id).values_list("cmpd_id__smiles", flat=True)]
        if len(cmps) == 0:
            return get_empty_image()
        fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024) for x in cmps]
        sims = DataStructs.BulkTanimotoSimilarity(AllChem.GetMorganFingerprintAsBitVect(tmpmol, 2, nBits=1024), fps)
        ind = max(enumerate(sims), key=lambda x: x[1])[0]
        mycmp = cmps[ind]
        dj_mols = Molecule.objects.filter(prot_id__target_id=target_id).filter(cmpd_id__smiles=Chem.MolToSmiles(mycmp, isomericSmiles=True))
        dj_mol = dj_mols[0]
    # Locate the target from this
    target_id = dj_mol.prot_id.target_id
    # Now iterate over the atoms in this molecule. Everyone starts of with a score of 100. Everyone is in a pharmacophore group. If it contributes to a positive pharmacophore group you add 10. If it contributes to a negative pharmacophore group you minus 10.
    my_mol = Chem.MolFromMolBlock(str(dj_mol.sdf_info))
    atms = my_mol.GetAtoms()
    weights = {a.GetIdx(): 0.1 for a in atms}
    # Find all the atoms matching the pharmacophore substructures for this type
    fdefName = os.path.join(os.path.join(os.path.split(sys.argv[0])[0], 'data/media'), 'BaseFeatures.fdef')
    if not os.path.isfile(fdefName):
        fdefName = os.path.join('data', 'media', 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    features = factory.GetFeaturesForMol(my_mol)
    # These maps indicate the properties to do with act, inact maps
    if int(maps) == 1:
        ## These are hard coded
        my_maps = MMPDiffMap.objects.filter(type__in=["act", "inact"], target_id=target_id,  activity_type=str(["IC50", "Ki"]))
        act_points = my_maps[0].actmappoint_set.all()
        inact_points = my_maps[1].actmappoint_set.all()
    # This map indicates pharmacophoric interests
    elif int(maps) == 3:
        my_points = PharmaPoint.objects.filter(mol_id__prot_id__target_id=target_id)
    if int(maps) == 1:
        for x in weights:
            weights[x] = 0.1
        # Now loop through the features
        for feat in features:
            fp = PharmaPoint()
            fp.x_com = feat.GetPos().x
            fp.y_com = feat.GetPos().y
            fp.z_com = feat.GetPos().z
            fp.type = feat.GetType()
            fp.ids = feat.GetAtomIds()
            # Now do the same again for the inactivity points
            for inact in inact_points:
                #
                if inact.mmp_comparsion_id.activity_change < 0.5 or inact.mmp_comparsion_id.num_change > 3:
                    continue
                pp = inact.pharma_id
                # Calculate the euclidean distance between this featur position and the other pharmapoistion
                if eucl_dist(pp, fp) > 4:
                    continue
                # If they are the same type
                if str(pp.uniq_id.smiles) == str(fp.type):
                    for x in fp.ids:
                        weights[x] -= exp(-eucl_dist(pp, fp)) * 300
    elif int(maps) == 3:
        # A map showing the pharmacophoric fit to the data -> i.e. conserved pharmacophoric features
        for x in weights:
            weights[x] = 0.1
        # Now loop through the features
        for feat in features:
            fp = PharmaPoint()
            fp.x_com = feat.GetPos().x
            fp.y_com = feat.GetPos().y
            fp.z_com = feat.GetPos().z
            fp.type = feat.GetType()
            fp.ids = feat.GetAtomIds()
            # Go through the activity points
            # Now loop through the features
            for pp in my_points:
                # Calculate the euclidean distance between this feature position and the other pharmapoistion
                if eucl_dist(pp, fp) > 4.0:
                    continue
                if str(pp.uniq_id.smiles) == str(fp.type):
                    for x in fp.ids:
                        weights[x] += exp(-eucl_dist(pp, fp)) * 100
    # Now do the explicity H's for choice 1
    if int(maps) == 1:
        # Add explicit H's
        hmy_mol = AllChem.AddHs(my_mol)
        AllChem.ConstrainedEmbed(hmy_mol, my_mol)
        my_mol = hmy_mol
        emy_mol = AllChem.EditableMol(my_mol)
        # And get this Conf
        conf = my_mol.GetConformer()
        # Loop through the explicit H's in the model to produce the extension list
        atoms = my_mol.GetAtoms()
        for atm in atoms:
            # Find the H's
            if atm.GetAtomicNum() == 1:
                # Now add this to the weight list
                ap = PharmaPoint()
                ap.id = atm.GetIdx()
                emy_mol.ReplaceAtom(ap.id, Chem.Atom(53))
                weights[ap.id] = 0.1
                ap.x_com = conf.GetAtomPosition(ap.id).x
                ap.y_com = conf.GetAtomPosition(ap.id).y
                ap.z_com = conf.GetAtomPosition(ap.id).z
                for act in act_points:
                    pp = act.pharma_id
                    weights[ap.id] += exp(-eucl_dist(pp, ap)) * 100
        my_mol = emy_mol.GetMol()
        Chem.SanitizeMol(my_mol)

    # Now compute 2D coords
    AllChem.Compute2DCoords(my_mol)
    output = StringIO.StringIO()
    try:
        if int(maps) == 1:
            out_map = make_similarity_maps(my_mol, weights, colorMap=cm.RdBu, alpha=0.00)
            out_map.savefig(output, format="PNG", bbox_inches='tight', dpi=35)
        else:
            out_map = make_similarity_maps(my_mol, weights, colorMap=cm.jet, alpha=0.00)
            out_map.savefig(output, format="PNG", bbox_inches='tight', dpi=35)
        contents = output.getvalue()
        return contents
    except ValueError:
        #Except if the values are all nought!
        return get_empty_image()


def get_list_inds(my_ls):
    """Function to find the indexes for points based on 1:4,5,7:9 type notation.
    Takes a string.
    Returns a list"""
    ls = my_ls.split(",")
    me = [range(int(x.split(":")[0]), int(x.split(":")[1]) + 1) for x in ls if ":" in x]
    import itertools
    me_out = list(itertools.chain(*me))
    me_out.extend([int(x) for x in ls if ":" not in x])
    return me_out
