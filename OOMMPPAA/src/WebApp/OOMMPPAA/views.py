from django.http import HttpResponse,HttpResponseRedirect
from django.shortcuts import render, get_object_or_404
from django.views.decorators.csrf import csrf_exempt
import tempfile
import WebApp.views as WV
from MMPMaker.functions import find_pharma_changes
from MMPMaker.models import MMPDiffMap
from IOhandle.models import Protein,Target,Project,Molecule
from rdkit import Chem
from rdkit.Chem import AllChem
import json
import os
import sys


def index(request):
    targets = Target.objects.exclude(title="DUMMY")
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets}
    return render(request, 'OOMMPPAA/index.html', context)


def fileupload(request):
    WV.arg_list = {}
    targets = Target.objects.exclude(title="DUMMY")
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets, "DATADIR": os.path.join(os.path.split(sys.argv[0])[0], "data") }
    return render(request, 'OOMMPPAA/fileupload.html', context)


def make_dialoge(my_type, message, out_f, element):
    make_out = '''<div class="alert alert-dismissable alert-''' + my_type + '''">
<button type="button" class="close" data-dismiss="alert">x</button>
''' + message + """</a>
</div>"""
    out_d = {"HTML": make_out, "SUCCESS": my_type, "FILENAME": out_f, "ELEMENT": element}
    return HttpResponse(json.dumps(out_d))


@csrf_exempt
def uploadme(request):
  # So now validate that the file is appropriate
    if 'activity' in request.FILES or 'append-activity' in request.FILES:#
        element = 'activity'
        try:
            my_file = request.FILES['activity']
        except KeyError:
            my_file = request.FILES['append-activity']
        # Check it's a CSV file
        file_text = my_file.read()
        out_f = tempfile.NamedTemporaryFile(mode='w', delete=False)
        out_f.write(file_text)
        header = file_text.split("\n")[0].rstrip()
        headings = header.split(",")
        print headings
        all_fields = ["smiles", "Activity", "ID", "operator"]
        req_fields = ["smiles", "Activity"]
        if len(headings) == 1:
            return make_dialoge("danger", "Error - not CSV file")
        # Check it has the required fields
        elif len([x for x in req_fields if x not in headings]) > 0:
            return make_dialoge("danger", "Error - missing required fields: " + ",".join([x for x in req_fields if x not in headings]), out_f.name, element)
        # Check it has the optional fields - warning if not
        elif len([x for x in all_fields if x not in headings]) > 0:
            return make_dialoge("warning", "Warning - missing required fields: " + ",".join([x for x in all_fields if x not in headings]), out_f.name, element)
        else:
            return make_dialoge("success", "Success - all fields in CSV file", out_f.name, element)
    elif 'molecules' in request.FILES or 'append-molecules' in request.FILES:
        element = 'molecules'
        try:
            my_file = request.FILES['molecules']
        except KeyError:
            my_file = request.FILES['append-molecules']
        out_f = tempfile.NamedTemporaryFile(mode='w', delete=False)
        out_f.write(my_file.read())
        out_f.close()
        # Check it's an SD file - warning if some molecules are not valid
        mols = AllChem.SDMolSupplier(out_f.name)
#        in_sd.close()
        err_counter = 0
        good_counter = 0

        for m in mols:
            if m is None:
                err_counter += 1
            else:
                good_counter +=1
        if good_counter == 0:
            return make_dialoge("danger", "Error - Not valid SD file - check input using RDKit", out_f.name, element)
        elif err_counter > 0:
            return make_dialoge("warning", "Warning - " + str(err_counter) + " molecule out of " + str(len(mols)) + " read incorrectly", out_f.name, element)
        else:
            return make_dialoge("success", "Success - SD file read correctly", out_f.name, element)
    elif 'protein' in request.FILES or 'append-protein' in request.FILES:
        element = 'protein'
        try:
            my_file = request.FILES['protein']
        except KeyError:
            my_file = request.FILES['append-protein']
        # No need to check if it's a valid protein - not required for processing
        out_f = tempfile.NamedTemporaryFile(mode='w', delete=False)
        my_text = my_file.read()
        out_f.write(my_text)
        out_f.close()
        if len(my_text.split("\n")) > 10:
            return make_dialoge("success", "Success - probably a real protein file", out_f.name, element)
        else:
            return make_dialoge("warning", "Warning - not a real protein file", out_f.name, element)
    else:
        return make_dialoge("", "", "", "")


def success(request):
    targets = Target.objects.exclude(title="DUMMY")
    # No projects at the moment as this hasn't been set up yet
    context = {"targets": targets}
    return render(request, 'OOMMPPAA/success.html', context)


def demoviewer(request, target_id):
    target = get_object_or_404(Target, pk=target_id)
    # Add in a filter for proteins with actual molecules
    targets = Target.objects.exclude(title="DUMMY")
    maps = MMPDiffMap.objects.filter(target_id=target_id).filter(activity_type=str("['IC50', 'Ki']"))
    context = {"target": target, "maps": maps, "targets": targets}
    # The class defining the list of ways of
    return render(request, 'OOMMPPAA/demoviewer.html', context)


def viewer(request, target_id):
    target = get_object_or_404(Target, pk=target_id)
    target.mol = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)[0].prot_id.code
    targets = Target.objects.exclude(title="DUMMY")
    mols = Molecule.objects.filter(prot_id__target_id=target_id).exclude(prot_id__code__startswith=target.title)
    # Add in a filter for proteins with actual molecules
    maps = MMPDiffMap.objects.filter(target_id=target_id).filter(activity_type=str("['IC50', 'Ki']"))
    context = {"mols":  mols, "target": target, "maps": maps, "targets": targets}
    # The class defining the list of ways of
    return render(request, 'OOMMPPAA/viewer.html', context)