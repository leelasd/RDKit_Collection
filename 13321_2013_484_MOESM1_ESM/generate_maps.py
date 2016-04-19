#
# Python script to generate similarity maps for four fingerprints:
# atom pairs (AP), Morgan2, CountMorgan2, and FeatMorgan2
# and two machine-learning methods: random forest and naive Bayes
#
# Example: ligands of dopamine D3 receptor extracted from ChEMBL
#
#
# Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import cPickle, gzip, numpy, copy, math
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, rdmolops
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import rdMolDescriptors
from matplotlib import cm
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import BernoulliNB

# helper functions
def getNormalizedWeights(weights):
    '''Normalizes a set of weight vectors'''
    for i in range(len(weights)):
        tmp = [math.fabs(j) for j in weights[i]]
        current_max = max(tmp)
        weights[i] = [j/current_max for j in weights[i]]
    return weights

def generateSimilarityMaps(mols, weights, fp):
    '''Generates a similarity map for a set of molecules and weights'''
    # colormap to use
    mycm = cm.PiYG
    # loop over molecules
    for i,m in enumerate(mols):
        fig = Draw.MolToMPL(m, coordScale=1.5, size=(250,250))
        # the values 0.02 and 0.01 can be adjusted for the size of the molecule
        x,y,z = Draw.calcAtomGaussians(m, 0.02, step=0.01, weights=weights[i])
        # use the maximum absolute peak as maximum scale
        maxscale = max(math.fabs(numpy.min(z)), math.fabs(numpy.max(z)))
        # this does the coloring
        fig.axes[0].imshow(z, cmap=mycm, interpolation='bilinear', origin='lower', extent=(0,1,0,1), vmin=-maxscale, vmax=maxscale)
        # this draws 10 contour lines
        # alternatively also the z values for the lines can be specified
        fig.axes[0].contour(x, y, z, 10, colors='k', alpha=0.5)
        # this writes the figure in a file
        fig.savefig('pics/mol'+str(i+1)+'_'+fp+'.png', bbox_inches='tight')

# hard coded variables
bit_size = 1024
radius = 2

# load reference compound and two test molecules
mols = []
for line in open('data/cmps.dat', 'r'):
    mols.append(Chem.MolFromSmiles(line))

# load training actives and inactives for ML methods
training = []
for line in open('data/training_cmps.dat', 'r'):
    line = line.rstrip().split()
    # line contains SMILES and active/inactive information
    training.append(line)

# precalculate fingerprints for reference compound
ref_morgan2 = AllChem.GetMorganFingerprintAsBitVect(mols[0], radius, bit_size)
ref_cmorgan2 = AllChem.GetMorganFingerprint(mols[0], radius)
ref_fmorgan2 = AllChem.GetMorganFingerprintAsBitVect(mols[0], radius, bit_size, useFeatures=True)
ref_ap = Pairs.GetAtomPairFingerprint(mols[0])

# precalculate fingerprints and bit information for test molecules
fps_morgan2 = []
fps_cmorgan2 = []
fps_fmorgan2 = []
fps_ap = []
info_morgan2 = []
info_cmorgan2 = []
info_fmorgan2 = []
num_mols = len(mols) - 1
mols = [mols[i+1] for i in range(num_mols)] # remove reference cmp from list
for m in mols:
    info = {}
    fps_morgan2.append(AllChem.GetMorganFingerprintAsBitVect(m, radius, bit_size, bitInfo=info))
    info_morgan2.append(info)
    info = {}
    fps_cmorgan2.append(AllChem.GetMorganFingerprint(m, radius, bitInfo=info))
    info_cmorgan2.append(info)
    info = {}
    fps_fmorgan2.append(AllChem.GetMorganFingerprintAsBitVect(m, radius, bit_size, useFeatures=True, bitInfo=info))
    info_fmorgan2.append(info)
    fps_ap.append(Pairs.GetAtomPairFingerprint(m))

### ATOM PAIRS
print "generate atom pairs similarity maps"
# calculate weights
mol_weights = []
for i,m in enumerate(mols):
    weights = []
    orig_simil = DataStructs.DiceSimilarity(ref_ap, fps_ap[i])
    matrix = rdmolops.GetDistanceMatrix(m)
    for at1 in range(m.GetNumAtoms()):
        new_fp = copy.deepcopy(fps_ap[i])
        for at2 in range(m.GetNumAtoms()):
            bit = Pairs.pyScorePair(m.GetAtomWithIdx(at1), m.GetAtomWithIdx(at2), matrix[at1][at2])
            new_fp[bit] -= 1
        new_simil = DataStructs.DiceSimilarity(ref_ap, new_fp)
        weights.append(orig_simil - new_simil)
    mol_weights.append(weights)
# normalization
mol_weights = getNormalizedWeights(mol_weights)
# draw similarity maps
generateSimilarityMaps(mols, mol_weights, 'ap')

### MORGAN2
print "generate morgan2 similarity maps"
# calculate weights
mol_weights = []
for i,m in enumerate(mols):
    weights = []
    orig_simil = DataStructs.DiceSimilarity(ref_morgan2, fps_morgan2[i])    
    # get bits for each atom
    bitmap = [~DataStructs.ExplicitBitVect(1024) for x in range(m.GetNumAtoms())]
    for bit, es in info_morgan2[i].iteritems():
        for at1, rad in es:
            if rad == 0: # for radius 0
                bitmap[at1][bit] = 0
            else: # for radii > 0
                env = Chem.FindAtomEnvironmentOfRadiusN(m, rad, at1)
                amap = {}
                submol = Chem.PathToSubmol(m, env, atomMap=amap)
                for at2 in amap.keys():
                    bitmap[at2][bit] = 0
    # loop over atoms
    for at1 in range(m.GetNumAtoms()):
        new_fp = fps_morgan2[i] & bitmap[at1]
        new_simil = DataStructs.DiceSimilarity(ref_morgan2, new_fp)
        weights.append(orig_simil-new_simil)
    mol_weights.append(weights)
# normalization
mol_weights = getNormalizedWeights(mol_weights)
# draw similarity maps
generateSimilarityMaps(mols, mol_weights, 'morgan2')

### COUNTMORGAN2
print "generate countmorgan2 similarity maps"
# calculate weights
mol_weights = []
for i,m in enumerate(mols):
    weights = []
    orig_simil = DataStructs.DiceSimilarity(ref_cmorgan2, fps_cmorgan2[i])
    # get bits for each atom
    bitmap = [[] for x in range(m.GetNumAtoms())]
    for bit, es in info_cmorgan2[i].iteritems():
        for at1, rad in es:
            if rad == 0: # for radius 0
                bitmap[at1].append(bit)
            else: # for radii > 0
                env = Chem.FindAtomEnvironmentOfRadiusN(m, rad, at1)
                amap = {}
                submol = Chem.PathToSubmol(m, env, atomMap=amap)
                for at2 in amap.keys():
                    bitmap[at2].append(bit)
    # loop over atoms
    for at1 in range(m.GetNumAtoms()):
        new_fp = copy.deepcopy(fps_cmorgan2[i])
        for bit in bitmap[at1]:
            new_fp[bit] -= 1
        new_simil = DataStructs.DiceSimilarity(ref_cmorgan2, new_fp)
        weights.append(orig_simil-new_simil)
    mol_weights.append(weights)
# normalization
mol_weights = getNormalizedWeights(mol_weights)
# draw similarity maps
generateSimilarityMaps(mols, mol_weights, 'cmorgan2')

### FEATMORGAN2
print "generate featmorgan2 similarity maps"
# calculate weights
mol_weights = []
for i,m in enumerate(mols):
    weights = []
    orig_simil = DataStructs.DiceSimilarity(ref_fmorgan2, fps_fmorgan2[i])    
    # get bits for each atom
    bitmap = [~DataStructs.ExplicitBitVect(1024) for x in range(m.GetNumAtoms())]
    for bit, es in info_fmorgan2[i].iteritems():
        for at1, rad in es:
            if rad == 0: # for radius 0
                bitmap[at1][bit] = 0
            else: # for radii > 0
                env = Chem.FindAtomEnvironmentOfRadiusN(m, rad, at1)
                amap = {}
                submol = Chem.PathToSubmol(m, env, atomMap=amap)
                for at2 in amap.keys():
                    bitmap[at2][bit] = 0
    # loop over atoms
    for at1 in range(m.GetNumAtoms()):
        new_fp = fps_fmorgan2[i] & bitmap[at1]
        new_simil = DataStructs.DiceSimilarity(ref_fmorgan2, new_fp)
        weights.append(orig_simil-new_simil)
    mol_weights.append(weights)
# normalization
mol_weights = getNormalizedWeights(mol_weights)
# draw similarity maps
generateSimilarityMaps(mols, mol_weights, 'fmorgan2')

### RANDOM FOREST WITH MORGAN2
print "generate random forest similarity maps"
# train random forest
training_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 2, 1024) for x,y in training]
training_labels = [y for x,y in training]
rf = RandomForestClassifier(n_estimators=100, max_depth=2, min_samples_split=2, min_samples_leaf=1, n_jobs=1)
rf.fit(training_fps, training_labels)
# calculate weights
mol_weights = []
for i,m in enumerate(mols):
    weights = []
    orig_pp = rf.predict_proba(fps_morgan2[i])[0][1]
    # get bits for each atom
    bitmap = [~DataStructs.ExplicitBitVect(1024) for x in range(m.GetNumAtoms())]
    for bit, es in info_morgan2[i].iteritems():
        for at1, rad in es:
            if rad == 0: # for radius 0
                bitmap[at1][bit] = 0
            else: # for radii > 0
                env = Chem.FindAtomEnvironmentOfRadiusN(m, rad, at1)
                amap = {}
                submol = Chem.PathToSubmol(m, env, atomMap=amap)
                for at2 in amap.keys():
                    bitmap[at2][bit] = 0
    # loop over atoms
    for at1 in range(m.GetNumAtoms()):
        new_fp = fps_morgan2[i] & bitmap[at1]
        new_pp = rf.predict_proba(new_fp)[0][1]
        weights.append(orig_pp-new_pp)
    mol_weights.append(weights)
# normalization
mol_weights = getNormalizedWeights(mol_weights)
# draw similarity maps
generateSimilarityMaps(mols, mol_weights, 'rf')

### NAIVE BAYES WITH MORGAN2
print "generate naive bayes similarity maps"
# train random forest
nb = BernoulliNB()
nb.fit(training_fps, training_labels)
# calculate weights
mol_weights = []
for i,m in enumerate(mols):
    weights = []
    orig_pp = nb.predict_log_proba(fps_morgan2[i])[0][1]
    # get bits for each atom
    bitmap = [~DataStructs.ExplicitBitVect(1024) for x in range(m.GetNumAtoms())]
    for bit, es in info_morgan2[i].iteritems():
        for at1, rad in es:
            if rad == 0: # for radius 0
                bitmap[at1][bit] = 0
            else: # for radii > 0
                env = Chem.FindAtomEnvironmentOfRadiusN(m, rad, at1)
                amap = {}
                submol = Chem.PathToSubmol(m, env, atomMap=amap)
                for at2 in amap.keys():
                    bitmap[at2][bit] = 0
    # loop over atoms
    for at1 in range(m.GetNumAtoms()):
        new_fp = fps_morgan2[i] & bitmap[at1]
        new_pp = nb.predict_log_proba(new_fp)[0][1]
        weights.append(orig_pp-new_pp)
    mol_weights.append(weights)
# normalization
mol_weights = getNormalizedWeights(mol_weights)
# draw similarity maps
generateSimilarityMaps(mols, mol_weights, 'nb')

