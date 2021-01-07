from __future__ import print_function
from __future__ import division

from pyLM.units import *

import h5py
import numpy as np
import os

filename_lm='morph_dend.lm'
filename_morph='CA1dend_small.h5'


NA = 6.022e23
##
##
##

print('Set a simulation space.')
latticeSpacing=nm(8) 
#sim=RDME.RDMESimulation(dimensions=micron(2.048,2.048,2.048), \
#                        spacing=latticeSpacing)

print('latticeSpacing: ', latticeSpacing)

print('Define regions')
print('Insert dendritic geometry.')
with h5py.File(filename_morph,'r') as r:
    dendrite                   = r['dendrite'][()]
    dendrite_not_mitochondrion = r['dendrite not mitochondrion'][()]
    PSD                        = r['PSD'][()]
    membrane_area              = r['membrane areas'][()]
    pitch          = r['unit length per voxel (um)'][()]
    num_voxel_dend = r['voxel num of dendrite not mitochondrion'][()]

volume    = np.sum(dendrite > 0)
#morpho    = buildAnyShape(sim, volume, domains, membrane_area, PSD)
print('num_voxel_dend: ', volume)
print('num_voxel_dend: ', num_voxel_dend)
volume_in_L  = num_voxel_dend * latticeSpacing * latticeSpacing * latticeSpacing * 1000
volume_in_m3 = num_voxel_dend * latticeSpacing * latticeSpacing * latticeSpacing
print('Volume in um3: ',volume_in_m3*(1e6)*(1e6)*(1e6))
print('Volume in L: ', volume_in_L )
print('Volume in fL: ', volume_in_L *1e15 )

#   6000 moleules per 100uM CaM and 0.1fL Spine

number_1umol = NA /(1e6)
number_1uM   = number_1umol * volume_in_L
print('number_per_1uM:', number_1uM)


print('Set molecules')
conc_Ca  = 100 # uM
conc_CaM = 100 # uM
conc_CB  = 120 # uM
conc_CN  = 0.5 # uM
num_Ca   = int(conc_Ca  * number_1uM)
num_CaM  = int(conc_CaM * number_1uM) 
num_CB   = int(conc_CB  * number_1uM)
num_CN   = int(conc_CN  * number_1uM)

print('num_Ca : ', num_Ca )
print('num_CaM: ', num_CaM)
print('num_CB : ', num_CB )
print('num_CN : ', num_CN )
