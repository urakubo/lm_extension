from __future__ import print_function
from __future__ import division

from pyLM import *
from pyLM.units import *
import lmUtils
from setMolecules import setMolecules
from pySTDLM import *
# from pySTDLM.PostProcessing import *
from pySTDLM.StandardReactionSystems import *
# from pySTDLM.StandardCells import *
import h5py
import numpy as np
import os

filename_lm='morph_dend.lm'
filename_morph='CA1dend_small.h5'

if os.path.exists(filename_lm):
    os.system('rm '+filename_lm)

ext  = 'default'
cyt  = 'cytoplasm'
psd  = 'psd'
domains = {ext: 0, cyt: 1, psd: 2}

NA = 6.022e23
##
##
##

print('Set a simulation space.')
latticeSpacing=nm(8) 
sim=RDME.RDMESimulation(dimensions=micron(2.048,2.048,2.048), \
                        spacing=latticeSpacing)


print('Define regions')
print('Insert dendritic geometry.')
with h5py.File(filename_morph,'r') as r:
    dendrite                   = r['dendrite_not_PSD'][()]
    dendrite_not_mitochondrion = r['dendrite not mitochondrion'][()]
    PSD                        = r['dendrite_and_PSD'][()]
    membrane_area              = r['membrane areas'][()]
    pitch          = r['unit length per voxel (um)'][()]
    num_voxel_dend = r['voxel num of dendrite not mitochondrion'][()]

padded    = (dendrite_not_mitochondrion > 0) * domains[cyt]
morpho    = lmUtils.buildAnyShape(sim, padded, domains, membrane_area)


print('Define molecules')
molecules = setMolecules(sim, cyt)

print('num_voxel_dend: ', num_voxel_dend)
volume_in_L  = num_voxel_dend * latticeSpacing * latticeSpacing * latticeSpacing * 1000
number_1umol = NA /(1e6)
number_1uM   = number_1umol * volume_in_L
print('number_per_1uM:', number_1uM)


print('Set molecules')
conc_Ca  = 100 # uM
conc_CaM = 100 # uM
conc_CB  = 30  # uM
num_Ca   = int(conc_Ca  * number_1uM)
num_CaM  = int(conc_CaM * number_1uM) 
num_CB   = int(conc_CB  * number_1uM)

print('num_Ca : ', num_Ca )
print('num_CaM: ', num_CaM)
print('num_CB : ', num_CB )


morpho.addCytosolicMolecules('Ca'  , num_Ca  , cyt) # Absolute number
morpho.addCytosolicMolecules('N0C0', num_CaM , cyt)
morpho.addCytosolicMolecules('CB'  , num_CB  , cyt)


morpho.addMembraneMolecules('PMCA' ,100) # 100 per um2

print('Set reactions')
molecules.setReactions(cyt)

print('Set up times')
sim.setTimestep(microsecond(1.0))
sim.setWriteInterval(0.1)
sim.setLatticeWriteInterval(0.1)
sim.setSimulationTime(1.0)


print('Save simulation setup.')
sim.save(filename_lm)

for molecular_name in molecules.cyt_mol:
    particleNum=sim.particleMap[molecular_name]
    print(molecular_name, str(particleNum))
    
for molecular_name in molecules.sur_mol:
    particleNum=sim.particleMap[molecular_name]
    print(molecular_name, str(particleNum))
   
# x_lsize = sim.lattice.getXSize()
# y_lsize = sim.lattice.getYSize()
# z_lsize = sim.lattice.getZSize()


# print('Run the simulation.')
# reps=1
# sim.run(filename_lm, "lm::rdme::MpdRdmeSolver", reps)

