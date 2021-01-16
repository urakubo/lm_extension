from __future__ import print_function
from __future__ import division

from pyLM import *
from pyLM.units import *
from pySTDLM import *
from pySTDLM.StandardReactionSystems import *
import h5py
import numpy as np
import os

from lib.lmUtils import buildAnyShape
from setMolecules import setMolecules



filename_morph='CA1_small.h5'
#filename_morph='CA1_ssmall.h5'

filename_lm='lms/_model.lm'

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


print('Define regions')
print('Insert dendritic geometry.')
with h5py.File(filename_morph,'r') as r:
    dendrite                   = r['dendrite'][()]
    dendrite_not_mitochondrion = r['dendrite not mitochondrion'][()]
    PSD                        = r['PSD'][()]
    membrane_area              = r['membrane areas in volume'][()]
    pitch          = r['unit length per voxel (um)'][()]
    num_voxel_dend = r['voxel num of dendrite not mitochondrion'][()]

volume    = (dendrite_not_mitochondrion > 0) * domains[cyt]

latticeSpacing=micron(pitch)
nx,ny,nz = dendrite.shape
# print("Input volume size (x,y,z): ", nx,ny,nz)
min_lattice_size = 32 
lx = np.ceil(1.0*nx/min_lattice_size)*min_lattice_size
ly = np.ceil(1.0*ny/min_lattice_size)*min_lattice_size
lz = np.ceil(1.0*nz/min_lattice_size)*min_lattice_size
# print("Lattice size (x,y,z): ", lx,ly,lz)
sim=RDME.RDMESimulation(dimensions=micron(lx*pitch,ly*pitch,lz*pitch), \
                        spacing=latticeSpacing)
morpho    = buildAnyShape(sim, volume, domains, membrane_area, PSD)

print('num_voxel_dend: ', num_voxel_dend)
volume_in_L  = num_voxel_dend * latticeSpacing * latticeSpacing * latticeSpacing * 1000
print('volume_in_fL:', volume_in_L * 1e15)
number_1umol = NA /(1e6)
V_Na       = number_1umol * volume_in_L
number_1uM = V_Na
print('number_per_1uM:', number_1uM)

print('Define molecules')
molecules = setMolecules(sim, cyt)
molecules.setDiffusion()
molecules.setReactions()


print('Set molecules')
conc_Ca  = 100 # uM
conc_CaM = 100 # uM
conc_CB  = 120 # uM
#conc_CN  = 0.5 # uM
conc_CN  = 5 # uM
num_Ca   = int(conc_Ca  * number_1uM)
num_CaM  = int(conc_CaM * number_1uM) 
num_CB   = int(conc_CB  * number_1uM)
num_CN   = int(conc_CN  * number_1uM)

print('num_Ca : ', num_Ca )
print('num_CaM: ', num_CaM)
print('num_CB : ', num_CB )
print('num_CN : ', num_CN )


morpho.addCytosolicMolecules('Ca'  , num_Ca  , cyt) # Absolute number
morpho.addCytosolicMolecules('N0C0', num_CaM , cyt)
morpho.addCytosolicMolecules('CB'  , num_CB  , cyt)
morpho.addCytosolicMolecules('CN'  , num_CN  , cyt)
#   6000 moleules per 100uM CaM and 0.1fL Spine

# 488e-12 [mol/m2] => (6.02e23) * 488e-12 / 1e12  [number/um2] 

PMCA  =  488e-12*NA/1e12 # number per um2
NCX   =  488e-12*NA/1e12 # number per um2
NMDAR =  500e-12*NA/1e12 # number per um2

morpho.addMembraneMolecules('PMCA' , PMCA)
morpho.addMembraneMolecules('NCX' , NCX)
morpho.addPSDMolecules('NR' , NMDAR)

print('Set up times')
sim.setTimestep(microsecond(3.0)) # 3.2 too big, 2.5 OK
sim.setWriteInterval(0.05)
sim.setLatticeWriteInterval(0.05)
sim.setSimulationTime(20.0)
#sim.setSimulationTime(2.0)


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

