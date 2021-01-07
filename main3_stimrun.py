from __future__ import print_function
from __future__ import division

import h5py
import numpy as np
import os, sys, shutil
import subprocess as s

filename_lm      = 'CA1_small_model.lm'
filename_prerun  = 'CA1_small_run_pre.lm'
filename_stimrun = 'CA1_small_run_stim.lm'

# filename_lm      = 'CA1_ssmall_model.lm'
# filename_prerun  = 'CA1_ssmall_run_pre.lm'
# filename_stimrun = 'CA1_ssmall_run_stim.lm'

if not os.path.isfile(filename_prerun):
    print('No lm file         : ', filename_lm)
    sys.exit()
if not os.path.isfile(filename_prerun):
    print('No prerun file     : ', filename_prerun)
    
if os.path.isfile(filename_stimrun):
    print('Stimrun file exists: ', filename_stimrun)
    sys.exit()
    # os.remove(filename_stimrun)

##
## Activate NMDAR (NR)
##
## NR == 27, NR_Glu == 28 NR_O == 29
##
## Lattice: uint8
## InitialSpeciesCounts: dtype=uint32
##
    

with h5py.File(filename_prerun,'r') as f:
    TimePoints = f['Simulations']['0000001']['Lattice'].keys()
    TimePoints.sort()
    print('Time ID: ', TimePoints[-1])
    Lattice = f['Simulations']['0000001']['Lattice'][TimePoints[-1]][()]
    SpeciesCount = f['Simulations']['0000001']['SpeciesCounts'][-1,:]

NR     = np.array(27).astype(np.uint8)
NR_Glu = np.array(28).astype(np.uint8)
NR_NUM = np.count_nonzero(Lattice == NR)

print('NR_NUM in Lattie : ', NR_NUM)
print('SpeciesCount[NR] : ', SpeciesCount[NR-1])
# print('SpeciesCount     : ', SpeciesCount)

Lattice[Lattice == NR] = NR_Glu
SpeciesCount[NR_Glu-1] = SpeciesCount[NR_Glu-1] + SpeciesCount[NR-1]
SpeciesCount[NR-1]     = 0
# sys.exit()

print('Create stimrun file: ', filename_stimrun)
shutil.copy(filename_lm, filename_stimrun)

with h5py.File(filename_stimrun,'a') as g:
    g['Model']['Diffusion']['Lattice'][()] = Lattice
    g['Model']['Reaction']['InitialSpeciesCounts'][()] = SpeciesCount


rep = '1'
# 0-2: repeats 0,1,2
# -sp: RDME
# -sl lm::rime::MpdRdmeSolver

com = ['lm','-r', rep, '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_stimrun]
#com = ['mpirun','-np','2','lm','-cr','1/2','-gr','1/2','-r', rep, '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_lm]
#com = ['lm','-cr','2','-g','0,1','-gr','2','-r', rep, '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_lm]
print(' '.join(com))
s.call(com)
