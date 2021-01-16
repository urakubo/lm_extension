from __future__ import print_function
from __future__ import division

import h5py
import numpy as np
import os, sys, shutil
import subprocess as s

filename_lm      = 'lms/_model.lm'
filename_prerun  = ['lms/_prerun_result.lm']
filename_stimrun_prefix = 'lms/_stimrun_'

# filename_lm      = 'CA1_ssmall_model.lm'
# filename_prerun  = 'CA1_ssmall_run_pre.lm'
# filename_stimrun = 'CA1_ssmall_run_stim.lm'

if not os.path.isfile(filename_lm):
    print('No model file         : ', filename_lm)
    sys.exit()
if not os.path.isfile(filename_prerun[0]):
    print('No prerun file     : ', filename_prerun)

##
## Activate NMDAR (NR)
##
## NR == 27, NR_Glu == 28 NR_O == 29
##
## Lattice: uint8
## InitialSpeciesCounts: dtype=uint32
##
    
f = h5py.File(filename_prerun[-1],'r')
mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')
S = {}
for i in range(len(mnames)):
    S[mnames[i]] = i+1
#test = f['Parameters'].attrs['maxTime']
f.close()

#print(type(test))
#sys.exit()
period = np.array('1.0  ',dtype=str)

for i in range(10):

    with h5py.File(filename_prerun[-1],'r') as f:
        TimePoints = f['Simulations']['0000001']['Lattice'].keys()
        TimePoints.sort()
        print('Time ID: ', TimePoints[-1])
        Lattice = f['Simulations']['0000001']['Lattice'][TimePoints[-1]][()]
        SpeciesCount = f['Simulations']['0000001']['SpeciesCounts'][-1,:]

        Lattice[Lattice == S['NR']] = S['NR_Glu']
        SpeciesCount[S['NR_Glu']-1] = SpeciesCount[S['NR_Glu']-1]\
                                      + SpeciesCount[S['NR']-1]
        SpeciesCount[S['NR']-1]     = 0
        id = '%02d' % i
        filename = filename_stimrun_prefix + id +'.lm'
        filename_prerun.append(filename)
        print('Create stimrun file: ', filename)
        shutil.copy(filename_lm, filename)

        with h5py.File(filename,'a') as g:
            g['Model']['Diffusion']['Lattice'][()] = Lattice
            g['Model']['Reaction']['InitialSpeciesCounts'][()] = SpeciesCount
            g['Parameters'].attrs['maxTime'] = period
            g['Model']['Reaction']['ReactionRateConstants'][78,0] = 8000.0
        com = ['lm','-r', '1', '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename]
        print(' '.join(com))
        s.call(com)
