from __future__ import print_function
from __future__ import division

import h5py
import numpy as np
import os, sys, shutil
import subprocess as s

filename_lm='CA1_small_model.lm'
filename_prerun  = 'CA1_small_run_pre.lm'

# filename_lm = 'CA1_ssmall_model.lm'
# filename_prerun = 'CA1_ssmall_run_pre.lm'

if os.path.isfile(filename_prerun):
    print('Prerun file exists : ', filename_prerun)
    sys.exit()

print('Create prerun file : ', filename_prerun)
shutil.copy(filename_lm, filename_prerun)


rep = '1'
# 0-2: repeats 0,1,2

# -sp: RDME
# -sl lm::rime::MpdRdmeSolver

com = ['lm','-r', rep, '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_prerun]
#com = ['mpirun','-np','2','lm','-cr','1/2','-gr','1/2','-r', rep, '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_lm]
#com = ['lm','-cr','2','-g','0,1','-gr','2','-r', rep, '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_lm]
print(' '.join(com))
s.call(com)
