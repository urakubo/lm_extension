from __future__ import print_function
from __future__ import division

import h5py
import numpy as np
import os
import subprocess as s

filename_lm='morph_dend.lm'
rep = '1'
# 0-2: repeats 0,1,2

# -sp: RDME
# -sl lm::rime::MpdRdmeSolver

com = ['lm','-r', rep, '-sp', '-sl','lm::rdme::MpdRdmeSolver','-f', filename_lm]
print(' '.join(com))
s.call(com)
