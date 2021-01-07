from __future__ import print_function
from __future__ import division

from pyLM import *
from pyLM.units import *
from pySTDLM import *
from pySTDLM.PostProcessing import *
from pySTDLM.StandardReactionSystems import *
# from pySTDLM.StandardCells import *
import h5py
import numpy as np
import os

from lib.lmUtils import buildAnyShape
from setMolecules import setMolecules


filename_lm='CA1_small_sim.lm'
filename_morph='CA1_small.h5'

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
