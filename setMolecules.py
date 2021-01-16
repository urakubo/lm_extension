from __future__ import print_function
from __future__ import division

from pyLM import *
from pyLM.units import *
import numpy as np
import os
NA = 6.022e23

class setMolecules:
    def __init__(self, sim, domain_name):
        self.sim    = sim
        self.domain_name = domain_name
        self.lx = sim.lattice.getXSize()
        self.ly = sim.lattice.getYSize()
        self.lz = sim.lattice.getZSize()

        print('lx, ly, lz:',  self.lx,  self.ly,  self.lz )
        
        # Cytosolic molecules
        self.cyt_mol  = ['Ca','N0C0','N0C1','N0C2','N1C0','N1C1','N1C2',\
                         'N2C0','N2C1','N2C2','CB','CBCa',\
                         'CN','N0C0_CN','N0C1_CN','N0C2_CN','N1C0_CN','N1C1_CN','N1C2_CN',\
                         'N2C0_CN','N2C1_CN', 'N2C2_CN']
        # Surface molecules
        self.sur_mol  = ['PMCA','PMCA_Ca', 'NCX', 'NCX_Ca','NR','NR_Glu','NR_O']

        self.sim.defineSpecies(self.cyt_mol)
        self.sim.defineSpecies(self.sur_mol)


    def setDiffusion(self):
        # Diffusion constants
        DCa          = 65  * 1e-12 ## Sould be specified in SI (233 um2/s)
        DCB          = 20  * 1e-12 ## Sould be specified in SI
        DCaM         = 7.6 * 1e-12 ## Sould be specified in SI
        DCN          = 0.12* 1e-12 #

        self.d_mol    = {'Ca'  :DCa  ,
                         'N0C0':DCaM ,
                         'N0C1':DCaM ,
                         'N0C2':DCaM ,
                         'N1C0':DCaM ,
                         'N1C1':DCaM ,
                         'N1C2':DCaM ,
                         'N2C0':DCaM ,
                         'N2C1':DCaM ,
                         'N2C2':DCaM ,
                         'CN':DCN ,                         
                         'N0C0_CN':DCN ,
                         'N0C1_CN':DCN ,
                         'N0C2_CN':DCN ,
                         'N1C0_CN':DCN ,
                         'N1C1_CN':DCN ,
                         'N1C2_CN':DCN ,
                         'N2C0_CN':DCN ,
                         'N2C1_CN':DCN ,
                         'N2C2_CN':DCN ,
                         'CB'  :DCB ,
                         'CBCa':DCB}
    
        # Warning set them after domain setting
        cyt = self.sim.modifyRegion( self.domain_name )
        for k, v in self.d_mol.items():
            cyt.setDiffusionRate(k, rate=v)


    def setReactions(self):

        Spacing = self.sim.latticeSpacing
        volume_in_L  = Spacing * Spacing * Spacing * 1000\
                       * self.lx * self.ly * self.lz
        V_Na         = NA * volume_in_L # (number per 1mol) * volume_in_L = NA*V
        on = 1e6 / V_Na # uM-1 => M-1
        # print('Spacing (m)', Spacing)
        # print('V_Na: ', V_Na)

        
        # Ca-CaM binding reactions
        kon_CB   = on * 75
        kof_CB   = 1  * 29.5

        # Ca-CaM binding reactions
        kon_TN   = on * 100
        kof_TN   = 1  * 2500
        kon_RN   = on * 150
        kof_RN   = 1  * 750
        kon_TC   = on * 4.0
        kof_TC   = 1  * 40.0
        kon_RC   = on * 10.0
        kof_RC   = 1  * 9.25

        # Ca-CaMCaN binding reactions
        kon_CN_TN   = on * 100
        kof_CN_TN   = 1  * 30  ####
        kon_CN_RN   = on * 150
        kof_CN_RN   = 1  * 12  ####
        kon_CN_TC   = on * 4.0
        kof_CN_TC   = 1  * 0.4 ####
        kon_CN_RC   = on * 10.0
        kof_CN_RC   = 1  * 0.6 ####

        # CaM-CaN binding reactions
        Mult = 1128.5
        kon_N0C0   = on  * 7.98  * 1e-8 * Mult
        kof_N0C0   = 1   * 3.19  * 1e-4 * Mult
        kon_N1C0   = on  * 6.65  * 1e-6 * Mult
        kof_N1C0   = 1   * 3.19  * 1e-4 * Mult
        kon_N0C1   = on  * 1.334 * 1e-5 * Mult
        kof_N0C1   = 1   * 3.19  * 1e-4 * Mult
        kon_N1C1   = on  * 6.65  * 1e-4 * Mult
        kof_N1C1   = 1   * 3.19  * 1e-4 * Mult
        kon_N2C0   = on  * 4.16  * 1e-4 * Mult
        kof_N2C0   = 1   * 3.19  * 1e-4 * Mult
        kon_N0C2   = on  * 1.23  * 1e-4 * Mult
        kof_N0C2   = 1   * 3.19  * 1e-4 * Mult
        kon_N2C1   = on  * 4.16  * 1e-2 * Mult
        kof_N2C1   = 1   * 3.19  * 1e-4 * Mult
        kon_N1C2   = on  * 1.02  * 1e-2 * Mult
        kof_N1C2   = 1   * 3.19  * 1e-4 * Mult
        kon_N2C2   = on  * 0.64  * Mult
        kof_N2C2   = 1   * 3.19  * 1e-4 * Mult

        # Ca pumps binding reactions
        kon_PMCA   = on  * 150
        kof_PMCA   = 1   * 15
        kct_PMCA   = 1   * 12
        klk_PMCA   = 1   * 4.3

        kon_NCX    = on  * 300
        kof_NCX    = 1   * 300
        kct_NCX    = 1   * 600
        klk_NCX    = 1   * 19.4

        # Reactant
        Ca='Ca'
        N0C0 = 'N0C0'
        N0C1 = 'N0C1'
        N0C2 = 'N0C2'
        N1C0 = 'N1C0'
        N1C1 = 'N1C1'
        N1C2 = 'N1C2'
        N2C0 = 'N2C0'
        N2C1 = 'N2C1'
        N2C2 = 'N2C2'
        CB   = 'CB'
        CBCa = 'CBCa'
        #
        CN      = 'CN'
        N0C0_CN = 'N0C0_CN'
        N0C1_CN = 'N0C1_CN'
        N0C2_CN = 'N0C2_CN'
        N1C0_CN = 'N1C0_CN'
        N1C1_CN = 'N1C1_CN'
        N1C2_CN = 'N1C2_CN'
        N2C0_CN = 'N2C0_CN'
        N2C1_CN = 'N2C1_CN'
        N2C2_CN = 'N2C2_CN'
        #
        PMCA    = 'PMCA'
        PMCA_Ca = 'PMCA_Ca'
        NCX     = 'NCX'
        NCX_Ca  = 'NCX_Ca'
        NR      = 'NR'
        NR_Glu  = 'NR_Glu'
        NR_O    = 'NR_O'

        # CB binding
        cyt = self.sim.modifyRegion( self.domain_name )
        cyt.addReaction(reactant=(Ca,CB), product=CBCa   , rate=kon_CB)
        cyt.addReaction(reactant= CBCa  , product=(Ca,CB), rate=kof_CB)

        # Ca-CaM binding
        cyt.addReaction(reactant=(Ca,N0C0), product= N1C0, rate=kon_TN)
        cyt.addReaction(reactant= N1C0, product=(Ca,N0C0), rate=kof_TN)
        cyt.addReaction(reactant=(Ca,N1C0), product= N2C0, rate=kon_RN)
        cyt.addReaction(reactant= N2C0, product=(Ca,N1C0), rate=kof_RN)
        
        cyt.addReaction(reactant=(Ca,N0C1), product= N1C1, rate=kon_TN)
        cyt.addReaction(reactant= N1C1, product=(Ca,N0C1), rate=kof_TN)
        cyt.addReaction(reactant=(Ca,N1C1), product= N2C1, rate=kon_RN)
        cyt.addReaction(reactant= N2C1, product=(Ca,N1C1), rate=kof_RN)

        cyt.addReaction(reactant=(Ca,N0C2), product= N1C2, rate=kon_TN)
        cyt.addReaction(reactant= N1C2, product=(Ca,N0C2), rate=kof_TN)
        cyt.addReaction(reactant=(Ca,N1C2), product= N2C2, rate=kon_RN)
        cyt.addReaction(reactant= N2C2, product=(Ca,N1C2), rate=kof_RN)

        cyt.addReaction(reactant=(Ca,N0C0), product= N0C1, rate=kon_TC)
        cyt.addReaction(reactant= N0C1, product=(Ca,N0C0), rate=kof_TC)
        cyt.addReaction(reactant=(Ca,N0C1), product= N0C2, rate=kon_RC)
        cyt.addReaction(reactant= N0C2, product=(Ca,N0C1), rate=kof_RC)

        cyt.addReaction(reactant=(Ca,N1C0), product= N1C1, rate=kon_TC)
        cyt.addReaction(reactant= N1C1, product=(Ca,N1C0), rate=kof_TC)
        cyt.addReaction(reactant=(Ca,N1C1), product= N1C2, rate=kon_RC)
        cyt.addReaction(reactant= N1C2, product=(Ca,N1C1), rate=kof_RC)

        cyt.addReaction(reactant=(Ca,N2C0), product= N2C1, rate=kon_TC)
        cyt.addReaction(reactant= N2C1, product=(Ca,N2C0), rate=kof_TC)
        cyt.addReaction(reactant=(Ca,N2C1), product= N2C2, rate=kon_RC)
        cyt.addReaction(reactant= N2C2, product=(Ca,N2C1), rate=kof_RC)
 
        # Ca-CaMCaN binding
        cyt.addReaction(reactant=(Ca,N0C0_CN), product=N1C0_CN, rate=kon_CN_TN)
        cyt.addReaction(reactant=N1C0_CN, product=(Ca,N0C0_CN), rate=kof_CN_TN)
        cyt.addReaction(reactant=(Ca,N1C0_CN), product=N2C0_CN, rate=kon_CN_RN)
        cyt.addReaction(reactant=N2C0_CN, product=(Ca,N1C0_CN), rate=kof_CN_RN)

        cyt.addReaction(reactant=(Ca,N0C1_CN), product=N1C1_CN, rate=kon_CN_TN)
        cyt.addReaction(reactant=N1C1_CN, product=(Ca,N0C1_CN), rate=kof_CN_TN)
        cyt.addReaction(reactant=(Ca,N1C1_CN), product=N2C1_CN, rate=kon_CN_RN)
        cyt.addReaction(reactant=N2C1_CN, product=(Ca,N1C1_CN), rate=kof_CN_RN)

        cyt.addReaction(reactant=(Ca,N0C2_CN), product=N1C2_CN, rate=kon_CN_TN)
        cyt.addReaction(reactant=N1C2_CN, product=(Ca,N0C2_CN), rate=kof_CN_TN)
        cyt.addReaction(reactant=(Ca,N1C2_CN), product=N2C2_CN, rate=kon_CN_RN)
        cyt.addReaction(reactant=N2C2_CN, product=(Ca,N1C2_CN), rate=kof_CN_RN)

        cyt.addReaction(reactant=(Ca,N0C0_CN), product=N0C1_CN, rate=kon_CN_TC)
        cyt.addReaction(reactant=N0C1_CN, product=(Ca,N0C0_CN), rate=kof_CN_TC)
        cyt.addReaction(reactant=(Ca,N0C1_CN), product=N0C2_CN, rate=kon_CN_RC)
        cyt.addReaction(reactant=N0C2_CN, product=(Ca,N0C1_CN), rate=kof_CN_RC)

        cyt.addReaction(reactant=(Ca,N1C0_CN), product=N1C1_CN, rate=kon_CN_TC)
        cyt.addReaction(reactant=N1C1_CN, product=(Ca,N1C0_CN), rate=kof_CN_TC)
        cyt.addReaction(reactant=(Ca,N1C1_CN), product=N1C2_CN, rate=kon_CN_RC)
        cyt.addReaction(reactant=N1C2_CN, product=(Ca,N1C1_CN), rate=kof_CN_RC)

        cyt.addReaction(reactant=(Ca,N2C0_CN), product=N2C1_CN, rate=kon_CN_TC)
        cyt.addReaction(reactant=N2C1_CN, product=(Ca,N2C0_CN), rate=kof_CN_TC)
        cyt.addReaction(reactant=(Ca,N2C1_CN), product=N2C2_CN, rate=kon_CN_RC)
        cyt.addReaction(reactant=N2C2_CN, product=(Ca,N2C1_CN), rate=kof_CN_RC)

        # CaCaM-CaN binding
        cyt.addReaction(reactant=(N0C0,CN), product=N0C0_CN,  rate=kon_N0C0)
        cyt.addReaction(reactant=N0C0_CN, product=(N0C0,CN),  rate=kof_N0C0)
        cyt.addReaction(reactant=(N1C0,CN), product=N1C0_CN,  rate=kon_N1C0)
        cyt.addReaction(reactant=N1C0_CN, product=(N1C0,CN),  rate=kof_N1C0)
        cyt.addReaction(reactant=(N2C0,CN), product=N2C0_CN,  rate=kon_N2C0)
        cyt.addReaction(reactant=N2C0_CN, product=(N2C0,CN),  rate=kof_N2C0)

        cyt.addReaction(reactant=(N0C1,CN), product=N0C1_CN,  rate=kon_N0C1)
        cyt.addReaction(reactant=N0C1_CN, product=(N0C1,CN),  rate=kof_N0C1)
        cyt.addReaction(reactant=(N1C1,CN), product=N1C1_CN,  rate=kon_N1C1)
        cyt.addReaction(reactant=N1C1_CN, product=(N1C1,CN),  rate=kof_N1C1)
        cyt.addReaction(reactant=(N2C1,CN), product=N2C1_CN,  rate=kon_N2C1)
        cyt.addReaction(reactant=N2C1_CN, product=(N2C1,CN),  rate=kof_N2C1)

        cyt.addReaction(reactant=(N0C2,CN), product=N0C2_CN,  rate=kon_N0C2)
        cyt.addReaction(reactant=N0C2_CN, product=(N0C2,CN),  rate=kof_N0C2)
        cyt.addReaction(reactant=(N1C2,CN), product=N1C2_CN,  rate=kon_N1C2)
        cyt.addReaction(reactant=N1C2_CN, product=(N1C2,CN),  rate=kof_N1C2)
        cyt.addReaction(reactant=(N2C2,CN), product=N2C2_CN,  rate=kon_N2C2)
        cyt.addReaction(reactant=N2C2_CN, product=(N2C2,CN),  rate=kof_N2C2)



        ##
        ## Ca pump and NMDAR
        ##
        flux = 4000
        cyt.addReaction(reactant=(PMCA,Ca), product= PMCA_Ca , rate=kon_PMCA)
        cyt.addReaction(reactant=PMCA_Ca  , product= PMCA    , rate=kct_PMCA)
        cyt.addReaction(reactant=PMCA     , product=(PMCA,Ca), rate=klk_PMCA)

        cyt.addReaction(reactant=(NCX,Ca), product=NCX_Ca , rate=kon_NCX)
        cyt.addReaction(reactant=NCX_Ca  , product=(NCX,Ca) , rate=kof_NCX)
        cyt.addReaction(reactant=NCX_Ca  , product=NCX      , rate=kct_NCX)
        cyt.addReaction(reactant=NCX     , product=(NCX,Ca) , rate=klk_NCX)

        cyt.addReaction(reactant=NR_Glu, product=NR        , rate=50  )
        cyt.addReaction(reactant=NR_Glu, product=NR_O      , rate=200 )
        cyt.addReaction(reactant=NR_O  , product=NR_Glu    , rate=50  )
        cyt.addReaction(reactant=NR_O  , product=(NR_O,Ca) , rate=flux)
