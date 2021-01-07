import numpy as np
import h5py

input_lm_files    = ['CA1_small_run_pre.lm']
input_lm_files    = ['CA1_ssmall_run_pre.lm']
#output_file = "test3"
fig_dir = 'figs'
## Offscreen rendering
# mlab.options.offscreen = True
S = {'Ca':1, 'N0C0':2, 'N0C1':3, 'N0C2':4, 'N1C0':5,\
           'N1C1':6,'N1C2':7, 'N2C0':8, 'N2C1':9, 'N2C2':10,\
           'CB':11, 'CBCa':12, 'CN':13, 'N0C0_CN':14, 'N0C1_CN':15,\
           'N0C2_CN':16, 'N1C0_CN':17, 'N1C1_CN':18, 'N1C2_CN':19,\
           'N2C0_CN':20, 'N2C1_CN':21, 'N2C2_CN':22, 'PMCA':23,\
           'PMCA_Ca':24, 'NCX':25, 'NCX_Ca':26, 'NR':27,
           'NR_Glu':28, 'NR_O': 29}

Timepoints = [0]
Numbers    = np.zeros((1,len(S)),int)
for i, lmfile in enumerate(input_lm_files):
    print('file :', lmfile)
    f = h5py.File(lmfile,'r')
    tmp = f['Simulations']['0000001']['LatticeTimes'][()]
    tmp = tmp + Timepoints[-1]
    tmp = tmp.tolist()
    Timepoints.extend(tmp[1:])
    #
    tmp = f['Simulations']['0000001']['SpeciesCounts'][()]
    Numbers = np.append(Numbers, tmp[1:,:], axis=0)
    if i == 0:
        print(S.keys())
        print('Initial numbers: ', tmp[0,:])
    #
    f.close()
    
#print('TimePoints: ', Timepoints)
#print('Numbers   : ', Numbers)

import matplotlib.pyplot as plt

Targ = 'Ca'
Targ = 'NMDAR'
Targ = 'Ca'
Targ = 'Ca'
Targ = 'CaM'
Targ = 'CaN2'
Targ = 'CaN'
Targ = 'Ca'
Targ = 'CaM'

fig = plt.figure(figsize=(6,4))
ax=fig.add_subplot(111)

if Targ == 'Ca':
    ax.plot(Timepoints, Numbers[:,S[Targ]-1], label=Targ)

elif Targ == 'CaN':
    CaN = Numbers[:,S['N0C0_CN']-1] + Numbers[:,S['N0C1_CN']-1] + Numbers[:,S['N0C2_CN']-1]\
          + Numbers[:,S['N1C0_CN']-1] + Numbers[:,S['N1C1_CN']-1] + Numbers[:,S['N1C2_CN']-1]\
          + Numbers[:,S['N2C0_CN']-1] + Numbers[:,S['N2C1_CN']-1] + Numbers[:,S['N2C2_CN']-1]
    ax.plot(Timepoints, CaN, label=Targ)

elif Targ == 'CaN2':
    CaNs = ['CN','N0C0_CN', 'N0C1_CN', 'N0C2_CN', 'N1C0_CN', 'N1C1_CN', 'N1C2_CN',\
            'N2C0_CN','N2C1_CN','N2C2_CN']
    for name in CaNs:
        ax.plot(Timepoints, Numbers[:,S[name]-1], label=name )

elif Targ == 'CaM':
    CaMs = ['N0C0','N0C1' ];
    CaMs = ['N0C2','N1C0', 'N1C1','N1C2','N2C0', 'N2C1', 'N2C2'];
    for cam in CaMs:
        ax.plot(Timepoints, Numbers[:,S[cam]-1], label=cam )

elif Targ == 'NMDAR':
    NRs = ['NR_Glu','NR_O'];
    for name in NRs:
        ax.plot(Timepoints, Numbers[:,S[name]-1], label=name )

ax.set_position([0.2,0.2,0.7,0.6])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.title(Targ)
plt.xlabel('Time (s)')
plt.ylabel('Number')
hans, labs = ax.get_legend_handles_labels()
ax.legend(handles=hans,labels=labs, frameon=False)
plt.savefig(fig_dir+'/Prof_'+Targ+'.pdf')
plt.savefig(fig_dir+'/Prof_'+Targ+'.png',dpi=150)
plt.show()

