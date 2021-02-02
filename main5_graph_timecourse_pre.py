import numpy as np
import h5py

input_lm_files    = ['_prerun_result.lm']
output_figfile_prefix = "Pre"
output_figfile_dir    = 'figs'
## Offscreen rendering
# mlab.options.offscreen = True


## Define molecules and volume
cyt = 1
NA  = 6.022e23
f = h5py.File(input_lm_files[0],'r')
data = f['Model']['Diffusion']['LatticeSites'][()]
num_voxels = np.count_nonzero(data == cyt)
Spacing = f['Model']['Diffusion'].attrs['latticeSpacing']
volume_in_L  = num_voxels * Spacing * Spacing * Spacing * 1000
mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')
S = {}
for i in range(len(mnames)):
    S[mnames[i]] = i+1
f.close()
##

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
print('Num voxels  : ', num_voxels)
print('Spacing     : ', Spacing)
print('Volume in fL: ', volume_in_L * 1e15)
uMs = Numbers / NA * 1e6 / volume_in_L
Numbers = uMs

import matplotlib.pyplot as plt

Targ = 'Ca'
Targ = 'NMDAR'
Targ = 'Ca'
Targ = 'Ca'
Targ = 'CaM'
Targ = 'CaN2'
Targ = 'CaN'
Targ = 'CaN'
#Targ = 'CaM'

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
    #CaMs = ['N0C2','N1C0', 'N1C1','N1C2','N2C0', 'N2C1', 'N2C2'];
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
#plt.ylabel('Number')
plt.ylabel('(uM)')
hans, labs = ax.get_legend_handles_labels()
ax.legend(handles=hans,labels=labs, frameon=False)
plt.savefig(output_figfile_dir+'/'+ output_figfile_prefix + '_' + Targ + '.pdf')
plt.savefig(output_figfile_dir+'/'+ output_figfile_prefix + '_' + Targ + '.png',dpi=150)
plt.show()

