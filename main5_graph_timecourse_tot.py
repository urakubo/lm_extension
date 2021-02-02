import numpy as np
import h5py
import os
import matplotlib.pyplot as plt

Targ = {'N0C2':['N0C2'], 'N1C0':['N1C0'], 'N1C1':['N1C1'],\
        'N1C2':['N1C2'], 'N2C0':['N2C0'], 'N2C1':['N2C1'], 'N2C2':['N2C2']};
Targ = {'CaM':['N0C2','N1C0', 'N1C1','N1C2','N2C0', 'N2C1', 'N2C2']};
Targ = {'NMDAR': ['NR_Glu','NR_O']};
Targ = {'Ca':['Ca']}
Targ = {'N0C0_CN':['N0C0_CN'], 'N0C1_CN':['N0C1_CN'], 'N0C2_CN':['N0C2_CN'],\
        'N1C0_CN':['N1C0_CN'], 'N1C1_CN':['N1C1_CN'], 'N1C2_CN':['N1C2_CN'],\
        'N2C0_CN':['N2C0_CN'], 'N2C1_CN':['N2C1_CN'], 'N2C2_CN':['N2C2_CN']}
Targ = {'CaN':['N0C0_CN', 'N0C1_CN', 'N0C2_CN', 'N1C0_CN', 'N1C1_CN', 'N1C2_CN',\
               'N2C0_CN','N2C1_CN','N2C2_CN']}

title = sorted(Targ.keys())[0]

input_lm_folder    = 'lms'
input_lm_files    = ['_prerun_result.lm',\
                     '_stimrun_00.lm',\
                     '_stimrun_01.lm',\
                     '_stimrun_02.lm',\
                     '_stimrun_03.lm',\
                     '_stimrun_04.lm',\
                     '_stimrun_05.lm',\
                     '_stimrun_06.lm',\
                     '_stimrun_07.lm',\
                     '_stimrun_08.lm',\
                     '_stimrun_09.lm',\
                     '_stimrun_10.lm',\
                     '_stimrun_11.lm',\
                     '_stimrun_12.lm',\
                     '_stimrun_13.lm',\
                     '_stimrun_14.lm',\
                     '_stimrun_15.lm',\
                     '_stimrun_16.lm',\
                     '_stimrun_17.lm',\
                     '_stimrun_18.lm',\
                     '_stimrun_19.lm',\
                     '_stimrun_20.lm',\
                     '_stimrun_21.lm',\
                     '_stimrun_22.lm',\
                     '_stimrun_23.lm',\
                     '_stimrun_24.lm',\
                     '_stimrun_25.lm']

output_figfile_prefix = "Stim"
output_figfile_dir    = 'figs'

## Define molecules and volume
cyt = 1
NA  = 6.022e23
f = h5py.File( input_lm_folder + os.sep + input_lm_files[0],'r')
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
    filename = input_lm_folder + os.sep + lmfile
    print('file :', filename)
    f = h5py.File(filename, 'r')
    tmp = f['Simulations']['0000001']['LatticeTimes'][()]
    tmp = tmp + Timepoints[-1]
    tmp = tmp.tolist()
    Timepoints.extend(tmp[1:])
    #
    tmp = f['Simulations']['0000001']['SpeciesCounts'][()]
    Numbers = np.append(Numbers, tmp[1:,:], axis=0)
    f.close()
    
#print('TimePoints: ', Timepoints)
#print('Numbers   : ', Numbers)

uMs = Numbers / NA * 1e6 / volume_in_L
Numbers = uMs


toffset = 20
t = np.array(Timepoints)-toffset

fig = plt.figure(figsize=(6,4))
ax=fig.add_subplot(111)

for k, v in Targ.items():
    conc = np.zeros_like( uMs[:,S[v[0]]-1] )
    for iv in v:
        conc += uMs[:,S[iv]-1]
    ax.plot(t, conc, label=k)

ax.set_xlim([-2.5,12.5])

ax.set_position([0.2,0.2,0.7,0.6])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.title(title)
plt.xlabel('Time (s)')
plt.ylabel('Conc (uM)')
#plt.ylabel('Number')
hans, labs = ax.get_legend_handles_labels()
ax.legend(handles=hans,labels=labs, frameon=False)
plt.savefig(output_figfile_dir+'/'+ output_figfile_prefix + '_' + title + '.pdf')
plt.savefig(output_figfile_dir+'/'+ output_figfile_prefix + '_' + title + '.png',dpi=150)
plt.show()


