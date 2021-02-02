import numpy as np
import h5py
import os
import pickle



Targs = ['N0C0_CN', 'N0C1_CN', 'N0C2_CN', 'N1C0_CN', 'N1C1_CN', 'N1C2_CN',\
            'N2C0_CN','N2C1_CN','N2C2_CN']
#Targs = ['Ca']
print('Targs[0]: ', Targs[0])

input_lm_files    = ['_prerun_result.lm','_stimrun_result.lm']
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

input_label_file = "lm_annot/labels.hdf5"
output_figfile_prefix = "Stim_Spine_"
output_figfile_dir    = 'figs'

## Offscreen rendering
# mlab.options.offscreen = True


## Decode molecular names and volume
cyt = 1
NA  = 6.022e23
f = h5py.File( input_lm_folder + os.sep + input_lm_files[0],'r')

data = f['Model']['Diffusion']['LatticeSites'][()]
num_voxels = np.count_nonzero(data == cyt)
Spacing = f['Model']['Diffusion'].attrs['latticeSpacing']
volume_in_L  = num_voxels * Spacing * Spacing * Spacing * 1000

## Decode molecular names
mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')
S = {}
for i in range(len(mnames)):
    S[mnames[i]] = i+1
f.close()


## Load spine labels    
with h5py.File(input_label_file,'r') as f:
    labels = f['dendrite'][()]
ids_spine, nums_spine_voxels = np.unique(labels, return_counts=True)
ids_spine         = ids_spine[1:-2]
nums_spine_voxels = nums_spine_voxels[1:-2]

# ids_spine         = ids_spine[::4]
# nums_spine_voxels = nums_spine_voxels[::4]

vols_spine_in_L   =  nums_spine_voxels * Spacing * Spacing * Spacing * 1000
print('ids_spine: ', ids_spine)

_labels = labels
nx,ny,nz = _labels.shape
min_lattice_size = 32 
lx = np.ceil(1.0*nx/min_lattice_size)*min_lattice_size
ly = np.ceil(1.0*ny/min_lattice_size)*min_lattice_size
lz = np.ceil(1.0*nz/min_lattice_size)*min_lattice_size
lx = lx.astype(np.int)
ly = ly.astype(np.int)
lz = lz.astype(np.int)
labels = np.zeros((lx,ly,lz), dtype=np.int)
labels[:nx,:ny,:nz] = _labels
labels = labels.flatten()


with open('lm_annot/list.pickle','rb') as f:
    list = pickle.load(f)
cols = []
for id_spine in ids_spine:
    c = [x for x in list['list'] if x['id'] == id_spine]
    r = c[0]['r']/256.0
    g = c[0]['g']/256.0
    b = c[0]['b']/256.0
    cols.append((r,g,b))


## Obtain timepoints
Timepoints = [0]
for i, lmfile in enumerate(input_lm_files):
    filename = input_lm_folder + os.sep + lmfile
    #print('file :', filename)
    hfile = h5py.File(filename, 'r')
    tmp = hfile['Simulations']['0000001']['LatticeTimes'][()]
    hfile.close
    tmp = tmp + Timepoints[-1]
    tmp = tmp.tolist()
    Timepoints.extend(tmp[1:])

Timepoints = Timepoints[1:]

# print('Timepoints: ',Timepoints)
    
## Obtain molecular concs of spines at each timepoint
num_molecules_spine = []
for i, lmfile in enumerate(input_lm_files):
    filename = input_lm_folder + os.sep + lmfile
    print('file :', filename)
    hfile = h5py.File(filename, 'r')
    Frames = [key for key in hfile['Simulations']['0000001']['Lattice'].keys()]
    Frames.sort()
    Frames.pop(0)
    for f in Frames:
        particles = hfile['Simulations']['0000001']['Lattice'][f][:,:,:,:]
        num_molecules_spine_time_i = []
        for id_spine in ids_spine:
            tmp_num = 0
            targ_spine_label = (labels == id_spine)
            for j in range(particles.shape[3]):
                p = particles[:,:,:,j].flatten()
                pp = p[targ_spine_label]
                for Targ in Targs:
                    tmp_num += np.count_nonzero( pp == S[Targ] )
            num_molecules_spine_time_i.append(tmp_num)
        print('num_molecules_spine_time_i: ', num_molecules_spine_time_i)
        num_molecules_spine.append(num_molecules_spine_time_i)
    hfile.close

num_molecules_spine = np.array(num_molecules_spine)
print()

uMs = num_molecules_spine / NA * 1e6 / vols_spine_in_L
toffset = 20
t = np.array(Timepoints[:])-toffset
np.savez('num_molecules_spine_'+Targs[0]+'.npz', t=t, num_molecules_spine=num_molecules_spine)

# tmp = np.load('num_molecules_spine.npz')
# num_molecules_spine = tmp['num_molecules_spine']

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6,4))
ax=fig.add_subplot(111)
#for i, id_spine in enumerate(ids_spine):
#    ax.plot(Timepoints, num_molecules_spine[:,i], label=str(id_spine), color=cols[i] )
for i, id_spine in enumerate(ids_spine):
    ax.plot(t, uMs[:,i], label=str(id_spine), color=np.random.rand(3,) )

ax.set_position([0.2,0.2,0.7,0.6])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim([-2.5,12.5])
plt.title(Targ)
plt.xlabel('Time (s)')
plt.ylabel('Conc (uM)')

#plt.ylabel('Number')
hans, labs = ax.get_legend_handles_labels()
ax.legend(handles=hans,labels=labs, frameon=False)
plt.savefig(output_figfile_dir+'/'+ output_figfile_prefix + Targs[0] + '.pdf')
plt.savefig(output_figfile_dir+'/'+ output_figfile_prefix + Targs[0] + '.png',dpi=150)
plt.show()
