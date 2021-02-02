import numpy as np
import h5py
import mcubes
from mayavi import mlab
from mayavi.api import OffScreenEngine
import trimesh

lm_files    = ['_prerun_result.lm',\
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
                     '_stimrun_25.lm',\
                     '_stimrun_26.lm',\
                     '_stimrun_27.lm',\
                     '_stimrun_28.lm',\
                     '_stimrun_29.lm',\
                     '_stimrun_30.lm',\
                     '_stimrun_31.lm']
input_lm_files = ['lms/'+ fname for fname in lm_files]
input_morpho_file = "CA1_small.h5"

#input_lm_files    = ['_prerun_result.lm','_stimrun_result.lm']

output_file = "test_"

## Offscreen rendering
mlab.options.offscreen = True

f = h5py.File(input_lm_files[0],'r')
mnames  = f['Parameters'].attrs['speciesNames'].decode().split(',')
S = {}
for i in range(len(mnames)):
    S[mnames[i]] = i+1
f.close()


with h5py.File(input_morpho_file,'r') as f:
    memb_v    = f['membrane vertices'][()]
    memb_f    = f['membrane faces'][()]
    PSD_ids   = f['PSD ids in membrane faces'][()]
    mito_v    = f['mitochondrion vertices'][()]
    mito_f    = f['mitochondrion faces'][()]
    pitch = f['unit length per voxel (um)'][()]

memb_v = memb_v / pitch
mito_v = mito_v / pitch



# Plot surface mesh
mlab.figure(bgcolor=(1.0,1.0,1.0), size=(400,740))
mlab.view(-180, 90, 700, [120.0, 120.0, 120.0 ])
mlab.triangular_mesh(mito_v[:,0] , mito_v[:,1] , mito_v[:,2] , mito_f, color=(1.0,1.0,0.6), opacity=0.6)
mlab.triangular_mesh(memb_v[:,0], memb_v[:,1], memb_v[:,2], memb_f, color=(0.7,0.7,0.7)  , opacity=0.3)
mlab.triangular_mesh(memb_v[:,0], memb_v[:,1], memb_v[:,2], memb_f[PSD_ids,:], color=(1,0,0)  , opacity=0.3)

target_molecules = ['N0C1','N0C2','N1C0', 'N1C1','N1C2','N2C0', 'N2C1', 'N2C2']
target_molecules = ['N0C1','N0C2','N1C0', 'N1C1','N1C2','N2C0', 'N2C1', 'N2C2']
#target_molecules = ['N0C2', 'N1C1','N1C2','N2C0', 'N2C1', 'N2C2']
target_molecules = ['N0C0_CN', 'N0C1_CN', 'N0C2_CN', 'N1C0_CN', 'N1C1_CN', 'N1C2_CN',\
            'N2C0_CN','N2C1_CN','N2C2_CN']
col=(0,1,0)

#target_molecules = ['Ca']
#col = (0,0,1)
toffset = 20

print('Plot.')
image_id    = 0
Time_offset = 0
for lmfile in input_lm_files:

    print('file :', lmfile)
    hfile = h5py.File(lmfile,'r')
    Timepoints = hfile['Simulations']['0000001']['LatticeTimes'][()]
    Timepoints = Timepoints + Time_offset
    Timepoints = Timepoints.tolist()
    Timepoints = Timepoints[1:]
    Time_offset = Timepoints[-1]
    
    Frames = list(hfile['Simulations']['0000001']['Lattice'].keys())
    Frames.sort()
    Frames.pop(0)
    print('Timepoints: ', Timepoints)
    
    for t, f in zip(Timepoints, Frames):
        particles = hfile['Simulations']['0000001']['Lattice'][f][:,:,:,:]
        Molecule = []
        NR   = []
        for i in range(particles.shape[3]):
            NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR']).tolist() )
            NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR_O']).tolist() )
            NR.extend( np.flatnonzero(particles[:,:,:,i] == S['NR_Glu']).tolist() )
        NR = np.unravel_index(NR, particles[:,:,:,0].shape )
        plot_points2 = mlab.points3d(NR[0], NR[1], NR[2], color=(1,0,0), scale_factor=2.0,line_width=0.1)
        for i in range(particles.shape[3]):
            for id in target_molecules:
                Molecule.extend( np.flatnonzero(particles[:,:,:,i] == S[id] ).tolist() )
        if  Molecule != []: 
            M  = np.unravel_index(Molecule, particles[:,:,:,0].shape )
            plot_points1 = mlab.points3d(M[0], M[1], M[2], color=col, scale_factor=2.0,line_width=0.1)
  
        text = mlab.text(0.1,0.05,"{0:.3f} s".format(t-toffset) ,color=(0,0,0), width=0.3)
        #print(dir(text.property.font_size))
        text.property.font_size = 10
        text.property.font_family = 'arial'
        mlab.view(-180.0-image_id*0.5, 90.0, 700, [120., 120., 120.])
        ffname = './pngs/'+ output_file +str(image_id).zfill(4)+ '.png'
        mlab.savefig(ffname)
        if  Molecule != []: 
            plot_points1.remove()
        plot_points2.remove()
        text.remove()
        # plot_points3.remove()
        image_id = image_id + 1
    hfile.close()

# ffmpeg -r 10 -i pngs/test_%04d.png -pix_fmt yuv420p output_CaN.mp4

