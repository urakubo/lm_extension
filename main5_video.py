import numpy as np
import h5py
import mcubes
from mayavi import mlab
from mayavi.api import OffScreenEngine
import trimesh

input_morpho_file = "CA1_small.h5"
input_lm_files    = ['CA1_small_run_pre.lm','CA1_small_run_stim.lm']

input_morpho_file = "CA1_small.h5"
input_lm_files    = ['_prerun_result.lm','_stimrun_result.lm']

output_file = "test_"

## Offscreen rendering
# mlab.options.offscreen = True

    
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
mlab.figure(bgcolor=(1.0,1.0,1.0), size=(400,741))
mlab.view(-180, 90, 700, [120.0, 120.0, 120.0 ])
mlab.triangular_mesh(mito_v[:,0] , mito_v[:,1] , mito_v[:,2] , mito_f, color=(1.0,1.0,0.6), opacity=0.6)
mlab.triangular_mesh(memb_v[:,0], memb_v[:,1], memb_v[:,2], memb_f, color=(0.7,0.7,0.7)  , opacity=0.3)
mlab.triangular_mesh(memb_v[:,0], memb_v[:,1], memb_v[:,2], memb_f[PSD_ids,:], color=(1,0,0)  , opacity=0.3)


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
        particles   = hfile['Simulations']['0000001']['Lattice'][f][:,:,:,:]
        Ca   = []
        NR   = []
        NR_Glu = []
        NR_O = []
        for i in range(particles.shape[3]):
            Ca.extend(     np.flatnonzero(particles[:,:,:,i] == 1 ).tolist() )
            NR.extend(     np.flatnonzero(particles[:,:,:,i] == 27).tolist() )
            NR_Glu.extend( np.flatnonzero(particles[:,:,:,i] == 28).tolist() )
            NR_O.extend(   np.flatnonzero(particles[:,:,:,i] == 29).tolist() )
            # PMCA    = np.nonzero(particles[:,:,:] == 13)
            # PMCA_Ca = np.nonzero(particles[:,:,:] == 14)
        Ca = np.unravel_index(Ca, particles[:,:,:,0].shape )
        NR = np.unravel_index(NR, particles[:,:,:,0].shape )
        plot_points1 = mlab.points3d(Ca[0], Ca[1], Ca[2], color=(0,0,1), scale_factor=2.0,line_width=0.1)  
        plot_points2 = mlab.points3d(NR[0], NR[1], NR[2], color=(1,0,0), scale_factor=2.0,line_width=0.1)  
        # plot_points2 = mlab.points3d(PMCA[0], PMCA[1], PMCA[2], color=(0,1,0), scale_factor=2.0,line_width=0.1)  
        # plot_points3 = mlab.points3d(PMCA_Ca[0], PMCA_Ca[1], PMCA_Ca[2], color=(0,1,0.5), scale_factor=2.0,line_width=0.1)  
        text = mlab.text(0.1,0.05,"{0:.2f} s".format(t) ,color=(0,0,0), width=0.3)
        #print(dir(text.property.font_size))
        text.property.font_size = 10
        text.property.font_family = 'arial'
        mlab.view(-180.0-image_id, 90.0, 700, [120., 120., 120.])
        ffname = './pngs/'+ output_file +str(image_id).zfill(4)+ '.png'
        mlab.savefig(ffname)
        plot_points1.remove()
        plot_points2.remove()
        text.remove()
        # plot_points3.remove()
        image_id = image_id + 1
    hfile.close()

# ffmpeg -r 10 -i pngs/test_%04d.png -pix_fmt yuv420p output.mp4



### ffmpeg -f image2 -r 10 -i pngs/test3%04d.png -profile:v main test3.mp4 -pass 2


