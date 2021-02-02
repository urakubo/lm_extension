import numpy as np
import h5py
import mcubes
from mayavi import mlab
#from itertools import count
from mayavi.api import OffScreenEngine
import trimesh

# input_model_file  = 'CA1_ssmall_model.lm'
# input_morpho_file = "CA1_ssmall.h5"

input_model_file = "lms/_model.lm"
input_morpho_file = "CA1_small.h5"
output_file = 'morph_mito.png'
output_file = 'morph_mito_psd.png'
output_file = 'morph_pmca.png'
output_file = 'morph_mito_nmda.png'


with h5py.File(input_model_file,'r') as f:
    particles = f['Model']['Diffusion']['Lattice'][:,:,:,:]

with h5py.File(input_morpho_file,'r') as f:
    memb_v    = f['membrane vertices'][()]
    memb_f    = f['membrane faces'][()]
    PSD_ids   = f['PSD ids in membrane faces'][()]
    mito_v    = f['mitochondrion vertices'][()]
    mito_f    = f['mitochondrion faces'][()]
    pitch = f['unit length per voxel (um)'][()]

memb_v = memb_v / pitch
mito_v = mito_v / pitch
print('pitch: ', pitch)

# Plot surface mesh
mlab.figure(bgcolor=(1.0,1.0,1.0), size=(400,741))
mlab.view(-180, 90, 700, [120.0, 120.0, 120.0 ])
#mlab.triangular_mesh(mito_v[:,0] , mito_v[:,1] , mito_v[:,2] , mito_f, color=(1.0,1.0,0.6), opacity=0.6)
mlab.triangular_mesh(memb_v[:,0], memb_v[:,1], memb_v[:,2], memb_f, color=(0.7,0.7,0.7)  , opacity=0.3)
mlab.triangular_mesh(memb_v[:,0], memb_v[:,1], memb_v[:,2], memb_f[PSD_ids,:], color=(1,0,0)  , opacity=0.3)

#Plot intracellular molecules
Ca   = []
NR   = []
PMCA = []
NCX  = []
for i in range(particles.shape[3]):
            Ca.extend(     np.flatnonzero(particles[:,:,:,i] == 1 ).tolist() )
            PMCA.extend( np.flatnonzero(particles[:,:,:,i] == 23).tolist() )
            NR.extend(     np.flatnonzero(particles[:,:,:,i] == 27).tolist() )
            NCX.extend(   np.flatnonzero(particles[:,:,:,i] == 25).tolist() )

Ca   = np.unravel_index(Ca, particles[:,:,:,0].shape )
NR   = np.unravel_index(NR, particles[:,:,:,0].shape )
PMCA = np.unravel_index(PMCA, particles[:,:,:,0].shape )
NCX  = np.unravel_index(NCX, particles[:,:,:,0].shape )


print('Num Ca    : ', Ca[0].shape[0])
print('Num PMCA  : ', PMCA[0].shape[0])
print('Num NMDAR : ', NR[0].shape[0])
# plot_points1 = mlab.points3d(Ca[0], Ca[1], Ca[2], color=(0,0,1), scale_factor=2.0,line_width=0.1)
plot_points2 = mlab.points3d(NR[0], NR[1], NR[2], color=(1,0,0), scale_factor=2.0,line_width=0.1)
# plot_points3 = mlab.points3d(PMCA[0], PMCA[1], PMCA[2], color=(0,1,0), scale_factor=2.0,line_width=0.1)

#
xvnum,yvnum,zvnum = particles[:,:,:,0].shape
Zoff = 0
xvnum = xvnum / 2
voxel_length = 50
length = voxel_length*pitch
text = "{0:.1f} um".format(length)
mlab.plot3d( [xvnum,xvnum],[0,voxel_length],[Zoff,Zoff],color=(0.7,0.7,0.7),tube_radius=2.5)
mlab.text3d( xvnum, voxel_length, Zoff-30, text, scale=10,color=(0.2,0.2,0.2))
# 0.02 um x 100

mlab.savefig(output_file)
mlab.show()

