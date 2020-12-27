import numpy as np
import h5py
import mcubes
from mayavi import mlab
from itertools import count
from mayavi.api import OffScreenEngine


input_file  = 'morph_dend.lm'
output_file = 'morph_dend.png'

# input_file  = 'morph_capsid.lm'
# output_file  = 'morph_capsid.png'

with h5py.File(input_file,'r') as f:
    morph = f['Model']['Diffusion']['LatticeSites'][()]
    particles = f['Model']['Diffusion']['Lattice'][:,:,:,0]
 
print('Morpho shape: ', morph.shape)
dendrite = (morph[:,:,:] == 1)
PSD    = (morph[:,:,:] == 2)

print('dendrite: ', np.sum(dendrite))
print('PSD     : ', np.sum(PSD))

vertices , faces  = mcubes.marching_cubes(dendrite, 0)
#pvertices, pfaces = mcubes.marching_cubes(PSD, 0)

print('Surface mesh was generated.')

# Plot surface mesh
mlab.figure(bgcolor=(1.0,1.0,1.0), size=(800,600))
mlab.view(-90, -90, 2, [0.50, 0.1277, 0 ],90)
mlab.triangular_mesh(vertices[:,0] , vertices[:,1] , vertices[:,2] , faces , color=(0.6,0.6,0.6), opacity=0.3)
#mlab.triangular_mesh(pvertices[:,0], pvertices[:,1], pvertices[:,2], pfaces, color=(1,0.5,0.5)  , opacity=0.6)

# Plot intracellular molecules

print('particles.shape: ', particles.shape)
Ca = np.nonzero(particles[:,:,:] == 1)
PCMA = np.nonzero(particles[:,:,:] == 13)
print('Num Ca   : ', Ca[0].shape[0])
print('Num PCMA : ', PCMA[0].shape[0])
#print('Num tot_particle     : ', np.sum(tot_particles))

plot_points = mlab.points3d(Ca[0], Ca[1], Ca[2], color=(0,0,1), scale_factor=2.0,line_width=0.1)  
plot_points = mlab.points3d(PCMA[0], PCMA[1], PCMA[2], color=(0,1,0), scale_factor=2.0,line_width=0.1)  

mlab.savefig(output_file)
# mlab.show()

