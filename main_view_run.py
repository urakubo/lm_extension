import numpy as np
import h5py
from mayavi import mlab
from itertools import count
from mayavi.api import OffScreenEngine
import mcubes

input_file  = "morph_dend.lm"
output_file = "test3"

## Offscreen rendering

# mlab.options.offscreen = True


with h5py.File(input_file,'r') as f:
    morph = f['Model']['Diffusion']['LatticeSites'][()]
    particles = f['Simulations']['0000001']['Lattice']['0000000000'][:,:,:,0]


 
print('Morpho shape: ', morph.shape)
dendrite = (morph[:,:,:] == 1)
PSD    = (morph[:,:,:] == 2)

print('dendrite: ', np.sum(dendrite))
print('PSD     : ', np.sum(PSD))



vertices , faces  = mcubes.marching_cubes(dendrite, 0)
#pvertices, pfaces = mcubes.marching_cubes(PSD, 0)

print('Surface mesh was generated.')

# Plot surface mesh
mlab.figure(bgcolor=(1.0,1.0,1.0), size=(1200,1000))
mlab.view(-90, -90, 2, [0.50, 0.1277, 0 ],90)
mlab.triangular_mesh(vertices[:,0] , vertices[:,1] , vertices[:,2] , faces , color=(0.6,0.6,0.6), opacity=0.3)
#mlab.triangular_mesh(pvertices[:,0], pvertices[:,1], pvertices[:,2], pfaces, color=(1,0.5,0.5)  , opacity=0.6)



loc = np.nonzero(particles[:,:,:] > 0)

f = h5py.File(input_file,'r')
# zoom = 0.5

print('Phase0')

for i in range(0,10):
    particles   = f['Simulations']['0000001']['Lattice'][str(i).zfill(10)][:,:,:,0]
    Ca = np.nonzero(particles[:,:,:] == 1)
    PCMA = np.nonzero(particles[:,:,:] == 13)
    plot_points1 = mlab.points3d(Ca[0], Ca[1], Ca[2], color=(0,0,1), scale_factor=2.0,line_width=0.1)  
    plot_points2 = mlab.points3d(PCMA[0], PCMA[1], PCMA[2], color=(0,1,0), scale_factor=2.0,line_width=0.1)  
    mlab.view(-180.0-i*2, 90.0, 700, [120., 120., 120.])
    ffname = './pngs/'+ output_file +str(i).zfill(4)+ '.png'
    mlab.savefig(ffname)
    plot_points1.remove()
    plot_points2.remove()
print('Phase1')

# ffmpeg -f image2 -r 10 -i pngs/test3%04d.png -vcodec copy -acodec copy test3.mov -pass 2


