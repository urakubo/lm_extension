import numpy as np
import h5py
import mcubes
from mayavi import mlab
#from itertools import count
from mayavi.api import OffScreenEngine
import trimesh
import pickle
from skimage import measure
from skimage import morphology

# input_model_file  = 'CA1_ssmall_model.lm'
# input_morpho_file = "CA1_ssmall.h5"

input_model_file = "lms/_model.lm"
input_morpho_file = "CA1_small.h5"
input_label_file = "lm_annot/labels.hdf5"
output_file = 'morph_domains.png'

with open('lm_annot/list.pickle','rb') as f:
    list = pickle.load(f)

with h5py.File(input_label_file,'r') as f:
    labels = f['dendrite'][()]

print('labels.shape: ', labels.shape)

ids = np.unique(labels)
ids = ids[1:-2]
print('ids: ', ids)
    
with h5py.File(input_morpho_file,'r') as f:
    dendrite     = f['dendrite'][()]
    memb_v    = f['membrane vertices'][()]
    memb_f    = f['membrane faces'][()]
    PSD_ids   = f['PSD ids in membrane faces'][()]
    mito_v    = f['mitochondrion vertices'][()]
    mito_f    = f['mitochondrion faces'][()]
    pitch = f['unit length per voxel (um)'][()]


print('dendrite.shape: ', dendrite.shape)

dendrite = dendrite.astype(np.bool) ^ (labels > 0).astype(np.bool)
memb_v, memb_f, normals, values = measure.marching_cubes(dendrite, 0.5, spacing=(1,1,1))

trimesh.constants.tol.merge = 1e-7
mesh = trimesh.Trimesh(vertices=memb_v, faces=memb_f)
mesh_smooth = trimesh.smoothing.filter_laplacian(mesh, iterations=5)
memb_v = mesh_smooth.vertices
memb_f = mesh_smooth.faces

#memb_v = memb_v / pitch
mito_v = mito_v / pitch
print('pitch: ', pitch)

# Plot surface mesh
mlab.figure(bgcolor=(1.0,1.0,1.0), size=(400,741))
mlab.view(-180, 90, 700, [120.0, 120.0, 120.0 ])
mlab.triangular_mesh(mito_v[:,0] , mito_v[:,1] , mito_v[:,2] , mito_f, color=(1.0,1.0,0.6), opacity=0.6)
mlab.triangular_mesh(memb_v[:,0], memb_v[:,1], memb_v[:,2], memb_f, color=(0.7,0.7,0.7)  , opacity=0.3)
#mlab.triangular_mesh(memb_v[:,0], memb_v[:,1], memb_v[:,2], memb_f[PSD_ids,:], color=(1,0,0)  , opacity=0.3)

for id in ids:
    l = (labels == id).astype(np.bool)
    v, f, normals, values = measure.marching_cubes(l, 0.5, spacing=(1,1,1))
    mesh = trimesh.Trimesh(vertices=v, faces=f)
    mesh_smooth = trimesh.smoothing.filter_laplacian(mesh, iterations=5)
    v = mesh_smooth.vertices
    f = mesh_smooth.faces
    c = [x for x in list['list'] if x['id'] == id]
    r = c[0]['r']/256.0
    g = c[0]['g']/256.0
    b = c[0]['b']/256.0
    mlab.triangular_mesh(v[:,0], v[:,1], v[:,2], f, color=(r,g,b)  , opacity=0.3)

    
#
xvnum,yvnum,zvnum = dendrite.shape
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

