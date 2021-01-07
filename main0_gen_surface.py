import h5py
import numpy as np
import os

from lib.calcSurfaceArea import calcSurfaceArea

pitch = 0.02; # 20nm
#filename_morph='CA1.h5'
#filename_morph='CA1_ssmall.h5'
filename_morph='CA1_small.h5'

print('Load geometry.')
delkeys = ['membrane areas', 'unit length per voxel (um)',\
           'Smoothing: vertices', 'Smoothing: faces', 'Smoothing: area per face',\
           'voxel num of dendrite not mitochondrion','dendrite not mitochondrion',\
           'membrane areas in volume',\
           'membrane vertices', 'membrane faces','membrane area per face', 'PSD ids in membrane faces',\
           'mitochondrion areas in volume', 'mitochondrion vertices', 'mitochondrion faces',\
           'mitochondrion area per face']
with h5py.File(filename_morph,'a') as f:
    print(f.keys())
    dendrite     = f['dendrite'][()]
    PSD          = f['PSD'][()]
    Mitochondria = f['Mitochondrion'][()]
    for key in delkeys:
        if key in f.keys():
            del f[key]

dendrite = dendrite.astype(np.bool)
PSD      = PSD.astype(np.bool)
Mitochondria = Mitochondria.astype(np.bool)
_Mitochondria = np.logical_not( Mitochondria )

memb_areas, memb_verts, memb_faces, memb_area_per_face, id_face_psd = calcSurfaceArea(pitch, dendrite, PSD = PSD)
mito_areas, mito_verts, mito_faces, mito_area_per_face = calcSurfaceArea(pitch, _Mitochondria)

dendrite_not_mitochondria = dendrite ^ Mitochondria
loc = np.where(dendrite_not_mitochondria > 0)[0].shape[0]

with h5py.File(filename_morph,'a') as w:
    w['unit length per voxel (um)'] = pitch
    w['membrane areas in volume']   = memb_areas  
    w['membrane vertices']      = memb_verts
    w['membrane faces']         = memb_faces
    w['membrane area per face'] = memb_area_per_face
    w['PSD ids in membrane faces']   = id_face_psd
    w['mitochondrion areas in volume'] = mito_areas
    w['mitochondrion vertices']      = mito_verts
    w['mitochondrion faces']         = mito_faces
    w['mitochondrion area per face'] = mito_area_per_face
    w['voxel num of dendrite not mitochondrion']  = loc
    w['dendrite not mitochondrion']  = dendrite_not_mitochondria

