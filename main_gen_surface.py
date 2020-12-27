import h5py
import numpy as np
import os


from calcSurfaceArea import calcSurfaceArea

pitch = 0.008; # 8nm each
filename_morph='CA1dend_small.h5'

print('Load geometry.')
with h5py.File(filename_morph,'a') as f:
    dendrite = f['dendrite'][()]
    Mitochondrion    = f['Mitochondrion'][()]
    del f['membrane areas']
    del f['unit length per voxel (um)']
    del f['Smoothing: vertices']
    del f['Smoothing: faces']
    del f['Smoothing: area per face']


dendrite = dendrite.astype(np.bool)
# PSD     = dendrite_and_PSD.astype(np.bool)

surface_areas, smooth_vertices, smooth_faces, smooth_area_per_face = calcSurfaceArea(dendrite, pitch)

dendrite_not_mitochondrion = dendrite ^ Mitochondrion
loc = np.where(dendrite_not_mitochondrion > 0)[0].shape[0]

with h5py.File(filename_morph,'a') as w:
    w['membrane areas']  = surface_areas
    w['unit length per voxel (um)'] = pitch
    w['Smoothing: vertices']  = smooth_vertices
    w['Smoothing: faces']    = smooth_faces
    w['Smoothing: area per face']  = smooth_area_per_face
    w['voxel num of dendrite not mitochondrion']  = loc
    w['dendrite not mitochondrion']  = dendrite_not_mitochondrion

