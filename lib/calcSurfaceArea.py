from skimage import measure
from skimage import morphology
import trimesh
import h5py
import numpy as np
import os

def calcSurfaceArea(volume, pitch, num_iterations = 15):
	#
	# @param volume			input volume			(numpy 3d array, bool preferred)
	# @param pitch			unit length per voxel		(float)
	# @param num_iterations		iter num for smoothing		(int)
	# @return surface_areas		surface areas in voxel space(numpy 3d array)
	# @return smooth_vertices	vertices of smoothing mesh	(numpy 3xX array)
	# @return smooth_faces		faces of smoothing mesh		(numpy 3xX array)
        # @return Smooth_area_per_face  Areas of faces			(numpy 1xX array)

	volume = volume.astype(np.bool)
	print('volume.shape: ', volume.shape)
	xvnum, yvnum, zvnum = volume.shape

	# Obtain the region of membrane
	cytosol = morphology.binary_erosion(volume, morphology.ball(1))
	memb    = volume ^ cytosol

	# Generate linked list
	memb_id = np.array(np.where(memb)).T

	# Marching cube
	v_march, f_march, normals, values = measure.marching_cubes(volume, 0.5, spacing=(1,1,1))
	v_march = v_march - 1

	# Smoothing
	trimesh.constants.tol.merge = 1e-7
	mesh = trimesh.Trimesh(vertices=v_march, faces=f_march)
	mesh_smooth = trimesh.smoothing.filter_laplacian(mesh, iterations=15)
	mesh.merge_vertices()
	mesh.remove_degenerate_faces()
	mesh.remove_duplicate_faces()
	v_smooth = mesh_smooth.vertices
	f_smooth = mesh_smooth.faces

	# Obrain the surface area of each trianglar face.
	f_vec_ac = v_smooth[f_smooth[:,0]] - v_smooth[f_smooth[:,2]]
	f_vec_bc = v_smooth[f_smooth[:,1]] - v_smooth[f_smooth[:,2]]
	f_outer  = np.cross(f_vec_ac, f_vec_bc) 
	f_areas  = np.linalg.norm(f_outer, axis=1) / 2
	# np.sum(f_areas) == measure.mesh_surface_area(v_smooth,f_smooth)

	# Obtain the location ids of trianglar faces in the voxel space
	face_loc  = ( v_smooth[f_smooth[:,0]]+v_smooth[f_smooth[:,1]]+v_smooth[f_smooth[:,2]] ) / 3.0
	face_voxel = np.round( face_loc ).astype(np.int)
	face_voxel = (face_voxel < 0) + (face_voxel >= 0) * face_voxel
	face_voxel[:,0] = (face_voxel[:,0] >= xvnum) * (xvnum-1) + (face_voxel[:,0] < xvnum) * face_voxel[:,0]
	face_voxel[:,1] = (face_voxel[:,1] >= yvnum) * (yvnum-1) + (face_voxel[:,1] < yvnum) * face_voxel[:,1]
	face_voxel[:,2] = (face_voxel[:,2] >= zvnum) * (zvnum-1) + (face_voxel[:,2] < zvnum) * face_voxel[:,2]

	# Set membrane area in the voxel space
	memb_area = np.zeros_like(volume, dtype=np.float)
	for i in range(f_areas.shape[0]):
		j = 1
		while (j < 20):
			loc = np.where(morphology.ball(j, dtype=np.bool))
			loc = np.array(loc).T - j + face_voxel[i,:]
			loc = loc[np.all(loc >= 0, axis=1),:]
			loc = loc[loc[:,0] < xvnum,:]
			loc = loc[loc[:,1] < yvnum,:]
			loc = loc[loc[:,2] < zvnum,:]
			tmp = memb[loc[:,0],loc[:,1],loc[:,2]]
			if np.sum(tmp) > 0:
				memb_area[loc[:,0],loc[:,1],loc[:,2]] += tmp.astype(np.float) / np.sum(tmp) * f_areas[i] # f_areas[i]
				break;
			else :
				j = j + 1
				# print(j, i)
				# print(tmp)
				# print(loc)

	surface_areas  = memb_area * pitch * pitch
	smooth_vertices = v_smooth * pitch
	smooth_faces    = f_smooth
	smooth_area_per_face = f_areas * pitch * pitch
	return surface_areas, smooth_vertices, smooth_faces, smooth_area_per_face






