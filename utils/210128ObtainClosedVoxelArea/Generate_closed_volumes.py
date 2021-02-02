##
import sys
import os, glob, pickle
from os import path, pardir
import json
import sqlite3
import h5py
import numpy as np
from Params import Params
from skimage import measure
from skimage import morphology
import trimesh
import gzip
import pymeshfix
import pyvista as pv

from Shared import Shared

main_dir = path.abspath(path.dirname(sys.argv[0]))
##
class GenerateClosedvolumes(Shared):
	def __init__(self):
		pass


	def AnalyzeAnnotFile(self, params):
		targ = Params()
		targ.SetUserInfoAnnotator(params['Empty Folder for Annotator'])

		surface_path  = targ.surfaces_whole_path
		paint_path    = targ.paint_path

		whole_mesh_filenames = glob.glob(os.path.join(surface_path, "*.stl"))
		ids_volumes  = {}

		ph = params['Pitch in X (um)']
		pw = params['Pitch in Y (um)']
		pz = params['Pitch in Z (um)']

		##
		with h5py.File(targ.volume_file, 'r') as f:		
			ids_volume = f['volume'][()]
		ids_volume = (ids_volume > 0).astype(np.int)
		new_labels = np.zeros_like(ids_volume).astype(np.int)
		##

		## Load surface meshes
		id = 1
		whole_mesh_filename = os.path.join(surface_path, str(id).zfill(10)+".stl")
		whole_mesh = trimesh.load( whole_mesh_filename )
		surf_vertices = whole_mesh.vertices
		surf_faces    = whole_mesh.faces

		surf_vertices[:,0] /= pw
		surf_vertices[:,1] /= ph
		surf_vertices[:,2] /= pz
		pitch = 1

		whole_mesh_name_wo_ext  = os.path.splitext(os.path.basename(whole_mesh_filename))[0]
		part_mesh_name_wildcard = os.path.normpath(os.path.join(paint_path, whole_mesh_name_wo_ext+"-*.pickle"))
		part_mesh_filenames = glob.glob(part_mesh_name_wildcard)

		## Check whether painted meshes
		if part_mesh_filenames == [] :
			return False

		ids = []
		for part_mesh_filename in part_mesh_filenames :
			with open(part_mesh_filename, 'rb') as file:
				data = pickle.load(file)
			closed_mesh = self.GetClosedTrimesh(surf_vertices, surf_faces, data['painted'])
			if closed_mesh.volume is None :
				continue
			###
			id = os.path.basename(part_mesh_filename) 
			id = os.path.splitext(id)[0]
			id = int( id.split('-')[1] )
			print('ID: ', id,', Volume:', closed_mesh.volume)
			ids.append(id)

			#filename = os.path.join(targ.paint_path, str(id).zfill(4)+'.stl')
			#closed_mesh.export(file_obj=filename)
			#continue
			###

			part_faces   = closed_mesh.faces
			part_verts   = closed_mesh.vertices
			unique_ids_verts = np.unique(np.ravel(part_faces))
			unique_verts = part_verts[unique_ids_verts]

			wmin = np.floor( np.min(unique_verts[:,0])).astype(int)
			hmin = np.floor( np.min(unique_verts[:,1])).astype(int)
			zmin = np.floor( np.min(unique_verts[:,2])).astype(int)

			wmax = np.floor( np.max(unique_verts[:,0])).astype(int)
			hmax = np.floor( np.max(unique_verts[:,1])).astype(int)
			zmax = np.floor( np.max(unique_verts[:,2])).astype(int)

			print('wmin, wmax, wdiff: ', wmin, wmax, wmax-wmin)
			print('hmin, hmax, hdiff: ', hmin, hmax, hmax-hmin)
			print('zmin, zmax, zdiff: ', zmin, zmax, zmax-zmin)

			## Trimesh
			# v = self.GetVolumeTrimesh(closed_mesh)
			## PyVista
			v = self.GetVolumePyVista(part_verts, part_faces, wmin, hmin, zmin, wmax, hmax, zmax)
			#print('dir(v)  : ', dir(v))
			#print('v.keys(): ', v.keys())
			#print('v[12]: ', v[12])
			wnum = v.shape[0]
			hnum = v.shape[1]
			znum = v.shape[2]

			wmin += 1
			hmin += 1
			zmin += 1

			print('wnum, hnum, znum : ', wnum, hnum, znum)
			new_labels[wmin:wmin+wnum , hmin:hmin+hnum, zmin:zmin+znum] += v.astype(np.int) * id


####
#### Dilution, clarification, etc
####

		new_labels_processed = np.zeros_like(ids_volume)
		for id in ids :
			
			## Pickup the target area
			print('Dilution, id: ', id)
			target_area = morphology.binary_dilation(new_labels == id, selem=morphology.ball(1), out=None).astype(np.int)
			
			## Pickup segmented areas
			labels_to_pickup_segmentation = morphology.label(ids_volume * (target_area == 0))
			us, counts   = np.unique(labels_to_pickup_segmentation, return_counts=True)
			#print('us                : ', us)
			#print('segmentation count: ', counts)
			segmented_us = us[counts < 30]
			#print('segmented us      : ', segmented_us)
			
			segments = np.zeros_like(ids_volume)
			for segmented_u in segmented_us:
				segments += (labels_to_pickup_segmentation == segmented_u).astype(np.int)

			## Merge target area with segmented areas if they are connected.
			target_plus_segment = target_area*ids_volume + segments
			labels_to_remove_segmented_target = morphology.label( target_plus_segment > 0 )
			u, counts   = np.unique(labels_to_remove_segmented_target, return_counts=True)
			labels_segmented_target_removed = (labels_to_remove_segmented_target == u[counts.argsort()[-2]]).astype(np.int)
			
			## Assign the id to (target area in cytosol) and (segmented area).
			new_labels_processed += labels_segmented_target_removed*id

		with h5py.File('labels.hdf5', 'w') as f:
			f.create_dataset('dendrite', data=new_labels_processed)

		with h5py.File('labeled_cytosol.hdf5', 'w') as f:
			f.create_dataset('dendrite', data=new_labels_processed + ids_volume)


		return ids_volume, new_labels



	def GetVolumeTrimesh(self, closed_mesh):
		v = closed_mesh.voxelized(pitch = pitch)
		print('v.matrix.shape: ', v.matrix.shape)
		wnum = v.matrix.shape[0]
		hnum = v.matrix.shape[1]
		znum = v.matrix.shape[2]
		return v.matrix

	def GetVolumePyVista(self, verts, faces, wmin, hmin, zmin, wmax, hmax, zmax):
		verts = np.array(verts)
		faces = np.array(faces)
		num = faces.shape[0]
		faces = np.hstack([np.ones([num,1]).astype(int)*3,faces])
		surf = pv.PolyData(np.array(verts), np.array(faces))
		ix = np.arange(wmin, wmax, 1)
		iy = np.arange(hmin, hmax, 1)
		iz = np.arange(zmin, zmax, 1)

		x, y, z = np.meshgrid(ix, iy, iz)
		grid = pv.StructuredGrid(x, y, z)
		ugrid = pv.UnstructuredGrid(grid)
		selection = ugrid.select_enclosed_points(surf, tolerance=0.0, check_surface=False)
		mask = selection.point_arrays['SelectedPoints'].view(np.bool)
		voxels = mask.reshape([iz.shape[0] , ix.shape[0], iy.shape[0] ])
		voxels = voxels.transpose((1, 2, 0))
		return voxels




if __name__ == "__main__":

	params = {}
	params['Hdf5 file containing segmentation volume'] = os.path.join(main_dir, 'CA1_small.h5')
	params['Container name']     = 'dendrite'
	params['Empty Folder for Annotator']  = os.path.join(main_dir, 'annot_lm')

	params['Pitch in X (um)'] = 0.02
	params['Pitch in Y (um)'] = 0.02
	params['Pitch in Z (um)'] = 0.02

	params['Downsampling factor in X'] = 1
	params['Downsampling factor in Y'] = 1
	params['Downsampling factor in Z'] = 1
	
	p = GenerateClosedvolumes()
#	p.GenerateAnnotFile(params)
	p.AnalyzeAnnotFile(params)



