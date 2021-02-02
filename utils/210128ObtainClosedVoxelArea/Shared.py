
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


class Shared:



	def SharedPreprocess(self, params, comm_title):

		print('Annotator folder is being generated for', comm_title)
		targ = Params()
		targ.SetUserInfoAnnotator(params['Empty Folder for Annotator'])

		if os.path.isdir(targ.surfaces_path) == False:
			os.makedirs(targ.surfaces_path, exist_ok=True)
			os.makedirs(targ.surfaces_whole_path, exist_ok=True)
		if os.path.isdir(targ.skeletons_path) == False:
			os.makedirs(targ.skeletons_path, exist_ok=True)
			os.makedirs(targ.skeletons_whole_path, exist_ok=True)
		if os.path.isdir(targ.volume_path) == False:
			os.makedirs(targ.volume_path, exist_ok=True)
		if os.path.isdir(targ.paint_path) == False:
			os.makedirs(targ.paint_path, exist_ok=True)

		return targ



	def SharedGenerateInfoFile(self, ids_volume, surfaces_segment_info_json_file):
		ids_nums = np.unique(ids_volume, return_counts=True)
		ids   = ids_nums[0]
		names = [str(id).zfill(10) for id in ids]
		sizes = ids_nums[1]
		colormap = np.random.randint(255, size=(ids.shape[0], 3), dtype='int')

		if ids[0] == 0:
			ids   = np.delete(ids, 0)
			names.pop(0)
			sizes = np.delete(sizes, 0)
			colormap = np.delete(colormap, 0, 0)

		ids      = ids.tolist()
		sizes    = sizes.tolist()
		colormap = colormap.tolist()

		print('Constainer shape: ', ids_volume.shape)
		print('IDs  : ', ids)
		print('names: ', names)
		print('sizes: ', sizes)
		print('cols : ', colormap)

		##
		keys = ['id', 'name', 'size']
		data_dict = [dict(zip(keys, valuerecord)) for valuerecord in zip(ids, names, sizes)]

		for i in range(len(data_dict)):
			col = {'confidence': 0, 'r': colormap[i][0], 'g': colormap[i][1],  'b': colormap[i][2],  'act': 0}
			data_dict[i].update(col)

		print('data_dict: ', data_dict)

		with open( surfaces_segment_info_json_file , 'w') as f:
			json.dump(data_dict, f, indent=2, ensure_ascii=False)

		return ids


	def SharedPostProcess(self, params, targ, ids_volume):
		##
		print("params['Pitch in X (um)']", params['Pitch in X (um)'])
		print("params['Downsampling factor in X']", params['Downsampling factor in X'])

		ph = params['Pitch in X (um)']
		pw = params['Pitch in Y (um)']
		pz = params['Pitch in Z (um)']
		ch = int(params['Downsampling factor in X'])
		cw = int(params['Downsampling factor in Y'])
		cz = int(params['Downsampling factor in Z'])
		ph *= ch
		pw *= cw
		pz *= cz
		wmax = ids_volume.shape[0]
		hmax = ids_volume.shape[1]
		zmax = ids_volume.shape[2]
		##
		with h5py.File(targ.volume_file, 'w') as f:		
			f.create_dataset('volume', data=ids_volume)
		##
		data_dict = {
		    	'boundingbox_voxel':{
		    		'x': hmax,
		    		'y': wmax,
		    		'z': zmax
		    		},
		    	'boundingbox_um':{
		    		'x': ph * hmax,
		    		'y': pw * wmax,
		    		'z': pz * zmax
		    		},
		    	'pitch_um':{
		    		'x': ph,
		    		'y': pw,
		    		'z': pz
		    		},
				}
		with open( targ.surfaces_volume_description_json_file , 'w') as f:
			json.dump(data_dict, f, indent=2, ensure_ascii=False)

		print('')
		print('Annotator folder created.')
		print('')

		return True


	def GetClosedTrimesh(self,v,f,data):

		unzipped_tri = gzip.decompress(data)
		sub_face_id = []
		for i in range( f.shape[0] ) :
			if (unzipped_tri[i*3:i*3+3] == b'\x01\x01\x01') :
				sub_face_id.append(f[i,:])

		part_mesh = pymeshfix.MeshFix(v, np.array(sub_face_id))
		part_mesh.repair()
		part_mesh.plot() # Visualization of cloased meshes

#		voxels = pv.voxelize(part_mesh, density=part_mesh.length/200)

#		p = pv.Plotter()
#		p.add_mesh(voxels, color=True, show_edges=True, opacity=0.5)
#		p.add_mesh(surface, color="lightblue", opacity=0.5)
#		p.show()


		closed_v = part_mesh.v # numpy np.float array
		closed_f = part_mesh.f # numpy np.int32 array

		closed_mesh = trimesh.Trimesh(vertices=closed_v, faces=closed_f)

		closed_mesh.merge_vertices()
		closed_mesh.remove_degenerate_faces()
		closed_mesh.remove_duplicate_faces()

		# print("Volume: ", closed_mesh.volume)
		# Numpy-stl (mesh) にも同ルーチン有
		return closed_mesh

	def GenerateAnnotFile(self, params):

		comm_title = ""
		targ = self.SharedPreprocess(params, comm_title)
		##
		with h5py.File(params['Hdf5 file containing segmentation volume'], 'r') as f:
			if params['Container name'] not in f.keys():
				print('No container: ', params['Container name'])
				return False
			ids_volume = f[ params['Container name'] ][()]

		ids = self.SharedGenerateInfoFile(ids_volume, targ.surfaces_segment_info_json_file)
		self.SharedPostProcess(params, targ, ids_volume)

		## Special surfaces
		ph = params['Pitch in X (um)']
		pw = params['Pitch in Y (um)']
		pz = params['Pitch in Z (um)']

		# Target
		id = 1
		print('ids: ',ids)
		print('id : ',id)
		# Marching cube

		v_march, f_march, normals, values = measure.marching_cubes_lewiner(ids_volume==id, 0.5, spacing=(1,1,1))
		v_march = v_march - 1

		# Scaling 
		v_march[:,0] *= pw
		v_march[:,1] *= ph
		v_march[:,2] *= pz
		mesh = trimesh.Trimesh(vertices=v_march, faces=f_march)
#		trimesh.constants.tol.merge = 1e-7
#		mesh_smooth = trimesh.smoothing.filter_laplacian(mesh, iterations=5)

		print('vertices:', mesh.vertices.shape)
		print('faces   :', mesh.faces.shape)
		filename = os.path.join(targ.surfaces_whole_path, str(id).zfill(10)+'.stl')
		mesh.export(file_obj=filename)

		return True

