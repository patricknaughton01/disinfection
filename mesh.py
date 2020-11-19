from mesh_to_sdf import mesh_to_voxels

import trimesh
import skimage
import skimage.measure

mesh = trimesh.load('lumps.off')

voxels = mesh_to_voxels(mesh, 32, pad=True)

vertices, faces, normals, _ = skimage.measure.marching_cubes(voxels, level=0.1)
mesh = trimesh.Trimesh(vertices=vertices, faces=faces, vertex_normals=normals)
mesh.show()
