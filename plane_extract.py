import klampt
import time
import merge_triangle_mesh
import numpy as np
import sys

from klampt import vis
from klampt.vis import colorize
from klampt import math
from klampt import io
from klampt.model.create import primitives

def build_triangle_adjacency(verts, inds):
    v_map = []
    for i in range(len(verts)):
        v_map.append([])
    for i in range(len(inds)):
        for j in range(len(inds[i])):
            v_map[inds[i][j]].append(i)
    t_neighbors = -np.ones((inds.shape), dtype=np.int64)
    for i in range(len(inds)):
        count = {}
        for j in range(len(inds[i])):
            for k, id in enumerate(v_map[inds[i][j]]):
                if id in count:
                    count[id] += 1
                else:
                    count[id] = 1
        ind = 0
        for j, (k,v) in enumerate(count.items()):
            if v == 2:
                t_neighbors[i][ind] = k
                ind += 1
            if ind == 3:
                break
    return t_neighbors

def get_rot_frame(z_world_target):
    z_world_target = np.array(z_world_target)
    if np.linalg.norm(z_world_target) == 0:
        print(f"Error, z_world_target {z_world_target} has norm 0")
        sys.exit()
    y_hat = np.cross(z_world_target, np.array([1, 0, 0]))
    if np.linalg.norm(y_hat) == 0:
        y_hat = np.cross(z_world_target, np.array([0, 1, 0]))
    y_hat = y_hat / np.linalg.norm(y_hat)
    x_hat = np.cross(y_hat, z_world_target)
    x_hat = x_hat / np.linalg.norm(x_hat)
    return np.array([x_hat, y_hat, z_world_target]).T

world = klampt.WorldModel()

obj = world.makeRigidObject("tm")
g = obj.geometry()
g.loadFile("meshes/rod.off")
offset = 0.1
g.setCollisionMargin(offset)
vg = g.convert("VolumeGrid", 0.1)
g.setCollisionMargin(0)
obj.appearance().setColor(0,0,0,0)
inflated_obj = world.makeRigidObject("itm")
inflated_obj.geometry().set(vg.convert("TriangleMesh"))
inflated_obj.appearance().setColor(0,1,1)
infl_tm = inflated_obj.geometry().getTriangleMesh()
v = np.array(infl_tm.vertices).reshape(-1, 3)
inds = np.array(infl_tm.indices).reshape(-1, 3)
v, inds = merge_triangle_mesh.dedup_triangle_mesh(v, inds)
# colorize.colorize(inflated_obj, np.random.rand((len(inds))),
#     colormap="jet", feature="faces")
# print("Number of triangles: ", len(inds))
# print("Python defects: ", np.sum(build_triangle_adjacency(v, inds) == -1))
# print(len(v))
# print(v)
# for vert in v:
#     primitives.sphere(0.01, vert, world=world).appearance().setColor(*np.random.rand((3)))
my_pf = merge_triangle_mesh.get_empty_plane_finder()
ret_val = merge_triangle_mesh.merge_triangle_mesh(my_pf, v, inds, 0.5)
hm_ret = merge_triangle_mesh.get_heightmaps(my_pf, 0.1, 0.1)
planes = np.array(ret_val[0])
triangles = ret_val[1]
colors = np.zeros((len(inds), 4))
colors[:, 3] = 1
for i, p in enumerate(planes):
    c = np.random.rand((3))
    # primitives.sphere(0.01, p[1], world=world).appearance().setColor(*c)
    # primitives.sphere(0.01, p[1]+p[0], world=world).appearance().setColor(*c)
    R = get_rot_frame(p[0]).T.flatten()
    # primitives.box(1, 1, 0.001, R=R, t=p[1],
    #     world=world).appearance().setColor(*c, 0.25)
    # primitives.box(0.001, 0.001, 1, R=R, t=p[1]
    #     + math.so3.apply(R, 0.5*np.array([0, 0, 1])),
    #     world=world).appearance().setColor(*c)
    n_map = hm_ret[0][i]
    wp_map = hm_ret[1][i]
    for j, row in enumerate(wp_map):
        for k, val in enumerate(row):
            R = get_rot_frame(n_map[j][k]).T.flatten()
            n_len = 0.1
            primitives.sphere(0.01, val, world=world).appearance().setColor(*c)
            primitives.box(0.001, 0.001, n_len, R=R, t=np.array(val)
                + np.array(math.so3.apply(R, (n_len/2)*np.array([0, 0, 1]))),
                world=world).appearance().setColor(*c)
    for tri in triangles[i]:
        colors[tri, :3] = c + 0.1*np.random.rand(3)
colorize.colorize(inflated_obj, colors,
    feature="faces")
print("Number of planes found: ", len(planes))
vis.add("world", world)
vis.debug(world)
