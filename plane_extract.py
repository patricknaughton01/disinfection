import klampt
import time
import merge_triangle_mesh
import numpy as np

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

world = klampt.WorldModel()

obj = world.makeRigidObject("tm")
g = obj.geometry()
g.loadFile("rod.off")
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
# colorize.colorize(inflated_obj, np.random.rand((len(inds))),
#     colormap="jet", feature="faces")
# print("Number of triangles: ", len(inds))
# print("Python defects: ", np.sum(build_triangle_adjacency(v, inds) == -1))
# print(len(v))
# print(v)
# for vert in v:
#     primitives.sphere(0.01, vert, world=world).appearance().setColor(*np.random.rand((3)))
ret_val = merge_triangle_mesh.merge_triangle_mesh(v, inds, 0.3)
planes = np.array(ret_val[0])
triangles = ret_val[1]
colors = np.zeros((len(inds), 4))
colors[:, 3] = 1
for i, p in enumerate(planes):
    c = np.random.rand((3))
    # primitives.sphere(0.01, p[1], world=world).appearance().setColor(*c)
    # primitives.sphere(0.01, p[1]+p[0], world=world).appearance().setColor(*c)
    z_hat = np.array([0,0,1])
    rot_axis = np.cross(z_hat, p[0])
    theta = np.linalg.norm(rot_axis)
    if theta != 0:
        rot_axis /= theta
        theta = np.arcsin(theta).item()
        if z_hat @ p[0] < 0:
            theta += np.pi
        R = math.so3.from_axis_angle((rot_axis, theta))
    else:
        R = np.eye(3).flatten()
    # primitives.box(1, 1, 0.001, R=R, t=p[1],
    #     world=world).appearance().setColor(*c, 0.25)
    # primitives.box(0.001, 0.001, 1, R=R, t=p[1] + math.so3.apply(R, 0.5*z_hat),
    #     world=world).appearance().setColor(*c)

    for tri in triangles[i]:
        colors[tri, :3] = c + 0.1*np.random.rand(3)
colorize.colorize(inflated_obj, colors,
    feature="faces")
print("Number of planes found: ", len(planes))
vis.add("world", world)
vis.debug(world)
