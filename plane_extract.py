import klampt
import time
import merge_triangle_mesh
import numpy as np
import sys
import random

from klampt import vis
from klampt.vis import colorize
from klampt import math
from klampt import io
from klampt.model.create import primitives
from disinfection import Planner, WipeSurface, Wiper

def main():
    world = klampt.WorldModel()

    surface = world.makeRigidObject("surface")
    surface.geometry().loadFile("meshes/lumps.off")
    offset = 0.1
    surface.geometry().setCollisionMargin(offset)
    vg = surface.geometry().convert("VolumeGrid", 0.1)
    surface.geometry().setCollisionMargin(0)
    surface.appearance().setColor(0,0,0,0)
    inflated_obj = world.makeRigidObject("itm")
    inflated_obj.geometry().set(vg.convert("TriangleMesh", 0.1))
    inflated_obj.appearance().setColor(0,1,1, 0.1)
    infl_tm = inflated_obj.geometry().getTriangleMesh()
    v = np.array(infl_tm.vertices).reshape(-1, 3)
    inds = np.array(infl_tm.indices).reshape(-1, 3)
    v, inds = merge_triangle_mesh.dedup_triangle_mesh(v, inds)
    print("Loaded triangle mesh and inflated mesh")

    my_pf = merge_triangle_mesh.get_empty_plane_finder()
    ret_val = merge_triangle_mesh.merge_triangle_mesh(my_pf, v, inds, 0.5)
    print("Merged planes")
    merge_triangle_mesh.build_heightmaps(my_pf, 0.1, 0.1)
    hm_ret = merge_triangle_mesh.get_heightmaps(my_pf)
    metadata = merge_triangle_mesh.get_heightmap_metadata(my_pf);
    print("Retreived heightmaps")
    planes = np.array(ret_val[0])
    triangles = ret_val[1]
    # colors = np.zeros((len(inds), 4))
    # colors[:, 3] = 0.1
    # for i, p in enumerate(planes):
    #     c = np.random.rand((3))
    #     # primitives.sphere(0.01, p[1], world=world).appearance().setColor(*c)
    #     # primitives.sphere(0.01, p[1]+p[0], world=world).appearance().setColor(*c)
    #     R = get_rot_frame(p[0]).T.flatten()
    #     # primitives.box(1, 1, 0.001, R=R, t=p[1],
    #     #     world=world).appearance().setColor(*c, 0.25)
    #     # primitives.box(0.001, 0.001, 1, R=R, t=p[1]
    #     #     + math.so3.apply(R, 0.5*np.array([0, 0, 1])),
    #     #     world=world).appearance().setColor(*c)
    #     n_map = hm_ret[0][i]
    #     wp_map = hm_ret[1][i]
    #     # for j, row in enumerate(wp_map):
    #     #     for k, val in enumerate(row):
    #     #         R = get_rot_frame(n_map[j][k]).T.flatten()
    #     #         n_len = 0.1
    #     #         primitives.sphere(0.01, val, world=world).appearance().setColor(*c)
    #     #         primitives.box(0.001, 0.001, n_len, R=R, t=np.array(val)
    #     #             + np.array(math.so3.apply(R, (n_len/2)*np.array([0, 0, 1]))),
    #     #             world=world).appearance().setColor(*c, 0.1)
    #     for tri in triangles[i]:
    #         colors[tri, :3] = c + 0.1*np.random.rand(3)
    # colorize.colorize(inflated_obj, colors,
    #     feature="faces")
    print("Number of planes found: ", len(planes))

    wiper_obj = world.makeRigidObject("wiper")
    wiper_obj.geometry().loadFile("meshes/wiper.off")
    wiper_obj.geometry().set(wiper_obj.geometry().convert("VolumeGrid"))
    wiper_obj.appearance().setColor(0.5,0.5,0.5,0.75)
    wiper_handle = world.makeRigidObject("wiper_handle")
    wiper_handle.geometry().loadFile("meshes/wiper_handle.off")
    wiper_handle.geometry().set(wiper_handle.geometry().convert("VolumeGrid"))
    wiper = Wiper(wiper_obj, wiper_handle, rows=10,
        cols=10, lam=0, gamma_0=1.0, beta_0=100)
    ws = WipeSurface("surface", surface, wiper)

    # Just optimize wiping the first patch
    plane_ind = 0
    print(metadata[plane_ind])
    min_x, max_x, min_y, max_y = *metadata[plane_ind]
    start = np.array((get_random(min_x, max_x), get_random(min_y, max_y)))
    end = np.array((get_random(min_x, max_x), get_random(min_y, max_y)))
    num_steps = 10
    step = np.linalg.norm(end - start) / num_steps

    vis.add("world", world)
    # vis.debug(world)
    vis.show()
    ind = 0
    dt = 0.1
    while vis.shown():
        if ind < num_steps:
            pt = start + ind * (end - start) / num_steps
            # TODO: Find pt in grid of sample points and interpolate - need to
            # know spacing and location of center.
            ind += 1
        time.sleep(dt)
    vis.show(False)


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


def get_random(lower, upper):
    assert(lower <= upper)
    return random.random() * (upper - lower) + lower


if __name__ == "__main__":
    main()
