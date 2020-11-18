import klampt
import time
import sys
import numpy as np
import cvxpy as cp
import faulthandler
# import matplotlib.pyplot as plt
import scipy as sp
import scipy.spatial as spatial
import scipy.sparse as sparse
import scipy.optimize as opt
import math as pmath
import numba
faulthandler.enable()

from klampt import vis
from klampt import math
from klampt import io
from klampt.vis import colorize
from klampt.model.create import primitives
from numba import jit
from multiprocessing import Process, Manager

world = klampt.WorldModel()
DEBUG = False
TIMING = False
DISPLAY = True

def main():
    global world
    obj = world.makeRigidObject("tm")
    oort = 1/(2**0.5)
    g = obj.geometry()
    g.loadFile("keys.off")
    # obj.setTransform([1, 0, 0, 0, -oort, oort, 0, -oort, -oort], [0,0.2,0])
    # vg_obj = world.makeRigidObject("vg")
    # vg_obj.geometry().loadFile("keys.off")
    # vg_obj.geometry().set(vg_obj.geometry().convert("VolumeGrid"))
    # vg_obj.geometry().setCollisionMargin(0.5)
    # vg_obj.appearance().setColor(0.3, 0.5, 0.5, 0.0)

    wiper_obj = world.makeRigidObject("wiper")
    wiper_obj.geometry().loadFile("wiper.off")
    wiper_obj.geometry().set(wiper_obj.geometry().convert("VolumeGrid"))
    wiper_obj.appearance().setColor(0.5,0.5,0.5,0.2)
    # wiper_obj.appearance().setColor(1,0,1,1)
    wiper_handle = world.makeRigidObject("wiper_handle")
    wiper_handle.geometry().loadFile("wiper_handle.off")
    wiper_handle.geometry().set(wiper_obj.geometry().convert("VolumeGrid"))

    dt = 0.1
    step = 0.01
    max = 1.0
    min = -0.1
    state = "left"

    sizes = [10, 30, 50]
    RUNS = 1
    if TIMING:
        RUNS = 11
    for i in range(1):#len(sizes)):
        wiper = Wiper(wiper_obj, wiper_handle, rows=sizes[i],
            cols=sizes[i], lam=0, gamma_0=1.0, beta_0=100)
        ws = WipeSurface("tm", obj, wiper)
        # R, t = wiper.getDesiredTransform([0,0,0], [0, -0.25, 1], [1, 0, 0], 0.1)
        # wiper.setTransform(R.T.flatten().tolist(), t.tolist())
        wiper.setTransform([1,0,0,0,1,0,0,0,1], [0.0,0.0,0.0])
        gamma, _ = wiper.eval_wipe_step(wiper.getTransform(), ([1,0,0,0,1,0,0,0,1], [0.0,0.01,0.0]), ws)
        if TIMING:
            start_time = time.monotonic()
            for j in range(RUNS):
                covered = ws.get_covered_triangles()
                t = np.random.rand(3) * 0.9
                t[2] = 0
                wiper.setTransform([1,0,0,0,1,0,0,0,1], t)
            end_time = time.monotonic()
            print(len(wiper.ray_t))
            print("\t\t".join(("Ray", "Opt", "Clean")))
            print("\t".join(("Mean", "Std", "Mean", "Std", "Mean", "std")))
            print(",".join((f"{np.mean(wiper.ray_t[1:]):.4g}",
                f"{np.std(wiper.ray_t[1:]):.4g}",
                f"{np.mean(wiper.opt_t[1:]):.4g}",
                f"{np.std(wiper.opt_t[1:]):.4g}",
                f"{np.mean(wiper.clean_t[1:]):.4g}",
                f"{np.std(wiper.clean_t[1:]):.4g}")))
            print("Average total time: {:.4g}".format( (end_time - start_time)/RUNS ))
    # ws.update_infection(ws.get_covered_triangles())
    # ws.update_infection(gamma)
    ws.update_colors()
    planner = Planner(ws, wiper)
    points = planner.gen_transforms(([1,0,0,0,1,0,0,0,1],[1,1,0]),
        ([1,0,0,0,1,0,0,0,1], [0.0,0.0,0.0]), 0.01)
    manager = Manager()
    shared_list = manager.list()
    for i in range(len(points)):
        shared_list.append(None)
    parallel_get_covered_triangles(points, shared_list, ws.obj, ws.verts,
        ws.inds, wiper.max_h, wiper.width, wiper.height, wiper.rows,
        wiper.cols, wiper.tot, wiper.dx, wiper.dy, wiper.Qx, wiper.Qy,
        wiper.id_top_points, wiper.lam, wiper.id_norm, ws.t_normals,
        ws.t_neighbors)
    for cover in shared_list:
        if cover is not None:
            for k in cover:
                ws.infection_level[k] = 0
    ws.update_colors()
    if DISPLAY:
        vis.add("world", world)
        # vis.getViewport().setTransform(([0, -1, 0,
        #     -oort, 0, -oort,
        #     oort, 0, -oort], [-1.5, 0.5, 2]))
        vis.show()
        # time.sleep(3)
        ind = 0
        print("Num interpolated points", len(points))
        while vis.shown():
            # if ind < len(points):
            #     wiper.setTransform(*points[ind])
            #     ind += 1
            #     time.sleep(0.1)
            # sim.simulate(dt)
            # R, t = wiper.obj.getTransform()
            # if state == "left":
            #     if t[1] < max:
            #         move_t = list(t[:])
            #         move_t[1] += step
            #         start = time.monotonic()
            #         wiper.wipe_step((R,t), (R, move_t), ws)
            #         print(time.monotonic() - start)
            #     else:
            #         state = "right"
            # elif state == "right":
            #     if t[1] > min:
            #         move_t = list(t[:])
            #         move_t[1] -= step
            #         wiper.wipe_step((R,t), (R,move_t), ws)
            #     else:
            #         state = "stop"
            time.sleep(dt)
        vis.show(False)
        vis.kill()


@jit(nopython=True)
def reduction(infection, gamma):
    reduction = 0
    for key in gamma:
        reduction += infection[key] - infection[key]*gamma[key]
    return reduction


@jit(nopython=True)
def combine_gamma(gamma1, gamma2):
    for key in gamma2:
        if key in gamma1:
            gamma1[key] *= gamma2[key]
        else:
            gamma1[key] = gamma2[key]
    return gamma1


class Planner:
    def __init__(self, surface, wiper, mode='online', plan_type='pfield'):
        self.surface = surface
        self.wiper = wiper
        self.plan_type = plan_type
        self.mode = mode
        self.wipes = []
        self.offsets = None
        self.cache = {}

    def gen_transforms(self, s_t, e_t, step_size):
        dist = math.se3.distance(e_t, s_t, Rweight=0.0)
        num_steps = int(dist / step_size)
        if num_steps <= 0:
            return [s_t, e_t]
        interpolants = []
        interpolator = math.se3.interpolator(s_t, e_t)
        for i in range(num_steps):
            t = interpolator(i / num_steps)
            interpolants.append(t)
        if dist % step_size != 0:
            interpolants.append(e_t)
        return interpolants

    def eval_wipe(self, x):
        """Evaluate a wipe whose parameters are encoded in the 5-vector x
        x[:2]:  start position of origin of wiper (x,y)
        x[2:4]: end position of origin of wiper(x,y)
        x[4]:  rotation of the wiper about z axis
        """
        print(x)
        key = tuple(x)
        if key in self.cache:
            return self.cache[key]
        orig_transform = self.wiper.getTransform()
        mag = 0.01
        t_0 = np.append(x[:2], [0])
        t_1 = np.append(x[2:4], [0])
        dist = np.linalg.norm(t_1 - t_0)
        c = pmath.cos(x[4])
        s = pmath.sin(x[4])
        R = [c, s, 0, -s, c, 0, 0, 0, 1]
        self.wiper.setTransform(R, t_0)
        move_vec = t_1 - t_0
        start_cover = None
        total_gamma = numba.typed.Dict.empty(numba.types.int64, numba.types.float64)
        while np.linalg.norm(move_vec) > mag:
            t_step = t_0 + mag * move_vec / np.linalg.norm(move_vec)
            gamma, start_cover = self.wiper.eval_wipe_step((R, t_0),
                (R, t_step), self.surface, start_cover=start_cover)
            total_gamma = combine_gamma(total_gamma, gamma)
            t_0 = t_step
            move_vec = t_1 - t_0
        # move_vec is smaller than mag now so we can just move to the endpoint
        gamma, _ = self.wiper.eval_wipe_step((R, t_0), (R, t_1), self.surface, start_cover=start_cover)
        total_gamma = combine_gamma(total_gamma, gamma)
        total_reduction = reduction(self.surface.infection_level, total_gamma)
        self.wiper.setTransform(*orig_transform)
        self.cache[key] = -total_reduction + 0.01 * dist
        return self.cache[key]

    def get_wipe(self):
        if self.plan_type == 'pfield':
            mag = 0.01
            vals = [[-mag, 0, mag] for i in range(3)]
            R, t = self.wiper.getTransform()
            t = np.array(t)
            if self.offsets is None:
                self.offsets = np.array(self.get_dirs(vals))
                norm = np.linalg.norm(self.offsets, axis=1, keepdims=True)
                for i, n_val in enumerate(norm):
                    if n_val == 0:
                        norm[i] = 1
                self.offsets = mag * (self.offsets / norm)
            d_vals = []
            for i in range(len(self.offsets)):
                d_vals.append(np.sum(self.surface.infection_level *
                    self.wiper.wipe_step((R, t),
                    (R, t + self.offsets[i, :]), self.surface)))
            min_d = None
            min_ind = -1
            for i, val in enumerate(d_vals):
                if min_d is None or val < min_d:
                    min_d = val
                    min_ind = i
            return ((R, t), (R, t + self.offsets[min_ind, :]))
        elif self.plan_type == 'greedy':
            res = opt.minimize(self.eval_wipe, np.zeros((5,)),
                method='Powell', options={'xtol':0.1, 'ftol':0.1})
            print(res.x)
            print(f'Iterations: {res.nit}')
            print(f'F Evals: {res.nfev}')
            return res.x

    def get_dirs(self, values, dim=0):
        if dim == len(values) - 1:
            # Base case
            combos = []
            for v in values[dim]:
                combos.append([v])
            return combos
        combos = []
        for v in values[dim]:
            n_combos = self.get_dirs(values, dim+1)
            for c in n_combos:
                c.insert(0, v)
                combos.append(c)
        return combos


class WipeSurface:
    def __init__(self, name, obj, wiper):
        self.name = name
        self.obj = obj
        self.tm = self.obj.geometry().getTriangleMesh()
        self.verts = np.array(self.tm.vertices).reshape(
            (len(self.tm.vertices)//3,3))
        self.inds = np.array(self.tm.indices,dtype=np.int32).reshape(
            (len(self.tm.indices)//3,3))
        R, t = self.obj.getTransform()
        R = np.array(R).reshape(3, -1).T
        t = np.array(t)
        self.verts = (R @ self.verts.T).T + t
        if DEBUG:
            for i in range(100):
                primitives.sphere(0.001, self.verts[i], world=world).appearance().setColor(0,1,1)
        self.num_triangles = len(self.inds)
        self.v_map = None
        self.t_neighbors = -np.ones((self.inds.shape), dtype=np.int64)
        self.build_triangle_adjacency()
        self.t_normals = None
        self.build_triangle_normals()
        self.infection_level = np.ones(self.num_triangles) + np.random.rand(self.num_triangles)/10
        self.infection_level[0] = 0
        self.obj.appearance().setColor(0.0, 0.0, 1.0)
        self.update_colors()
        self.covered = []
        self.hit_dist_thresh = 0
        self.epsilon = 1e-5
        self.visited_triangles = []

        self.wiper = wiper

        self.h = cp.Variable(self.wiper.tot)
        self.objective = cp.Minimize(cp.sum_squares(self.h / self.wiper.max_h)
            + self.wiper.lam * (cp.quad_form(self.h, self.wiper.Qy)
            + cp.quad_form(self.h, self.wiper.Qx)))
        self.min_h = cp.Parameter((self.wiper.tot,), nonneg=True)
        self.constraints = [self.h <= self.wiper.max_h, self.h >= self.min_h]
        self.prob = cp.Problem(self.objective, self.constraints)

    def build_triangle_adjacency(self):
        self.v_map = []
        for i in range(len(self.verts)):
            self.v_map.append([])
        for i in range(self.num_triangles):
            for j in range(len(self.inds[i])):
                self.v_map[self.inds[i][j]].append(i)
        for i in range(self.num_triangles):
            count = {}
            for j in range(len(self.inds[i])):
                for k, id in enumerate(self.v_map[self.inds[i][j]]):
                    if id in count:
                        count[id] += 1
                    else:
                        count[id] = 1
            ind = 0
            for j, (k,v) in enumerate(count.items()):
                if v == 2:
                    self.t_neighbors[i][ind] = k
                    ind += 1
                if ind == 3:
                    break

    def build_triangle_normals(self):
        self.t_normals = np.empty((self.num_triangles, 3), dtype=np.float64)
        for i in range(self.num_triangles):
            a = self.verts[self.inds[i][0], :]
            b = self.verts[self.inds[i][1], :]
            c = self.verts[self.inds[i][2], :]
            v1 = b - a
            v2 = c - a
            cx = np.cross(v1, v2)
            self.t_normals[i, :] = cx / np.linalg.norm(cx)

    def get_covered_triangles(self):
        global world, DEBUG
        """Get the indices of the covered triangles by a wipe represented by
        the passed in wiper.
        """
        ray_s = time.monotonic()
        min_h = np.zeros(self.wiper.tot)
        # Store index of the triangle that this point raycasts to
        h_t_correspondence = -np.ones(self.wiper.tot, dtype=np.long)
        hit_flag = False
        for i in range(self.wiper.tot):
            start_pt = self.wiper.top_points[i][:3]
            hit, pt = self.obj.geometry().rayCast_ext(
                start_pt, self.wiper.norm)
            if hit >= 0:
                hit_flag = True
                pt = np.array(pt)
                if DEBUG:
                    primitives.sphere(0.001, pt, world=world).appearance().setColor(0,1,0)
                min_h[i] = self.wiper.max_h - np.linalg.norm(start_pt - pt)
                if min_h[i] >= 0:
                    h_t_correspondence[i] = hit
        if not hit_flag:
            return numba.typed.Dict.empty(numba.types.int64, numba.types.float64)
        min_h = np.maximum(min_h, 0)
        ray_e = time.monotonic()
        # Solve opt problem to find contacted triangles
        opt_s = time.monotonic()
        self.min_h.value = min_h
        result = self.prob.solve()
        h_val = self.h.value
        opt_e = time.monotonic()

        if DEBUG:
            for i in range(self.wiper.tot):
                pt = self.wiper.top_points[i][:3] + (self.wiper.max_h - h_val[i]) * self.wiper.norm
                primitives.sphere(0.002, pt, world=world).appearance().setColor(1,0,1)

        clean_s = time.monotonic()
        contact = np.abs(h_val - min_h) < 1e-3
        covered_triangles = numba.typed.Dict.empty(numba.types.int64,
            numba.types.float64)
        visited = numba.typed.Dict.empty(numba.types.int64, numba.types.boolean)
        for i in range(self.wiper.tot):
            tind = h_t_correspondence[i]
            if tind > -1 and visited.get(tind, False) == False:
                covered_triangles.update(interpolate_contact(self.verts,
                    self.inds, tind, visited, contact, h_val,
                    self.wiper.max_h, self.wiper.width,
                    self.wiper.height, self.wiper.rows, self.wiper.cols,
                    self.wiper.dx, self.wiper.dy, self.wiper.lam,
                    self.wiper.norm, self.wiper.H_i, self.t_normals,
                    self.t_neighbors
                ))
        clean_e = time.monotonic()
        if TIMING:
            self.wiper.ray_t.append(ray_e - ray_s)
            self.wiper.opt_t.append(opt_e - opt_s)
            self.wiper.clean_t.append(clean_e - clean_s)
        self.covered = covered_triangles
        return self.covered

    def clear_covered(self):
        self.covered = np.zeros(self.num_triangles)

    def update_infection(self, gamma):
        """Update the infection_level of the triangles in this mesh.

        Parameters
        ----------------
        gamma: np.ndarray (self.num_triangles,)
            for each triangle, the proportion of infection remaining
        """
        update_infection(self.infection_level, gamma)

    def update_colors(self):
        colorize.colorize(self.obj, self.infection_level,
            colormap="jet", feature="faces")


def parallel_get_covered_triangles(transforms, coverages, ws_obj, verts, inds,
    max_h, width, height, rows, cols, tot, dx, dy, Qx, Qy, id_top_points, lam,
    id_norm, normals, t_neighbors, offset=0, chunk=32):
    n = len(transforms)
    print(n, offset)
    if n > chunk:
        proc = Process(target=parallel_get_covered_triangles, args=[
            transforms[:n//2], coverages, ws_obj, verts, inds,
            max_h, width, height, rows, cols, tot, dx, dy, Qx, Qy,
            id_top_points, lam, id_norm, normals, t_neighbors, offset, chunk])
        proc.start()
        parallel_get_covered_triangles(transforms[n//2:], coverages, ws_obj,
            verts, inds, max_h, width, height, rows, cols, tot, dx, dy,
            Qx, Qy, id_top_points, lam, id_norm, normals, t_neighbors,
            offset=offset+n//2, chunk=chunk)
        proc.join()
    else:
        for t_ind, t in enumerate(transforms):
            R, p = t
            H = np.array(math.se3.homogeneous((R,p)))
            H_i = np.array(math.se3.homogeneous(math.se3.inv((R,p))))
            R_mat = H[:3,:3]
            norm = R_mat @ id_norm
            top_points = (H @ id_top_points.T).T

            min_h = np.zeros(tot)
            h_t_correspondence = -np.ones(tot, dtype=np.long)
            hit_flag = False
            for i in range(tot):
                start_pt = top_points[i][:3]
                hit, pt = ws_obj.geometry().rayCast_ext(
                    start_pt, norm)
                if hit >= 0:
                    hit_flag = True
                    pt = np.array(pt)
                    min_h[i] = max_h - np.linalg.norm(start_pt - pt)
                    if min_h[i] >= 0:
                        h_t_correspondence[i] = hit
            if not hit_flag:
                coverages[offset+t_ind] = {}
            h = cp.Variable(tot)
            objective = cp.Minimize(cp.sum_squares(h / max_h)
                + lam * (cp.quad_form(h, Qy)
                + cp.quad_form(h, Qx)))
            min_h = np.maximum(min_h, 0)
            constraints = [h <= max_h, h >= min_h]
            prob = cp.Problem(objective, constraints)
            result = prob.solve()
            h_val = h.value

            contact = np.abs(h_val - min_h) < 1e-3
            covered_triangles = numba.typed.Dict.empty(numba.types.int64,
                numba.types.float64)
            visited = numba.typed.Dict.empty(numba.types.int64,
                numba.types.boolean)
            for i in range(tot):
                tind = h_t_correspondence[i]
                if tind > -1 and visited.get(tind, False) == False:
                    covered_triangles.update(interpolate_contact(verts,
                        inds, tind, visited, contact, h_val, max_h, width,
                        height, rows, cols, dx, dy, lam, norm, H_i, normals,
                        t_neighbors
                    ))
            coverages[offset+t_ind] = dict(covered_triangles)


@jit(nopython=True)
def interpolate_contact(verts, inds, triangle, visited, grid, heights, max_h,
    width, height, rows, cols, dx, dy, lam, norm, H_i, normals, t_neighbors,
):
    open = [triangle]
    covered = numba.typed.Dict.empty(numba.types.int64, numba.types.float64)
    while len(open) > 0:
        triangle = open.pop()
        visited[triangle] = True
        centroid = np.zeros(4)
        sum_coords = np.zeros(3)
        for i in range(3):
            sum_coords += verts[inds[triangle, i], :]
        # centroid[:3] = np.mean(verts[inds[triangle, :], :], 0)
        centroid[:3] = sum_coords / 3
        centroid[3] = 1
        proj_c = H_i @ centroid
        if (0 < proj_c[0] < width and 0 < proj_c[1] < height
            and 0 <= proj_c[2] <= max_h
        ):
            # Bilinearly interpolate between the four h values around this pt
            # g's formatted as [col, row] or [x, y]
            #  g2 x  x g4
            #  g1 x  x g3
            g1 = [proj_c[0] // dx, proj_c[1] // dy]
            g2 = [g1[0], g1[1] + 1]
            g3 = [g1[0] + 1, g1[1]]
            g4 = [g3[0], g3[1] + 1]
            gs = [g1,g2,g3,g4]
            # How far is p from g1 in x and y normalized by distances between
            # adjacent grid points
            mx = (proj_c[0] % dx) / dx
            my = (proj_c[1] % dy) / dy
            vs = [grid[int(g[0] + cols * g[1])] for g in gs]
            i1 = (1 - mx) * vs[0] + mx * vs[2]
            i2 = (1 - mx) * vs[1] + mx * vs[3]
            # Exponential fall off for points far from the heights
            var = 1e-5/(1e-8 + lam)
            avg_h = np.mean(np.array([heights[int(g[0] + cols * g[1])] for g in gs]))
            coverage = (((1 - my) * i1 + my * i2)
                * np.exp(-0.5 * (proj_c[2] - avg_h)**2 / var)
                * max(0, (-norm.T @ normals[triangle])))
            if coverage > 0:
                covered[triangle] = coverage
            for i in range(len(t_neighbors[triangle])):
                if visited.get(t_neighbors[triangle, i], False) == False:
                    open.append(t_neighbors[triangle, i])
    return covered


def update_infection(infection, gamma):
    for key in gamma:
        infection[key] *= gamma[key]


class Wiper:
    def __init__(self, obj, handle, gamma_0=0.5, beta_0=1.0, rows=10,
        cols=10, lam=0
    ):
        self.obj = obj
        self.handle = handle
        self.gamma_0 = gamma_0
        self.beta_0 = beta_0
        self.id_norm = np.array([0,0,-1], dtype=np.float64)
        self.norm = self.id_norm.copy()
        self.H = np.eye(4, dtype=np.float64)
        self.H_i = np.eye(4, dtype=np.float64)
        self.R = np.eye(3, dtype=np.float64)
        self.t = np.zeros(3, dtype=np.float64)
        # Need at least 2 rows and cols to get four corners of the wiper
        self.rows = max(2, rows)
        self.cols = max(2, cols)
        self.tot = self.rows * self.cols
        self.opt_compression = 0.1
        self.max_h = 0.2
        self.width = 0.1
        self.height = 0.1
        self.handle_offset = np.array([
            [0,0,-1,self.width],
            [0,1,0,0],
            [1,0,0,self.max_h],
            [0,0,0,1]], dtype=np.float64)
        self.dx = self.width / (self.cols - 1)
        self.dy = self.height / (self.rows - 1)
        # lam is (1 - 2v) / (2(1 - v)) where v is Poisson's ratio
        self.lam = lam
        self.Qy = None
        self.Qx = None
        self.init_Qy()
        self.init_Qx()
        self.id_top_points = np.zeros((self.tot, 4), dtype=np.float64)
        self.top_points = self.id_top_points.copy()
        self.init_top_points()
        self.ray_t = []
        self.opt_t = []
        self.clean_t = []

    def init_Qy(self):
        diag = 2 * np.ones(self.tot)
        for i in range(self.cols):
            diag[i] = 1
            diag[-i-1] = 1
        self.Qy = sparse.diags([diag], [0], shape=(self.tot, self.tot))
        self.Qy += sparse.diags([-1, -1], [-self.cols, self.cols],
            shape=(self.tot, self.tot))
        self.Qy *= 1/self.dy**2

    def init_Qx(self):
        diag = 2 * np.ones(self.tot)
        for i in range(0, self.tot, self.cols):
            diag[i] = 1
            if i > 0:
                diag[i - 1] = 1
        diag[-1] = 1
        off_diag = -np.ones(self.tot-1)
        for i in range(0, self.tot, self.cols):
            if i > 0:
                off_diag[i-1] = 0
        self.Qx = sparse.diags([diag], [0], shape=(self.tot, self.tot))
        self.Qx += sparse.diags([off_diag, off_diag], [-1, 1],
            shape=(self.tot, self.tot))
        self.Qx *= 1/self.dx**2

    def init_top_points(self):
        for i in range(self.rows):
            for j in range(self.cols):
                ind = i * self.cols + j
                self.id_top_points[ind, :] = np.array([
                    j * self.dx,
                    i * self.dy,
                    self.max_h,
                    1.0
                ])

    def wipe_step(self, start, end, ws):
        """
        """
        # Grab the second element (translation),
        # s = time.monotonic()
        move_vec = np.array(klampt.math.se3.error(end, start)[3:])
        self.setTransform(*start)
        start_cover = ws.get_covered_triangles()

        self.setTransform(*end)
        end_cover = ws.get_covered_triangles()

        avg_cover = (start_cover + end_cover) / 2
        # print(f'Coverage: {time.monotonic() - s}')
        # s = time.monotonic()
        dists = np.zeros(ws.num_triangles)
        for i, c in enumerate(avg_cover):
            if c > 0:
                v = move_vec
                n = ws.t_normals[i, :]
                dists[i] = np.linalg.norm(v - ((v.T @ n)/(n.T @ n)) * n)
        # print(f'Dists: {time.monotonic() - s}')
        return self.gamma(1.0, avg_cover, 1.0, dists)
        # ws.update_infection(self.gamma(1.0, avg_cover, 1.0, dists))
        # ws.update_colors()

    def eval_wipe_step(self, start, end, ws, start_cover=None):
        """
        """
        # Grab the second element (translation)
        move_vec = np.array(klampt.math.se3.error(end, start)[3:])
        if start_cover is None:
            self.setTransform(*start)
            start_cover = ws.get_covered_triangles()

        self.setTransform(*end)
        end_cover = ws.get_covered_triangles()
        dists, avg_cover = compute_dists(start_cover, end_cover,
            move_vec, ws.t_normals)
        return compute_gamma(self.gamma_0, self.beta_0, avg_cover, dists), end_cover

    def setTransform(self, R, t):
        """Update transform of the wiper volume, update the relevant
        geometry data.
        """
        self.obj.setTransform(R, t)
        self.H = np.array(math.se3.homogeneous((R,t)))
        self.H_i = np.array(math.se3.homogeneous(math.se3.inv((R,t))))
        self.R = self.H[:3,:3]
        self.t = self.H[:3, 3]
        handle_h = self.H @ self.handle_offset
        self.handle.setTransform(handle_h[:3, :3].T.flatten(), handle_h[:3, 3])
        self.norm = self.R @ self.id_norm
        self.top_points = (self.H @ self.id_top_points.T).T

    def getTransform(self):
        return self.obj.getTransform()

    def getDesiredTransform(self, pt, z_d, move_vec, theta):
        pt = np.array(pt)
        z_d = np.array(z_d)
        move_vec = np.array(move_vec)
        y_d = np.cross(z_d, move_vec)
        n_y = np.linalg.norm(y_d)
        if n_y == 0:
            # Move doesn't make sense, moving in/out of surface does nothing
            return None, None
        y_d = y_d / np.linalg.norm(y_d)
        z_d = z_d / np.linalg.norm(z_d)
        x_d = np.cross(y_d, z_d)
        x_d = x_d / np.linalg.norm(x_d)
        R_d = np.array([x_d, y_d, z_d]).T
        c = pmath.cos(theta)
        s = pmath.sin(theta)
        rot_z = np.array([
            [c,-s,0],
            [s,c,0],
            [0,0,1]])
        R_d = R_d @ rot_z
        offset = R_d @ (self.opt_compression * self.id_norm)
        t = pt + offset
        return R_d, t

    def gamma(self, s, f, w, d):
        """Compute the proportion of infection remaining for a wipe at speed s
        using force f where the wiper is w dirty. d is the distance the wiper
        moved over each triangle.

        Parameters
        ----------------
        """
        return compute_gamma(self.gamma_0, self.beta_0, f, d)

    def s_function(self, s):
        return 1

    def f_function(self, f):
        return 1

    def w_function(self, w):
        return 1


@jit(nopython=True)
def compute_dists(start_dict, end_dict, v, n):
    avg_dict = numba.typed.Dict.empty(numba.types.int64, numba.types.float64)
    dists = numba.typed.Dict.empty(numba.types.int64, numba.types.float64)
    for key in start_dict:
        if key in end_dict:
            avg_dict[key] = (start_dict[key] + end_dict[key]) / 2
        else:
            avg_dict[key] = start_dict[key] / 2
    for key in end_dict:
        if key not in start_dict:
            avg_dict[key] = end_dict[key] / 2
    for key in avg_dict:
        dists[key] = np.linalg.norm(v
            - ((n[key, :] @ v) / (n[key, :] @ n[key, :])) * n[key, :])
    return dists, avg_dict


@jit(nopython=True)
def compute_gamma(gamma_0, beta_0, cover, dists):
    gamma = numba.typed.Dict.empty(numba.types.int64, numba.types.float64)
    for key in cover:
        gamma[key] = (1 - gamma_0 * cover[key]) ** (beta_0 * dists[key])
    return gamma


if __name__ == "__main__":
    main()
