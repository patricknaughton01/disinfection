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
faulthandler.enable()

from klampt import vis
from klampt import math
from klampt import io
from klampt.vis import colorize
from klampt.model.create import primitives

world = klampt.WorldModel()
DEBUG = False
TIMING = False
DISPLAY = True

def main():
    global world
    obj = world.makeRigidObject("tm")
    oort = 1/(2**0.5)
    # obj.setTransform([1, 0, 0, 0, -oort, oort, 0, -oort, -oort], [0,0.2,0])
    g = obj.geometry()
    g.loadFile("keys.off")
    ws = WipeSurface("tm", obj)

    wiper_obj = world.makeRigidObject("wiper")
    wiper_obj.geometry().loadFile("wiper.off")
    wiper_obj.geometry().set(wiper_obj.geometry().convert("VolumeGrid"))
    wiper_obj.appearance().setColor(0.5,0.5,0.5,0.2)

    dt = 0.1
    step = 0.01
    max = 1.0
    min = -0.1
    state = "left"

    sizes = [30, 30, 50]
    RUNS = 1
    if TIMING:
        RUNS = 11
    for i in range(1):#len(sizes)):
        wiper = Wiper(wiper_obj, rows=sizes[i], cols=sizes[i], lam=0.5)
        wiper.setTransform([1,0,0,0,1,0,0,0,1], [0.01,0.01,0.0])
        # wiper.setTransform([1,0,0,0,1,0,0,0,1], [0.0,0.0,-0.0])
        # wiper.setTransform([0.8678192,  0.0000000, -0.4968801,
        #                0.0000000,  1.0000000,  0.0000000,
        #                0.4968801,  0.0000000,  0.8678192 ],
        #     [0.01,0.01,0.1])

        covered = ws.get_covered_triangles(wiper)
        if TIMING:
            start_time = time.monotonic()
            for j in range(RUNS):
                covered = ws.get_covered_triangles(wiper)
                t = np.random.rand(3) * 0.9
                t[2] = -0.05
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
    ws.update_infection(1-covered)
    ws.update_colors()
    # wiper = Wiper(wiper_obj, rows=10, cols=10, lam=100, beta_0=10, gamma_0=0.8)
    # wiper.setTransform([1,0,0,0,1,0,0,0,1], [0.01,0.01,0.05])
    if DISPLAY:
        vis.add("world", world)
        # vis.getViewport().setTransform(([0, -1, 0,
        #     -oort, 0, -oort,
        #     oort, 0, -oort], [-1.5, 0.5, 2]))
        vis.show()
        # time.sleep(3)
        while vis.shown():
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


class Planner:
    def __init__(self, surface, wiper, mode='online', plan_type='pfield'):
        self.surface = surface
        self.wiper = wiper
        self.plan_type = plan_type
        self.mode = mode
        self.wipes = []
        self.offsets = None

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
    def __init__(self, name, obj):
        self.name = name
        self.obj = obj
        self.tm = self.obj.geometry().getTriangleMesh()
        self.verts = np.array(self.tm.vertices).reshape(
            (len(self.tm.vertices)//3,3))
        self.inds = np.array(self.tm.indices,dtype=np.int32).reshape(
            (len(self.tm.indices)//3,3))
        self.num_triangles = len(self.inds)
        self.v_map = None
        self.t_neighbors = None
        self.build_triangle_adjacency()
        self.t_normals = None
        self.build_triangle_normals()
        self.infection_level = np.ones(self.num_triangles) + np.random.rand(self.num_triangles)/10
        self.infection_level[0] = 0
        self.obj.appearance().setColor(0.0, 0.0, 1.0)
        self.update_colors()
        self.covered = np.zeros(self.num_triangles)
        self.hit_dist_thresh = 0
        self.epsilon = 1e-5
        self.visited_triangles = []

    def build_triangle_adjacency(self):
        self.v_map = []
        for i in range(len(self.verts)):
            self.v_map.append([])
        for i in range(self.num_triangles):
            for j in range(len(self.inds[i])):
                self.v_map[self.inds[i][j]].append(i)
        self.t_neighbors = -np.ones((self.inds.shape), dtype=np.int64)
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
        self.t_normals = np.empty((self.num_triangles, 3))
        sum = 0
        for i in range(self.num_triangles):
            a = self.verts[self.inds[i][0], :]
            b = self.verts[self.inds[i][1], :]
            c = self.verts[self.inds[i][2], :]
            v1 = b - a
            v2 = c - a
            cx = np.cross(v1, v2)
            sum += ((a[0] + b[0] + c[0])/3) * cx[0] * np.linalg.norm(cx)/2
            self.t_normals[i, :] = cx / np.linalg.norm(cx)
        if sum < 0:
            self.t_normals = -self.t_normals
        # for i in range(0, self.num_triangles, self.num_triangles//1000):
        #     pt = np.mean(self.verts[self.inds[i, :], :], axis=0) + self.t_normals[i, :]/100
        #     primitives.sphere(0.002, pt, world=world).appearance().setColor(1,0,1)

    def get_covered_triangles(self, wiper):
        global world, DEBUG
        """Get the indices of the covered triangles by a wipe represented by
        the passed in wiper. A triangle is considered covered if all of
        its vertices are inside the rigidObject.

        REMARK: Implementation/speedup: If the SUW is much larger than the
        wiper, it makes sense to store this as a list of covered triangles
        rather than an array with an element for each triangle denoting whether
        it was covered.
        """
        ray_s = time.monotonic()
        min_h = np.zeros(wiper.tot)
        # Store index of the triangle that this point raycasts to
        h_t_correspondence = -np.ones(wiper.tot, dtype=np.long)
        for i in range(wiper.tot):
            start_pt = wiper.top_points[i][:3]
            hit, pt = self.obj.geometry().rayCast_ext(
                start_pt, wiper.norm)
            if hit >= 0:
                pt = np.array(pt)
                if DEBUG:
                    primitives.sphere(0.001, pt, world=world).appearance().setColor(0,1,0)
                min_h[i] = wiper.max_h - np.linalg.norm(start_pt - pt)
                if min_h[i] >= 0:
                    h_t_correspondence[i] = hit
        min_h = np.maximum(min_h, 0)
        ray_e = time.monotonic()
        # Solve opt problem to find contacted triangles
        opt_s = time.monotonic()
        h = cp.Variable(wiper.tot)
        objective = cp.Minimize(cp.sum_squares(h / wiper.max_h)
            + wiper.lam * (cp.quad_form(h, wiper.Qy)
            + cp.quad_form(h, wiper.Qx)))
        constraints = [h <= wiper.max_h, h >= min_h]
        prob = cp.Problem(objective, constraints)
        result = prob.solve()
        h_val = h.value
        opt_e = time.monotonic()

        if DEBUG:
            for i in range(wiper.tot):
                pt = wiper.top_points[i][:3] + (wiper.max_h - h_val[i]) * wiper.norm
                primitives.sphere(0.002, pt, world=world).appearance().setColor(1,0,1)

        clean_s = time.monotonic()
        contact = np.abs(h_val - min_h) < 1e-3
        covered_triangles = np.zeros(self.num_triangles)
        visited = np.zeros(self.num_triangles)
        for i in range(wiper.tot):
            tind = h_t_correspondence[i]
            if tind > -1:
                self.interpolate_contact(tind, wiper, visited, contact, h_val,
                    covered_triangles)
        clean_e = time.monotonic()
        if TIMING:
            wiper.ray_t.append(ray_e - ray_s)
            wiper.opt_t.append(opt_e - opt_s)
            wiper.clean_t.append(clean_e - clean_s)
        self.covered = covered_triangles
        return self.covered

    def interpolate_contact(self, triangle, wiper, visited, grid, heights,
        cover
    ):
        if visited[triangle]:
            return
        visited[triangle] = 1
        centroid = np.zeros(4)
        centroid[:3] = np.mean(self.verts[self.inds[triangle, :], :], axis=0)
        centroid[3] = 1
        hit, pt = self.obj.geometry().rayCast(centroid[:3], -wiper.norm)
        if (hit and 1e-5 < np.linalg.norm(np.array(pt) - centroid[:3])
            < wiper.max_h
        ):
            # Line of sight to wiper is blocked --> can't be in contact
            return
        proj_c = wiper.H_i @ centroid
        if (0 < proj_c[0] < wiper.width and 0 < proj_c[1] < wiper.height
            and 0 <= proj_c[2] <= wiper.max_h
        ):
            # Bilinearly interpolate between the four h values around this pt
            # g's formatted as [col, row] or [x, y]
            #  g2 x  x g4
            #  g1 x  x g3
            g1 = [proj_c[0] // wiper.dx, proj_c[1] // wiper.dy]
            g2 = [g1[0], g1[1] + 1]
            g3 = [g1[0] + 1, g1[1]]
            g4 = [g3[0], g3[1] + 1]
            gs = [g1,g2,g3,g4]
            # How far is p from g1 in x and y normalized by distances between
            # adjacent grid points
            mx = (proj_c[0] % wiper.dx) / wiper.dx
            my = (proj_c[1] % wiper.dy) / wiper.dy
            vs = [grid[wiper.flatten(g)] for g in gs]
            i1 = (1 - mx) * vs[0] + mx * vs[2]
            i2 = (1 - mx) * vs[1] + mx * vs[3]
            # Exponential fall off for points far from the heights
            var = 2e-5/(1e-8 + wiper.lam)
            avg_h = np.mean([heights[wiper.flatten(g)] for g in gs])
            cover[triangle] = (((1 - my) * i1 + my * i2)
                * np.exp(-0.5 * (proj_c[2] - avg_h)**2 / var)
                * max(0, (-wiper.norm.T @ self.t_normals[triangle])))
        else:   # Outside the confines of the wiper
            return
        for i in range(len(self.t_neighbors[triangle])):
            self.interpolate_contact(self.t_neighbors[triangle][i], wiper,
                visited, grid, heights, cover)

    def clear_covered(self):
        self.covered = np.zeros(self.num_triangles)

    def update_infection(self, gamma):
        """Update the infection_level of the triangles in this mesh.

        Parameters
        ----------------
        gamma: np.ndarray (self.num_triangles,)
            for each triangle, the proportion of infection remaining
        """
        self.infection_level *= gamma

    def update_colors(self):
        colorize.colorize(self.obj, self.infection_level,
            colormap="jet", feature="faces")


class Wiper:
    def __init__(self, obj, gamma_0=0.5, beta_0=1.0, rows=10,
        cols=10, lam=0
    ):
        self.obj = obj
        self.gamma_0 = gamma_0
        self.beta_0 = beta_0
        self.id_norm = np.array([0,0,-1])
        self.norm = self.id_norm.copy()
        self.H = np.eye(4)
        self.H_i = np.eye(4)
        self.R = np.eye(3)
        self.t = np.zeros(3)
        # Need at least 2 rows and cols to get four corners of the wiper
        self.rows = max(2, rows)
        self.cols = max(2, cols)
        self.tot = self.rows * self.cols
        self.max_h = 0.2
        self.width = 0.1
        self.height = 0.1
        self.dx = self.width / (self.cols - 1)
        self.dy = self.height / (self.rows - 1)
        # lam is (1 - 2v) / (2(1 - v)) where v is Poisson's ratio
        self.lam = lam
        self.Qy = None
        self.Qx = None
        self.init_Qy()
        self.init_Qx()
        self.id_top_points = np.zeros((self.tot, 4))
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
        start_cover = ws.get_covered_triangles(self)

        self.setTransform(*end)
        end_cover = ws.get_covered_triangles(self)

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

    def setTransform(self, R, t):
        """Update transform of the wiper volume, update the relevant
        geometry data.
        """
        self.obj.setTransform(R, t)
        self.H = np.array(math.se3.homogeneous((R,t)))
        self.H_i = np.array(math.se3.homogeneous(math.se3.inv((R,t))))
        self.R = self.H[:3,:3]
        self.t = self.H[:3, 3]
        self.norm = self.R @ self.id_norm
        for i in range(self.tot):
            self.top_points[i, :] = (self.H @ self.id_top_points[i, :])

    def getTransform(self):
        return self.obj.getTransform()

    def get_dists(self, ws, move_vec, covered_triangles):
        """Compute distance from vertex to edge of wiper
        """
        dists = np.zeros(ws.num_triangles)
        for i, c in enumerate(covered_triangles):
            if c > 0:
                cen = np.mean(self.verts[self.inds[i, :], :], axis=0)
                hit, pt = self.obj.geometry().rayCast(cen, -move_vec)
                if hit:
                    v = np.array(pt) - cen
                    n = ws.t_normals[i, :]
                    dists[i] = np.linalg.norm(v - ((v.T @ n)/(n.T @ n)) * n)
        return dists

    def flatten(self, pt):
        """Return flat index into a grid with self.rows rows and self.cols
        columns for a point pt in format [col, row]
        """
        return int(pt[0] + self.cols * pt[1])

    def gamma(self, s, f, w, d):
        """Compute the proportion of infection remaining for a wipe at speed s
        using force f where the wiper is w dirty. d is the distance the wiper
        moved over each triangle.

        Parameters
        ----------------
        s: float
        f: float
        w: float
        d: np.ndarray
        """
        return (self.s_function(s)
            * (1 - self.gamma_0 * f) * self.w_function(w)) ** (self.beta_0 * d)

    def s_function(self, s):
        return 1

    def f_function(self, f):
        return 1

    def w_function(self, w):
        return 1

if __name__ == "__main__":
    main()
