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
TIMING = True
DISPLAY = True

def main():
    global world
    obj = world.makeRigidObject("tm")
    g = obj.geometry()
    g.loadFile("lumps.off")
    ws = WipeSurface("tm", obj)
    oort = 1/(2**0.5)
    #ws.obj.setTransform([1, 0, 0, 0, -oort, oort, 0, -oort, -oort], [0,0.2,0])

    wiper_obj = world.makeRigidObject("wiper")
    wiper_obj.geometry().loadFile("wiper.off")
    wiper_obj.geometry().set(wiper_obj.geometry().convert("VolumeGrid"))
    wiper_obj.appearance().setColor(0.5,0.5,0.5,0.2)


    # sim = klampt.Simulator(world)
    # sim.setGravity((0,0,0))
    dt = 0.1
    step = 0.01
    max = 1.0
    min = -0.1
    state = "left"

    # ws_body = sim.body(world.rigidObject(obj.getID()))
    # ws_body.enableDynamics(False)
    #
    # wiper_body = sim.body(world.rigidObject(wiper_obj.getID()))
    # wiper_body.enableDynamics(False)
    # wiper_body.setVelocity([0,0,0], [0,0.0,0])

    sizes = [10, 30, 50]
    RUNS = 1
    if TIMING:
        RUNS = 11
    for i in range(len(sizes)):
        wiper = Wiper(wiper_obj, rows=sizes[i], cols=sizes[i], lam=100)
        #wiper.setTransform([1,0,0,0,1,0,0,0,1], [0.0,0.0,-0.0])
        wiper.setTransform([1,0,0,0,1,0,0,0,1], [0.01,0.01,-0.05])
        # wiper.setTransform([0.8678192,  0.0000000, -0.4968801,
        #                0.0000000,  1.0000000,  0.0000000,
        #                0.4968801,  0.0000000,  0.8678192 ],
        #     [0.01,0.01,0.1])

        # covered = ws.get_covered_triangles(wiper)[0]
        start_time = time.monotonic()
        for j in range(RUNS):
            covered = ws.get_covered_triangles(wiper)[0]
            t = np.random.rand(3) * 0.9
            t[2] = -0.05
            wiper.obj.setTransform([1,0,0,0,1,0,0,0,1], t)
        end_time = time.monotonic()
        if TIMING:
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
    #wiper.wipe_step(wiper.obj.getTransform(), ([1,0,0,0,1,0,0,0,1], [0.0,0.0,0.95]), ws)
    wiper.setTransform([1,0,0,0,1,0,0,0,1], [-0.1,-0.1,0.05])
    #wiper_obj.setTransform([1,0,0,0,1,0,0,0,1], [0.2,0,0])
    if DISPLAY:
        vis.add("world", world)
        vis.show()
        while vis.shown():
            # sim.simulate(dt)
            # R, t = wiper.obj.getTransform()
            # if state == "left":
            #     if t[1] < max:
            #         move_t = list(t[:])
            #         move_t[1] += step
            #         wiper.wipe_step((R,t), (R, move_t), ws)
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
        # st = time.monotonic()
        self.build_triangle_adjacency()
        # et = time.monotonic()
        # print(f"Building triangle adjacency took {et - st} seconds")
        self.infection_level = np.ones(self.num_triangles)
        #self.infection_level[0] = 0
        self.obj.appearance().setColor(0.0, 0.0, 1.0)
        self.update_colors()
        self.covered = np.zeros(self.num_triangles)
        self.hit_dist_thresh = 0
        self.epsilon = 1e-5

    def build_triangle_adjacency(self):
        # vs = time.monotonic()
        self.v_map = []
        for i in range(len(self.verts)):
            self.v_map.append([])
        for i in range(self.num_triangles):
            for j in range(len(self.inds[i])):
                self.v_map[self.inds[i][j]].append(i)
        # ve = time.monotonic()
        # print(f"Vertices: {ve - vs}")
        # ts = time.monotonic()
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
            # assert ind==3, f"Only found {ind} neighbors for triangle {i}"
        # te = time.monotonic()
        # print(f"Triangles: {te - ts}")

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
        objective = cp.Minimize(cp.sum_squares(h)
            + wiper.lam * (cp.quad_form(h, wiper.Q)))
        constraints = [h <= wiper.max_h, h >= min_h]
        prob = cp.Problem(objective, constraints)
        result = prob.solve()
        h_val = h.value
        opt_e = time.monotonic()

        if DEBUG:
            for i in range(wiper.tot):
                if h_val[i] - min_h[i] < 1e-3:
                    pt = wiper.top_points[i][:3] + (wiper.max_h - h_val[i]) * wiper.norm
                    primitives.sphere(0.002, pt, world=world).appearance().setColor(1,0,1)

        clean_s = time.monotonic()
        contact = np.abs(h_val - min_h) < 1e-3
        covered_vertices = np.zeros(len(self.verts))
        covered_triangles = np.zeros(self.num_triangles)
        visited = np.zeros(self.num_triangles)
        for i in range(wiper.tot):
            tind = h_t_correspondence[i]
            if tind > -1:
                self.interpolate_contact(tind, wiper, visited, contact,
                    covered_triangles)
        clean_e = time.monotonic()
        if TIMING:
            wiper.ray_t.append(ray_e - ray_s)
            wiper.opt_t.append(opt_e - opt_s)
            wiper.clean_t.append(clean_e - clean_s)
        self.covered = covered_triangles
        return self.covered, covered_vertices

    def interpolate_contact(self, triangle, wiper, visited, grid, cover):
        if visited[triangle]:
            return
        visited[triangle] = 1
        centroid = np.zeros(4)
        centroid[:3] = np.mean(self.verts[self.inds[triangle]], axis=0)
        centroid[3] = 1
        hit, pt = self.obj.geometry().rayCast(centroid[:3], -wiper.norm)
        if hit and np.linalg.norm(np.array(pt) - centroid[:3]) > 1e-8:
            # Line of sight to wiper is blocked --> can't be in contact
            return
        proj_c = wiper.H_i @ centroid
        if 0 <= proj_c[0] <= wiper.width and 0 <= proj_c[1] <= wiper.height:
            # Bilinearly interpolate between the four h values around this pt
            dx = wiper.width / max(1, wiper.cols - 1)
            dy = wiper.height / max(1, wiper.rows - 1)
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
            vs = [grid[wiper.flatten(g)] for g in gs]
            i1 = (1 - mx) * vs[0] + mx * vs[2]
            i2 = (1 - mx) * vs[1] + mx * vs[3]
            cover[triangle] = (1 - my) * i1 + my * i2
        else:   # Outside the confines of the wiper
            return
        for i in range(len(self.t_neighbors[triangle])):
            self.interpolate_contact(self.t_neighbors[triangle][i], wiper,
                visited, grid, cover)

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
    def __init__(self, obj, compliance=0.01, gamma_0=0.5, beta_0=1.0, rows=10,
        cols=10, lam=0
    ):
        self.obj = obj
        self.compliance = compliance
        self.gamma_0 = gamma_0
        self.beta_0 = beta_0
        self.id_norm = np.array([0,0,-1])
        self.norm = self.id_norm.copy()
        self.H = np.eye(4)
        self.H_i = np.eye(4)
        self.R = np.eye(3)
        self.t = np.zeros(3)
        self.rows = rows
        self.cols = cols
        self.tot = self.rows * self.cols
        self.Q = None
        self.init_Q()
        # plt.spy(self.Q)
        # plt.show()
        self.max_h = 0.2
        self.width = 0.1
        self.height = 0.1
        self.lam = lam
        self.id_top_points = np.zeros((self.tot, 4))
        self.top_points = self.id_top_points.copy()
        self.init_top_points()
        self.ray_t = []
        self.opt_t = []
        self.clean_t = []

    def init_Q(self):
        # Border terms
        self.Q = sparse.diags([-1,-1], [-self.cols, self.cols],
            shape=(self.tot, self.tot))
        # Handle beginnings and endings of rows not having left/right neighbors
        off_diag = -np.ones(self.tot-1)
        for i in range(0, self.tot, self.cols):
            if i > 0:
                off_diag[i-1] = 0
        self.Q += sparse.diags([off_diag, off_diag], [-1, 1], shape=(self.tot, self.tot))
        # Squared terms
        self.Q += sparse.diags([4], [0], shape=(self.tot, self.tot))
        # Edges
        for i in range(self.cols):
            # Top row
            self.Q[i, i] = 3
            # Bottom row
            ind = (self.rows - 1) * self.cols + i
            self.Q[ind, ind] = 3
        for i in range(self.rows):
            # First col
            ind = (i * self.cols)
            self.Q[ind, ind] = 3
            # Last col
            ind = (i * self.cols) + self.cols - 1
            self.Q[ind, ind] = 3
        # Corners
        self.Q[0, 0] = 2                            # TL
        self.Q[self.cols - 1, self.cols - 1] = 2    # TR
        ind = (self.rows - 1) * self.cols           # BL
        self.Q[ind, ind] = 2
        ind = self.rows * self.cols - 1             # BR
        self.Q[ind, ind] = 2

    def init_top_points(self):
        for i in range(self.rows):
            for j in range(self.cols):
                ind = i * self.cols + j
                dr = max(1, self.rows - 1)
                dc = max(1, self.cols - 1)
                self.id_top_points[ind, :] = np.array([
                    self.width * j / dc,
                    self.height * i / dr,
                    self.max_h,
                    1.0
                ])

    def wipe_step(self, start, end, ws):
        """
        """
        # Grab the second element (translation),
        move_vec = np.array(klampt.math.se3.error(end, start)[3:])
        self.obj.setTransform(*start)
        start_cover, start_cover_verts = ws.get_covered_triangles(self)
        start_vert_dists = self.get_vert_dists(ws, move_vec, start_cover_verts)
        ws.clear_covered()

        self.obj.setTransform(*end)
        end_cover, end_cover_verts = ws.get_covered_triangles(self)
        end_vert_dists = self.get_vert_dists(ws, -move_vec, end_cover_verts)
        ws.clear_covered()

        both_cover = np.logical_and(start_cover, end_cover)
        d = np.zeros(ws.num_triangles)
        d += both_cover * math.se3.distance(start, end, Rweight=0.0)

        # (S U E) - B
        # u_minus_i = np.logical_and(np.logical_or(start_cover, end_cover),
        #     np.logical_not(both_cover))
        # u_minus_i_d = start_vert_dists + end_vert_dists
        # u_minus_i_d = u_minus_i_d[u_minus_i]
        # Abstract this into function? Or can I just deal with the union?
        # d[u_minus_i] = np.mean()
        for i, c in enumerate(start_cover):
            if c and not both_cover[i]:
                avg_dist = np.mean(start_vert_dists[ws.inds[i, :]])
                d[i] = avg_dist
        for i, c in enumerate(end_cover):
            if c and not both_cover[i]:
                avg_dist = np.mean(end_vert_dists[ws.inds[i, :]])
                d[i] = avg_dist
        ws.update_infection(self.gamma(1.0, 1.0, 1.0, d))
        ws.update_colors()

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

    def get_vert_dists(self, ws, move_vec, covered_vertices):
        """Compute distance from vertex to edge of wiper
        """
        vert_dists = np.zeros(ws.verts.shape[0])
        for i, c in enumerate(covered_vertices):
            if c:
                hit, pt = self.obj.geometry().rayCast(ws.verts[i], -move_vec)
                if hit:
                    vert_dists = np.linalg.norm(
                        np.array(pt) - ws.verts[i, :])
        return vert_dists

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
        return (self.gamma_0 * self.s_function(s)
            * self.f_function(f) * self.w_function(w)) ** (self.beta_0 * d)

    def s_function(self, s):
        return 1

    def f_function(self, f):
        return 1

    def w_function(self, w):
        return 1

if __name__ == "__main__":
    main()
    vis.kill()
