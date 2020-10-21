import klampt
import time
import sys
import numpy as np
import cvxpy as cp
import faulthandler
import matplotlib.pyplot as plt
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
DISPLAY = False

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
    #wiper_obj.setTransform([1,0,0,0,1,0,0,0,1], [0.0,0.0,-0.0])
    wiper_obj.setTransform([1,0,0,0,1,0,0,0,1], [0.01,0.01,-0.05])

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

    sizes = [30, 30, 50]
    RUNS = 1
    if TIMING:
        RUNS = 11
    for i in range(1):#len(sizes)):
        wiper = Wiper(wiper_obj, rows=sizes[i], cols=sizes[i], lam=300)
        start_time = time.monotonic()
        for j in range(RUNS):
            covered = ws.get_covered_triangles(wiper)[0]
            t = np.random.rand(3) * 0.9
            t[2] = -0.05
            wiper.obj.setTransform([1,0,0,0,1,0,0,0,1], t)
        end_time = time.monotonic()
        if TIMING:
            print(len(wiper.ray_t))
            print("\t\t".join(("Ray", "Opt")))
            print("\t".join(("Mean", "Std", "Mean", "Std")))
            print("\t".join((f"{np.mean(wiper.ray_t[1:]):.4g}",
                f"{np.std(wiper.ray_t[1:]):.4g}",
                f"{np.mean(wiper.opt_t[1:]):.4g}",
                f"{np.std(wiper.opt_t[1:]):.4g}")))
            print("Average total time: {:.4g}".format( (end_time - start_time)/RUNS ))
    ws.update_infection(np.logical_not(covered).astype(np.float32))
    ws.update_colors()
    #wiper.wipe_step(wiper.obj.getTransform(), ([1,0,0,0,1,0,0,0,1], [0.0,0.0,0.95]), ws)
    wiper_obj.setTransform([1,0,0,0,1,0,0,0,1], [-0.1,-0.1,0.05])
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
        # self.build_triangle_adjacency()
        self.infection_level = np.ones(self.num_triangles)
        #self.infection_level[0] = 0
        self.obj.appearance().setColor(0.0, 0.0, 1.0)
        self.update_colors()
        self.covered = np.zeros(self.num_triangles)
        self.hit_dist_thresh = 0
        self.epsilon = 1e-5

    def build_triangle_adjacency(self):
        self.v_map = {}
        for i in range(len(self.verts)):
            self.v_map[i] = np.nonzero(np.sum(self.inds == i, axis=1))
        self.t_neighbors = -np.ones((self.inds.shape))
        for i in range(self.num_triangles):
            count = {}
            for j in range(len(self.inds[i])):
                for k, id in enumerate(self.v_map[j]):
                    if id in count:
                        count[id] += 1
                    else:
                        count[id] = 1
            print(count)
            break

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
        _, t = wiper.obj.getTransform()
        for i in range(wiper.tot):
            start_pt = np.array(t) + wiper.top_point_offsets[i]
            hit, pt = self.obj.geometry().rayCast_ext(
                start_pt, wiper.norm)
            if hit >= 0:
                pt = np.array(pt)
                if DEBUG:
                    primitives.sphere(0.001, pt, world=world).appearance().setColor(0,1,0)
                min_h[i] = np.max(wiper.max_h - np.linalg.norm(start_pt - pt), 0)
                h_t_correspondence[i] = hit
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
                pt = np.array(t) + wiper.top_point_offsets[i] + (wiper.max_h - h_val[i]) * wiper.norm
                primitives.sphere(0.002, pt, world=world).appearance().setColor(1,0,1)

        clean_s = time.monotonic()
        contact = np.abs(h_val - min_h) < 1e-3
        covered_vertices = np.zeros(len(self.verts))
        for i in range(wiper.tot):
            tind = h_t_correspondence[i]
            if contact[i] and tind > -1:
                covered_vertices[self.inds[tind, :]] = 1
        covered_triangles = np.sum(covered_vertices[self.inds], axis=1) >= 3
        clean_e = time.monotonic()
        if TIMING:
            wiper.ray_t.append(ray_e - ray_s)
            wiper.opt_t.append(opt_e - opt_s)
        self.covered = covered_triangles
        return self.covered, covered_vertices

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
        # TODO: This needs to change based on the configuration of the gripper,
        # for now, it's always pointing down.
        # Outward facing normal of the surface of the wiper in the EE frame.
        self.norm = np.array([0,0,-1])
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
        self.top_point_offsets = np.zeros((self.tot, 3))
        for i in range(self.rows):
            for j in range(self.cols):
                ind = i * self.cols + j
                dr = max(1, self.rows - 1)
                dc = max(1, self.cols - 1)
                self.top_point_offsets[ind, :] = np.array([
                    self.height* i / dr,
                    self.width * j / dc,
                    self.max_h
                ])
        self.ray_t = []
        self.opt_t = []

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

    def get_vert_dists(self, ws, move_vec, covered_vertices):
        vert_dists = np.zeros(ws.verts.shape[0])
        for i, c in enumerate(covered_vertices):
            if c:
                hit, pt = self.obj.geometry().rayCast(ws.verts[i], -move_vec)
                if hit:
                    vert_dists = np.linalg.norm(
                        np.array(pt) - ws.verts[i, :])
        return vert_dists

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
