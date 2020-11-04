import disinfection
import klampt
import time
import numpy as np

from klampt import vis

world = klampt.WorldModel()


def main():
    global world
    obj = world.makeRigidObject("tm")
    oort = 1/(2**0.5)
    # obj.setTransform([1, 0, 0, 0, -oort, oort, 0, -oort, -oort], [0,0.2,0])
    g = obj.geometry()
    g.loadFile("keys.off")
    ws = disinfection.WipeSurface("tm", obj)

    wiper_obj = world.makeRigidObject("wiper")
    wiper_obj.geometry().loadFile("wiper.off")
    wiper_obj.geometry().set(wiper_obj.geometry().convert("VolumeGrid"))
    wiper_obj.appearance().setColor(0.5,0.5,0.5,0.2)

    dt = 0.1
    size = 10
    max_wipe = 3000
    wipe_ind = 0
    wiper = disinfection.Wiper(wiper_obj, rows=size, cols=size, lam=0.5,
        gamma_0=0.8, beta_0=10)

    planner = disinfection.Planner(ws, wiper)
    planner.get_wipe()

    vis.add("world", world)
    vis.getViewport().setTransform(([0, -1, 0, -1, 0, 0, 0, 0, -1],
        [0.5, 0.5, 2.8]))
    vis.show()
    p_times = []
    s = time.monotonic()
    while vis.shown():
        if wipe_ind < max_wipe:
            print(f"{wipe_ind}/{max_wipe-1}")
            print(wiper.getTransform()[1])
            ps = time.monotonic()
            best_wipe = planner.get_wipe()
            p_times.append(time.monotonic() - ps)
            gamma = wiper.wipe_step(best_wipe[0], best_wipe[1], ws)
            ws.update_infection(gamma)
            ws.update_colors()
            wipe_ind += 1
        elif wipe_ind == max_wipe:
            e = time.monotonic()
            total_time = e - s
            print(f"Average planning time: {np.mean(p_times)}")
            print(f"Std. dev planning time: {np.std(p_times)}")
            print(f"Total time: {total_time}")
            print(f"Average time per iteration: {total_time/max_wipe}")
            wipe_ind += 1
        time.sleep(dt)
    vis.show(False)
    vis.kill()


if __name__ == "__main__":
    main()
