import subprocess
import sys
from typing import Dict

def run_sim(setting_file: str) -> None:

    command = f'bash -c "source activate root; conda activate fpm; fpm @runmp -- sphere.toml"'
    subprocess.run(command, check=True, shell=True)

def make_settings(user_dict: Dict, file: str) -> None:
    source = {"name": "uniform", "nphotons": 1000000, "directions": "-z"}
    grid = {"nxg":200, "nyg":200, "nzg":200, "xmax":1.0, "ymax":1.0, "zmax":1.0}
    geometry = {"geom_name": "sphere_scene", "num_spheres":10}
    output = {"fluence":f"sphere_scene_{geometry['num_spheres']}.nrrd", "render_geom": "false"}
    simulation = {"iseed": 123456789, "tev": "false"}

    defaults = {"source": source, "grid": grid, "geometry": geometry, "output": output, "simulation": simulation}
    # check all user keys are valid
    for user_key in user_dict:
        for key in defaults:
            if user_key in defaults[key]:
                defaults[key][user_key] = user_dict[user_key]
                if user_key == "num_spheres":
                    defaults["output"]["fluence"] = f"sphere_scene_{user_dict[user_key]}.nrrd"
    
    # write out setting file
    with open("res/" + file, "w") as f:
        for key in defaults:
            f.write(f"[{key}]\n")
            for subkey in defaults[key]:
                output = defaults[key][subkey]
                if type(output) == str:
                    if output == "false" or output == "true":
                        f.write(f"{subkey}={output}\n")
                    else:
                        f.write(f"{subkey}='{output}'\n")
                else:
                    f.write(f"{subkey}={output}\n")
            f.write("\n")


def squared(x, c):
    return c*x**2


def power(x, c):
    return x**c


def linear(x, m, c):
    return m*x+c


def get_r2(xdata, ydata, f):
    res = ydata - f(xdata, *popt)
    ss_res = np.sum(res**2)
    ss_tot = np.sum((ydata - np.mean(ydata))**2)
    return 1. - (ss_res / ss_tot)


if __name__ == "__main__":
    import time
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.optimize import curve_fit

    times = [[100, 30.874497413635254], [200, 67.9496521949768], [300, 119.79060006141663], [400, 184.2490565776825], [500, 276.18548679351807], [600, 302.8419463634491], [700, 560.928364276886], [800, 496.095983505249], [900, 644.8481585979462], [1000, 739.9614040851593], [1100, 807.1722948551178], [1200, 898.0271279811859], [1300, 882.4315783977509], [1400, 1139.8273587226868], [1500, 1354.6854293346405], [1600, 1380.4796850681305], [1800, 1793.0182490348816], [1900, 1926.0771131515503], [4000,9586.09]]

    # x, y = zip(*times)
    # x = np.array(x)
    # y = np.array(y) - y[0]
    # plt.scatter(x, y)
    # plt.xlabel("num of spheres")
    # plt.ylabel("time/s")

    # xdata = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1800, 1900, 4000]
    # xx = np.linspace(100, 5000, 1900)

    # funcs = [linear, squared, power]
    # for func in funcs:
    #     popt, pconv = curve_fit(func, xdata, y)
    #     plt.plot(xx, func(xx, *popt), label=func.__name__ + f", r^2={get_r2(x, y, func):.2f}")
    #     print(func.__name__,*popt)
    # plt.legend()

    # plt.show()

    num_spheres = [100, 500, 1000, 1500, 2000, 4000]
    for num_sphere in num_spheres:
        user_dict = {"num_spheres": num_sphere}
        make_settings(user_dict, "sphere.toml")
        print(f"running: {num_sphere}")
        start = time.time()
        run_sim("")
        end = time.time()
        times.append([num_sphere, end - start])
with open("times.dat", "w") as file:
    for time in times:
        file.write(f"{time[0]}, {time[1]:.2f}" + "\n")
