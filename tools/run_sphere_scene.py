import subprocess
import sys
from typing import Dict


def run_sim(setting_file: str) -> None:

    command = f'bash -c "source activate root; conda activate fpm; fpm @runmp -- sphere.toml"'
    subprocess.run(command, check=True, shell=True)


def make_settings(user_dict: Dict, file: str) -> None:
    source = {"name": "uniform", "nphotons": 1000000, "directions": "-z"}
    grid = {"nxg": 200, "nyg": 200, "nzg": 200, "xmax": 1.0, "ymax": 1.0, "zmax": 1.0}
    geometry = {"geom_name": "sphere_scene", "num_spheres": 10}
    output = {"fluence": f"sphere_scene_{geometry['num_spheres']}.nrrd", "render_geom": "false"}
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


    for num_sphere in num_spheres:
        user_dict = {"num_spheres": num_sphere}
        make_settings(user_dict, "sphere.toml")
        start = time.time()
        run_sim("")
        end = time.time()
        times.append([num_sphere, end - start])
with open("times.dat", "w") as file:
    for time in times:
        file.write(f"{time[0]}, {time[1]:.2f}" + "\n")
