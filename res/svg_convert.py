from pysvg.parser import parse  # type: ignore[import]
from pysvg.shape import Path  # type: ignore[import]

import matplotlib.pyplot as plt  # type: ignore[import]
from svg.path import parse_path  # type: ignore[import]
from svg.path.path import Move, Close, Line  # type: ignore[import]

from collections import OrderedDict
import os
import re

import numpy as np

SMALL_SIZE = 8+15
MEDIUM_SIZE = 10+18
BIGGER_SIZE = 14+18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


_types = {"float": np.float32, "double": np.float64,
          "int32": np.int32, "int64": np.int64}


def _get_value_type(key, value):

    if key in ["dimension", "space dimension"]:
        return int(value)
    elif key in ["endian", "encoding"]:
        return value
    elif key in ["type"]:
        return _types[value]
    elif key in ["sizes"]:
        return [int(x) for x in value.split()]
    # else:
    #     raise NotImplementedError


def read_header(file):

    it = iter(file)
    magic_line = next(it)
    if hasattr(magic_line, "decode"):
        need_decode = True
        magic_line = magic_line.decode("ascii", "ignore")  # type: ignore[assignment]
        if not magic_line.startswith("NRRD"):  # type: ignore[arg-type]
            raise NotImplementedError

    header = OrderedDict()

    for line in it:
        if need_decode:
            line = line.decode("ascii", "ignore")  # type: ignore[assignment]

        line = line.rstrip()
        if line.startswith("#"):  # type: ignore[arg-type]
            continue
        elif line == "":
            break

        key, value = re.split(r":=?", line, 1)  # type: ignore[type-var]
        key, value = key.strip(), value.strip()  # type: ignore[attr-defined]

        value = _get_value_type(key, value)

        header[key] = value

    return header


def read_data(file, header):

    size = np.array(header["sizes"])
    dtype = header["type"]

    total_data_points = size.prod(dtype=np.int64)
    dtype_size = np.dtype(dtype).itemsize
    file.seek(-1 * dtype_size * total_data_points, os.SEEK_END)
    data = np.fromfile(file, dtype=dtype, sep="")
    data = data.reshape((size))
    return data


def read_nrrd(file):
    with open(file, "rb") as fh:
        hdr = read_header(fh)
        data = read_data(fh, hdr)
    return data, hdr


eps = 1e-5


def write_segment(x1, x2, y1, y2, i, fh):
    if x1 == x2:
        x1 += eps
    if y1 == y2:
        y1 += eps

    fh.write(f"a = vector({x1:.7f}d0, {y1:.7f}d0, 0.0)\n")
    fh.write(f"b = vector({x2:.7f}d0, {y2:.7f}d0, 0.0)\n")
    fh.write(f"seg({i}) = segment(a, b, mus, mua, hgg, n, layer)\n")
    fh.write(f"ex({i}) = extrude(seg({i}), .5)\n\n")
    i += 1
    return i


fig, axs = plt.subplots(1, 2, figsize=(19.2, 10.8))

with open("svg.f90", "w") as fh:

    maxx = 299.15545999999995
    maxy = 368.92027
    minx = 194.75158
    miny = 197.11304

    x = []
    y = []
    k = 1
    svg = parse('crest-simple.svg')
    paths = svg.getAllElements()[2].getElementsByType(Path)
    for path in paths:
        # dicts = a[i].getAttributes()
        # string = dicts["style"]
        # pos = string.find("stroke-width")
        # string = string[pos:]
        # endpos = string.find(";")
        # width = string[13:endpos]

        segments = parse_path(path.get_d())._segments
        for i in segments:
            if isinstance(i, Line):
                start = i.start
                end = i.end
                x1 = ((start.real - minx) / maxx) - 0.5
                x2 = ((end.real - minx) / maxx) - 0.5
                y1 = ((start.imag - miny) / maxy) - 0.5
                y2 = ((end.imag - miny) / maxy) - 0.5
                x.append(x1)
                x.append(x2)
                y.append(y1)
                y.append(y2)
                # k = write_segment(x1, x2, y1, y2, k, fh)
            elif isinstance(i, Move):
                axs[1].plot(x, y)
                x = []
                y = []
            elif isinstance(i, Close):
                start = i.start
                end = i.end
                x1 = ((start.real - minx) / maxx) - 0.5
                x2 = ((end.real - minx) / maxx) - 0.5
                y1 = ((start.imag - miny) / maxy) - 0.5
                y2 = ((end.imag - miny) / maxy) - 0.5
                x.append(x1)
                x.append(x2)
                y.append(y1)
                y.append(y2)

    # plt.show()


ves, _ = read_nrrd("../data/jmean/logo.nrrd")
img = ves[200, :, :]
print(img.shape)
axs[0].imshow(img[::-1, :], extent=[-1., 1., -1., 1.], aspect="auto")
axs[0].set_xlim([-0.5638903504803154, 0.5638903504803154])
axs[0].set_ylim([-0.5638903504803154, 0.5638903504803154])
axs[0].set_xlabel("Distance/cm")
axs[0].set_ylabel("Distance/cm")
axs[1].set_xlabel("Distance/cm")

plt.savefig("crest-sdf-svg.png", dpi=100, bbox_inches="tight")
# plt.show()