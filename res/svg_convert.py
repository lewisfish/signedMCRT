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

eps = 1e-5


def write_segment(x1, x2, y1, y2, i, fh):
    if x1 == x2:
        x1 += eps
    if y1 == y2:
        y1 += eps

    fh.write(f"a = vector({x1:.7f}_wp, {y1:.7f}_wp, 0.0_wp)\n")
    fh.write(f"b = vector({x2:.7f}_wp, {y2:.7f}_wp, 0.0_wp)\n")
    fh.write(f"seg({i}) = segment(a, b, opt(2), layer)\n")
    fh.write(f"array({i}) = extrude(seg({i}), .5_wp)\n\n")
    i += 1
    return i


fig, axs = plt.subplots(1, 1, figsize=(19.2, 10.8))

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
                axs.plot(x, y)
                k = write_segment(x1, x2, y1, y2, k, fh)
            elif isinstance(i, Move):
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
    # axs.invert_xaxis()
    axs.invert_yaxis()
    plt.show()

# plt.savefig("crest-sdf-svg.png", dpi=100, bbox_inches="tight")
# plt.show()