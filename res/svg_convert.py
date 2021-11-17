from pysvg.parser import parse
from pysvg.shape import Line, Path

import matplotlib.pyplot as plt
from svg.path import parse_path
from svg.path.path import Move, Close, Line as Liner


def write_segment(x1, x2, y1, y2, i, fh):
    fh.write(f"a = vector({x1:.4f}d0, {y1:.4f}d0, 0.0)\n")
    fh.write(f"b = vector({x2:.4f}d0, {y2:.4f}d0, 0.0)\n")
    fh.write(f"seg({i}) = segment(a, b, mus, mua, hgg, n, layer)\n\n")
    i += 1
    return i


with open("test.txt", "w") as fh:

    maxx = 299.15545999999995
    maxy = 368.92027
    minx = 194.75158
    miny = 197.11304

    x = []
    y = []
    k = 1
    s = parse('crest-simple.svg')
    a = s.getAllElements()[2].getElementsByType(Path)
    for i in range(len(a)):
        b = parse_path(a[i].get_d())
        c = b._segments
        for i in c:
            if isinstance(i, Liner):
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
                k = write_segment(x1, x2, y1, y2, k, fh)

            elif isinstance(i, Move):
                plt.plot(x, y)
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

    plt.show()
print(max(x), max(y))