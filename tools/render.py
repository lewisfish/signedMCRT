from argparse import ArgumentParser
from collections import OrderedDict
import os
import re
import struct

import numpy as np
from skimage import measure  # type: ignore[import]

_types = {"float": np.float32, "double": np.float64,
          "int32": np.int32, "int64": np.int64,
          "uchar": np.ubyte}


def read_nrrd(file):
    with open(file, "rb") as fh:
        hdr = read_header(fh)
        data = read_data(fh, hdr)
    return data, hdr


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


def _get_value_type(key, value):

    if key in ["dimension", "space dimension"]:
        return int(value)
    elif key in ["endian", "encoding"]:
        return value
    elif key in ["type"]:
        return _types[value]
    elif key in ["sizes"]:
        return [int(x) for x in value.split()]
    else:
        pass
        # raise NotImplementedError


def read_data(file, header):

    size = np.array(header["sizes"])
    dtype = header["type"]

    total_data_points = size.prod(dtype=np.int64)
    dtype_size = np.dtype(dtype).itemsize
    file.seek(-1 * dtype_size * total_data_points, os.SEEK_END)
    data = np.fromfile(file, dtype=dtype, sep="")
    data = data.reshape((size))
    return data


def write_binary_stl(path, points):
    n = len(points) // 3
    dtype = np.dtype([
        ('normal', ('<f', 3)),
        ('points', ('<f', 9)),
        ('attr', '<H'),
    ])
    a = np.zeros(n, dtype=dtype)
    a['points'] = np.array(points, dtype='float32').reshape((-1, 9))
    with open(path, 'wb') as fp:
        fp.write(b'\x00' * 80)
        fp.write(struct.pack('<I', n))
        fp.write(a.tobytes())


def marching_cubes(volume, level=0):
    verts, faces, _, _ = measure.marching_cubes(volume, level)
    return verts[faces].reshape((-1, 3))


parser = ArgumentParser()
parser.add_argument("-f", "--file", type=str, default="model.dat",
                    help="Name of file to be plotted.")
parser.add_argument("-n", "--ndim", type=int, default=200,
                    help="Gives the dimensions of the data cube to be examined.")
parser.add_argument("-l", "--levels", type=np.float32, default=0.0,
                    help="Isolevel to render a contour on.")

args = parser.parse_args()

if ".nrrd" in args.file:
    data, hdr = read_nrrd(args.file)
else:
    data = np.fromfile(args.file, sep="", dtype=np.float64)
    data = data.reshape(args.ndim, args.ndim, args.ndim)

points = marching_cubes(data, level=args.levels)
points -= np.amin(points)
points /= np.amax(points)
if ".nrrd" in args.file:
    write_binary_stl(f"{args.file[:-5]}.stl", points)
else:
    write_binary_stl(f"{args.file[:-4]}.stl", points)
