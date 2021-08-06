from argparse import ArgumentParser
import struct

from skimage import measure  # type: ignore[import]
import numpy as np  # type: ignore[import]


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
parser.add_argument("-f", "--file", type=str,
                    help="Name of file to be plotted.")
parser.add_argument("-n", "--ndim", type=int, default=200,
                    help="Gives the dimensions of the data cube to be examined.")
parser.add_argument("-l", "--levels", type=np.float32, default=0.0,
                    help="Isolevel to render a contour on.")

args = parser.parse_args()


data = np.fromfile(args.file, sep="", dtype=np.float64)
data = data.reshape(args.ndim, args.ndim, args.ndim)

points = marching_cubes(data, level=args.levels)
points -= np.amin(points)
points /= np.amax(points)
write_binary_stl("out.stl", points)
