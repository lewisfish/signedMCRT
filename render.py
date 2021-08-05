from skimage import measure
import numpy as np
import struct


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


file = "model.dat"
data = np.fromfile(file, sep="", dtype=np.float64)
data = data.reshape(401, 401, 401)

points = marching_cubes(data)
points -= np.amin(points)
points /= np.amax(points)
write_binary_stl("out.stl", points)
