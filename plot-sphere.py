import numpy as np
import matplotlib.pyplot as plt


sdffile = "data/jmean/jmean.dat"
voxfile = "../voxelMCRT/data/jmean/jmean.raw"

fig, axs = plt.subplots(1, 2)
axs = axs.ravel()

data = np.fromfile(sdffile, sep="", dtype=np.float64)
data = data.reshape((200, 200, 200), order="F")
circle1 = plt.Circle((100, 100), 50, facecolor=None, edgecolor="r", fill=False)
circle2 = plt.Circle((100, 100), 50, facecolor=None, edgecolor="r", fill=False)

axs[0].imshow(data[:, 100, ::-1].T, vmin=0, vmax=10)
axs[0].add_patch(circle1)
axs[0].set_title("SDF")
axs[0].set_xlabel("Distance/cm")
axs[0].set_ylabel("Distance/cm")


data = np.fromfile(voxfile, sep="", dtype=np.float64)
data = data.reshape((200, 200, 200), order="F")
axs[1].imshow(data[:, 100, ::-1].T, vmin=0, vmax=10)
axs[1].add_patch(circle2)
axs[1].set_title("Voxel")
axs[1].set_xlabel("Distance/cm")
axs[1].set_ylabel("Distance/cm")

plt.show()
