import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def explode(data):
    size = np.array(data.shape)*2
    data_e = np.zeros(size - 1, dtype=data.dtype)
    data_e[::2, ::2, ::2] = data
    return data_e

if __name__ == "__main__":
    # load voxel file
    subject = 'Julie'
    date = '20231012'
    vol = np.load(f"{subject}_{date}_merge_vol_20mm.npy")
    vol = vol.T

    # shringe voxel size --> save computation time
    vol = vol[::4, ::4, ::12]
    cmap = ['red', 'salmon', 'sienna', 'silver',
            'tan', 'white', 'violet', 'wheat', 'yellow']
    vol = explode(vol)
    
    # plot 3D voxel
    colors = np.empty(list(vol.shape) + [4], dtype=np.float32)
    alpha = 0.5
    colors[vol == 1] = [1, 0, 0, alpha]
    colors[vol == 2] = [0, 1, 0, alpha]
    colors[vol == 3] = [0, 0, 1, alpha]
    colors[vol == 4] = [1, 1, 0, alpha]
    colors[vol == 5] = [1, 0, 1, alpha]
    colors[vol == 6] = [0, 1, 1, 0.1]
    colors[vol == 7] = [1, 1, 1, 1]
    colors[vol == 8] = [0, 0, 0, 1]
    colors[vol == 9] = [0.5, 0.5, 0.5, 1]
    edgecolor = [1, 1, 1, 0]

    x, y, z = np.indices(np.array(vol.shape) + 1).astype(float) // 2
    x[0::2, :, :] += 0.05
    y[:, 0::2, :] += 0.05
    z[:, :, 0::2] += 0.05
    x[1::2, :, :] += 0.95
    y[:, 1::2, :] += 0.95
    z[:, :, 1::2] += 0.95
    
    ax = plt.figure().add_subplot(projection='3d')
    ax.voxels(x, y, z, vol, facecolors=colors, edgecolor=edgecolor)
    ax.view_init(25, 128)
    plt.savefig(os.path.join("pic", subject, "20mmvol3d.png"), dpi=300, format='png', bbox_inches='tight')
    plt.show()