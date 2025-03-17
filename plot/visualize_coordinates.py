#!/usr/bin/python

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D


# Plot appearance settings.
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "STIXGeneral"
})

# Parse arguments.
def parse_args():
    parser = argparse.ArgumentParser(
        description = "Visualize particles in 3D space from a data file"
    )
    parser.add_argument(
        'filename',
        type = str,
        help = "Path to the data file with x, y, z coordinates"
    )
    return parser.parse_args()

# Plot the data from the data file.
def plot_particles(filename):
    data = np.loadtxt(filename)

    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    fig = plt.figure()
    ax = fig.add_subplot(
        111,
        projection = '3d'
    )

    ax.scatter(
        x, y, z,
        c = 'r',
        marker = 'o'
    )

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Particle Visualization')

    output_filename = os.path.splitext(filename)[0] + '.pdf'
    plt.savefig(output_filename, format = 'pdf')
    plt.close()

    print(f"Particle coordinates plot saved as {output_filename}")

if __name__ == '__main__':
    args = parse_args()
    plot_particles(args.filename)
