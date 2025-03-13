#!/usr/bin/python

import os
import argparse
import numpy as np
import imageio.v2 as imageio
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Plot appearance settings.
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "STIXGeneral"
})

# Parse arguments.
def parse_args():
    parser = argparse.ArgumentParser(
        description = "Visualize particle trajectories over time"
    )
    parser.add_argument(
        'filename',
        type = str,
        help = "Path to the data file with time, x, y, z coordinates"
    )
    parser.add_argument(
        '-o', '--output',
        type = str,
        default = 'output.gif',
        help = "Output file name (default: output.gif)"
    )
    return parser.parse_args()

# Create a GIF of the particle trajectories over time.
def create_gif(filename, output_filename):
    data = np.loadtxt(filename)
    times = np.unique(data[:, 0])

    images = []

    for time in times:
        # Extract the particles' positions at the current time.
        time_data = data[data[:, 0] == time]
        x = time_data[:, 1]
        y = time_data[:, 2]
        z = time_data[:, 3]

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
        ax.set_title(f'Time: {time:.4f}')

        # Save the current plot to a temporary file
        temp_file = 'temp_plot.png'
        plt.savefig(temp_file)
        plt.close()

        print(f'Generated plot for time {time:.4f}')

        # Read the saved plot as an image and append it to the images list.
        image = imageio.imread(temp_file)
        images.append(image)

        os.remove(temp_file)

    # Create the GIF from the list of images.
    # Adjust duration for timing between frames.
    imageio.mimsave(output_filename, images, duration=0.5)

    print(f"GIF saved as {output_filename}")

if __name__ == '__main__':
    args = parse_args()
    create_gif(args.filename, args.output)
