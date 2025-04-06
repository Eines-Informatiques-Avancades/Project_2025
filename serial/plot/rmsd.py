#!/usr/bin/python3

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt


# Plot appearance settings.
plt.style.use('./mplstyle/science.mplstyle')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "STIXGeneral"
})

# Parse arguments.
def parse_args():
    parser = argparse.ArgumentParser(
        description = "Plot RMSD from a data file."
    )
    parser.add_argument(
        'filename',
        type = str,
        help = "Path to the data file with the RMSD data."
    )
    return parser.parse_args()

# Plot the data from the data file.
def plot_rmsd(filename):
    data = np.loadtxt(filename)
    frames = data[:, 0]
    rmsd_values = data[:, 1]

    plt.plot(
        frames, rmsd_values,
        #  marker = 'o', linestyle = '-',
        color = 'r', label = 'RMSD'
    )

    plt.xlabel(r'$t$')
    plt.ylabel(r'$RMSD(Å)$')
    plt.title(r'Root mean square deviation')

    output_filename = os.path.splitext(filename)[0] + '.pdf'
    plt.savefig(output_filename, format = 'pdf')
    plt.close()

    print(f"RMSD plot saved as {output_filename}")


if __name__ == '__main__':
    args = parse_args()
    plot_rmsd(args.filename)
