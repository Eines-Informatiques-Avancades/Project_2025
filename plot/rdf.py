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
        description = "Plot RDF from a data file."
    )
    parser.add_argument(
        'filename',
        type = str,
        help = "Path to the data file with the RDF data."
    )
    return parser.parse_args()

# Plot the data from the data file.
def plot_rdf(filename):
    data = np.loadtxt(filename)
    r = data[:, 0]
    g_r = data[:, 1]

    plt.plot(r, g_r, label = 'RDF', color = 'b')

    plt.xlabel(r'$r(Å)$')
    plt.ylabel(r'$g(r)$')
    plt.title(r'Radial distribution function')

    output_filename = os.path.splitext(filename)[0] + '.pdf'
    plt.savefig(output_filename, format = 'pdf')
    plt.close()

    print(f"RDF plot saved as {output_filename}.")


if __name__ == '__main__':
    args = parse_args()
    plot_rdf(args.filename)
