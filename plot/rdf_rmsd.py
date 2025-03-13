#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D


def plot_rdf(filename="rdf_data.txt"):
    data = np.loadtxt(filename)
    r = data[:, 0]
    g_r = data[:, 1]

    plt.figure(figsize=(8,6))
    plt.plot(r, g_r, label="RDF", color='b')
    plt.xlabel("Distance r")
    plt.ylabel("g(r)")
    plt.title("Radial Distribution Function (RDF)")
    plt.legend()
    plt.grid()

    rdf_output = "rdf_plot.png"
    plt.savefig(rdf_output, dpi=300)
    plt.show()
    print(f"RDF plot saved as {rdf_output}")


def plot_rmsd(filename="rmsd_data.txt"):
    data = np.loadtxt(filename)
    frames = data[:, 0]
    rmsd_values = data[:, 1]

    plt.figure(figsize=(8,6))
    plt.plot(frames, rmsd_values, marker='o', linestyle='-', color='r', label="RMSD")
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD Over Time")
    plt.legend()
    plt.grid()

    rmsd_output = "rmsd_plot.png"
    plt.savefig(rmsd_output, dpi=300)
    plt.show()
    print(f"RMSD plot saved as {rmsd_output}")


plot_rdf()
plot_rmsd()
