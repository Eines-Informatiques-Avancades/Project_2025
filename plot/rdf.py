#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

def plot_rdf(filename="rdf_data.txt"):
    data = np.loadtxt(filename)
    r = data[:, 0]
    g_r = data[:, 1]

    plt.figure(figsize=(8,6))
    plt.plot(r, g_r, label="RDF", color='b')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$g(r)$')
    plt.title(r'Radial distribution function')
    plt.legend()
    plt.grid()

    rdf_output = "rdf_plot.pdf"
    plt.savefig(rdf_output, format = 'pdf')
    plt.show()
    print(f"RDF plots saved successfully.")


plot_rdf()

