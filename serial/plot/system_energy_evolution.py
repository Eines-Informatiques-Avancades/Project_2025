#!/usr/bin/python3

import os
import numpy as np
import matplotlib.pyplot as plt

# Plot appearance settings.
plt.style.use('./mplstyle/science.mplstyle')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "STIXGeneral"
})

# Function to plot the energy data
def plot_energy(data_file, column, label, color):
    data = np.loadtxt(data_file)

    plt.plot(data[:, 0], data[:, column], label = label, color = color)

    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')

    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15),
        fontsize=10,
        ncol = 3
    )


#
# Plot the kinetic, potential and total energies in a single plot.
#

energy_inputfiles = [
    ('../output/thermodynamics.dat', 1, 'U', '#00B945'),
    ('../output/thermodynamics.dat', 2, 'K', '#0C5DA5'),
    ('../output/thermodynamics.dat', 3, 'E', '#FF2C00')
]

for data_file, column, label, color in energy_inputfiles:
    plot_energy(data_file, column, label, color)

output_file = os.path.join(os.path.dirname(data_file), 'energy_evolution.pdf')

plt.savefig(output_file, format = 'pdf')
plt.close()

print(f"Energy plot saved as {output_file}")
