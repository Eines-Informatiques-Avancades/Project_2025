#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

# Plot appearance settings.
plt.style.use('./mplstyle/science.mplstyle')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "STIXGeneral"
})

# Function to plot the energy data
def plot_energy(data_file, label, color):
    data = np.loadtxt(data_file)

    plt.plot(data[:, 0], data[:, 1], label = label, color = color)

    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')
    #  plt.legend(loc = 'lower right', fontsize = 10)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fontsize=10, ncol = 3)


#
# Plot the kinetic, potential and total energies in a single plot.
#

output_file = 'energy_evolution.pdf'

energy_inputfiles = [
    ('../output/lj_potential.dat', 'U', '#00B945'),
    ('../output/kinetic_energy.dat', 'K', '#0C5DA5'),
    ('../output/total_energy.dat', 'E', '#FF2C00')
]

for data_file, label, color in energy_inputfiles:
    plot_energy(data_file, label, color)

plt.savefig(output_file, format = 'pdf')
plt.close()

print(f"Energy plot saved as {output_file}")
