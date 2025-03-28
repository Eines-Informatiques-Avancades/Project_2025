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

# Function to plot the temperature data
def plot_temperature(data_file, label, color):
    data = np.loadtxt(data_file)

    plt.plot(data[:, 0], data[:, 4], label = label, color = color)

    plt.xlabel(r'$t$')
    plt.ylabel(r'$T_{inst}$')


input_file = '../output/thermodynamics.dat'
output_file = os.path.join(
    os.path.dirname(input_file),
    'temperature_evolution.pdf'
)

plot_temperature(input_file, 'T', '#0C5DA5')

plt.savefig(output_file, format = 'pdf')
plt.close()

print(f"Temperature plot saved as {output_file}")
