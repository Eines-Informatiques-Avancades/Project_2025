#!/usr/bin/python3

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

    plt.plot(data[:, 0], data[:, 1], label = label, color = color)

    plt.xlabel(r'$t$')
    plt.ylabel(r'$T_{inst}$')


output_file = 'temperature_evolution.pdf'

plot_temperature('../temperature_inst.dat', 'T', '#0C5DA5')

plt.savefig(output_file, format = 'pdf')
plt.close()

print("Temperature plot saved successfully.")
