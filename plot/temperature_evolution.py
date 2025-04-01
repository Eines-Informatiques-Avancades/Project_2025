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

# Parameters used to simulate different noble gases
atom_dictionary={'atom': ['Ar', 'Kr', 'Xe'],
                 'epsilon': [143.224, 201.2, 282.29],
                 'sigma': [3.3527, 3.5709, 3.8924]}

# Function to read input_parameters.in in order to determine atom type to plot
def read_line_from_file(file_name, line_number):
    with open(file_name, 'r') as file:
        for line_no, line in enumerate(file):
            if line_no == line_number:
                file.close()
                return line
        else:
            file.close()
            raise ValueError('line %s does not exist in file %s' % (line_number, file_name))

line = read_line_from_file('../input_parameters.in', 1)

# Function that sets parameters according to atom_type chosen per last simmulation
def set_atom_parameters(atom_dictionary, line):
    for i in range(len(atom_dictionary['atom'])):
        if atom_dictionary['atom'][i] == line[:2]:
            atom = atom_dictionary['atom'][i]
            epsilon = atom_dictionary['epsilon'][i]
            sigma = atom_dictionary['sigma'][i]
            return atom, epsilon, sigma
atom = set_atom_parameters(atom_dictionary,line)[0]
epsilon = set_atom_parameters(atom_dictionary,line)[1]
sigma = set_atom_parameters(atom_dictionary,line)[2]

# Function to plot the temperature data
def plot_temperature(data_file, atom, epsilon, label, color):
    data = np.loadtxt(data_file)

    # Initialise the subplot function using number of rows and columns
    figure, axis = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, gridspec_kw=None)
    figure.subplots_adjust(hspace=0.5)
    
    # Plot for adimensional data
    axis[0].plot(data[:, 0], data[:, 4], label = label, color = color)
    axis[0].set_title("Adimensional simmulation")
    
    # Plot using chosen atom parameters
    axis[1].plot(data[:, 0], data[:, 4]*epsilon, label = label, color = color)
    axis[1].set_title(atom + " simulation")
    plt.xlabel(r'$t$')
    plt.ylabel(r'$T_{inst}$(K)')


input_file = '../output/thermodynamics.dat'
output_file = os.path.join(
    os.path.dirname(input_file),
    'temperature_evolution.pdf'
)

plot_temperature(input_file, atom, epsilon, 'T', '#0C5DA5')

plt.savefig(output_file, format = 'pdf')
plt.close()

print(f"Temperature plot saved as {output_file}")

