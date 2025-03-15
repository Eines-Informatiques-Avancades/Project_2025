import numpy as np
import matplotlib.pyplot as plt

def plot_rmsd(filename="rmsd_data.txt"):
    data = np.loadtxt(filename)
    frames = data[:, 0]
    rmsd_values = data[:, 1]

    plt.figure(figsize=(8,6))
    plt.plot(frames, rmsd_values, marker='o', linestyle='-', color='r', label="RMSD")
    plt.xlabel(r'$Time$')
    plt.ylabel(r'$RMSD(Å)$')
    plt.title(r'Root mean square deviation')
    plt.legend()
    plt.grid()

    rmsd_output = "rmsd_plot.pdf"
    plt.savefig(rmsd_output, format = 'pdf')
    plt.show()
    print(f"RMSD plot saved as {rmsd_output}")

plot_rmsd()