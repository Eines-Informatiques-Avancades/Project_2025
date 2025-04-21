import numpy as np
import matplotlib.pyplot as plt

# Kr LJ potential real unit parameters
sigma = 3.6274  # Å
epsilon_by_kB = 162.58  # K

# T* reduced temperature
T_star = 1.0
T_real = T_star * epsilon_by_kB  # T = T* eps / Kb


m = ( 83.798 / (6.022e23) )/ 1000  # 1 Kr atom mass in Kg
kB = 1.380649e-23  # Boltzmann constant in J/K
epsilon = epsilon_by_kB * kB  # eps in J
sigma_m = sigma * 1e-10  # sigma in m

tau = sigma_m * np.sqrt(m / epsilon)  # tau (in s) = sigma * sqrt(m / epsilon)
tau_ps = tau * 1e12  # convert it to ps

# load rmsd datas in reduced units
rmsd_data = np.loadtxt('rmsd.dat')
time_reduced = rmsd_data[:, 0]
rmsd_reduced = rmsd_data[:, 1]

# convert it to real units
rmsd_real = rmsd_reduced * sigma  # RMSD in Å
time_real = time_reduced * tau_ps  # ps


# Plot RMSD in real units
plt.figure()
plt.plot(time_real, rmsd_real)
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (Å)")
plt.title(f"RMSD of Kr at {T_real:.2f} K")
plt.grid(True)
plt.tight_layout()
plt.savefig("rmsd_real_units.png")

# load rdf.dat
rdf_data = np.loadtxt('rdf.dat')
r_reduced = rdf_data[:, 0]
g_r = rdf_data[:, 1]  # rdf is a relative quantity, reduced units does not affect it.
r_real = r_reduced * sigma  # interatomic distance in Å


# Plot RDF in real units
plt.figure()
plt.plot(r_real, g_r)
plt.xlabel("r (Å)")
plt.ylabel("g(r)")
plt.title(f"Radial Distribution Function of Kr at {T_real:.2f} K")
plt.grid(True)
plt.tight_layout()
plt.savefig("rdf_real_units.png")

plt.show()

