import numpy as np
import matplotlib.pyplot as plt

# Kr LJ potential real unit parameters
sigma = 3.6274  # Å
epsilon_by_kB = 162.58  # K

kB = 1.380649e-23  # J/K
m = (83.798 / 6.022e23) / 1000  # 1 Kr atom mass in Kg

epsilon = epsilon_by_kB * kB  # eps in J
tau_ps = sigma * 1e-10 * np.sqrt(m / epsilon) * 1e12  # tau (in ps) = sigma * sqrt(m / epsilon)
energy_conv = epsilon * (6.022e23) / 1000  # epsilon in KJ/mol


# load thermodynamics.dat
data = np.loadtxt("thermodynamics.dat", comments="#")
time = data[:, 0]
lj_pot = data[:, 1]
ke = data[:, 2]
etot = data[:, 3]
temp = data[:, 4]

# convert it to real units
time_ps = time * tau_ps
pe_kJmol = lj_pot * energy_conv
ke_kJmol = ke * energy_conv
etot_kJmol = etot * energy_conv
temp_K = temp * epsilon_by_kB

# Instantaneous temperature variation plot
plt.figure()
plt.plot(time_ps, temp_K, lw=1)
plt.xlabel("Time (ps)")
plt.ylabel("Temperature (K)")
plt.title("Instantaneous temperature variation")
plt.grid(True)
plt.tight_layout()
plt.savefig("temperature_vs_time.png")


# Kinetic, potential, total energies plot
plt.figure()
plt.plot(time_ps, pe_kJmol, label="U", lw=1)
plt.plot(time_ps, ke_kJmol, label="KE", lw=1)
plt.plot(time_ps, etot_kJmol, label="E = U + KE", lw=1)
plt.xlabel("Time (ps)")
plt.ylabel("Energy (kJ/mol)")
plt.title("Energies evolution")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("energy_components_vs_time.png")

plt.show()
