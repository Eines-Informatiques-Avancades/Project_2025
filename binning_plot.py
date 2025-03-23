import pandas as pd
import matplotlib.pyplot as plt

# Leer el archivo binning_columns.out, ignorando las líneas que empiezan con #
df = pd.read_csv('binning_columns.out', comment='#', 
                 delim_whitespace=True, 
                 header=None, 
                 names=['ColName', 'BinSize', 'Mean', 'Variance', 'StdDev'])

# Agrupar los datos por el nombre de la columna (primera columna)
groups = df.groupby('ColName')
n_groups = len(groups)

# Crear una figura con un subplot por cada grupo
fig, axs = plt.subplots(n_groups, 1, figsize=(8, 4 * n_groups), sharex=True)

# Si solo hay un grupo, convertir axs a una lista para poder iterar
if n_groups == 1:
    axs = [axs]

# Para cada grupo, generar el plot correspondiente
for ax, (col_name, group) in zip(axs, groups):
    ax.plot(group['BinSize'], group['StdDev'], marker='o', linestyle='-')
    ax.set_title(f'Binning para {col_name}')
    ax.set_xlabel('  ')
    ax.set_ylabel('Desviación Estándar (StdDev)')
    ax.grid(True)

plt.tight_layout()
plt.show()
