import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

# 读取 RDF 数据并保存图像
def plot_rdf(filename="rdf_data.txt"):
    data = np.loadtxt(filename)
    r = data[:, 0]
    g_r = data[:, 1]

    plt.figure(figsize=(8,6))
    plt.plot(r, g_r, label="RDF", color='b')
    plt.xlabel("Distance r")
    plt.ylabel("g(r)")
    plt.title("Radial Distribution Function (RDF)")
    plt.legend()
    plt.grid()

    rdf_output = "rdf_plot.png"
    plt.savefig(rdf_output, dpi=300)
    plt.show()
    print(f"RDF plot saved as {rdf_output}")

# 读取 RMSD 数据并保存图像
def plot_rmsd(filename="rmsd_data.txt"):
    data = np.loadtxt(filename)
    frames = data[:, 0]
    rmsd_values = data[:, 1]

    plt.figure(figsize=(8,6))
    plt.plot(frames, rmsd_values, marker='o', linestyle='-', color='r', label="RMSD")
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD Over Time")
    plt.legend()
    plt.grid()

    rmsd_output = "rmsd_plot.png"
    plt.savefig(rmsd_output, dpi=300)
    plt.show()
    print(f"RMSD plot saved as {rmsd_output}")

# 读取 traj.xyz 并保存粒子轨迹图像
def plot_trajectories(filename="traj.xyz"):
    data = np.loadtxt(filename)

    frame_ids = data[:, 0].astype(int)  # 提取帧编号
    x = data[:, 1]
    y = data[:, 2]
    z = data[:, 3]

    unique_frames = np.unique(frame_ids)  # 获取所有帧编号
    num_atoms = len(frame_ids) // len(unique_frames)  # 计算原子数

    # 重新整理数据：每个原子单独存储
    trajectories = {i: {'x': [], 'y': [], 'z': []} for i in range(num_atoms)}

    for i in range(len(frame_ids)):
        atom_id = i % num_atoms  # 根据原子编号分组
        trajectories[atom_id]['x'].append(x[i])
        trajectories[atom_id]['y'].append(y[i])
        trajectories[atom_id]['z'].append(z[i])

    # 绘制 3D 轨迹
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    for atom_id, traj in trajectories.items():
        ax.plot(traj['x'], traj['y'], traj['z'], marker='o', label=f'Atom {atom_id+1}')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Particle Trajectories')
    ax.legend()

    traj_output = os.path.splitext(filename)[0] + '_trajectory.png'
    plt.savefig(traj_output, dpi=300)
    plt.show()

    print(f"Trajectory plot saved as {traj_output}")

# 调用所有函数
plot_rdf()
plot_rmsd()
plot_trajectories("traj.xyz")
