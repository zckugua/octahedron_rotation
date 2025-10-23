import matplotlib
matplotlib.use("Agg")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def cluster_z_layers(z_values, tolerance=0.5):
    z_values_sorted = sorted(z_values)
    clustered = []
    current_group = []

    for z in z_values_sorted:
        if not current_group:
            current_group.append(z)
        elif abs(z - current_group[-1]) <= tolerance:
            current_group.append(z)
        else:
            clustered.append(np.mean(current_group))
            current_group = [z]

    if current_group:
        clustered.append(np.mean(current_group))

    return clustered

def plot_rotation_heatmaps_binned(data_file="rotation_angles_xy_with_coords.dat",
                                   output_dir="rotation_heatmaps",
                                   bin_width=1.0,
                                   z_tolerance=0.5):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(data_file, delim_whitespace=True)

    # ✅ 进行 z 层聚类
    z_values = df["z"].unique()
    z_layers = cluster_z_layers(z_values, tolerance=z_tolerance)

    for z in z_layers:
        layer = df[np.abs(df["z"] - z) <= z_tolerance].copy()

        layer.loc[:, "x_bin"] = (layer["x"] / bin_width).round().astype(int)
        layer.loc[:, "y_bin"] = (layer["y"] / bin_width).round().astype(int)

        x_sorted = sorted(layer["x_bin"].unique())
        y_sorted = sorted(layer["y_bin"].unique())
        x_map = {val: i for i, val in enumerate(x_sorted)}
        y_map = {val: i for i, val in enumerate(y_sorted)}

        grid = np.full((len(y_sorted), len(x_sorted)), np.nan)

        for _, row in layer.iterrows():
            xi = x_map[row["x_bin"]]
            yi = y_map[row["y_bin"]]
            grid[yi, xi] = row["rotation_angle_deg"]

        plt.figure(figsize=(6, 6))
        vmax = np.nanmax(np.abs(grid))
        plt.imshow(grid, cmap="seismic", origin="lower", vmin=-vmax, vmax=vmax)
        plt.colorbar(label="Rotation Angle (deg)")
        plt.title(f"Si-O Octahedral Rotation (z ≈ {z:.2f})")
        plt.xticks(range(len(x_sorted)), x_sorted)
        plt.yticks(range(len(y_sorted)), y_sorted)
        plt.grid(False)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/rotation_heatmap_z{z:.2f}.png")
        plt.close()

    print(f"[✔] 热图已保存至目录：{output_dir}/")

if __name__ == "__main__":
    plot_rotation_heatmaps_binned()

