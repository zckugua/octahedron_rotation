import numpy as np
import pandas as pd
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder

def calculate_si_octahedral_rotation(input_file, output_file="rotation_angles_xy_with_coords.dat", cutoff=2.5):
    # Step 1: Load structure
    pipeline = import_file(input_file, input_format="vasp")
    data = pipeline.compute()

    positions = data.particles.positions
    type_property = data.particles['Particle Type']
    particle_types = data.particles.particle_types

    Si_id = particle_types.type_by_name('Si').id
    O_id = particle_types.type_by_name('O').id

    finder = CutoffNeighborFinder(cutoff, data)
    results = []

    # Step 2: Loop over all Si atoms
    for i in range(len(positions)):
        if type_property[i] != Si_id:
            continue

        si_pos = positions[i]
        o_neighbors = []

        for neigh in finder.find(i):
            j = neigh.index
            if type_property[j] == O_id:
                o_neighbors.append(j)

        if len(o_neighbors) != 6:
            continue

        # Step 3: Select 4 O atoms closest in z
        o_pos = np.array([positions[j] for j in o_neighbors])
        dz = np.abs(o_pos[:, 2] - si_pos[2])
        top4_indices = np.argsort(dz)[:4]
        o_pos_xy = o_pos[top4_indices]

        # Step 4: Compute signed angles
        vecs = o_pos_xy - si_pos
        vecs[:, 2] = 0
        vecs_norm = vecs[:, :2] / np.linalg.norm(vecs[:, :2], axis=1)[:, None]

        ideal_dirs = np.array([
            [1, 0],
            [-1, 0],
            [0, 1],
            [0, -1]
        ])

        individual_angles = []

        for v in vecs_norm:
            best_angle = None
            min_diff = float('inf')
            for ref in ideal_dirs:
                dot = np.dot(v, ref)
                angle = np.arccos(np.clip(dot, -1.0, 1.0))
                cross = np.cross(ref, v)
                signed_angle = np.rad2deg(np.sign(cross) * angle)
                if angle < min_diff:
                    min_diff = angle
                    best_angle = signed_angle
            individual_angles.append(best_angle)

        avg_angle = np.mean(individual_angles)

        results.append((i, si_pos[0], si_pos[1], si_pos[2], round(si_pos[2], 1), avg_angle))

    # Step 5: Output to .dat file
    with open(output_file, 'w') as f:
        f.write("Si_index\tx\ty\tz\tz_rounded\trotation_angle_deg\n")
        for row in results:
            f.write("\t".join(f"{val:.6f}" if isinstance(val, float) else str(val) for val in row) + "\n")

    print(f"[✔] 保存成功：{output_file}")

if __name__ == "__main__":
    calculate_si_octahedral_rotation("POSCAR_fixed")

