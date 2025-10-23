import numpy as np
from ovito.io import import_file
from ovito.data import CutoffNeighborFinder
from scipy.spatial.transform import Rotation as R

def compute_rotation_matrix_kabsch(P, Q):
    P_cent = P - P.mean(axis=0)
    Q_cent = Q - Q.mean(axis=0)
    H = P_cent.T @ Q_cent
    U, S, Vt = np.linalg.svd(H)
    R_kabsch = Vt.T @ U.T
    if np.linalg.det(R_kabsch) < 0:
        Vt[-1, :] *= -1
        R_kabsch = Vt.T @ U.T
    return R_kabsch

def extract_euler_angles_from_rotation_matrix(R_mat):
    r = R.from_matrix(R_mat)
    return r.as_euler('xyz', degrees=True)

def unwrap_relative_positions(center, positions, cell_lengths):
    delta = positions - center
    delta -= np.round(delta / cell_lengths) * cell_lengths
    return delta  # 返回相对坐标（而不是绝对坐标）

def build_si_o_neighbor_map(pipeline, cutoff=2.2):
    data = pipeline.compute()
    positions = data.particles.positions
    types = data.particles['Particle Type']
    pt = data.particles.particle_types
    si_id = pt.type_by_name('Si').id
    o_id = pt.type_by_name('O').id
    finder = CutoffNeighborFinder(cutoff, data)

    si_o_map = {}
    for i in range(len(positions)):
        if types[i] != si_id:
            continue
        neighbors = [n.index for n in finder.find(i) if types[n.index] == o_id]
        if len(neighbors) == 6:
            si_o_map[i] = sorted(neighbors)
    return si_o_map

def calculate_rotation_angles_kabsch_relative(ideal_file, deformed_file, output_file="rotation_angles_kabsch.dat", cutoff=2.2):
    ideal_pipeline = import_file(ideal_file, input_format="vasp")
    deformed_pipeline = import_file(deformed_file, input_format="vasp")
    ideal_data = ideal_pipeline.compute()
    deformed_data = deformed_pipeline.compute()

    ideal_positions = ideal_data.particles.positions
    deformed_positions = deformed_data.particles.positions
    cell_lengths = np.diag(deformed_data.cell[:3, :3])  # 假设正交晶格

    si_o_map = build_si_o_neighbor_map(ideal_pipeline, cutoff)

    results = []

    for si_index, o_indices in si_o_map.items():
        si_pos_ideal = ideal_positions[si_index]
        si_pos_deformed = deformed_positions[si_index]

        o_pos_ideal = np.array([ideal_positions[j] for j in o_indices])
        o_pos_deformed = np.array([deformed_positions[j] for j in o_indices])

        # ✅ 使用相对坐标 + PBC 展开处理
        relative_o_ideal = unwrap_relative_positions(si_pos_ideal, o_pos_ideal, cell_lengths)
        relative_o_deformed = unwrap_relative_positions(si_pos_deformed, o_pos_deformed, cell_lengths)

        # ✅ 用相对坐标计算旋转
        R_mat = compute_rotation_matrix_kabsch(relative_o_ideal, relative_o_deformed)
        angles = extract_euler_angles_from_rotation_matrix(R_mat)

        results.append((si_index, si_pos_deformed[0], si_pos_deformed[1], si_pos_deformed[2], *angles))

    # 保存结果
    with open(output_file, 'w') as f:
        f.write("Si_index\tx\ty\tz\trot_x_deg\trot_y_deg\trot_z_deg\n")
        for row in results:
            f.write("\t".join(f"{val:.6f}" if isinstance(val, float) else str(val) for val in row) + "\n")

    print(f"[✔] 已保存旋转角结果：{output_file}")

if __name__ == "__main__":
    calculate_rotation_angles_kabsch_relative("POSCAR", "POSCAR_fixed")

