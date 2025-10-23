from pymatgen.core import Structure
import numpy as np

# 读取初始POSCAR和POSCAR-10000
structure_ref = Structure.from_file("POSCAR")            # 参考结构
structure_md = Structure.from_file("POSCAR_AVG")       # MD中提取的结构

# 检查是否原子数量一致
if len(structure_ref) != len(structure_md):
    raise ValueError("两个结构中的原子数不一致！")

# 创建新的坐标列表
new_coords = []

for site_ref, site_md in zip(structure_ref, structure_md):
    diff = site_md.frac_coords - site_ref.frac_coords
    # 对于大于0.5的坐标偏移，减去1
    corrected_coords = site_md.frac_coords - np.where(diff > 0.5, 1.0, 0.0)
    new_coords.append(corrected_coords)

# 创建新的结构对象
structure_new = structure_md.copy()
structure_new.remove_sites(range(len(structure_new)))  # 清除旧结构
for i, site in enumerate(structure_md):
    structure_new.append(site.specie, new_coords[i], coords_are_cartesian=False)

# 写出新的POSCAR文件
structure_new.to(fmt="poscar", filename="POSCAR_fixed")
print("已生成修正后的结构文件：POSCAR_fixed")

