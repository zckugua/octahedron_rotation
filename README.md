## 项目概述
用于分析 Si–O 八面体在分子动力学/结构数据中的旋转行为，包含：
- 从 `XDATCAR` 计算各原子的平均分数坐标并写出 `POSCAR_AVG`
- 对平均结构进行坐标修正与（可选）晶胞标准化
- 基于 OVITO 邻居搜索计算八面体在 xy 平面上的平均旋转角
- 使用 Kabsch 算法计算理想与形变八面体之间的三轴欧拉角
- 绘制按 z 分层的旋转角热图

## 目录结构
- `get_avg_pos_butch.py`：MPI 并行读取 `POSCAR` 与 `XDATCAR`，计算每个原子的平均位移/平均坐标，输出 `POSCAR_AVG`。
- `fix_poscar.py`：以 `POSCAR` 与 `POSCAR_AVG` 对照，修正超过 0.5 的分数坐标偏移，输出 `POSCAR_fixed`。
- `get_primitive_cell.py`：基于 `POSCAR_AVG` 写出 `POSCAR_primitive` 与 `POSCAR_conventional`（可选）。
- `cal_rotation_angle.py`：使用 OVITO 搜索 Si 的 6 个 O 邻居，取 z 方向最近的 4 个 O，计算 xy 平面的有符号平均旋转角，输出 `rotation_angles_xy_with_coords.dat`。
- `cal_rotation_angle_kabsch.py`：以理想结构（`POSCAR`）与形变结构（`POSCAR_fixed`）对应八面体点集，Kabsch 求旋转矩阵并提取 `xyz` 欧拉角，输出 `rotation_angles_kabsch.dat`。
- `plot_rotation_heatmap.py`：读取 `rotation_angles_xy_with_coords.dat`，按 z 近似分层聚类并绘制热图到 `rotation_heatmaps/`。

## 依赖
- 运行环境：Python 3.8+
- 核心库：
  - 数值/绘图：`numpy`, `pandas`, `matplotlib`, `tqdm`
  - I/O/并行：`psutil`, `mpi4py`
  - 结构处理：`pymatgen`
  - 分析/邻居搜索：`ovito`
  - 旋转与欧拉角：`scipy`

安装示例：
```bash
pip install numpy pandas matplotlib tqdm psutil mpi4py pymatgen ovito scipy
```

## 基本流程
1) 计算平均结构（需要 `POSCAR`、`XDATCAR` 在当前目录）
```bash
# 建议并行运行（按需调整核数）
mpirun -n 8 python get_avg_pos_butch.py
# 产物：POSCAR_AVG
```

2) 修正分数坐标的周期跳变
```bash
python fix_poscar.py
# 产物：POSCAR_fixed
```

3)（可选）导出原胞/常规胞
```bash
python get_primitive_cell.py
# 产物：POSCAR_primitive, POSCAR_conventional
```

4) 计算 xy 平面有符号平均旋转角（基于 4 个最邻近 z 的 O）
```bash
python cal_rotation_angle.py
# 默认读取：POSCAR_fixed
# 产物：rotation_angles_xy_with_coords.dat
```

5) 绘制分层热图
```bash
python plot_rotation_heatmap.py
# 读入：rotation_angles_xy_with_coords.dat
# 产物目录：rotation_heatmaps/
```

6)（可选）Kabsch 三轴欧拉角（理想 vs 形变）
```bash
python cal_rotation_angle_kabsch.py
# 读取：POSCAR（理想），POSCAR_fixed（形变）
# 产物：rotation_angles_kabsch.dat
```

## 输入/输出格式
- `POSCAR`/`XDATCAR`：VASP 标准文件（`get_avg_pos_butch.py` 假设分数坐标 Direct 格式）
- `POSCAR_AVG`/`POSCAR_fixed`：VASP 结构文件
- `rotation_angles_xy_with_coords.dat`：
  - 列：`Si_index  x  y  z  z_rounded  rotation_angle_deg`
- `rotation_angles_kabsch.dat`：
  - 列：`Si_index  x  y  z  rot_x_deg  rot_y_deg  rot_z_deg`
- 热图输出：`rotation_heatmaps/rotation_heatmap_z{z:.2f}.png`

## 可调参数与注意事项
- `get_avg_pos_butch.py`
  - `chunk_size`：流式读取 XDATCAR 的块大小，默认 1000；支持 MPI 并行，打印内存消耗。
- `cal_rotation_angle.py`
  - `cutoff`（默认 2.5 Å）：Si–O 邻居判定半径；需保证能找到 6 个 O。
  - 仅使用与 Si 在 z 方向最近的 4 个 O 计算 xy 面旋转。
- `cal_rotation_angle_kabsch.py`
  - `cutoff`（默认 2.2 Å）构建 Si–O 映射；
  - 使用相对坐标并做 PBC 展开；提取欧拉角顺序为 `xyz`（度数）。
  - 代码中对晶格长度的获取假设正交晶格：`cell_lengths = diag(cell[:3,:3])`，非正交请自行修改。
- `plot_rotation_heatmap.py`
  - `bin_width`（默认 1.0）：xy 方向网格化宽度；
  - `z_tolerance`（默认 0.5）：z 分层聚类容差。
