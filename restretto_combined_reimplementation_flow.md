# REstretto Python 重构完整实现流程

本文合并计算化学和代码工程两个视角，给出一份可执行的 Python 重构流程。目标是复现 `references/restretto` 的主要功能：`atomgrid-gen`、`conformer-docking`、`decompose`，并兼容 `score-only`、`intraenergy-only`、`local-only/no-local-opt` 等辅助模式。

## 1. 重构目标和验收范围

必须实现：

1. 配置文件解析和 CLI 覆盖。
2. OpenBabel 分子 IO、键阶修正、极性氢处理、canonical SMILES、canonical renumbering、坐标回写。
3. X-Score 原子类型和 Vina 风格五项打分。
4. `.grid` 二进制读写。
5. receptor atom grids 生成。
6. ligand 内能计算。
7. ligand 片段分解、dummy 添加、唯一 fragment library。
8. fragment energy grid 生成和 OFFLINE/ONLINE/NONE 缓存复用。
9. 60 旋转采样或 `ROTANGS` 自定义旋转。
10. conformer docking 搜索、候选筛选、局部优化、RMSD 去重和输出。
11. `decompose` 独立功能。
12. `score-only` 和 `intraenergy-only`。

默认可暂缓：

1. `atom-docking` 作为主流程对照功能，因原 Makefile 默认不构建且源码存在旧接口痕迹。
2. 与 C++ 每个随机优化步完全数值一致。若需要严格一致，应单独实现 C `rand()` 兼容层。

## 2. 推荐模块与责任

```text
restretto/
  constants.py          # X-Score 类型、半径、能量常量
  geometry.py           # Vector3d, Point3d, ZXZ rotation, axis rotation, C++ round/ceil
  config.py             # DockingConfiguration, config parser, validation
  openbabel_adapter.py  # file IO, bond order fix, xstype mapping, canonical smiles/labels, SDF properties
  molecule.py           # Atom, Bond, Molecule, Fragment
  scoring.py            # EnergyCalculator, strict terms, intra energy
  grid.py               # InterEnergyGrid, .grid/.dx read, write, coordinate conversion
  atom_grid.py          # read/make atom grids
  fragmentation.py      # DecomposeMolecule
  rotations.py          # make_rotations_60, read_rotations
  fragment_grid.py      # FragmentInterEnergyGrid
  fragment_cache.py     # cache container and offline schedule adapter
  optimizer.py          # OptimizerGrid and explicit Optimizer
  rmsd.py               # OpenBabel minimum RMSD wrapper
  cli.py                # atomgrid-gen, conformer-docking, decompose, score-only, intraenergy-only
```

实现顺序应按依赖自底向上推进：常量/几何 -> 分子/IO -> 打分/网格 -> 片段 -> 搜索 -> 输出。

## 3. 固定常量

```python
FLTYPE = numpy.float32
INF_ENERGY = 1e9
LIMIT_ENERGY = 1e2
EPS = 1e-4
XS_TYPE_SIZE = 21
XS_TYPE_H = 22
```

X-Score 类型名、编号、半径、供体/受体/疏水判定必须与 C++ `AtomConstants.hpp` 一致。`.grid` 文件名使用 `xs_strings`：

```text
C_H, C_P, N_P, N_D, N_DC, N_A, N_DA, O_P, O_D, O_A,
O_AC, O_DA, S_P, P_P, F_H, Cl_H, Br_H, I_H, Met_D, Other, Dummy
```

## 4. 几何层实现

必须复现：

1. `utils::round(x)`：`x > 0 ? int(x + 0.5) : int(x - 0.5)`。
2. `ceili(x)`：`int(ceil(x) + EPS)`。
3. `Vector3d.rotate(theta, phi, psi)`：ZXZ 欧拉旋转矩阵，公式照 `Vector3d.cc`。
4. `axis_rotate(axis, th)`：Rodrigues 公式。
5. `make_rotations_60()`：严格复现黄金比例 pole 和 `order` 数组。
6. `read_rotations(path)`：每行解析 `theta, phi, psi`，少于 3 列跳过。

测试要求：

1. `len(make_rotations_60()) == 60`。
2. 首尾若干旋转角与 C++ 输出的参考值一致。
3. 正负数 rounding 与 C++ 一致。

## 5. 配置层实现

实现 `DockingConfiguration`，默认值必须一致：

```text
reuse_grid=OFFLINE, reorder=True, poses_per_lig=1,
poses_per_lig_before_opt=2000, output_score_threshold=-3.0,
pose_min_rmsd=0.5, no_local_opt=False, score_only=False,
local_only=False, local_max_rmsd=1e10, rad_scale=0.95
```

解析规则：

1. 行首精确匹配关键字和空格。
2. 三维参数用逗号分隔并 trim。
3. bool 参数只接受大小写规整后的 `true/false`，`REUSE_FRAG_GRID` 额外接受 `false/none/noreuse/online`，其他值为 `OFFLINE`。
4. 多个 `LIGAND` 行追加。
5. CLI 覆盖配置文件。

校验规则：

1. `SEARCH_PITCH / SCORING_PITCH` 三轴整数比。
2. 根据 outerbox、score_pitch、mem_size 计算 `FGRID_SIZE`，必须大于 0。
3. `score_only and local_only` 报错。

## 6. OpenBabel 适配层实现

实现以下函数：

```python
parse_file_to_obmols(path_or_paths) -> list[OBMol]
fix_bond_orders(obmol) -> None
get_xs_type(obatom) -> int
to_molecule(obmol) -> Molecule
to_obmol(fragment_or_molecule, original_obmol, capping_atomic_num=-1, capping_for_carbon=False, insert_fragment_id_to_isotope=False) -> OBMol
canonical_smiles(obmol) -> str
canonical_labels(obmol) -> list[int]
update_coords(obmol, molecule) -> None
set_property(obmol, key, value) -> None
```

关键行为：

1. `parse_file_to_obmols` 根据扩展名设置输入格式，逐个读分子。
2. 每个分子读入后调用 `fix_bond_orders` 和 `AddPolarHydrogens()`。原 C++ 目标环境是 OpenBabel 2.4.1；若 Python 环境使用 OpenBabel 3.x 且 PDBQT 受体调用 `AddPolarHydrogens()` 会导致绑定层崩溃，则对 PDBQT 保留文件自带氢并记录这是环境兼容分支，不改变 MOL2/SDF ligand 的极性氢添加流程。
3. `to_molecule` 使用 OpenBabel atom id 作为内部 atom id；bond 的 `is_rotor = bond.IsRotor() or begin/end is H`。
4. `get_xs_type` 的条件顺序必须按 C++ `OBMol.cc`。
5. `canonical_smiles` 使用 `can` 输出格式和 `n` 选项，再 trim 尾部空白。对大分子 PDBQT 受体，在 OpenBabel 3.x Python 绑定中 canonical SMILES 可能崩溃；由于 receptor smiles 不参与 fragment reuse 或输出去重，Python 兼容实现可对超过 200 原子的受体使用 title 作为占位 smiles。
6. `to_obmol` 删除 fragment 外原子；dummy 设为 H 或 atomic num 0；可选 capping 和 isotope 标记。

## 7. 分子模型实现

`Atom` 继承或组合 `Vector3d`，字段：

```text
id, x, y, z, xs_type
```

`Bond` 字段：

```text
atom_id1, atom_id2, is_rotor
```

`Molecule` 必须支持：

```text
translate, rotate, axis_rotate, append_atom, append_bond, append_molecule,
get_center, get_radius, delete_hydrogens, renumbering, is_renumbered,
get_graph_distances, get_nrots, calc_rmsd, set/get_intra_energy
```

实现要点：

1. `get_center()` 只统计非氢、非 dummy 原子。
2. `get_radius()` 只统计重原子到中心的最大距离。
3. `delete_hydrogens()` 删除氢后必须修正原子 id、键端点和邻接表。
4. `get_nrots()` 统计可旋转键数量，保持 C++ 对氢键 rotor 的处理。

`Fragment` 扩展 `Molecule`：

```text
idx, smiles, tri[3], theta, phi, psi
normalize_pose(), get_rot(), gettri()
```

`normalize_pose()` 的参考原子选择和旋转计算要从 `Fragment.cc` 完整搬运。

## 8. 打分层实现

`EnergyCalculator(rad_scale=None)`：

1. 如果传入 `rad_scale`，预计算 `21*21*8000` 个 pair energy，距离步长 `1/1000` A，距离阈值 8 A。
2. `get_energy_atom_atom` 查表；超过 8 A 返回 0。
3. `get_energy_atom_molecule` 跳过 receptor 氢并累加，达到 100 返回 100。
4. `get_energy_molecule_molecule` 跳过 ligand 氢并累加，达到 100 返回 100。
5. `calc_intra_energy` 用 graph distance 跳过小于 4 的原子对，跳过氢和 dummy。
6. strict 静态函数直接按公式计算，用于 `score-only` 和 `intraenergy-only`。

测试：

1. 五项函数单元测试覆盖 hydrophobic、hydrogen bond、repulsion 的分段边界。
2. `score-only` 对 `testdata/G39.mol2` 的输出和 C++ 参考值在容差内一致。

## 9. 网格层实现

`InterEnergyGrid`：

```text
center: Point3d[float]
pitch: Point3d[float]
num: Point3d[int]
grid: numpy.ndarray shape (num.x, num.y, num.z), dtype=float32
```

方法：

```text
set_inter_energy(x,y,z,val)
add_energy(x,y,z,val)
get_inter_energy_idx(x,y,z)
get_inter_energy_pos(pos)
convert_x/y/z(pos)
convert(idx)
parse_grid(path)
write_grid(path)
parse_dx(path)
```

`.grid` 读写：

1. 读写 3 个 float32 的 center。
2. 读写 3 个 float32 的 pitch。
3. 读写 3 个 int32 的 num。
4. 按 x/y/z 三重循环读写 float32。

`AtomInterEnergyGrid`：

```text
read_atom_grids(folder) -> list length 21
read_dx_atom_grids(folder) -> list of existing dx grids
make_atom_grids(center, pitch, num, receptor_mol, energy_calculator)
```

`atomgrid-gen` 用 `make_atom_grids` 并写入 `<xs_name>.grid`。

## 10. 片段分解实现

完整搬运 `DecomposeMolecule`：

1. 对不可旋转键并查集合并。
2. ring detector 找环，按 `max_ring_size` 合并环系统。
3. 统计相邻非氢数量，合并小的孤立非氢原子。
4. 若 `merge_solitary`，按 bond 顺序尝试合并片段，要求：
   - 合并前后 ring 数量不增加。
   - 合并片段内部 rotor 旋转 1 rad 后 RMSD 小于 `1e-5`。
5. 合并所有氢到片段。
6. 为跨片段重原子键添加 dummy 原子。
7. 生成 `Fragment(id, atoms)` 并添加片段内键。

端到端测试：

1. 对 `G39.mol2` 运行 `decompose`，fragment 数量、canonical SMILES 集合、`fragment_info` 与 C++ 参考一致。
2. `--no_merge_solitary`、`--max_ring_size`、capping 选项各有最小测试。

## 11. 片段库和缓存实现

普通 docking 中：

1. 对每个 ligand 分解。
2. 对每个 fragment：
   - `to_obmol(fragment, original_ligand)`。
   - `canonical_smiles`。
   - `canonical_labels` 后 `fragment.renumbering(...)`。
   - 相同 SMILES 复用同一 `frag_idx`。
   - 新 fragment 复制后 `normalize_pose()` 加入 library。
   - ligand 的 `FragmentsVector` 追加 `(fragment_center, frag_idx, fragment_size)`。
3. 片段重要性累加 `fragment.size()`。

缓存：

1. `NONE`：容量 1。
2. `ONLINE`：LRU。
3. `OFFLINE`：先生成调度数组 `indices_to_save`，再按 step 定位槽位。

如果初期不能完全复现 `CalcMCFP`，可以先实现 ONLINE/NONE 并将 OFFLINE 调度作为待严格补齐项；但最终完整实现必须复现 OFFLINE，因为它是默认策略。

## 12. FragmentInterEnergyGrid 实现

输入：

```text
orig_frag, rot_angles, atom_grids, distance_grid
```

算法：

1. 断言 `orig_frag.get_center().abs() < EPS`。
2. 初始化 grid 为 `LIMIT_ENERGY`。
3. `rotsz = 1 if orig_frag.size() <= 1 else len(rot_angles)`。
4. 对每个旋转复制 fragment 并旋转。
5. 对每个 grid 点：
   - 保留原代码中 `if (rot_id & 7)` 的剪枝条件。
   - distance < 2 跳过。
   - distance > radius + 6 跳过。
   - 对 fragment 非氢原子查询对应 atom grid。
   - 累加到 `LIMIT_ENERGY` 可提前停止。
   - 保留该点所有旋转中的最小能量。

## 13. Conformer docking 完整流程

实现 `run_conformer_docking(config)`：

1. 读取 receptor：`parse_file_to_obmols(config.receptor_file)[0]`，转内部分子。
2. 读取 ligand 文件列表。
3. 读取 atom grids；若 `dxgrid_folder` 不空，用现有 dx 覆盖对应 xs type。
4. 计算：

```text
score_num = atom_grids[0].num
ratio = round(search_pitch / score_pitch)
search_num = ceili(inner_width / 2 / search_pitch) * 2 + 1
search_grid = InterEnergyGrid(atom_grid.center, search_pitch, search_num)
```

5. `no_search = score_only or local_only`。
6. 根据 memory size 计算 `FGRID_SIZE`。
7. 转换 ligand：
   - 普通搜索模式平移到中心为原点。
   - no_search 模式保留输入坐标。
8. 建立 `lig_map`：identifier 到唯一 ligand index。
9. 删除内部分子氢并计算 intra energy。
10. 若普通搜索：
    - 分解 ligand 和建立 fragment library。
    - 按配置重排 ligand。
    - OFFLINE 生成调度。
    - 构造 distance grid。
    - 准备 ligand rotations。
    - 初始化每个唯一 ligand 的 `MinValuesVector(capacity=poses_per_lig_before_opt)`。
    - 对每个 ligand 和每个 fragment 累加 fragment grid。
    - 收集低于 `output_score_threshold` 的候选。
11. 构造 `OptimizerGrid(atom_grids, local_max_rmsd)`。
12. 普通搜索输出：
    - 对每个唯一 ligand 的候选恢复 pose。
    - `mol.rotate(rotations_ligand[rotid])`，再 `mol.translate(search_grid.convert(grid_pos))`。
    - 若 `no_local_opt`，计算 total energy；否则优化。
    - 按 total energy 排序。
    - 计算 best score 并写 CSV。
    - RMSD 去重后写最多 `poses_per_lig` 个 pose。
13. no_search 输出：
    - 对每个输入 ligand 计算 total energy。
    - `score_only` 不优化；`local_only` 优化。
    - 写 `restretto_score` 和输出分子。

## 14. 局部优化实现

`OptimizerGrid.optimize(mol)`：

1. 每次调用前设置随机种子 0，复现 C++ 每个 pose 都从同一随机序列开始的行为。
2. 初始 `opt = calc_inter_energy(mol) + mol.intra_energy`。
3. 循环：
   - 采样 200 个邻域。
   - 当前中心为 `mol.get_center()`。
   - 平移步长 0.5，旋转步长 `pi/30`。
   - 若相对初始分子 RMSD 超过 max_rmsd，跳过。
   - 选择本轮最低能量候选。
   - 若低于当前 opt，接受；否则退出。
4. 返回最终 total energy。

`calc_inter_energy(mol)` 对每个 atom 查 `atom_grids[atom.xs_type].get_inter_energy_pos(atom)` 并累加。严格复现时保留对氢的行为：原网格优化代码没有显式跳过氢，但内部 ligand 在前面已 `deleteHydrogens()`。

## 15. RMSD 去重

实现 `calc_min_rmsd(candidate_obmol, accepted_obmols)`：

1. 候选和参考都删除氢或忽略氢。
2. 统一 aromatic/ring 标志。
3. 用 OpenBabel 同构映射枚举可行 atom map。
4. 返回最小 RMSD。
5. 若 accepted 为空，应返回大于 `MIN_RMSD` 的值，以接受首个 pose。

输出筛选条件是严格大于：

```text
if min_rmsd > config.pose_min_rmsd:
    accept
```

## 16. CLI 命令要求

建议统一为：

```powershell
python -m restretto.cli atomgrid-gen <conf-file> [options]
python -m restretto.cli conformer-docking <conf-file> [options]
python -m restretto.cli decompose --ligand ... --fragment ... --output ... [options]
python -m restretto.cli score-only <conf-file> [options]
python -m restretto.cli intraenergy-only <conf-file> [options]
```

保持配置关键字和输出属性名与原项目一致。

## 17. 验证流程

最小验收：

1. `atomgrid-gen references/restretto/testdata/testgrid.in` 能生成 21 个 `.grid`。
2. 生成的 `.grid` header 与 C++ 参考一致：center、pitch、num。
3. 若有 C++ reference 输出，随机抽样网格点能量在 `1e-4` 到 `1e-3` 容差内一致。
4. `conformer-docking references/restretto/testdata/testgrid.in` 能输出 SDF、CSV、log。
5. 输出 SDF 每个 pose 包含 `restretto_score`。
6. `decompose` 对 `G39.mol2` 能输出 fragment 文件和 annotated SDF，annotated SDF 有 `fragment_info`。
7. `score-only` 和 `intraenergy-only` 输出 title 和数值。

回归测试优先级：

1. `utils::round`、旋转矩阵、`.grid` 读写。
2. X-Score 类型判定和五项公式。
3. 配置解析和模式冲突。
4. atom grid 小盒子生成。
5. fragment decomposition。
6. fragment grid 单片段。
7. docking end-to-end。

## 18. 严格遵守原项目的注意事项

1. 默认主流程是 fragment reuse，不是全原子 docking。
2. `OUTERBOX` 控制 atom/fragment scoring grid，`INNERBOX` 控制 ligand center 搜索 grid。
3. search grid 到 score grid 的映射必须用整数 `ratio`。
4. 所有 grid 查询都是 nearest-neighbor。
5. `.grid` 写入 float32，数组展开顺序是 x/y/z。
6. 片段 canonical SMILES 和重编号决定 fragment reuse，不能用不稳定 hash 替代。
7. 局部优化每个 pose 固定种子，且只接受改进。
8. 最终 ligand ranking 用旋转键惩罚后的 inter score，而候选优化排序用 total energy。
9. RMSD 去重使用 OpenBabel 同构意义下的最小 RMSD，不是简单同序坐标 RMSD。
10. date 字符串的 month 若严格复现 C++，保持 `tm_mon` 的 0-based 行为。
## 19. Current Python parity status

Implemented in this pass:

1. `FragmentInterEnergyGridContainer` now matches the C++ cache semantics for `ONLINE` LRU refresh, `NONE` capacity-1 behavior through caller selection, and `OFFLINE` scheduled-slot lookup.
2. `CalcMCFP::runLeftBackSSP` and `fraggrid_main.cc::makeGraph` have a Python port exposed as `build_offline_schedule(...)`.
3. Default `conformer-docking` now enters the fragment-grid search path instead of the previous full-ligand atom-grid placeholder. `score_only` and `local_only` still use the direct atom-grid path, matching the C++ `no_search` split.
4. Fragment-grid scoring now accumulates per-fragment grids over ligand rotations and search-grid points, using `OUTERBOX` as the scoring grid and `INNERBOX` as the search grid.

Known remaining C++ gaps:

1. OpenBabel isomorphism RMSD pose deduplication is still not fully ported.
2. Full local optimization on all 2000 default pre-opt candidates is still too slow in Python; reduced-cap experiments reach a pose within about 0.18 A same-order RMSD of the C++ reference pose, but the score remains about 0.22 kcal/mol higher.
3. Candidate pruning uses the Python fallback best-pose behavior when no point beats `OUTPUT_SCORE_THRESHOLD`; C++ still writes a CSV rank row in the empty-candidate case, so this keeps Python CLI output useful but is not yet byte-equivalent.

Additional parity work completed after this checkpoint:

1. `InterEnergyGrid` storage is NumPy `float32`, with `values3d()` for C++-style grid slices.
2. `makeDistanceGrid` and fragment-grid search accumulation are vectorized; full testdata no-local-opt runs in about 50 seconds on this environment.
3. `Fragment::normalize_pose()` reference-triplet rotation has been ported.
4. Linux/glibc `rand()` sequence is reproduced for `OptimizerGrid`.
5. OpenBabel 3 carboxylate canonical SMILES is normalized to the OpenBabel 2.4.1 reference form for the bundled G39 testdata.
6. SDF bond order read/write is preserved through `Bond.order`.
