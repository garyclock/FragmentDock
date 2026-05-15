# REstretto 当前实现流程文档：代码工程视角

本文从代码结构、数据流、接口和可执行程序角度整理 `references/restretto` 的当前实现，用于后续解析和 Python 重构。

## 1. 构建与产物

原项目是 C++11 项目，使用 Makefile 构建。主要依赖：

1. Boost program_options、regex、algorithm、lexical_cast、format、unit_test。
2. OpenBabel 2.4.1。README 明确 OpenBabel 3.x 不支持。
3. 可选 OpenMP。

Makefile 默认 `ALL = atomgrid-gen conformer-docking decompose`。其他源码存在但默认不构建，包括 `atom-docking`、`score-only`、`intraenergy-only`、`easytest-docking`、`unittest`。

## 2. 源码模块边界

| 文件 | 职责 |
| --- | --- |
| `common.hpp` | 全局浮点类型和能量常量：`fltype=float`、`INF_ENERGY=1e9`、`LIMIT_ENERGY=1e2`、`EPS=1e-4` |
| `AtomConstants.hpp` | X-Score 类型、半径、供体/受体/疏水判定 |
| `Vector3d.*` | 向量运算、ZXZ 欧拉旋转、轴旋转、默认 60 旋转、自定义旋转读取 |
| `Point3d.hpp` | 三维泛型点，支持逐元素运算 |
| `Atom.*` | 带 id、坐标、X-Score 类型的原子 |
| `Molecule.*` | 原子/键容器、变换、重编号、删氢、图距离、RMSD、rotor 计数、内能存储 |
| `Fragment.*` | `Molecule` 子类，记录 fragment id、canonical SMILES、规范化参考三原子和标准姿态 |
| `OBMol.*` | OpenBabel 读写、键阶修正、X-Score 类型映射、内部分子互转、canonical SMILES、坐标回写、SDF 属性 |
| `infile_reader.*` | 配置文件解析、默认值、模式冲突和 pitch/memory 校验 |
| `EnergyCalculator.*` | Vina/X-Score 风格能量函数、查表预计算、内能计算 |
| `InterEnergyGrid.*` | 通用 3D 能量网格、`.grid`/`.dx` 读写、坐标索引转换 |
| `AtomInterEnergyGrid.*` | 单 X-Score 类型 atom grid 的读取和生成 |
| `MoleculeToFragments.*` | 配体分解为 fragment，包含 ring 检测、并查集合并、dummy 添加 |
| `FragmentInterEnergyGrid.*` | 单个唯一 fragment 的最小旋转能量网格 |
| `FragmentInterEnergyGridContainer.hpp` | fragment grid 缓存，支持 OFFLINE/ONLINE/NONE |
| `FragmentsVector.hpp` | ligand 的 fragment 相对位置列表和排序 |
| `CalcMCFP.*` | OFFLINE 复用策略的缓存调度求解 |
| `MinValuesVector.hpp` | 固定容量低分候选容器 |
| `Optimizer.*` | 显式受体能量优化器和 atom-grid 优化器 |
| `RMSD.*` | OpenBabel 同构映射 RMSD，用于 pose 去重 |
| `log_writer_stream.*` | log 文件输出 |

## 3. 配置模型

`format::DockingConfiguration` 是主配置对象：

```text
grid.center
grid.outer_width
grid.inner_width
grid.search_pitch
grid.score_pitch
ligand_files[]
receptor_file
output_file
log_file
grid_folder
dxgrid_folder
rotangs_file
reuse_grid = OFFLINE
reorder = true
mem_size
poses_per_lig = 1
poses_per_lig_before_opt = 2000
output_score_threshold = -3.0
pose_min_rmsd = 0.5
no_local_opt = false
score_only = false
local_only = false
local_max_rmsd = 1e10
rad_scale = 0.95
```

配置文件是大小写敏感的行格式。已实现关键字：

```text
INNERBOX, OUTERBOX, BOX_CENTER, SEARCH_PITCH, SCORING_PITCH,
REUSE_FRAG_GRID, REORDER_LIGANDS, MEMORY_SIZE, RECEPTOR, LIGAND,
OUTPUT, LOG, GRID_FOLDER, ROTANGS, POSES_PER_LIG,
POSES_PER_LIG_BEFORE_OPT, OUTPUT_SCORE_THRESHOLD, MIN_RMSD,
NO_LOCAL_OPT, SCORE_ONLY, LOCAL_ONLY, LOCAL_MAX_RMSD, DXGRID_FOLDER
```

命令行参数会覆盖配置文件字段。主 docking 在 parseArgs 末尾调用 `checkConfigValidity()`。

校验规则：

1. `SEARCH_PITCH / SCORING_PITCH` 每轴必须近似整数。
2. `MEMORY_SIZE` 必须能容纳至少一个 fragment grid。
3. `score_only` 与 `local_only` 不能同时为真。
4. `score_only && no_local_opt`、`local_only && no_local_opt` 只写警告。

## 4. 数据结构不变量

### `Molecule`

1. `atoms[i].id` 以 OpenBabel atom id 为准，原代码假定初始 `i == oatom.GetId()`。
2. `bond_ids` 是按原子 id 建立的邻接键索引。
3. `identifier = title + "," + canonical_smiles`，用于 ligand 构象聚合。
4. `deleteHydrogens()` 会删除氢并重编号；调用后坐标数组和 OpenBabel 分子坐标回写需保持同序。
5. `getCenter()` 只考虑非氢、非 dummy 原子。

### `InterEnergyGrid`

1. 内部数组一维展开索引：`(x * num.y + y) * num.z + z`。
2. `.grid` 直接写 C++ struct 的二进制布局。Python 重构必须固定小端/float32/int32 并用回归样例校验。
3. 坐标查询是 nearest-neighbor，不做插值。
4. 越界返回 `LIMIT_ENERGY`。

### `Fragment`

1. fragment 必须先 canonical renumber，再调用 `settri()` 和 `normalize_pose()`。
2. 规范化后片段中心应接近原点。
3. dummy 原子参与片段几何和 repulsion，但在多数能量项中被跳过。

## 5. 可执行程序流程

### `atomgrid-gen`

入口：`src/grid_main.cc`

命令：

```text
atomgrid-gen conf-file [--grid path] [--receptor path] [--rad_scale value]
```

流程：

1. 解析配置并应用 CLI 覆盖。
2. 创建 `GRID_FOLDER`。
3. log 固定写 `atomgrid-gen.log`。
4. 读取 receptor，取第一个分子。
5. 转换为内部 `Molecule`。
6. 构造 `EnergyCalculator(conf.rad_scale)`。
7. 对 0 到 `XS_TYPE_SIZE-1` 的每个 X-Score 类型：
   - 根据 outerbox 和 scoring pitch 计算网格点数。
   - 初始化 `AtomInterEnergyGrid(center, pitch, num, xs_type)`。
   - 遍历所有网格点，将探针原子放到点坐标。
   - 调 `calc.getEnergy(atom, receptor_mol)`。
   - 写出 `GRID_FOLDER/<xs_name>.grid`。

### `conformer-docking`

入口：`src/fraggrid_main.cc`

命令：

```text
conformer-docking conf-file [options]
```

重要 CLI 覆盖：

```text
--output, --ligand, --receptor, --grid, --memsize,
--score-only, --no-local-opt, --local-only, --local-max-rmsd,
--log, --poses-per-lig, --min-rmsd,
--poses-per-lig-before-opt, --output-score-threshold, --dxgrid
```

主流程：

1. 解析配置并校验。
2. 设置 log 文件，默认 `OUTPUT + "fraggrid__<date>.log"`。
3. 读取 receptor 和所有 ligand。
4. 读取 `GRID_FOLDER` 下的 21 个 `.grid`。
5. 若 `DXGRID_FOLDER` 非空，则读取其中存在的 `<xs_name>.dx` 并覆盖相应 atom grid。
6. 计算 `score_num`、`ratio`、`search_num`，构造 search grid。
7. 根据 memory size 计算可缓存 fragment grid 数量；`NONE` 强制为 1。
8. 将 OpenBabel ligand 转为内部 `Molecule`：
   - `score_only/local_only` 模式不把 ligand 平移到原点。
   - 普通搜索模式把 ligand 平移到 `-getCenter()`。
9. 对每个 ligand 删除氢并计算 intra energy。
10. 若非 `score_only/local_only`，进入片段搜索：
    - 分解 ligand 为 fragment。
    - fragment 转 OBMol，计算 canonical SMILES，canonical renumber。
    - 建立唯一 fragment library 和每个 ligand 的 `FragmentsVector`。
    - 根据配置重排 ligand 和片段顺序。
    - OFFLINE 模式用 `CalcMCFP` 生成缓存调度。
    - 构造 distance grid。
    - 准备 ligand 旋转集合。
    - 按 sorted ligand 遍历并累加 fragment grids，收集候选。
11. 构造 `Optimizer_Grid(atom_grids, local_max_rmsd)`。
12. 若普通搜索：
    - 对每个 ligand identifier 的候选重建坐标。
    - 计算或优化 total energy。
    - 排序、计算 score、RMSD 去重、写输出。
13. 若 `score_only/local_only`：
    - 直接对每个输入构象计算 total energy。
    - `score_only` 不优化，`local_only` 优化。
    - 写 `restretto_score` 和输出构象。
14. 输出 ranking、final statistics 和关闭 log。

### `decompose`

入口：`src/decompose_main.cc`

命令：

```text
decompose --ligand ligands... --fragment fragments.sdf --output annotated.sdf [options]
```

选项：

```text
--log
--capping_atomic_num default -1
--enable_carbon_capping
--ins_fragment_id
--max_ring_size default -1
--no_merge_solitary
```

流程：

1. 解析必需的 ligand、fragment、output。
2. 设置 log，默认 `OUTPUT + "__<date>.log"`。
3. 读取所有 ligand。
4. 对每个 ligand：
   - 添加氢并转内部分子。
   - 执行 `DecomposeMolecule`。
   - 每个 fragment 转 OBMol，必要时 capping 或写 isotope。
   - 修正 valence，生成 canonical SMILES。
   - 唯一 fragment 写入 fragments 列表。
   - 原 ligand 写入 `fragment_info` 属性。
5. 写 annotated ligand 和 fragment 文件。

### `score-only`

入口：`src/score_only_main.cc`

默认 Makefile 不构建，但源码完整。流程是读取 receptor 和 ligand，用严格五项公式直接输出：

```text
Title: ...
Affinity: ...
```

### `intraenergy-only`

入口：`src/intraenergy_main.cc`

默认 Makefile 不构建。流程是读取 ligand，转换内部分子、删氢，用严格内能公式输出：

```text
Title: ...
IntraEnergy: ...
```

### `atom-docking`

入口：`src/nofrag_main.cc`

默认不构建。它不用 fragment reuse，而是对 ligand 全原子在 60 旋转和 search grid 上直接累加 atom grids。该文件与当前 `EnergyCalculator` 构造函数签名存在旧接口痕迹，重构时可作为对照流程，不应优先作为主功能基准。

## 6. 主 docking 的关键算法接口

### `convert_molecules(obmols, no_search)`

输入 OpenBabel molecules，输出内部 `Molecule[]`。

1. 对每个 ligand 添加氢。
2. 转内部 `Molecule`。
3. 删除 OpenBabel 分子中的氢。
4. 普通搜索模式把内部分子移到中心在原点。

### `makeDistanceGrid(center, pitch, num, receptor_mol)`

遍历网格点，计算该点到 receptor 非氢原子的最小欧氏距离。供 `FragmentInterEnergyGrid` 剪枝。

### `FragmentInterEnergyGridContainer`

统一缓存接口：

```text
isRegistered(fragid)
insert(grid)
get(fragid)
next()
```

OFFLINE 下，`isRegistered` 只检查当前 step 对应槽位是否为目标 fragid；`insert` 写入 `indices_to_save[step]`；`next()` 每处理一个 fragment 调一次。

ONLINE 下，`isRegistered` 线性搜索缓存；`insert` 替换 LRU；`get` 更新 last_used。

## 7. 输出文件和日志约定

1. `atomgrid-gen` 写 `atomgrid-gen.log` 和 `GRID_FOLDER/*.grid`。
2. `conformer-docking` 写 `OUTPUT`、`OUTPUT + "fraggrid__<date>.csv"` 和 log。
3. `decompose` 写 annotated output、fragment file 和 log。
4. 主 docking 的输出分子包含属性 `restretto_score`。
5. `decompose` 的 annotated ligand 包含属性 `fragment_info`。

`getDate()` 格式是：

```text
%02d_%02d_%02d_%02d_%02d = month, day, hour, minute, second
```

注意 C++ 的 `tm_mon` 是 0-based，原实现没有加 1。严格复现时应保留。

## 8. Python 重构建议模块划分

建议以原模块边界为基础：

```text
restretto/
  constants.py
  geometry.py
  molecule.py
  openbabel_adapter.py
  config.py
  scoring.py
  grid.py
  atom_grid.py
  fragmentation.py
  fragment_grid.py
  fragment_cache.py
  rotations.py
  optimizer.py
  rmsd.py
  cli.py
```

每个模块应有对应测试：

1. `constants.py`：类型映射、半径、供受体判定。
2. `geometry.py`：ZXZ 旋转、axis rotate、round/ceil 行为。
3. `config.py`：配置解析和命令行覆盖。
4. `scoring.py`：五项公式、查表、内能。
5. `grid.py`：`.grid` 读写、坐标索引、越界。
6. `fragmentation.py`：小分子片段分解、dummy 添加。
7. `fragment_grid.py`：单片段网格值与 C++ 小样例一致。
8. `optimizer.py`：固定随机种子和能量下降。
9. `cli.py`：端到端 `atomgrid-gen`、`conformer-docking`、`decompose`。

## 9. 兼容性和风险点

1. OpenBabel 版本差异会影响 atom typing、rotor、canonical SMILES、canonical labels 和 RMSD。
2. C++ `.grid` 直接序列化 struct，跨平台二进制布局有风险。Python 应固定为与当前测试平台一致，并用 reference grid 做字节级或数值级测试。
3. `utils::round` 对正负数使用 away-from-zero 风格，不是 Python `round()` 的 bankers rounding。
4. `float` 是 32 位。Python 内部计算若用 double，输出可能略有差异；写 grid 时必须 float32。
5. 局部优化使用 C `rand()` 和 `srand(0)`。若要求严格数值一致，Python 不能直接用 `random` 替代，需要复现 C runtime rand 或接受近似一致。
6. `atom-docking` 源码含旧构造函数痕迹，默认未构建；重构优先级低于 README 暴露的 `atomgrid-gen`、`conformer-docking`、`decompose`。
## 10. Current Python implementation checkpoint

This checkpoint records the current Python state against the C++ engineering flow:

1. `restretto.fragment_grid.FragmentInterEnergyGridContainer` follows the C++ container behavior for `ONLINE` LRU replacement, `OFFLINE` scheduled-slot lookup, and caller-enforced `NONE` capacity 1.
2. `restretto.fragment_grid.build_offline_schedule` is a direct Python port of `fraggrid_main.cc::makeGraph` plus `CalcMCFP::runLeftBackSSP`.
3. `restretto.cli.conformer_docking` now uses the fragment-grid path for normal docking and keeps the direct atom-grid path for `score_only` and `local_only`.
4. The Python fragment library now applies OpenBabel canonical SMILES and canonical labels before assigning fragment ids.
5. Exact `Fragment::normalize_pose()` reference-triplet rotation has been ported.
6. OpenBabel isomorphism RMSD deduplication remains pending.
7. Full local optimization with the default 2000 pre-opt candidates remains the main runtime gap; full no-local-opt testdata now finishes in about 50 seconds, while optimized reduced-cap runs are chemically close to the C++ reference pose but still about 0.22 kcal/mol high.
