# REstretto 当前实现流程与模块级接口说明

本文档面向后续解析、复现与重构，目标是把 `references/restretto` 当前源码中的主要行为、数据流、模块边界和依赖关系整理成一份可以直接照着实现的说明。

本文档覆盖的范围是当前仓库中的主要可执行程序与公共模块，不覆盖历史分支或论文层面的抽象实现。

## 1. 项目定位

REstretto 是一个面向虚拟筛选的蛋白-配体对接工具。当前实现的核心思想是：

- 先对受体预计算原子类型能量网格。
- 再把配体拆成 fragment。
- 通过 fragment grid 重用减少重复计算。
- 在搜索后做局部优化和 RMSD 去重。
- 输出构象文件和分数。

当前源码中的主要可执行文件如下：

| 可执行文件 | 作用 |
| --- | --- |
| `atomgrid-gen` | 生成受体原子相互作用能量网格 |
| `conformer-docking` | 主 docking 流程，包含片段分解、fragment reuse、搜索、局部优化、输出 |
| `decompose` | 将配体分解为 fragment，并写回 fragment 信息 |
| `atom-docking` | 不使用 fragment reuse 的原子级对照流程 |
| `score-only` | 只做打分，不做搜索 |
| `intraenergy-only` | 只算配体内能 |
| `easytest-docking` | 简化调试输出 |
| `unittest` | Boost 单元测试 |

## 2. 依赖清单

### 2.1 构建依赖

| 依赖 | 说明 |
| --- | --- |
| C++11 编译器 | `Makefile` 默认使用 `g++`，也可通过 `CXX` 覆盖 |
| GNU Make | 使用 `make all` 构建 |
| Boost headers/libs | 代码与链接都依赖 Boost |
| OpenBabel 2.4.x | README 明确要求 OpenBabel 2.4.1，OpenBabel 3.x 不支持 |
| zlib | `STATIC=Y` 时链接参数包含 `-lz` |
| OpenMP（可选） | `OpenMP=Y` 时启用 `-fopenmp` |
| Unix 风格工具链 | `mkdir`、`rm` 等命令被 `Makefile` 使用 |

### 2.2 必须链接/使用的 Boost 组件

| 组件 | 使用位置 |
| --- | --- |
| `boost_program_options` | 命令行参数解析 |
| `boost_regex` | `Makefile` 中的通用链接依赖 |
| `boost_algorithm` | 字符串处理、split、trim、join |
| `boost_lexical_cast` | 配置文件数值解析 |
| `boost_format` | 日期和日志字符串格式化 |
| `boost_iostreams` | 日志流封装 |
| `boost_unit_test_framework` | 单元测试 |

### 2.3 OpenBabel 依赖

| 依赖 | 说明 |
| --- | --- |
| `openbabel/mol.h` | 分子对象 `OBMol` |
| `openbabel/obconversion.h` | 文件读写与格式转换 |
| `openbabel/griddata.h` | OpenDX 网格读取 |
| `openbabel/canon.h` | canonical smiles / 重编号 |
| `openbabel/babelconfig.h`、`op.h`、`graphsym.h`、`query.h` | RMSD、同构、结构标准化支持 |

### 2.4 运行依赖

| 依赖 | 说明 |
| --- | --- |
| 输入 receptor 文件 | 通常为 `.pdb` 或 OpenBabel 可读格式 |
| 输入 ligand 文件 | `.mol2`、`.sdf` 等 OpenBabel 可读格式 |
| 预计算的 atom grid 文件 | `GRID_FOLDER/<xs_type>.grid`，或可选 `.dx` 覆盖 |
| 片段重用配置 | `REUSE_FRAG_GRID`、`MEMORY_SIZE`、`REORDER_LIGANDS` |

### 2.5 测试依赖

| 依赖 | 说明 |
| --- | --- |
| Boost Unit Test Framework | `unittest` 目标使用 |
| `testdata/` | 提供最小可运行样例 |
| `conformer-docking`、`atomgrid-gen` 生成物 | 主流程需要先生成网格再跑 docking |

### 2.6 独立环境建议

当前仓库提供了 devcontainer 入口：

- `.devcontainer/Dockerfile` 基于 `kyanagis/restretto:1.1`
- `.devcontainer/devcontainer.json` 指定 VS Code 容器开发环境

如果要在独立虚拟环境里复现，推荐策略是：

| 方案 | 说明 |
| --- | --- |
| Docker / devcontainer | 最稳妥，直接复用项目原始依赖栈 |
| 自建 Linux 容器 | 安装相同版本 OpenBabel 2.4.x 与 Boost |
| 本地系统环境 | 仅适合已经具备兼容 OpenBabel 2.4.x 的场景 |

## 3. 目录与模块边界

| 目录/文件 | 模块职责 |
| --- | --- |
| `src/common.hpp` | 通用类型、能量上下限、日志头文件 |
| `src/AtomConstants.hpp` | X-Score 原子类型、半径、供体/受体判定 |
| `src/Vector3d.*` | 3D 向量运算、旋转、轴旋转 |
| `src/Point3d.hpp` | 三维离散/连续坐标模板 |
| `src/Atom.*` | 带类型的原子对象 |
| `src/Molecule.*` | 原子、键、图距离、变换、RMSD 辅助 |
| `src/Fragment.*` | fragment 姿态标准化与中心化 |
| `src/OBMol.*` | OpenBabel 与自定义分子对象互转 |
| `src/infile_reader.*` | 配置文件解析与合法性校验 |
| `src/InterEnergyGrid.*` | 3D 能量网格抽象，支持 `.grid`/`.dx` |
| `src/AtomInterEnergyGrid.*` | 单原子类型的受体网格 |
| `src/FragmentInterEnergyGrid.*` | 单 fragment 的能量网格 |
| `src/FragmentInterEnergyGridContainer.*` | fragment grid 缓存与重用策略 |
| `src/MoleculeToFragments.*` | 分子分解成 fragment |
| `src/EnergyCalculator.*` | X-Score 风格打分函数 |
| `src/Optimizer.*` | 局部优化 |
| `src/RMSD.*` | 结构去重与最小 RMSD |
| `src/log_writer_stream.*` | 日志输出 |

## 4. 配置模型

### 4.1 配置文件语义

配置文件按行读取，使用关键字 + 空格 + 参数的形式。当前实现支持的主要字段如下：

| 字段 | 含义 |
| --- | --- |
| `INNERBOX X, Y, Z` | 配体中心搜索范围 |
| `OUTERBOX X, Y, Z` | 原子网格覆盖范围 |
| `BOX_CENTER X, Y, Z` | 网格中心 |
| `SEARCH_PITCH X, Y, Z` | 搜索步长 |
| `SCORING_PITCH X, Y, Z` | 评分网格步长 |
| `MEMORY_SIZE N` | fragment grid 缓存容量，单位 MB |
| `RECEPTOR path` | 受体文件路径 |
| `LIGAND path` | 配体文件路径，可多行 |
| `OUTPUT path` | 输出文件路径 |
| `GRID_FOLDER path` | `.grid` 文件目录 |
| `DXGRID_FOLDER path` | 可选 OpenDX 网格目录，用于覆盖 `.grid` |
| `ROTANGS path` | 自定义旋转角文件 |
| `REUSE_FRAG_GRID` | `OFFLINE`、`ONLINE` 或 `NONE` |
| `REORDER_LIGANDS` | 是否按片段重要度重排 |
| `POSES_PER_LIG` | 每个 ligand 最终输出构象数 |
| `POSES_PER_LIG_BEFORE_OPT` | 局部优化前保留候选数 |
| `OUTPUT_SCORE_THRESHOLD` | 搜索阶段的候选分数阈值 |
| `MIN_RMSD` | 输出姿态之间的最小 RMSD |
| `NO_LOCAL_OPT` | 是否跳过局部优化 |
| `SCORE_ONLY` | 是否只评分不搜索 |
| `LOCAL_ONLY` | 是否只做局部搜索 |
| `LOCAL_MAX_RMSD` | 局部优化允许的最大 RMSD |

### 4.2 校验规则

`DockingConfiguration::checkConfigValidity()` 强制检查以下约束：

- `SEARCH_PITCH / SCORING_PITCH` 每个维度都必须是整数比。
- `MEMORY_SIZE` 必须足够容纳至少一个 fragment grid。
- `SCORE_ONLY` 与 `LOCAL_ONLY` 不能同时为真。
- `NO_LOCAL_OPT` 与 `SCORE_ONLY`、`LOCAL_ONLY` 的组合会触发警告语义，但仍保持现有实现行为。

## 5. 模块级接口说明

### 5.1 `format::DockingConfiguration`

| 成员 | 类型 | 说明 |
| --- | --- | --- |
| `grid` | `SearchGrid` | 包含 center、outer/inner box、search/scoring pitch |
| `ligand_files` | `std::vector<std::string>` | 配体文件列表 |
| `receptor_file` | `std::string` | 受体文件 |
| `output_file` | `std::string` | 输出文件 |
| `log_file` | `std::string` | 日志文件 |
| `grid_folder` | `std::string` | atom grid 目录 |
| `dxgrid_folder` | `std::string` | 可选 OpenDX grid 目录 |
| `rotangs_file` | `std::string` | 自定义旋转角文件 |
| `reuse_grid` | `ReuseStrategy` | `OFFLINE` / `ONLINE` / `NONE` |
| `reorder` | `bool` | 是否重排 ligand |
| `mem_size` | `int64_t` | fragment grid 缓存大小 MB |
| `poses_per_lig` | `int64_t` | 最终输出构象数 |
| `poses_per_lig_before_opt` | `int64_t` | 优化前候选数 |
| `output_score_threshold` | `fltype` | 候选阈值 |
| `pose_min_rmsd` | `fltype` | 输出去重阈值 |
| `no_local_opt` | `bool` | 跳过局部优化 |
| `score_only` | `bool` | 仅评分 |
| `local_only` | `bool` | 仅局部搜索 |
| `local_max_rmsd` | `fltype` | 局部搜索上限 |
| `rad_scale` | `fltype` | 能量预计算半径缩放 |

### 5.2 `format::ParseInFile(const char*)`

输入：配置文件路径。  
输出：填充后的 `DockingConfiguration`。  
行为：

- 逐行解析。
- 识别固定关键字。
- 覆盖命令行传入的字段。
- 最终调用 `checkConfigValidity()`。

### 5.3 `format::ParseFileToOBMol(...)`

输入：

- 单个文件路径
- 文件路径数组
- 输入流 + 格式字符串

输出：`std::vector<OpenBabel::OBMol>`。  
行为：

- 按格式读取所有分子。
- 修正部分 bond order。
- 添加 polar hydrogens。

### 5.4 `format::toFragmentMol(const OBMol&)`

输入：OpenBabel 分子。  
输出：`fragdock::Molecule`。  
行为：

- 为每个原子建立 `Atom(id, pos, xs_type)`。
- 通过 OpenBabel 键对象建立 `Bond`。
- rotor 键由 `IsRotor()` 或氢键接触推断。
- 保存 title 和 canonical smiles。

### 5.5 `fragdock::Molecule`

核心接口：

| 接口 | 作用 |
| --- | --- |
| `translate()` | 整体平移 |
| `rotate()` | 欧拉角旋转 |
| `axisRotate()` | 围绕任意轴旋转 |
| `append()` | 添加原子、键或子分子 |
| `getCenter()` | 计算非氢/非 dummy 原子中心 |
| `getRadius()` | 计算重原子半径 |
| `deleteHydrogens()` | 删除氢并重编号 |
| `getGraphDistances()` | 计算图最短距离 |
| `getNrots()` | 计算可旋转键数 |
| `calcRMSD()` | 结构 RMSD，带图一致性检查 |
| `setIntraEnergy()` / `getIntraEnergy()` | 保存配体内能 |

### 5.6 `fragdock::Fragment`

核心接口：

| 接口 | 作用 |
| --- | --- |
| `normalize_pose()` | 把 fragment 平移到原点并旋转到标准姿态 |
| `getRot()` | 返回标准化旋转参数 |
| `gettri()` | 返回用于定义姿态的三原子索引 |
| `setIdx()` / `getIdx()` | fragment 库编号 |
| `setSmiles()` / `getSmiles()` | 归一化后唯一标识 |

规范化逻辑：

- 单原子 fragment 只记录一个参考点。
- 多原子 fragment 选择三原子构型定义标准旋转。
- `normalize_pose()` 会先平移到中心，再计算并施加标准旋转。

### 5.7 `fragdock::EnergyCalculator`

核心接口：

| 接口 | 作用 |
| --- | --- |
| `EnergyCalculator(rad_scale)` | 预计算 atom pair 能量表 |
| `getEnergy(atom, atom)` | 原子对相互作用能 |
| `getEnergy(atom, molecule)` | 单原子对受体能量 |
| `getEnergy(molecule, molecule)` | 分子-分子能量 |
| `calcIntraEnergy(molecule)` | 配体内能 |
| `gauss1/gauss2/repulsion/hydrophobic/hydrogenBond` | 具体项 |
| `getEnergy_strict()` | 直接公式计算，不走查表 |
| `getIntraEnergy_strict()` | 直接公式计算配体内能 |

打分公式保持当前源码常数：

- `score = -0.035579 * gauss1 - 0.005156 * gauss2 + 0.840245 * repulsion - 0.035069 * hydrophobic - 0.587439 * hydrogenBond`
- `final_score = score / (1 + 0.05846 * nrots)`

### 5.8 `fragdock::InterEnergyGrid`

核心接口：

| 接口 | 作用 |
| --- | --- |
| `setInterEnergy()` | 写网格值 |
| `addEnergy()` | 累加网格值 |
| `getInterEnergy()` | 读网格值 |
| `convertX/Y/Z()` | 坐标到索引 |
| `convert()` | 索引到坐标 |
| `parseGrid()` | 读 `.grid` |
| `parseDx()` | 读 `.dx` |
| `writeFile()` | 写 `.grid` |

约定：

- 网格越界时返回 `LIMIT_ENERGY`。
- 网格存储顺序为 x-y-z 线性展开。

### 5.9 `fragdock::AtomInterEnergyGrid`

核心接口：

| 接口 | 作用 |
| --- | --- |
| `AtomInterEnergyGrid(center, pitch, num, xs_type)` | 单原子类型网格容器 |
| `readAtomGrids(folder)` | 读取所有 `.grid` |
| `readDxAtomGrids(folder)` | 读取所有 `.dx` |
| `makeAtomGrids(...)` | 直接构建网格 |

### 5.10 `fragdock::FragmentInterEnergyGrid`

核心接口：

| 接口 | 作用 |
| --- | --- |
| `FragmentInterEnergyGrid(orig_frag, rot_angles, atom_grids, distance_grid)` | 构建 fragment 能量网格 |
| `getGrid()` | 获取内部网格 |

构建逻辑：

- 对每个 fragment 旋转姿态做能量累加。
- 使用 `distance_grid` 做碰撞与远离剪枝。
- 只保留所有姿态中的最小能量。

### 5.11 `fragdock::FragmentInterEnergyGridContainer`

核心接口：

| 接口 | 作用 |
| --- | --- |
| `insert(grid)` | 插入 fragment grid |
| `isRegistered(fragid)` | 判断缓存命中 |
| `get(fragid)` | 读取 grid |
| `next()` | 推进调度步 |

策略：

- `ONLINE`：LRU 替换。
- `OFFLINE`：按预先求出的顺序存放。
- `NONE`：只保留一个 grid。

### 5.12 `fragdock::DecomposeMolecule`

输入：`Molecule`。  
输出：`std::vector<Fragment>`。  
行为：

- 先根据 rotor 键和 ring 结构划分连通块。
- 处理 solitary atom 合并。
- 给 fragment 边界处加入 dummy 原子。
- 为每个 fragment 保留内部键关系。

### 5.13 `fragdock::Optimizer` / `Optimizer_Grid`

核心接口：

| 接口 | 作用 |
| --- | --- |
| `optimize(mol, ec)` | 基于显式能量函数局部优化 |
| `optimize(mol)` | 基于 atom grid 的局部优化 |
| `calcTotalEnergy(mol)` | `inter + intra` |

局部优化行为：

- 固定随机种子 `srand(0)`。
- 每轮随机采样 200 个邻域候选。
- 平移步长约 `0.5`。
- 旋转步长约 `pi/30`。
- RMSD 超过 `local_max_rmsd` 的候选直接跳过。
- 只接受更优解，直到收敛。

### 5.14 `OpenBabel::calc_minRMSD`

输入：

- 一个输出分子
- 一组参考分子

输出：最小 RMSD。  
行为：

- 去氢。
- 统一芳香性和 ring 标志。
- 以骨架同构映射计算最小 RMSD。

用于：

- docking 输出姿态去重
- 保证最终输出构象互相有足够差异

## 6. 端到端流程

### 6.1 标准流程

1. 先运行 `atomgrid-gen` 生成 receptor atom grid。
2. 再运行 `conformer-docking` 进行 fragment 搜索、评分、局部优化和输出。
3. 如果需要分析输入分子的 fragment 组成，再运行 `decompose`。

### 6.2 `atomgrid-gen` 的输入输出

输入：

- 配置文件
- receptor

输出：

- `GRID_FOLDER/<xs_type>.grid`
- 日志文件

### 6.3 `conformer-docking` 的输入输出

输入：

- 配置文件
- receptor
- ligands
- `GRID_FOLDER` 中的 atom grids

输出：

- `OUTPUT`
- `OUTPUT + fraggrid__<timestamp>.csv`
- 日志文件

### 6.4 `decompose` 的输入输出

输入：

- ligand 集合
- fragment 文件

输出：

- annotated ligand 文件
- fragments 文件
- 日志文件

## 7. 复现实现时必须保留的行为

| 约束 | 原因 |
| --- | --- |
| OpenBabel 2.4.x 兼容 | 解析、canonical smiles、RMSD、grid 读取都依赖旧 API |
| 原子类型枚举保持一致 | 能量函数和网格命名依赖它 |
| `.grid` 二进制格式保持一致 | 方便直接复用已有网格文件 |
| `SEARCH_PITCH / SCORING_PITCH` 整数比规则不变 | 影响搜索索引映射 |
| fragment reuse 的三种策略不变 | 影响主性能特征 |
| `srand(0)` 和局部搜索参数不变 | 保证结果可复现 |
| RMSD 标准化逻辑不变 | 影响输出去重 |
| affinity 公式不变 | 影响结果对齐 |

## 8. 建议的重构拆分顺序

1. 配置解析与文件 IO。
2. 原子类型、`Atom`、`Molecule`、`Fragment`。
3. `EnergyCalculator` 与 `.grid` 格式。
4. OpenBabel 适配层。
5. `MoleculeToFragments` 分解器。
6. `FragmentInterEnergyGrid` 与缓存容器。
7. 主 docking 搜索流程。
8. 局部优化与 RMSD 去重。
9. `decompose`、`score-only`、`intraenergy-only`、`atom-docking` 等辅助可执行文件。

## 9. 最小可验证清单

- 能用相同依赖编译出 `atomgrid-gen` 和 `conformer-docking`。
- 能用 `testdata/testgrid.in` 生成 `.grid`。
- 能用生成的 `.grid` 对 `testdata/G39.mol2` 执行 docking。
- 能运行 `unittest` 并通过 `EnergyCalculator` 的基本测试。
- 能运行 `decompose` 生成 fragment 输出。
