# REstretto 当前实现流程文档：计算化学视角

本文从计算化学逻辑解释 `references/restretto` 的当前实现。目标是让 Python 重构时能复现原项目的主要科学假设、打分项、搜索空间、片段复用和输出判据。

## 1. 工具定位

REstretto 是面向虚拟筛选的蛋白-配体对接工具。核心思想不是直接对每个配体的每个原子位姿重复计算受体相互作用，而是：

1. 对受体预先生成 X-Score 原子类型能量网格。
2. 将配体按可旋转键和环系统分解为片段。
3. 为唯一片段预计算片段-受体能量网格。
4. 在配体平移/旋转搜索时，通过片段网格平移叠加估计整配体能量。
5. 对候选位姿做局部优化、旋转键惩罚校正、RMSD 去重，输出构象和分数。

## 2. 分子输入与化学标准化

原实现依赖 OpenBabel 2.4.x。重构时如果使用 Python，OpenBabel/RDKit 适配层必须显式复现以下行为：

1. 按文件扩展名读取 receptor 和 ligand。
2. 每个输入分子经过 `fixBondOrders` 修正：
   - carboxyl oxygen 相关键阶设为 `2 + formal_charge`。
   - arginine `Ng+` 氮相关键在 formal charge 为 `+1` 时设为 2，为 0 时设为 1。
   - 修正后清除隐式价态和氢添加感知标记。
3. 读取后添加极性氢 `AddPolarHydrogens()`。
4. ligand 转换为内部分子对象前，docking 主流程会再调用 `AddHydrogens()`，转换后再从 OpenBabel 分子删除氢。
5. 内部 `Molecule` 用 OpenBabel canonical SMILES 和 title 形成 identifier，用于同一配体构象去重聚合。

## 3. X-Score 原子类型

每个原子必须映射到原实现的 21 个网格类型或额外氢类型：

| 类型 | 编号 | 含义 |
| --- | ---: | --- |
| `C_H` | 0 | 疏水碳 |
| `C_P` | 1 | 极性碳 |
| `N_P` | 2 | 非供体/受体氮 |
| `N_D` | 3 | 氢键供体氮 |
| `N_DC` | 4 | 带正电、受体判定氮 |
| `N_A` | 5 | 氢键受体氮 |
| `N_DA` | 6 | 供受体氮 |
| `O_P` | 7 | 非供体/受体氧 |
| `O_D` | 8 | 供体氧 |
| `O_A` | 9 | 中性受体氧 |
| `O_AC` | 10 | 负电受体氧 |
| `O_DA` | 11 | 供受体氧 |
| `S_P` | 12 | 硫 |
| `P_P` | 13 | 磷 |
| `F_H` | 14 | 疏水氟 |
| `Cl_H` | 15 | 疏水氯 |
| `Br_H` | 16 | 疏水溴 |
| `I_H` | 17 | 疏水碘 |
| `Met_D` | 18 | 金属，按供体处理 |
| `Other` | 19 | 其他 |
| `Dummy` | 20 | 片段边界 dummy |
| `H` | 22 | 氢，仅内部判定使用，不生成网格 |

范德华半径必须保持：

```text
C_H 1.9, C_P 1.9, N_* 1.8, O_* 1.7, S_P 2.0, P_P 2.1,
F_H 1.5, Cl_H 1.8, Br_H 2.0, I_H 2.2, Met_D 1.2,
Other 2.0, Dummy 1.5
```

疏水类型是 `C_H/F_H/Cl_H/Br_H/I_H`。受体类型是 `N_A/N_DA/O_A/O_AC/O_DA`。供体类型是 `N_D/N_DC/N_DA/O_D/O_DA/Met_D`。

## 4. 相互作用能量

原实现使用 AutoDock Vina 风格的 5 项打分，距离截断为 8 A，能量上限常量 `LIMIT_ENERGY = 100`。

令：

```text
r = 两原子欧氏距离
d = r - (radius(type1) + radius(type2)) * rad_scale
```

在预计算查表中，`rad_scale` 由构造 `EnergyCalculator(rad_scale)` 传入。`atomgrid-gen` 默认使用配置中的 `rad_scale = 0.95`，主 docking 片段内部和候选优化使用 `EnergyCalculator(1.0)` 或网格值。

单对原子能量：

```text
gauss1      = exp(-(2*d)^2)
gauss2      = exp(-((d - 3.0) * 0.5)^2)
repulsion   = d > 0 ? 0 : d^2
hydrophobic = 0, unless both hydrophobic
              d <= 0.5 -> 1
              0.5 < d < 1.5 -> 1.5 - d
              d >= 1.5 -> 0
hydrogen    = 0, unless donor/acceptor pair
              d <= -0.7 -> 1
              -0.7 < d < 0 -> -1.428571*d
              d >= 0 -> 0

pair_score = -0.035579*gauss1
             -0.005156*gauss2
             +0.840245*repulsion
             -0.035069*hydrophobic
             -0.587439*hydrogen
```

氢和 dummy 不参与 gauss、hydrophobic、hydrogen 项。repulsion 对 dummy 仍保留。受体氢在原子-受体求和中跳过；配体氢在配体-受体求和中跳过。

最终 affinity 校正：

```text
final_score = inter_energy / (1 + 0.05846 * ligand_rotatable_bond_count)
```

其中 `inter_energy = total_energy - ligand_intra_energy`。主流程用总能量排序，但写入排名和 SDF 属性时使用上述 `final_score`。

## 5. 配体内能

配体内能用同一对原子相互作用模型计算，但只统计：

1. 非氢、非 dummy 原子。
2. 图距离至少为 4 的原子对。
3. 若累加能量达到 `LIMIT_ENERGY`，快速返回 `LIMIT_ENERGY`。

图距离来自分子键图的最短路径。该内能在 docking 搜索网格初始化时作为基线加入每个候选位姿。

## 6. 原子能量网格

`atomgrid-gen` 为 21 个 X-Score 网格类型分别生成 `<xs_type>.grid`。

网格定义：

```text
center = BOX_CENTER
pitch  = SCORING_PITCH
num    = ceil(OUTERBOX / 2 / SCORING_PITCH) * 2 + 1
```

每个网格点放置一个该 X-Score 类型的探针原子，对 receptor 所有非氢原子求相互作用能。若能量和达到 `LIMIT_ENERGY`，该点能量截断为 100。

`.grid` 是二进制格式：

1. `Point3d<float> center`
2. `Point3d<float> pitch`
3. `Point3d<int> num`
4. 按 `x -> y -> z` 三重循环写入 `float` 能量值。

坐标到索引采用最近网格点：

```text
idx_x = round((x - center_x) / pitch_x + (num_x - 1) / 2)
```

越界能量返回 `LIMIT_ENERGY`。

## 7. 配体片段分解

`DecomposeMolecule` 的目标是把刚性或近似刚性的部分作为片段，减少重复计算。

分解规则：

1. 初始时，所有不可旋转键连接的原子合并到同一并查集。
2. 检测环；若 `max_ring_size == -1` 或环大小不超过阈值，则环内原子合并。
3. 对非氢孤立小片段做合并：
   - 单原子且相邻重原子数不超过 2 时，倾向并入相邻较大片段。
   - 若两个片段均为小片段，也可以合并。
4. 若 `merge_solitary=true`，继续尝试沿键合并孤立片段；合并必须不生成新环，并且合并后内部可旋转键旋转 1 rad 后 RMSD 小于 `1e-5`。
5. 所有氢最终合并到其相连片段。
6. 跨片段重原子键会在每个片段边界处加入 dummy 原子，dummy 坐标取相邻片段真实原子的坐标，用于定义片段姿态和排斥项。
7. 片段内保留原有键关系。

## 8. 片段规范化与唯一化

每个片段转回 OpenBabel 分子后：

1. 计算 canonical SMILES。
2. 根据 OpenBabel canonical labels 对内部片段原子重编号。
3. 相同 SMILES 的片段共享同一个 `frag_idx`。
4. 第一次遇到的唯一片段会复制一份，执行 `normalize_pose()` 后加入片段库。

`normalize_pose()` 要求片段中心接近原点，并用 `settri()` 选择的 1 到 3 个参考原子定义标准姿态。片段中心、半径和标准姿态必须和原逻辑一致，否则片段网格无法复用。

## 9. 旋转采样

默认旋转集合为 `makeRotations60()`，返回 60 个 ZXZ 欧拉角。算法基于二十面体方向构造：

1. 由黄金比例构造两个 pole。
2. 按固定 `order = {0,0,1,1,1,2,1,2,2,2,3,-1}` 迭代。
3. 每个 order 生成 5 个旋转，共 60 个。

也可用 `ROTANGS` 文件读取自定义旋转，每行是 `theta, phi, psi`。

## 10. 片段能量网格

片段网格把一个规范化片段放到每个受体网格点，并在 60 个片段旋转中取最小相互作用能。

构建流程：

1. 输入规范化片段、旋转集合、21 个 atom grids、distance grid。
2. 初始片段网格每点设为 `LIMIT_ENERGY`。
3. 单原子片段只计算 1 个旋转；多原子片段计算 60 个旋转。
4. distance grid 存储每个网格点到最近非氢受体原子的距离，用于剪枝：
   - 距离小于 2 A 视为碰撞，跳过。
   - 距离大于 `fragment_radius + 6 A` 视为过远，跳过。
   - 注意原代码仅对部分旋转分支执行该剪枝条件，重构若追求严格一致要保留这个行为。
5. 对片段中每个非氢原子，将片段原子相对坐标加到当前网格点坐标，查询对应 X-Score atom grid。
6. 原子能量累加，达到 `LIMIT_ENERGY` 后停止。
7. 若当前旋转能量更低，更新该网格点。

## 11. 整配体搜索

搜索空间：

```text
search_pitch = SEARCH_PITCH
search_num   = ceil(INNERBOX / 2 / SEARCH_PITCH) * 2 + 1
score_num    = atom_grid.num
ratio        = round(SEARCH_PITCH / SCORING_PITCH)
```

配置校验要求 `SEARCH_PITCH / SCORING_PITCH` 在三个方向上都是整数比。search grid 和 score grid 使用同一 center。

对每个 ligand：

1. ligand 内部坐标先平移到质心/几何中心为原点。
2. 分解出的片段记录其相对 ligand center 的向量。
3. 对每个 ligand 旋转，把每个片段相对向量旋转后除以 scoring pitch 并四舍五入，得到片段网格偏移。
4. 对每个 ligand 旋转和每个 search grid 点，初始化 score grid 为 ligand intra energy。
5. 对 ligand 每个片段：
   - 按 `frag_idx` 从缓存获取或生成片段网格。
   - 将片段网格按片段相对偏移叠加到 ligand score grid。
6. 收集所有 `score < OUTPUT_SCORE_THRESHOLD` 的候选位姿，按每个 ligand identifier 保留 `POSES_PER_LIG_BEFORE_OPT` 个最低分候选。

## 12. 片段网格复用策略

`REUSE_FRAG_GRID` 有三种策略：

1. `OFFLINE`：默认。先根据所有 ligand 的片段出现顺序构图，用 `CalcMCFP` 求解给定缓存容量下每一步应保存哪个 grid。
2. `ONLINE`：LRU。缓存未命中时替换最久未使用的片段网格。
3. `NONE`：缓存大小强制为 1，相当于每次只保留当前片段。

若 `REORDER_LIGANDS=true`，先按片段重要性重排：

1. 片段重要性累加 `fragment.size()`。
2. ligand 内部片段按重要性排序。
3. ligand 按其片段序列排序，增加复用机会。

## 13. 局部优化

`Optimizer_Grid` 是主 docking 使用的优化器。

1. 固定随机种子 `srand(0)`，保证可重复。
2. 初值是候选位姿的 `atom_grid inter energy + ligand intra energy`。
3. 每轮采样 200 个邻域候选：
   - 平移增量每轴均匀采样 `[-0.5, 0.5]` A。
   - 欧拉角增量每轴均匀采样 `[-pi/30, pi/30]`。
   - 围绕当前分子中心旋转，再平移。
4. 若候选相对初始位姿 RMSD 超过 `LOCAL_MAX_RMSD`，跳过。
5. 只接受能量更低的邻域最优候选；没有改进时停止。

`NO_LOCAL_OPT=true` 时跳过优化，仅重新计算候选总能量并排序。`SCORE_ONLY=true` 时不做搜索，只对输入构象打分。`LOCAL_ONLY=true` 时从输入构象开始只做局部优化。

## 14. 输出与去重

主 docking 输出：

1. `OUTPUT`：SDF/MOL 等 OpenBabel 可写格式，坐标更新为 docked pose。
2. 每个输出分子写入属性 `restretto_score`。
3. `OUTPUT + "fraggrid__<date>.csv"`：每个 ligand identifier 的最佳分。
4. log 文件：默认 `OUTPUT + "fraggrid__<date>.log"`。

每个 ligand identifier 的候选处理：

1. 对候选位姿按优化后的 total energy 升序排序。
2. 计算 best intra、best inter 和 best score。
3. 从低能到高能选择最多 `POSES_PER_LIG` 个输出姿态。
4. 用 OpenBabel 同构映射计算候选与已接受姿态的最小 RMSD。
5. 只有 `min_rmsd > MIN_RMSD` 才接受该姿态。

## 15. 辅助功能的化学语义

### `score-only`

对输入配体构象直接计算受体-配体五项能量，输出 title 和 affinity。它使用严格公式逐项求和，不依赖网格。

### `intraenergy-only`

对输入配体构象计算图距离至少为 4 的配体内能，输出 title 和 intra energy。

### `decompose`

独立输出片段库和带片段信息的 ligand：

1. 输入 ligand 文件。
2. 添加氢并转换为内部分子。
3. 按片段规则分解。
4. 每个片段转回 OpenBabel 分子，可选择 capping、碳 capping 和用 isotope 写入 fragment id。
5. canonical SMILES 去重，唯一片段写入 fragment 文件。
6. 原 ligand 添加 `fragment_info` 属性，值是片段 SMILES 逗号拼接。

## 16. Python 重构必须保留的科学约束

1. 原子类型、半径、供体/受体/疏水判定必须一致。
2. 五项打分公式、权重、截断距离和旋转键惩罚必须一致。
3. 原子网格的中心、点数、pitch、二进制布局和最近邻查询必须一致。
4. 配体内能只统计图距离至少为 4 的非氢/非 dummy 原子对。
5. 片段分解必须保留环、rotor、dummy、孤立片段合并规则。
6. canonical SMILES 和 canonical renumbering 要稳定，否则 fragment reuse 行为会改变。
7. 60 个默认旋转、ZXZ 旋转定义、坐标平移中心必须一致。
8. 候选筛选阈值、优化前候选数、局部优化随机种子、RMSD 去重阈值必须一致。
9. `score-only`、`local-only`、`no-local-opt` 的模式差异必须保留。
