# REstretto `--full-rotation` 6D 片段格点执行流程

本文档用于在**不影响当前默认行为**的前提下，为当前 `restretto` 项目增加一个新的对接选项 `--full-rotation`。

目标不是重写现有流程，而是在保留原始 3D 逻辑的同时，增加一条可选的 6D 片段格点对接路径：

- 默认情况下仍走原来的 3D 平移 + 现有打分逻辑
- 开启 `--full-rotation` 时，切换到 3D 平移 + 3D 旋转的离散 6D 搜索
- 旋转角采样参考 `references/restretto` 的旋转逻辑
- 不做旋转插值，直接取最近邻角度对应的能量/分数
- 当多个角度距离相同时，优先选择打分更高/能量更低的角度

## 1. 当前问题确认

当前项目的默认实现存在以下事实：

1. `InterEnergyGrid` 只表示 `x/y/z` 三维平移网格。
2. `conformer-docking` 入口只对输入构象做现有打分，没有 6D 片段网格搜索。
3. `decompose` 只把分子切成 fragment，没有构建“旋转维度上的片段格点”。
4. 当前片段相关逻辑本质上是“fragment 分解 + 3D 平移搜索/打分”，不是完整的 6D 搜索。

因此，新增 `--full-rotation` 是合理且独立的扩展，不应修改默认路径。

## 2. 设计原则

### 2.1 默认路径保持不变

未指定 `--full-rotation` 时：

- 继续使用当前实现
- 继续保持当前输出结构、阈值、排序和测试行为
- 不改变已有 CLI 含义

### 2.2 新功能独立成分支

`--full-rotation` 只影响 docking 入口的搜索策略：

- 复用现有的分子解析、fragment 分解、能量函数、输出逻辑
- 新增一个 6D 搜索分支
- 新分支不反向侵入旧分支

### 2.3 旋转采用最简单离散方案

旋转维度先不插值，直接做最近邻匹配：

- 先建立离散旋转采样集合
- 对任意连续旋转输入，映射到最近的旋转 bin
- 若两个候选角度的距离相同：
  - 优先选择更优分数
  - 如果是能量场，则优先选择更低能量

## 3. 新增能力清单

### 3.1 CLI 新选项

在 `conformer-docking` 对应入口增加：

- `--full-rotation`

语义：

- `false` 或未传入：走原逻辑
- `true`：走 6D 片段格点对接逻辑

### 3.2 6D 片段格点

新增一个逻辑上的 6D 搜索空间：

- 3 个平移维度：`x, y, z`
- 3 个旋转维度：`rot_x, rot_y, rot_z` 或等价旋转参数

注意：

- 不要求把 6D 做成连续插值网格
- 不要求改动现有 `.grid` 二进制布局
- 只要逻辑上支持 `(translation, rotation)` 的离散搜索即可

### 3.3 最近邻旋转选择

对旋转角的处理方式：

1. 为 fragment 或 ligand 生成一组离散初始旋转角。
2. 对任意候选姿态，找到最近的角度 bin。
3. 直接使用该 bin 的能量/打分。
4. 不对角度做线性插值、球面插值或三线性旋转插值。

## 4. 建议实现结构

建议在当前项目中增加以下新模块或等价结构：

- `rotation.py`
  - 负责初始旋转采样、最近邻旋转 bin 选择、tie-break 规则
- `fragment_grid_6d.py`
  - 负责 6D 片段格点的索引、查询与缓存
- `docking_full_rotation.py`
  - 负责 `--full-rotation` 的完整对接流程

现有模块保持用途：

- `grid.py`：保留 3D 平移格点
- `energy.py`：保留当前能量函数
- `model.py`：保留 Atom/Molecule/Fragment 的几何基础
- `decompose.py`：保留 fragment 分解
- `cli.py`：增加参数分支，不改默认行为

## 5. 执行流程

以下流程是 `--full-rotation` 模式的建议执行顺序。

### Step 1: 读取配置与输入

1. 解析配置文件。
2. 解析 receptor 与 ligands。
3. 读取现有 3D atom grids。
4. 进行配置合法性校验。
5. 若用户没有传 `--full-rotation`，直接返回当前默认流程。

### Step 2: 进入 full-rotation 分支

如果用户传入 `--full-rotation`：

1. 对每个 ligand 执行 fragment 分解。
2. 保留每个 fragment 的原子坐标、键关系、rotor 关系和中心位置。
3. 建立 fragment 库，按 canonical 标识去重。
4. 为 fragment 预生成离散旋转采样集合。

### Step 3: 初始角度采样

这里参考 `references/restretto` 的做法，但只保留“离散初始采样”的思想。

推荐流程：

1. 先定义一组标准初始旋转角集合。
2. 可以沿用参考项目中的 60 个固定旋转，或者从其逻辑上生成等价集合。
3. 每个 fragment 都以这些初始角作为候选旋转起点。
4. 如果用户提供自定义角度文件，则优先读取自定义集合。

这一步的要求是：

- 采样必须稳定
- 顺序必须可复现
- 对同一输入，候选旋转集合必须一致

### Step 4: 构建 6D 逻辑格点

1. 对每个 fragment 建立 `(x, y, z, rot_idx)` 搜索空间。
2. 平移维度仍沿用当前 `INNERBOX / SEARCH_PITCH / SCORING_PITCH` 的离散化思路。
3. 旋转维度使用离散旋转 bin。
4. 每个 6D 点对应一个明确的姿态和能量值。
5. 若某个旋转 bin 已经可由最近邻角度表示，则只记录这个 bin，不做插值。

### Step 5: 计算 fragment 在各旋转 bin 下的分数

1. 对每个 fragment。
2. 对每个旋转 bin。
3. 将 fragment 旋转到该 bin 对应的姿态。
4. 将 fragment 平移到每个候选 3D 网格点。
5. 计算该姿态在 receptor 上的能量。
6. 记录到 6D 搜索表中。

这里的关键是：

- 旋转维度不插值
- 角度就近取 bin
- 若出现多个候选角度与目标角度距离相同，按更低能量优先

### Step 6: 生成候选构象并直接打分

这一阶段不再强制在 6D 格点上做额外优化，而是参考 `reference/restretto` 的实现思路：

1. 先用旋转采样生成一批离散候选姿态。
2. 将这些姿态映射到 fragment 的当前位置与朝向。
3. 对每个候选 ligand 构象，直接调用现有能量函数或等价的片段能量汇总逻辑进行打分。
4. 不做 6D 网格上的连续优化，也不做旋转插值。
5. 直接按打分结果对候选构象排序。
6. 只保留前 `poses_per_lig_before_opt` 个候选，或者按现有规则保留对应数量的高分姿态。

这一步的重点是：

- `full-rotation` 用于扩展候选构象空间
- 最终打分仍然是“对候选 ligand 构象直接评估”
- 不是必须做“6D 网格上的全局最优搜索”
- 这与参考实现中“生成候选 -> 评分 -> 排序 -> 可选局部优化”的结构保持一致

### Step 7: tie-break 规则

当多个角度距离相同时，必须按以下顺序处理：

1. 先比较能量/分数。
2. 分数更优者优先。
3. 若分数完全相同，再用固定顺序打破平局。

建议固定顺序：

- 先比较 `rotation bin` 索引
- 再比较平移网格索引
- 再比较输入 ligand 的原始顺序

这样可以确保结果确定性。

### Step 8: 局部优化与最终输出

1. 对排序后的最佳候选做局部优化。
2. 如果 `NO_LOCAL_OPT` 开启，则跳过优化。
3. 优化后的姿态回写到原始 ligand。
4. 使用现有输出逻辑写 SDF/CSV。
5. 如果启用了姿态去重，沿用现有 RMSD 规则或其 6D 扩展版本。

## 6. 具体接口建议

### 6.1 `cli.py`

在 `conformer-docking` 的参数中增加：

- `--full-rotation`

建议逻辑：

- 如果 `args.full_rotation` 为 `False`，调用当前实现
- 如果 `args.full_rotation` 为 `True`，调用新 6D 分支

### 6.2 `config.py`

可增加配置字段：

- `full_rotation: bool = False`

以及必要时的旋转采样参数：

- `rotation_bins`
- `rotation_step`
- `rotation_source`

如果你想尽量少改配置，也可以只把 `--full-rotation` 作为 CLI 选项，不落入配置文件。

### 6.3 `rotation.py`

建议暴露的接口：

- `make_initial_rotations(reference_mode=True)`
- `nearest_rotation_bin(rotation, bins)`
- `rotation_distance(a, b)`
- `choose_better_rotation(candidate_a, candidate_b)`

### 6.4 `fragment_grid_6d.py`

建议暴露的接口：

- `build_fragment_6d_grid(fragment, receptor, atom_grids, rotations)`
- `score_pose(fragment, translation_idx, rotation_idx)`
- `query_nearest(fragment, position, rotation)`

## 7. 测试计划

### 7.1 默认模式回归

必须确认：

- 不传 `--full-rotation` 时，现有命令行为不变
- `atomgrid-gen` 不受影响
- `score-only` 不受影响
- `intraenergy-only` 不受影响
- `decompose` 不受影响

### 7.2 6D 模式测试

建议新增测试：

1. `--full-rotation` 参数能被 CLI 正确识别。
2. 最近邻旋转 bin 能返回确定性结果。
3. 同距离 tie-break 会优先选择更低能量/更高分数。
4. 6D 模式下至少可以完成一次完整的 docking 流程。
5. 6D 模式不会破坏默认 3D 模式。

### 7.3 约束测试

建议加以下约束测试：

- 旋转不插值
- 同距离平局按能量优先
- 默认不开启 `--full-rotation`
- 旧输出格式仍可读

## 8. 实施顺序

建议按这个顺序做：

1. 增加 `--full-rotation` 参数，但先不启用新逻辑。
2. 提取旋转采样与 nearest-bin 规则。
3. 实现 6D 片段格点容器。
4. 接入 full-rotation 分支。
5. 加测试。
6. 验证默认路径回归无变化。

## 9. 成功标准

当以下条件全部满足时，认为任务完成：

- 默认对接流程保持原样
- `--full-rotation` 可用
- 6D 片段格点对接可跑通
- 旋转采用最近邻，不插值
- 同距离情况下按更优分数/更低能量优先
- 相关测试通过

## 10. 备注

这个扩展的核心不是“把 3D 网格变成真正连续 6D 张量”，而是：

- 保留现有 3D 平移格点
- 给 fragment 增加离散旋转维度
- 在 6D 采样空间里做直接搜索

这样可以最大限度降低对当前项目的侵入，同时满足你要的“6D 片段格点直接对接和打分”的选项化扩展。
