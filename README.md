# HeraGrid 计算程序

从 `../grid` 提取的正式计算代码，用于生成 HERA DIS 的 PDF 无关系数权重，再结合 PDF 和 alpha_s 快速计算截面。未复制测试代码、smoke 配置、历史 `.hgrid`、日志、缓存或拟合工程接入补丁；CMake 不包含测试目标。原 `grid` 目录未修改。

## 文件和依赖

- `src/main.cpp`：命令行入口、LHAPDF 回调、生成和精度校验流程。
- `src/config.cpp`：读取 YAML、数据集和 cuts。
- `src/builder.cpp`、`include/heragrid/builder.hpp`：直接积分、系数权重生成、并行和断点恢复。
- `src/grid.cpp`、`include/heragrid/grid.hpp`：grid 读写、截面定义和快速卷积。
- `configs/hera.yaml`：完整 HERA I+II 与 combined charm 的生产配置。
- `CMakeLists.txt`、`agent-env.sh`：构建和项目临时目录设置。

需要 CMake >= 3.19、C++17 编译器、Boost 头文件、yaml-cpp、LHAPDF，以及验证/计算所用的 LHAPDF PDF 集。独立程序不需要 ROOT、APPLgrid 或 Fortran。

**此目录仍依赖外部 CTEQDisPac 源码和实验数据，不是把整个物理后端和数据打包进来了。** 本机默认来源为 `/Data4/pyridineking/GlobalFitting/newglobalfitting`。迁移时分别修改 CMake 的 `NEWGLOBALFITTING_DIR` 和 YAML 的 `DataPath`、`Datasets`。

## 编译与运行

```bash
cd /Data4/pyridineking/Tension/heragrid
source ./agent-env.sh
cmake -S . -B .tmp/build -DCMAKE_BUILD_TYPE=Release
cmake --build .tmp/build -j4
.tmp/build/hera-grid --help
```

指定其他后端源码位置：

```bash
cmake -S . -B .tmp/build -DCMAKE_BUILD_TYPE=Release \
  -DNEWGLOBALFITTING_DIR=/path/to/newglobalfitting
cmake --build .tmp/build -j4
```

所指定目录下必须有 `CTEQDisPac/` 及 CMake 列出的 8 个后端源文件/头文件。构建只修改 `.tmp/build/cteq/` 中的副本：增加积分精度设置接口，并修正原基础电弱通道循环的数组越界，不修改外部源码。`backend_sources.sha256` 记录后端源文件指纹。

每次开启新 shell 都先 `source ./agent-env.sh`。临时脚本、缓存、构建和校验输出放在本目录 `.tmp/`，该目录已加入忽略规则。

### 1. 生成权重

```bash
.tmp/build/hera-grid generate configs/hera.yaml
```

生产配置输出：

```text
grid/hera_NLO.hgrid
grid/hera_NNLO.hgrid
.tmp/production/hera.work/       # 每个点的断点、校验 CSV 和摘要
```

`grid/` 输出目录由程序创建。程序先生成权重，再与直接积分比较；每个阶数只有校验通过才保存最终文件。这里的校验是生成程序的正式功能，不是单元测试代码。

### 2. 查看和计算

```bash
.tmp/build/hera-grid inspect grid/hera_NNLO.hgrid
.tmp/build/hera-grid evaluate grid/hera_NNLO.hgrid \
  --pdf CT18NNLO --member 0 --output .tmp/predictions_nnlo.csv
.tmp/build/hera-grid evaluate grid/hera_NLO.hgrid \
  --pdf CT18NLO --member 0 --output .tmp/predictions_nlo.csv
```

`--pdf` 必填，指定已安装的 LHAPDF 集；`--member` 默认 0，指定该集合的成员编号。PDF 和 running alpha_s 均取自该成员。换 PDF/member 不用重新生成权重。PDF 的 x/Q 范围必须覆盖所有卷积节点，程序不自动外推。

`evaluate` 默认使用文件保存的阶数；可在 NNLO 文件上加 `--order 1` 只取 NLO，不能用 NLO 文件计算 NNLO。CSV 包含点编号、数据集、原数据行号（从 0 开始，不含表头）、过程、x、Q2、y 和 prediction。默认 prediction 为约化截面。省略 `--output` 时写到 `<grid路径>.predictions.csv`；建议显式指定 `.tmp/` 下的路径。CSV 会覆盖同名文件。

已有原目录结果时也可以直接读，无须复制或重算：

```bash
.tmp/build/hera-grid evaluate ../grid/grid/hera_NNLO.hgrid \
  --pdf CT18NNLO --output .tmp/existing_predictions.csv
```

### 3. 更换 PDF 后检查精度

```bash
.tmp/build/hera-grid validate grid/hera_NNLO.hgrid \
  --pdf CT18NNLO --member 0 --rtol 0.001 --atol 1e-8 \
  --output .tmp/validation.csv
```

通过判据是 `abs(grid - direct) <= atol + rtol * abs(direct)`。另输出 `.tmp/validation.csv.summary.txt`。`--rtol` 默认 0.001，`--atol` 默认 1e-8；`--reference-rtol` 可指定直接积分的相对容差（大于 0、不超过 0.01）。默认参考积分使用生成容差与 `1e-5` 中较小者，递归深度至少 10。`validate` 使用文件中的阶数及所有点，并要求后端源指纹一致。省略 `--output` 时写到 `<grid路径>.check.csv`。

退出码：0 成功；1 配置、文件或运行错误；2 精度校验未通过。

## config 怎么设置

直接编辑 `configs/hera.yaml`。**YAML 中所有相对路径相对于该 YAML 所在目录；命令行路径相对于当前工作目录。** 字段名区分大小写；未知顶层字段会报错。以下“默认”指省略字段时的程序默认值，不一定是生产配置的取值。

### 基本配置

| 字段 | 默认/必填 | 含义和设置方法 |
| --- | --- | --- |
| `Backend` | `CTEQDisPac` | 目前只支持此后端。 |
| `Order` / `Orders` | 必填其一 | 单阶写 `Order: NLO` 或 `Order: NNLO`；两阶写 `Orders: [NLO, NNLO]`。也接受 1/2。不能同时写两个字段。两阶一起生成时只构建最高阶权重，再截取低阶。 |
| `Output` | 必填 | 输出文件前缀，不是完整文件名。例如 `../grid/hera` 自动加 `_NLO.hgrid` / `_NNLO.hgrid`，即输出到本项目 `grid/`。 |
| `WorkDirectory` | `Output` 加 `.work` | 断点和校验报告目录。建议保持 `../.tmp/production/hera.work`。 |
| `XNodes` | 64 | 每个实验点的卷积 x 节点数，范围 8–512。生产配置为 80；增大通常改善插值精度，但增加生成时间、文件大小和卷积开销。 |
| `Threads` | 1 | 生成阶段的工作线程数，范围 1–256；生产配置为 8。按可用 CPU 和内存设置，不代表校验和 evaluate 也使用这么多线程。 |
| `AllowNLOChargedCurrent` | false | NNLO 配置含 CC 点时必须显式设 true，接受后端 CC 仍为 NLO 的限制；这不会增加 CC NNLO 项。 |

### 质量、电弱参数和数值精度

| 字段 | 默认 | 含义和约束 |
| --- | --- | --- |
| `Masses.Charm` | 1.3 | charm 质量，GeV，必须大于源程序 strange 质量 0.2。 |
| `Masses.Bottom` | 4.75 | bottom 质量，GeV，必须大于 Charm。 |
| `Masses.Top` | 1e99 | top 质量，GeV，必须有限且大于 Bottom；极大默认值用于使 top 在相关能区退耦，不是物理 top 质量。 |
| `Electroweak.MW` | 80.835 | W 质量，GeV；正有限数。仅用于本程序非约化 CC 截面的归一化，不会统一替换后端内部所有电弱参数。 |
| `Electroweak.GF` | 1.16638e-5 | Fermi 常数，GeV^-2；正有限数，同样用于非约化 CC 归一化。 |
| `Integration.RelativeTolerance` | 1e-4 | 生成权重时后端积分的相对容差，要求 `0 < 值 <= 0.01`；越小通常越慢。与最终插值校验容差是两个不同概念。 |
| `Integration.MaxDepth` | 8 | 自适应积分最大递归深度，0–25；增大允许更深的积分细分。 |
| `Envelope.SmallXPower` | -0.3 | 对 xf 插值使用的展平函数 `x^a * (1-x)^b` 中的 a；有限且绝对值不超过 2。 |
| `Envelope.LargeXPower` | 1.0 | 上述 b，范围 0–5；用于改善大 x 插值。这个函数不固定最终 PDF。 |
| `Validation.PDF` | CT18NNLO | 生成后独立校验用的 LHAPDF 集，必须已安装。它不参与生成系数权重；两个输出阶数都用这个集校验。 |
| `Validation.Member` | 0 | 校验 PDF 成员编号，非负且须在该 PDF 集范围内。 |
| `Validation.RelativeTolerance` | 0.001 | 最终 grid/direct 比较的相对容差，正有限数。 |
| `Validation.AbsoluteTolerance` | 1e-8 | 最终比较的绝对容差，非负有限数，避免直接结果接近 0 时仅相对误差失效。 |

正式预测应选择与计算阶数一致的 PDF 演化和 alpha_s。提高精度时先看失败点的 CSV，尝试将 `XNodes` 从 80 增到 120；必要时减小积分容差或增加深度，并使用新的输出前缀和断点目录。

### 输入已有数据集

```yaml
DataPath: /Data4/pyridineking/GlobalFitting/newglobalfitting/share/datasets
Datasets:
  - /Data4/pyridineking/GlobalFitting/newglobalfitting/share/DefaultConfigs/HeraIplusII.yaml
  - /Data4/pyridineking/GlobalFitting/newglobalfitting/share/DefaultConfigs/HeraCombinedCharm.yaml
```

`Datasets` 是数据集描述 YAML 的列表，不是 `.dat` 文件列表；配置了它就必须有 `DataPath`。数据文件实际位置为 `DataPath / SubDir / DataFiles.Filename`。

数据集描述中的关键字段：

| 字段 | 含义 |
| --- | --- |
| `Name` | 保存到 grid 的数据集标识。 |
| `Class` | 必须是 `HeraData` 或 `HeraDataNLO`；生成阶数仍由顶层 `Order(s)` 决定。 |
| `Theory` | `Internal` 为 inclusive，`Internal_c` 为 charm，`Internal_b` 为 bottom。 |
| `SubDir` | 数据文件在 DataPath 下的子目录。 |
| `DataFiles[].Filename` | 数据文件名。 |
| `DataFiles[].Metadata.Type` | `NCep`、`NCem`、`CCep`、`CCem`。NC/CC 分别为中性流/带电流，ep/em 分别为正电子/电子与质子散射。 |
| `DataFiles[].Metadata.reduced` | 整数；0 表示非约化，非零表示约化。NC 只支持约化截面。 |
| `Cuts[]` | 用 `VarName`、`Min`、`Max` 筛选数据列，包含上下界；不支持非空 `Formula`，应提前筛选。 |

`.dat` 采用原工程的 `NData`、`NColumn`、`ColumnNames` 表头，必须包含 `x`、`Q2`、`y` 列。程序读取运动学和 cuts，不把实验中心值/误差写入理论权重。最终点数取决于数据和 cuts。

### 直接指定计算点

也可以用 `Points` 替代 `Datasets`，此时无需 `DataPath`；若两者都给出，先读数据集，再追加 Points，不自动去重。例如另存为 `configs/custom.yaml`：

```yaml
Backend: CTEQDisPac
Order: NLO
Output: ../grid/custom
WorkDirectory: ../.tmp/custom.work
XNodes: 80
Threads: 4
Masses: {Charm: 1.3, Bottom: 4.75, Top: 1.0e99}
Integration: {RelativeTolerance: 1.0e-4, MaxDepth: 8}
Points:
  - {Dataset: my_nc, Process: NCep, Flavor: 0, Reduced: true, x: 0.001, Q2: 12, y: 0.118}
Validation:
  PDF: CT18NLO
  Member: 0
  RelativeTolerance: 0.001
  AbsoluteTolerance: 1.0e-8
```

运行 `.tmp/build/hera-grid generate configs/custom.yaml`。

| Points 内字段 | 含义 |
| --- | --- |
| `Dataset` | 点所属的数据集名，默认 `inline`。 |
| `Process` | 必填，支持上表四种过程。 |
| `Flavor` | 默认 0：inclusive；14：charm；15：bottom。重味标签仅支持 NC。 |
| `Reduced` | 默认 true，计算约化截面；false 只支持 CC，沿用后端对应的 pb/GeV² 归一化。 |
| `x` | 必填，Bjorken x，要求 `0 < x < 1`。 |
| `Q2` | 必填，Q²，单位 GeV²，正有限数；不是 Q。 |
| `y` | 必填，非弹性度，范围 [0, 1]。需自行保证 x、Q²、y 与碰撞能量相容，程序不会根据能量补算 y。 |

至少需要一个通过 cuts 的数据点或一个显式 Point。

## 恢复、输出与物理限制

- 中断后以相同配置执行相同命令，可复用 `WorkDirectory` 中已完成的逐点断点。断点会核对后端指纹、运动学、质量及数值设置。改变这些设置时使用新的 WorkDirectory 和 Output。
- 最终 `.hgrid` 不覆盖已有文件。若两阶中已有一阶保存，重跑时应更换 Output；不要让两个生成进程共用 Output 或 WorkDirectory。
- `generate CONFIG --limit N` 只取前 N 个点，会产生不完整 grid；应使用独立输出/断点路径，不可将它绑定到完整数据集。
- 文件是本程序的 `.hgrid` 格式，不是 APPLgrid ROOT 格式。`NNLO` 沿用原后端物理内容：NC 的 NNLO 增量来自现有 photon 部分，CC 仍为 NLO；每个点记录实际阶数。
- 当前固定 `muF = muR = Q`。更改质量、运动学、截面定义、尺度处方或后端物理实现需要重新生成。更换 PDF/member/alpha_s 不需重建，但应重新检查插值精度。
- 核心卷积为 `sum W[p,a,j] * alpha_s(Q)^p * [x_j f_a(x_j,Q)]`。权重已包含相应系数函数常数，乘的是 alpha_s 本身，不要额外除以 `2*pi`。

在自己的 C++ 程序中只读取和卷积时，使用 `include/heragrid/grid.hpp` 和 `src/grid.cpp`，或链接 `.tmp/build/libhera_grid.a`；不必链接生成器/CTEQDisPac/yaml-cpp。`Grid::load()` 读文件，`Grid::convolve(xf, alpha)` 计算所有点。`xf(x,Q)` 返回 13 个 `x*f` 值，槽位 0–12 对应 -6–6，槽位 6 为胶子（LHAPDF PDG 21）；`alpha(Q)` 返回 running alpha_s，输入 Q 的单位为 GeV。具体 LHAPDF 回调参考 `src/main.cpp`。
