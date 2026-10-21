# SSA-Simulator：Gillespie 随机反应模拟器

基于 Gillespie 随机模拟算法（Stochastic Simulation Algorithm, SSA）的化学反应离散事件模拟器。
适合分子数较少、涨落不可忽略的体系：酶动力学、自由基链反应、快速平衡机理等。

## 特性

- **组合数倾向**：$a = k \times \prod_i C(n_i, |\nu_i|)$，多分子同种反应物严格按组合数计算，不足时倾向自动归零
- **可逆反应拆分为独立正逆反应**，符合化学惯例
- **物理时间推进**：$\Delta t = -\ln(u)/a_0$（指数分布事件间隔）
- **三种终止条件**：最大事件步数 / 最大模拟时间 / 滑动窗口稳态检测（相对均值变化 $\leq \varepsilon$）
- **两种 CSV 表头风格**：`反应物X` 前缀、`stoich_X` 前缀或无前缀直写物种名，均兼容
- **详细轨迹输出**：每步记录时间、触发反应、全部物种数量
- **Origin 风格可视化脚本**：论文级配图（Arial、四边框、刻度朝内、300 dpi）

## 项目结构

```
ssa_project/
├── README.md                 ← 本文件
├── requirements.txt          ← 绘图依赖（模拟器本身零第三方依赖）
├── ssa_simulator.py          ← 模拟器主程序
├── ssa_plot.py               ← Origin 风格可视化
├── examples/
│   ├── michaelis_menten.csv  ← 例1：酶动力学（E+S⇌ES→E+P）
│   ├── no_oxidation.csv      ← 例2：NO 氧化预平衡机理
│   └── methane_chlorination.csv ← 例3：甲烷氯代自由基链（24 基元反应）
└── outputs/                  ← 运行示例生成的轨迹与图（可再生）
```

## 安装

- Python ≥ 3.8
- 仅运行模拟器：**零第三方依赖**（纯标准库）
- 绘图脚本：`pip install numpy pandas matplotlib`

## 快速开始

```bash
# Michaelis-Menten 酶动力学
python ssa_simulator.py --csv examples/michaelis_menten.csv \
    --init "E=30,S=300,ES=0,P=0" --seed 7 \
    --max-steps 500000 --window 2000 --rel-eps 0.001 \
    --out outputs/mm_trajectory.csv

# 绘图（Origin 风格，300 dpi）
python ssa_plot.py outputs/mm_trajectory.csv --species E,S,ES,P \
    --smooth 20 --title "Michaelis-Menten" --out outputs/mm_fig.png
```

## 输入格式

CSV 第一行表头，两种风格等价：

| 风格 | 示例 |
|---|---|
| 无前缀（推荐） | `反应ID,反应名称,NO,N2O2,O2,NO2,反应速率` |
| 带前缀 | `反应ID,反应名称,反应物A,反应物B,反应物C,反应物D,反应速率` |

- **ID / 名称列**：`反应ID`/`id`/`ID`，`反应名称`/`name`（别名可混用）
- **速率列**：`反应速率`/`rate`/`k`
- **化学计量列**：表头即物种名；负数=反应物，正数=产物，0 可省略
- 同种自由基偶联写 `Cl:-2`（倾向自动用 $C(n,2)$），交叉偶联写 `X:-1, Y:-1`
- 可逆反应拆成两行，各自独立速率常数

## 命令行参数

### ssa_simulator.py

| 参数 | 默认 | 说明 |
|---|---|---|
| `--csv` | 必填 | 反应表路径 |
| `--init` | 必填 | 初始分子数，如 `"A=100,B=50,C=0"`，未列出的物种默认为 0 |
| `--max-steps` | 100000 | 最大事件步数 |
| `--max-time` | 1e9 | 最大模拟时间 |
| `--window` | 200 | 稳态检测滑动窗口（事件步） |
| `--rel-eps` | 0.01 | 窗口均值相对变化阈值；设 0 关闭稳态检测 |
| `--seed` | 随机 | 随机种子（复现实验必需） |
| `--out` | trajectory.csv | 轨迹输出路径 |

### ssa_plot.py

| 参数 | 默认 | 说明 |
|---|---|---|
| 位置参数 | — | 轨迹 CSV |
| `--x` | time | 横轴 `time` 或 `step` |
| `--species` | 全部 | 逗号分隔的物种子集 |
| `--ylog` | 关 | 对数纵轴 |
| `--smooth` | 1 | 滑动平均窗口（去涨落看趋势） |
| `--style` | line | `line` 或 `both`（线+符号，配 `--marker-every`） |
| `--hist` | 关 | 附"各反应触发次数"横向条形图 |
| `--hist-top` | 12 | 条形图只显示前 N 个反应，其余并入 others |
| `--figsize` | 120,80 | 图尺寸（毫米） |
| `--dpi` | 300 | 导出分辨率，支持 png/pdf/tiff/svg |
| `--out` | 自动 | 输出文件名 |

## 输出格式

轨迹 CSV 每行一个事件：

```
step,time,reaction,reaction_name,<物种1>,<物种2>,...
1,0.0021,1,E+S->ES,29,299,1,0
...
```

`step` 事件序号，`time` 模拟时间，`reaction` 触发的反应 ID，`reaction_name` 反应名，其余为反应后各物种分子数。

## 示例与预期结果

### 例 1：Michaelis-Menten（examples/michaelis_menten.csv）

$k_1=0.1$，$k_2=k_3=1.0$，E=30, S=300。
预期：SSA 与确定性 ODE 吻合；前 ~50 步出现 ES 爆发相（酶饱和），随后零级平台 $v \approx v_{\max}$。
本仓库 seed=7 运行：1232 步，S 300→0，P 0→300。

### 例 2：NO 氧化预平衡（examples/no_oxidation.csv）

2NO⇌N₂O₂（快，$k_1=10$, $k_2=100$），N₂O₂+O₂→2NO₂（慢，$k_3=0.01$），NO=100, O₂=50。
预期：总包速率 ~ $k_{\mathrm{eff}}\cdot C(\mathrm{NO},2)\cdot\mathrm{O_2}$，$k_{\mathrm{eff}} = k_1 k_3 / k_2$；N₂O₂ 涨落呈亚泊松分布
（$\mathrm{std} \approx 0.77\times\sqrt{\text{均值}}$），与二聚化 CME 精确平稳分布定量吻合——这同时验证了组合数倾向的正确性。

### 例 3：甲烷氯代（examples/methane_chlorination.csv）

24 个基元反应：引发（Cl₂→2Cl·, $k=0.02$）+ 传递 8 步（逐级氯代）+ 终止 15 步（全自由基偶联），CH₄=200, Cl₂=200。
预期（seed=3）：412 步；引发 6 次、传递 400 次 → 平均链长 ≈ 67；
产物 CH₄:68, CH₃Cl:70, CH₂Cl₂:39, CHCl₃:10, CCl₄:3，HCl 205；
C/H/Cl 原子守恒精确成立（建议对自建网络一律做守恒校验）。

## 原理速览

1. 倾向：$a_j = k_j \times \prod_i C(n_i, |\nu_i|)$（$\nu_i<0$ 的反应物；$n_i<|\nu_i|$ 时组合数为 0 → 不反应）
2. 总倾向 $a_0 = \sum_j a_j$；按 $a_j/a_0$ 加权随机选反应
3. 时间推进 $\Delta t = -\ln(u)/a_0$，$u\in(0,1)$ 均匀随机
4. 终止：步数/时间上限、总倾向归零（反应物耗尽）、或滑动窗口均值稳定

## 已知限制

- 倾向选择为线性扫描 O(R)，反应数很多时可换别名法/二叉堆
- 耗尽长尾（低计数时事件间隔极大）是 SSA 固有特性，可用 tau-leaping 加速
- 稳态判据对极缓慢漂移可能误判，重要运行请核对终止原因
- `--ylog` 时若数据含 0 会触发 matplotlib 警告（属正常）

## 许可

MIT。欢迎修改与再分发。
