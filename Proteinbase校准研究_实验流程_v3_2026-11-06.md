# Proteinbase Boltz‑2 ipTM 外部校准与跨靶点异质性

## 短篇研究报告实验流程 v3（主表已构建版）

版本日期：2026-11-06  
最近修订：2026-11-17（RQ4E 与 RQ4-C 已进入正常英文论文 v2；RQ5 随机森林模型复杂度敏感性已完成交付核对）  
版本状态：主分析单元、主 assay、主结局和严格结局已冻结；RQ1／RQ2、RQ3／RQ3‑E、RQ4、探索性 RQ5 及最小描述性 RQ6 均已完成；RQ4E 截距可解释误差诊断与 RQ4-C 共享预算 pilot 已纳入正常英文论文 v2；RQ5-RF 为新增事后探索性补充，正式计算与交付核对已完成，但尚未写入 v2 稿件、补充材料或 PDF。当前接续状态以第 17 节为准。  
原始数据：Proteinbase 2026-11-01  
输入文件：`data/proteinbase_latest.csv`  
文件大小：51,867,337 bytes  
SHA-256：`8eac08bfc931409e0703eda904cdd1066c9d4a2ceb46af43763f94c74dabb9cc`

> v3 为独立新文件，v1、v2 保留为方案演进记录，不作修改。
>
> 本版不再继续讨论选题，而是将已经核实的口径转化为可执行主表和后续分析协议。主表之外的元数据注释可以版本化补充；主 assay、主结局或严格结局如需改变，必须另起版本，不得静默覆盖。

---

## 0. 研究问题与短篇定位

### 0.1 核心问题

在六个具有不同阳性率和设计来源的靶点中：

> **原始 Boltz‑2 ipTM 能否直接解释为实验结合成功概率；其区分能力和校准关系如何随靶点及来源上下文变化；不使用新靶点标签时，统一概率映射能否迁移？**

### 0.2 四个主分析问题

- **RQ1：区分能力。** 各靶点内 ipTM 的分布、ROC-AUC 和 average precision（PR 曲线摘要）如何？
- **RQ2：原始校准与异质性。** 原始 ipTM 的可靠性、校准截距、校准斜率和 Brier 表现如何，靶点间差异多大？
- **RQ3：迁移与适配。** U 层零标签 LOTO 能否迁移；允许基率截距更新（I）及截距／斜率共同适配（F）后，折外表现如何变化？
- **RQ4：实验决策价值。** 固定阈值和固定预算下的命中率、富集和跨靶点资源分配表现如何？

探索性 RQ5／RQ6 分别评价界面指标增量和外部队列内的跨预测器描述性一致性，定义见第 16 节，不改变上述四个主分析问题。

### 0.3 定位边界

本研究是回顾性的外部校准与内部—外部混合验证，不是：

- Overath 完整 AF3 ipSAE＋结构过滤 pipeline 的验证；
- Boltz‑1 → Boltz‑2 的配对模型漂移研究；
- 靶点蛋白特征的因果机制研究；
- 实验时间上晚于 Overath 的前瞻性队列研究；
- 多指标性能排行榜；
- 通用 minibinder 成功概率或通用阈值的建立。

当前主队列包含全部设计类别，因此标题和摘要不得将 2,084 个配对整体称为 de novo minibinders。Miniprotein-only 和可确认 de novo 来源是范围敏感性分析。

---

## 1. 已锁定的核心决定

| 项目 | v3 冻结规则 |
|---|---|
| 主队列 | 六个靶点，全部设计类别，N=2,084 |
| 分析单元 | 设计 × 靶点，每个配对只保留一个主 assay |
| AAH-I | SPR 为主 assay；BLI 仅作确认性结果，不进入主模型 |
| 其他靶点 | 使用唯一观测 assay；出现多 assay 时脚本直接失败，不临时选规则 |
| 主结局 | 主 assay 内精确去重后的多数 binding，包含 Weak |
| 严格结局 | 主结局阳性的子集；最高 strength 为 Medium／Strong 才记 1 |
| strength 聚合 | 主 assay 内精确去重后取最高观察强度 |
| AAH-I 敏感性 | 排除 AAH-I，五靶点 N=1,951、主阳性 729 |
| 主模型 | 靶点特异截距和斜率，部分汇聚；不加设计层随机截距 |
| 相关性处理 | 序列簇划分＋campaign／source 分块敏感性 |
| I／F 分解 | 描述性保留；预测性改善声明必须有折外评价 |
| assay 结论 | 不作 assay 稳健性声明，只报告构成、选择和混杂 |
| Bootstrap | 至少 2,000 次有效重采样；报告重采样层级和 Monte Carlo 误差 |
| 汇总权重 | 配对加权和靶点等权分开报告 |
| 重叠报告 | 全库精确重叠 345（8.1%）；ipTM 主队列精确重叠 0；保留镜像差异说明 |
| 序列相似性 | 报告最大 identity 分布，不只报告 90% 二值结果 |

---

## 2. 主表和可执行文件

### 2.1 已生成文件

| 文件 | 用途 |
|---|---|
| `analysis/build_pair_master_v3.py` | 从原始 CSV 重建主表；输入哈希或关键质量门不符时失败 |
| `analysis/v3_outputs/pair_master_v3.csv` | 后续所有 RQ1–RQ4 的核心输入，一行一个主配对 |
| `analysis/v3_outputs/assay_unit_audit_v3.csv` | 所有 assay 单元的标签审计，包括 AAH-I 的 BLI 确认记录 |
| `analysis/v3_outputs/reconciliation_v3.json` | 机器可读的总体、分靶点和质量门核对结果 |
| `analysis/v3_outputs/reconciliation_v3.md` | 可直接阅读的口径核对表 |

重建命令：

```bash
python analysis/build_pair_master_v3.py
```

不得手工编辑输出 CSV。输入哈希变化、主配对数变化、AAH-I 规则破坏或关键样本量不一致时，脚本应停止，而不是继续生成“看起来可用”的表。

### 2.2 主表包含的核心字段

| 字段组 | 主要内容 |
|---|---|
| 标识 | `pair_id`、`design_id`、`target`、`source_row` |
| 设计与来源 | `sequence`、`author`、`designMethod`、`design_class` |
| assay | `primary_assay`、选择规则、全部观测 assay、逐 assay 标签审计 JSON |
| 预测指标 | 原始 `boltz2_iptm`、ipTM 记录数 |
| 主结局 | 原始／去重记录数、阳性数、多数标签、标签不一致标记 |
| 标签敏感性 | 全阳、任一阳、仅在标签一致时保留的派生结局 |
| strength | 原始／去重记录数、观察类别、最高强度 |
| 严格结局 | `strict_outcome`、缺失原因、从主阳性降级标记 |
| 冲突审计 | 主结局阴性但出现 Medium／Strong 的标记 |
| AAH-I 确认 | 是否存在 BLI、BLI 去重后结局 |

### 2.3 主 assay 规则

**AAH-I：**

- 133 个设计均有 SPR，因此 SPR 提供完整覆盖；
- 18 个设计另有 BLI，且 18 个 BLI 全为阳性；
- 这 18 个设计的 SPR 也全部为阳性；
- 其余 115 个 SPR-only 设计中只有 1 个阳性；
- 主分析使用 SPR 的 133／19，BLI 只作确认性描述和选择结构说明。

BLI 是不同的实验测量方法，但这 18 条不是额外的独立设计证据，也不能用于估计 BLI 的校准或验证 assay 稳健性。

**其他五个靶点：**当前均为单一 BLI，主 assay 为 BLI。若重建时发现多 assay，脚本必须失败并要求新版本规则。

---

## 3. 结局定义与冻结后的实际数字

### 3.1 主结局

对每一个设计 × 靶点 × 主 assay：

1. 保留原始 binding 记录；
2. 按 `design × target × assay × measuredAt × binding value` 折叠精确重复；
3. 以去重后阳性比例 ≥0.5 定义 `primary_outcome=1`，否则为 0；
4. 不跨 assay 做多数表决；
5. 多数标签是操作性结局，不自动声称每个时间戳都是独立实验重复。

主表同时保留去重前多数、全阳、任一阳和排除不一致标签所需字段。

### 3.2 严格结局

严格结局必须始终是主结局阳性的子集：

| 主结局与 strength 状态 | `strict_outcome` |
|---|---:|
| 主结局阴性 | 0 |
| 主结局阳性，最高 strength=Medium | 1 |
| 主结局阳性，最高 strength=Strong | 1 |
| 主结局阳性，最高 strength=None | 0 |
| 主结局阳性，最高 strength=Weak | 0 |
| 主结局阳性，但主 assay 无 strength | 缺失，不强行记 0 |

主结局阴性但出现 Medium／Strong 时，严格结局仍不得反向写成 1，必须由 `primary_negative_but_medium_or_strong` 单独标记。当前主队列为 0 条此类冲突。

“取最高强度”是本研究冻结的操作性聚合规则，不等于把多个测量解释为独立重复验证。

### 3.3 冻结后的主表核对

| 靶点 | 主 assay | 配对 | 主阳性 | 严格阳性 | Weak | Medium | Strong |
|---|---|---:|---:|---:|---:|---:|---:|
| Nipah glycoprotein G | BLI | 1,028 | 101 | 100 | 1 | 22 | 78 |
| PD-L1 | BLI | 277 | 198 | 159 | 39 | 87 | 72 |
| VEGF110 | BLI | 268 | 158 | 155 | 3 | 55 | 100 |
| SARS-CoV-2 RBD | BLI | 268 | 226 | 204 | 22 | 136 | 68 |
| AAH-I | SPR | 133 | 19 | 19 | 0 | 10 | 9 |
| TREM2 | BLI | 110 | 46 | 41 | 5 | 19 | 22 |
| **合计** | — | **2,084** | **748** | **678** | **70** | **329** | **349** |

补充核对：

- 主表行数：2,084；
- 唯一设计数：2,084；
- 跨两个以上主靶点的设计：0；
- 主结局阳性率：35.89%；
- 严格结局可用：2,084，缺失 0；
- 严格结局阳性率：32.53%；
- 从主阳性降级：70；
- 五靶点敏感性：1,951 配对、729 主阳性、659 严格阳性；
- AAH-I 使用 SPR 后，与旧跨 assay 聚合标签不同的配对：0；
- 主 assay 内精确去重导致主标签改变的配对：0；
- 各靶点内部 ipTM 无重复值；跨靶点合并时有 2 个配对共享一个数值，不构成靶点内 top-k 并列问题。

这些数字是**主表构建结果**，不是正式统计模型结果。

---

## 4. RQ1：靶点内区分能力

### 4.1 主要输出

按靶点分别报告：

- ROC-AUC；
- average precision（AP，PR 曲线摘要），并同时列出阳性率基线；不得与阳性率、其他数值积分口径的 PR-AUC 或跨靶点 pooled AP 混淆；
- 阳性、阴性 ipTM 分布；
- 有效设计数、阳性数和阴性数；
- ipTM 唯一值数和最大并列数；
- 序列簇 bootstrap 区间。

主文以逐靶点森林图和分布图为主，不把 pooled AUC 作为核心结论。

### 4.2 权重与 AAH-I 的特殊处理

- 配对加权汇总回答“当前所有送测配对中的平均表现”；
- 靶点等权汇总回答“六个靶点上的典型表现”；
- 两者不得混为同一个平均数；
- Nipah 的大样本不得在靶点等权结果中占主导；
- AAH-I 单独成节报告，不让其高 AUC 主导总体区分能力叙述；
- 五靶点排除 AAH-I 的结果是区分能力的必要对照。

AAH-I 的 19 个阳性中，18 个同时有 BLI 确认，只有 1 个为 SPR-only 阳性。此前描述性高 AUC 可能受到该结构影响，不能在没有对照的情况下外推为普遍高区分能力。

### 4.3 全局再校准与 AUC

全局单调再校准不改变同一靶点内排名，因此不能用合并折外 AUC 的变化声称靶点内排序能力提高。

---

## 5. RQ2：原始校准与靶点异质性

### 5.1 原始分数评价

将 `s = boltz2_iptm` 直接作为待检验的候选概率，报告：

- Brier score；
- log loss；
- 校准截距和斜率；
- 可靠性曲线；
- ECE 仅作次要描述。

不得先拟合再校准映射，再把拟合后的结果称为原始 ipTM 校准。

### 5.2 校准诊断模型

在各靶点中拟合：

\[
\operatorname{logit}\{P(Y=1)\}=\alpha+\beta x,
\qquad x=\operatorname{logit}\{\operatorname{clip}(s)\}.
\]

- 原始概率理想校准的参照为 \(\alpha=0,\beta=1\)；
- \(\alpha,\beta\) 联合描述整体高估或低估及分数压缩／扩张，不能把联合拟合的 \(\alpha\) 单独解释为固定斜率后的 calibration-in-the-large；
- \(\beta\) 表示分数 logit 与结局 logit 的比例关系；
- 校准斜率不等于 AUC；
- 同一数据上的诊断拟合不是折外预测性能。

### 5.3 可靠性图

采用共同固定区间：[0,0.2)、[0.2,0.4)、[0.4,0.6)、[0.6,0.8)、[0.8,0.9)、[0.9,0.95)、[0.95,1]。

每箱报告配对数、平均 ipTM、阳性数、观测阳性率和区间。稀疏箱只作低精度展示，不独立下结论。

### 5.4 层级模型

主模型为：

\[
\operatorname{logit}\{P(Y_i=1)\}=\alpha_{t[i]}+\beta_{t[i]}x_i.
\]

- 对 \(\alpha_t,\beta_t\) 做部分汇聚；
- 报告总体均值、靶点间异质性、收缩前后估计和不确定性；
- 使用弱信息先验或明确的正则化；
- 随机斜率不可识别时，按预设规则简化为随机截距＋共同斜率；
- 模型简化依据收敛和可识别性，不依据结果是否符合预期。

### 5.5 为什么不加设计层随机截距

当前主队列有 2,084 个唯一设计，没有设计跨主靶点重复；采用 AAH-I SPR 主 assay 后，每个设计只保留一条主记录。

因此，设计随机截距没有先前假设的设计内重复结构可用于处理相关性，不能替代：

- 序列近重复聚类；
- campaign／source 分块；
- 靶点层异质性建模。

若未来扩展队列并引入真正的多靶点重复设计，必须重新评估，而不是沿用这一结论。

### 5.6 Brier 参照

\[
BS=n^{-1}\sum_i(\hat p_i-y_i)^2,
\qquad BSS=1-\frac{BS_{model}}{BS_{reference}}.
\]

靶点内描述使用当前靶点恒定阳性率作为诊断参照；零标签迁移按 RQ3 分开报告训练折参照与测试靶点基率参照。

---

## 6. RQ3：零标签迁移与适配分解

> <span style="color:#1f4e79">2026-11-12 修订：主问题仍是 U 层零标签迁移。I／F 是信息递增的适配诊断，不得取代 U 层主结论；加入 Pinc 后能否改善迁移属于探索性 RQ5，不回写本节的主假设。</span>

### 6.1 分析集、源模型与训练权重

主结局为多数 binding（含 Weak），严格 Medium／Strong 结局为敏感性分析。六靶点分析保留 AAH-I，排除 AAH-I 的五靶点分析为强制对照；AAH-I 的 19 个阳性意味着高方差，而不是“不能拟合”。

令 \(x=\operatorname{logit}(\operatorname{clip}(\mathrm{ipTM}))\)，截断规则与 RQ2 冻结配置一致。每个留一靶点（LOTO）折中，仅在源靶点 \(j\neq t\) 上拟合靶点特异截距和共同斜率：

\[
\operatorname{logit}(p_{ij})=\alpha_j+\beta x_{ij}.
\]

- **靶点等权为主分析**：每个源靶点的总权重相等，靶点内配对等权，避免 Nipah 的大样本量主导源模型。
- **配对数加权为敏感性分析**：每个配对等权，回答当前送测配对加权下的平均表现。
- 源模型不是忽略靶点的单一 pooled logistic，也不把 RQ2 的描述性随机斜率层级模型直接当作已验证预测器。

### 6.2 U：未适配的零标签 LOTO

每次完全留出一个靶点，留出靶点不参与拟合、调参、特征变换、分箱或阈值选择。若源靶点配对与留出靶点属于同一序列簇，先从训练侧剔除并报告数量，禁止同一序列簇跨训练／测试。

主分析将源靶点截距等权平均，把共同斜率和平均截距原样用于新靶点：

\[
\bar\alpha_{-t}=\frac{1}{m_t}\sum_{j\neq t}\hat\alpha_j,
\qquad
\eta^{U}_{-t,i}=\bar\alpha_{-t}+\hat\beta_{-t}x_{ti},
\qquad
\hat p^{U}_{ti}=\operatorname{expit}(\eta^{U}_{-t,i}),
\]

其中 \(m_t\) 为该折源靶点数。配对加权敏感性使用与该权重口径一致的源模型及截距加权聚合，不能与主分析混报。

逐留出靶点报告 Brier、log loss、校准截距／斜率、可靠性曲线和两种 BSS。U 层只检验“不使用新靶点标签时绝对概率能否迁移”，不单独证明靶点内排序有效。

### 6.3 I：仅更新截距的基率适配

保留 U 的线性预测值为 offset，只更新新靶点截距：

\[
\eta^{I}_{t,i}=a_t+\eta^{U}_{-t,i},
\qquad
\hat p^{I}_{ti}=\operatorname{expit}(\eta^{I}_{t,i}).
\]

若目标基率为 \(\pi_t\)，截距必须满足下式，而不能简单写成 \(a_t=\operatorname{logit}(\pi_t)\)：

\[
\frac{1}{n_t}\sum_{i=1}^{n_t}
\operatorname{expit}(a_t+\eta^{U}_{-t,i})=\pi_t.
\]

- 回顾性分析若使用留出靶点全体标签计算的实际基率，必须标为 **`oracle base-rate scenario`**；它只回答“若已知真实基率，仅修正整体水平有多少帮助”，不是可部署的零标签结果。
- 预测性 I 的基率只能来自适配训练折或有记录的外部预实验／独立信息，不能使用测试折实际基率。外部基率须同时报告来源和样本量。
- I 层主要以新靶点恒定基率预测为参照，判断 ipTM 映射在基率之外是否仍有增量信息。

### 6.4 F：截距和斜率共同适配

描述性 F 定义为：

\[
\eta^{F}_{t,i}=a_t+\gamma_t\eta^{U}_{-t,i}.
\]

它允许新靶点的有效截距和斜率共同变化。使用新靶点全样本拟合的两参数结果仅称**描述性适配诊断**，不得称为折外预测改善，也不以 Brier 必然改善为前提。

只有预测性 F 可以支持“适配后预测更好”：

- 在新靶点内部按序列簇做分组 K 折，尽量保持各折结局构成；每折只用适配训练部分估计 \(a_t,\gamma_t\)，在折外评价。
- 使用预设正则化或部分汇聚控制小样本不稳定；截距不惩罚，斜率适配向“不改变 U 斜率”的方向收缩。不得根据全样本结果选择是否正则化。
- **拟冻结执行默认，待方案终审后写入运行配置**：序列簇 5 折、随机种子 `20261112`；正则化强度只能在适配训练折内部通过嵌套验证选择。K 值、正则化形式和候选网格在正式运行前冻结，不能看折外结果后调整。
- 若训练折无法形成两类结局或模型不可识别，按预设规则记为不可评价；信息不足时结论写“不确定”，不静默利用测试折补救。

### 6.5 非线性形态：仅限预设敏感性

RQ3 主分析保留上述线性 logistic 形式。非线性只增加一种低自由度敏感性：对 \(x=\operatorname{logit}(\operatorname{clip}(\mathrm{ipTM}))\) 使用 **3 个结点的限制性立方样条 logistic**。

结点分位数须在运行前冻结，结点位置只由各训练折数据确定，再原样用于对应测试折；不得用留出靶点选择结点或函数形态。该分析用于检验“线性假设是否遗漏明显形态”，不是多算法排行榜。随机森林等灵活 ML 仅属于 RQ5 的可选探索基准，不进入 RQ3 主分析。

### 6.6 两种 BSS 参照与逐层判定

- `BSS-transport`：以训练靶点恒定阳性率预测留出靶点为参照，对应零标签迁移场景；训练基率的汇总权重与源模型口径一致。
- `BSS-target-oracle`：以留出靶点自身恒定阳性率为参照，仅用于已知基率后的诊断，不得称为实际可获得的零标签基线。

两种 BSS 必须与原始 Brier、log loss 和校准曲线同时报告。判定规则如下：

| 层级 | 主要证据 | 判定边界 |
|---|---|---|
| U | 逐折 `BSS-transport`，辅以 `BSS-target-oracle`、Brier、log loss 和曲线 | 平均 `BSS-transport` 不改善，判为当前零标签映射未解决迁移；改善但逐靶点有正有负，或相对靶点基率仍无增量，判为有限／部分迁移并保留靶点依赖 |
| I | 回顾性情形看 `BSS-target-oracle`；预测性情形看折外 Brier 相对训练折基率预测的改善 | 不改善说明基率之外没有稳定增量；回顾性 oracle I 改善不能单独证明部署有效 |
| F | 新靶点内部折外 Brier 相对 U／预测性 I 的改善，log loss 和曲线作支持 | 样本内拟合或校准参数变化不能证明预测性适配成功 |

跨层统一解释：平均不改善＝未解决；平均改善但逐靶点方向不一致＝部分有效、仍有靶点依赖；指标冲突＝结论混合；有效样本或事件不足＝不确定。**不得采用“BSS＞0 且校准斜率置信区间包含 1”作为成功标准，斜率区间包含 1 不是良好校准的充分证据。**

逐 LOTO 折报告结果和折间范围。靶点内序列簇 bootstrap 可描述该折的条件精度，但不得把六个靶点的配对 bootstrap 宣称为对未来靶点普遍不确定性的覆盖。最终实验决策价值由 RQ4 判断，不能由 RQ3 的校准参数单独宣布“问题解决”。

若 U／I／F 不理想，RQ3 如实报告；预设的后续问题是“加入 Pinc 能否产生折外增量”，由 RQ5 在相同信息层级和划分下回答，不能事后把 RQ3 主暴露改为组合指标。

---

## <span style="color:#1f4e79">7. RQ4：实验筛选决策价值（2026-11-14 修订）</span>

> <span style="color:#1f4e79">完整冻结协议见 `analysis/RQ4_decision_value_frozen_protocol_2026-11-14.md`。RQ4-A 是预设主要决策分析；RQ4-B 是描述性 operating-point 分析；RQ4-C 是看到 RQ3／RQ3-E 后新增的探索性回顾性模拟。本轮不执行缺乏实验成本模型支持的 net benefit。</span>

### <span style="color:#1f4e79">7.1 RQ4-A：靶点内固定预算</span>

<span style="color:#1f4e79">选择各靶点原始 ipTM 前 5%、10%、20%，报告实际入选数、hits、precision、recall、enrichment、相对随机的命中增益和单位 hit 成本，并进行 2,000 次序列簇 bootstrap。绝对预算 k=10、20、50 是结果链形成后新增的探索性敏感性。单调再校准不改变同一靶点内 top-k，不声称校准提高本层命中率。</span>

### <span style="color:#1f4e79">7.2 RQ4-B：原始 ipTM 固定 operating points</span>

<span style="color:#1f4e79">保留 0.80、0.90、0.95，逐靶点报告入选数、precision、recall、富集和区间。它们不是成功概率阈值；空选集合必须保留，不搜索或推荐跨靶点通用截点。</span>

### <span style="color:#1f4e79">7.3 RQ4-C：共享预算与本地 pilot 成本</span>

<span style="color:#1f4e79">使用 RQ3 对应留出靶点的折外 U 作为 pilot 适配起点。U 与 ipTM 只在靶点内必然同序；不同靶点的 LOTO 源模型不同，不能拼成一个统一部署映射，因此共享预算基线只展示 raw ipTM。比较 raw 共享排序、靶点等额预算、随机分配、pilot-I 和 pilot-F。预设 B=120、300，m=10、20、40，只执行总 pilot 不超过 B 的组合并完整报告。</span>

<span style="color:#1f4e79">所有策略使用相同总预算；pilot 消耗计入成本，pilot hits 计入总命中，已测配对不得重复选择。pilot 按序列簇抽代表，后续排除同簇；主分析 2,000 次重复，五靶点和严格结局各 1,000 次，固定种子。小 pilot 使用预冻结的 Jeffreys 平滑和固定 ridge，明确区别于 RQ3 的大训练折适配。RQ4-C 只是多靶点竞争共享预算的回顾性模拟，不是实际部署或未来靶点泛化证据。</span>

---

## 8. 序列、来源和不确定性

### 8.1 序列近重复

- 精确序列重叠与相似序列重叠分开报告；
- 主聚类规则：90% identity、双向覆盖度 ≥80%；
- 严格敏感性：80% identity、覆盖度要求不变；
- 同一序列簇不得跨训练／测试折；
- 报告每条 Proteinbase 序列对 Overath 的最大 identity 分布，而非只给一个阈值结果；
- identity 阈值是操作性泄漏控制规则，不是“同一生物设计”的自然定律。

### 8.2 来源块

主文报告每个靶点的作者数、designMethod 数、campaign 或代理定义、最大来源占比和 assay 构成。

必要敏感性分析：

- 排除 AAH-I；
- 将 DeepMind 来源的 PD-L1、VEGF110、SARS-CoV-2 RBD 作为共享来源块排除或留出；
- 在 Nipah、TREM2 内，若 campaign 可识别，进行留一 campaign／提交来源分析；
- 独立块过少时，以影响分析和区间展示为主，不声称完成充分外部泛化。

### 8.3 区间估计

- 单上下文指标优先按序列簇 bootstrap；
- 有足够 campaign 时按 campaign／submission 块重采样；
- 独立块过少时不依赖普通聚类渐近推断；
- 至少 2,000 次有效重采样，固定随机种子；
- 同时报告配对加权和靶点等权结果；
- 层级模型区间必须披露少靶点、少来源块带来的不确定性。

---

## 9. 必做敏感性分析

| 编号 | 分析 | 目的 |
|---|---|---|
| S1 | 主结局 → 严格结局 | 检查 Weak 纳入对结论的影响 |
| S2 | 六靶点 → 排除 AAH-I 的五靶点 | 检查特殊阳性结构和单一来源影响 |
| S3 | <span style="color:#1f4e79">全阳、任一阳、排除标签不一致配对仅作标签审计；不再分别平行建模</span> | <span style="color:#1f4e79">仅 1／2,084 配对不一致，报告第 16.1 节规则表即可；严格结局仍执行 S1</span> |
| S4 | Miniprotein-only | 检查全类别结果能否外推至 minibinder 范围 |
| S5 | 可确认 de novo 来源 | 检查 de novo 范围；无可靠注释时不强行执行 |
| S6 | 排除／整体留出 DeepMind 来源块 | 检查共享设计流程主导 |
| S7 | 90%／80% identity 及最大 identity 分布 | 检查近重复和划分规则 |
| S8 | 留一靶点、留一 campaign／来源块 | 检查单一上下文主导 |
| S9 | 配对加权与靶点等权 | 检查大样本靶点主导 |
| S10 | 标记 ipTM 曾参与生成或筛选的设计 | 检查选择偏倚；元数据不足时列为限制 |

“仅 BLI”不再被称为 assay 稳健性分析，因为它会同时排除 SPR 主 assay 的 AAH-I 并改变靶点组成。

---

## 10. 重叠与时间独立性措辞

当前锁定报告口径：

- 全库与 Overath HF 镜像精确序列可匹配设计：345／4,280（8.1%）；
- 重叠设计—靶点配对：342／4,943（6.9%）；
- 当前 Boltz‑2 ipTM 主队列精确序列重叠：0；
- HF 镜像为 3,669 条，原文为 3,766 个设计，97 条差异尚需归档解释；
- 未证明镜像为严格子集前，不将 8.1%／6.9% 称为全数据集数学下限；
- 完成近重复分析前，只能称“精确序列不重叠”；
- `measuredAt` 不用于声称实验时间独立或前瞻性。

---

## 11. 靶点蛋白特征

主模型不加入蛋白长度、等电点、无序比例、疏水性等靶点层特征，也不进行自动特征筛选。

原因是这些特征在靶点内不变，其效应主要由六个靶点及更少的独立来源块识别，不能将 2,084 个配对误当成 2,084 个独立的靶点层观测。

正文最多保留一句限制说明，详细特征表和任何探索性结果放补充材料。不得据此作蛋白特征因果归因。

---

## 12. 主文结果组织与当前完成状态

> <span style="color:#1f4e79">2026-11-15 补记：RQ1–RQ6 已统一汇入 `analysis/paper_results_RQ1_RQ6_draft_2026-11-15.md`，旧 RQ1–RQ5 稿保留。RQ5 使用 2,074 个共同 complete-case，不能与原始 2,084 行结果作非同样本优劣比较；RQ6 仅为 Anthropic 外部队列内 720 个共同设计的同配对描述性补充。</span>

| 输出 | 内容与当前状态 |
|---|---|
| Fig 1 | 主表构建流程、六靶点来源与 assay 构成、Overath 重叠口径；**待制作** |
| Fig 2 | 各靶点 ipTM 分布、ROC-AUC／average precision；AAH-I 单独标注；**已完成**：`analysis/rq12_outputs/RQ1_target_discrimination.png`、`analysis/rq12_outputs/RQ1_ipTM_distributions.png` |
| Fig 3 | 原始可靠性曲线、校准截距／斜率及层级收缩；**已完成**：`analysis/rq12_outputs/RQ2_raw_reliability.png`、`analysis/rq12_outputs/RQ2_partial_pooling.png` |
| Fig 4 | LOTO 未适配迁移与 U／I／F 分解；<span style="color:#1f4e79">**RQ3 已完成**</span> |
| Fig 5 | <span style="color:#1f4e79">RQ3-E 外部 U、低事件压力测试与探索性 I／F；**已完成**，证据等级分开标注</span> |
| Fig 6 | <span style="color:#1f4e79">RQ4 固定预算、原始 operating points 和共享预算 pilot；**已完成并通过验收**</span> |
| Fig 7 | <span style="color:#1f4e79">RQ5 同样本 Pinc／ipTM 诊断、U 层组合增量与 I／F 适配；**已完成并通过主执行线验收**，四张结果图见 RQ5 正式运行目录</span> |
| Fig 8 | <span style="color:#1f4e79">RQ6 分数偏移、排序一致性、逐靶点 AUC 与原始 Brier；**已完成并通过主执行线验收**，整合稿图 17／18 使用高对比度副本</span> |
| Table 1 | 主表样本量、主／严格阳性、来源、assay 和 ipTM 覆盖；**核心数据已冻结，发表版待汇总** |
| Table 2 | Brier、BSS、校准参数、区分指标及区间；<span style="color:#1f4e79">**RQ1–RQ6 数值已具备，整合草稿已排版**</span> |
| Supplement | 严格结局和五靶点对照**已完成**；<span style="color:#1f4e79">RQ3／RQ3-E／RQ4／RQ5 的序列簇、权重、低事件、预算、pilot 及 complete-case 敏感性已具备；RQ6 两种结局及逐预测器区间已具备</span>；RQ5 保留全部源配对加权和逐靶点不利方向 |

### <span style="color:#1f4e79">12.1 RQ1／RQ2 已验收的核心结果</span>

- <span style="color:#1f4e79">六靶点等权 ROC-AUC 为 **0.675（95% CI 0.642–0.705）**，排除 AAH-I 后为 **0.620（0.582–0.655）**。AAH-I 的 AUC 为 0.950，必须单独解释，不能让其主导总体结论。</span>
- <span style="color:#1f4e79">主结局下，六靶点等权原始 Brier 为 **0.416**、log loss 为 **1.222**、BSS 为 **−1.417**；配对加权分别为 **0.474、1.386、−2.397**。严格结局和五靶点对照均未把原始 ipTM 变成已校准概率。</span>
- <span style="color:#1f4e79">主结局层级模型的靶点间标准差为 **τ(α)=3.021、τ(β)=0.824**；排除 AAH-I 后分别为 **1.820、0.167**。这支持上下文异质性，但不等于已证明新靶点的预测分布或通用映射。</span>

### <span style="color:#1f4e79">12.2 RQ2 的方法口径</span>

<span style="color:#1f4e79">RQ2 使用的是**同一 logistic 校准模型族的两层估计**：①逐靶点未汇聚 logistic；②随机截距＋随机斜率的层级部分汇聚 logistic。两者都作用于 \(x=\operatorname{logit}(\operatorname{clip}(\mathrm{ipTM}))\)，用于描述逐靶点校准与跨靶点收缩，不是多个竞争校准算法的排行榜。</span>

<span style="color:#1f4e79">四个分析集均收敛且信息矩阵正定；每靶点、每结局均有 2,000 次有效序列簇 bootstrap。自适应 Gauss–Hermite 9／11 节点的最大靶点估计差为 \(2.71\times10^{-6}\)。自适应积分是数值实现，bootstrap 是不确定性估计，均不得写成额外的校准算法。</span>

### <span style="color:#1f4e79">12.3 可复现结果索引</span>

- <span style="color:#1f4e79">主报告：`analysis/rq12_outputs/RQ1_RQ2_results_2026-11-09.md`。</span>
- <span style="color:#1f4e79">逐靶点及汇总表：`rq1_target_metrics.csv`、`rq2_target_metrics.csv`、`rq12_summary_metrics.csv`。</span>
- <span style="color:#1f4e79">层级模型表：`rq2_hierarchical_targets.csv`、`rq2_hierarchical_hyperparameters.csv`。</span>
- <span style="color:#1f4e79">运行与验收：`rq12_run_manifest.json`。以上 CSV／JSON 均位于 `analysis/rq12_outputs/`。</span>

主文结论首先基于逐靶点结果和五靶点对照，不使用单一 pooled 指标掩盖来源和基率差异。

---

## 13. 执行顺序

### A. 主表冻结——已完成

- [x] 输入哈希核验；
- [x] 六靶点主队列构建；
- [x] AAH-I SPR 主 assay；
- [x] 主结局和严格结局；
- [x] 全阳、任一阳和排除不一致标签所需字段；
- [x] 主表 2,084／748／678 核对；
- [x] 五靶点 1,951／729／659 核对；
- [x] AAH-I SPR／BLI 结构核对；
- [x] 严格阳性为主阳性子集；
- [x] 主表与旧 ipTM、跨 assay 标签对账。

### B. RQ1 与 RQ2——已完成并通过数值验收

- [x] 固定分析配置、序列簇和随机种子；
- [x] 生成逐靶点分布、ROC-AUC、average precision（PR 曲线摘要）；
- [x] 生成原始可靠性曲线和 Brier；
- [x] 拟合层级校准模型并完成收敛诊断；
- [x] 同步输出五靶点和严格结局对照。

> <span style="color:#1f4e79">2026-11-09 验收记录：主队列 2,084／748／678，五靶点 1,951／729／659；每靶点、每结局均保留 2,000 次有效序列簇 bootstrap。四个层级模型均保留随机截距＋随机斜率并收敛，信息矩阵正定。普通积分未通过稳定性检查，正式结果改用经独立密集网格核验的自适应积分；9／11 节点的最大靶点估计差为 2.71×10⁻⁶。完整数值表、图件、失败版本归档和解释边界见 `analysis/rq12_outputs/RQ1_RQ2_results_2026-11-09.md` 与 `analysis/rq12_outputs/rq12_run_manifest.json`。</span>

> <span style="color:#1f4e79">2026-11-15 当前主线：**RQ1／RQ2（已完成并验收） → RQ3（已完成并验收） → RQ3-E（U 确认性、I／F 事后探索性，已完成并验收） → RQ4（已完成并验收） → RQ5（探索性增量，已完成并验收） → RQ6（最小描述性补充，已完成并验收）**。Pinc 批量计算与结果性分析均归入 RQ5，不占独立 RQ 编号。RQ6 仅提供 Anthropic 外部队列内的同配对跨预测器证据，不承担完整跨模型迁移验证。</span>

### <span style="color:#1f4e79">C. RQ3 与 RQ3-E——已完成并验收</span>

- [x] <span style="color:#1f4e79">终审并冻结 K 值、随机种子、正则化形式／候选网格和样条结点分位数；</span>
- [x] <span style="color:#1f4e79">完成靶点等权 U 层 LOTO、配对加权敏感性及逐折指标；</span>
- [x] <span style="color:#1f4e79">完成 oracle I、描述性 F 和序列簇折外预测性 I／F；</span>
- [x] <span style="color:#1f4e79">完成 RCS、严格结局、五靶点对照及 F 层实现修正后的完整重跑；</span>
- [x] <span style="color:#1f4e79">完成 RQ3-E 外部零标签 U、低事件压力测试与事后探索性外部 I／F。</span>

<span style="color:#1f4e79">正式结果见 `analysis/rq3_outputs/RQ3_results_2026-11-13.md`、`analysis/rq3e_outputs/run_20261114_164950/RQ3E_confirmatory_720_results.md` 和 `analysis/rq3e_if_outputs/run_20261114_174403/RQ3E_IF_results.md`。RQ3-E 的 U 与 I／F 不合并为同一确认性证据等级。</span>

### <span style="color:#1f4e79">D. RQ4 与稳健性——已完成并通过独立验收</span>

- [x] <span style="color:#1f4e79">固定阈值、比例预算和绝对预算；</span>
- [x] <span style="color:#1f4e79">共享预算、pilot 成本、五靶点及严格结局对照；</span>
- [x] <span style="color:#1f4e79">Miniprotein-only；de novo 因缺少可靠统一注释而不强行执行；</span>
- [x] <span style="color:#1f4e79">完成 23 项独立数值验收、三张主图和结论边界审查。</span>

<span style="color:#1f4e79">结果：`analysis/rq4_outputs/run_20261114_231925/RQ4_decision_results.md`；验收：`analysis/audit_outputs/RQ4_acceptance_run_20261114_231925/ACCEPTANCE_RESULT.md`。空选 recall=0 的文字口径澄清单独留痕，不改变任何选择或命中结果。</span>

### <span style="color:#1f4e79">E. RQ5——探索性增量分析已完成并通过主执行线验收</span>

- [x] <span style="color:#1f4e79">完成官方实现核验、批量 Pinc、108 个深度 QC 独立复跑及扩展表验收；</span>
- [x] <span style="color:#1f4e79">在 2,074 个 complete-case 配对上完成 Pinc 与 ipTM 的可比 RQ1／RQ2 基础诊断和 2,000 次配对序列簇 bootstrap；</span>
- [x] <span style="color:#1f4e79">在相同样本、冻结 folds 和外测行上评价 Pinc 单指标及 ipTM＋Pinc 组合的 U／I／F 折外增量，完成五靶点、严格结局和源配对加权敏感性；</span>
- [x] <span style="color:#1f4e79">完成交付方 15 项单元测试、18 项独立验证和主执行线 22 项定向验收；</span>
- [ ] 其他界面指标和灵活 ML 仅作后续探索，逐项固定共同 complete-case 并报告覆盖。

<span style="color:#1f4e79">结果：`analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/HANDOFF_RESULT.md`；验收：`analysis/audit_outputs/RQ5_acceptance_run_20261115_111720/ACCEPTANCE_RESULT.md`。Pinc 未显示稳定区分优势；组合在靶点等权 U 层改善部分平均概率误差，但未同步改善 AUC／AP，Brier／BSS 优势对源权重敏感，F 相对 I 的额外收益很小且不稳定。不得宣称 RQ5 已解决通用迁移。</span>

### <span style="color:#1f4e79">F. RQ6 与 C04——最小描述性分析已完成并通过主执行线验收</span>

- [x] <span style="color:#1f4e79">C04／Anthropic 数据下载、字段／结局／CRO、序列重叠和预测器审计均已完成；RQ3-E 与 RQ6 已分别完成外部迁移评价和队列内跨预测器描述。</span>
- [x] <span style="color:#1f4e79">主 Proteinbase 队列的跨预测器 RQ6 取消：936 条可连接 `esmfold2_iptm` 全为 0，其他分数不是现成的同配对复杂物 ipTM；不从 PAE 重新转换。</span>
- [x] <span style="color:#1f4e79">正式分析使用 Anthropic 8 靶点／720 个共同 `uuid`，10 个预测器共享同一 Adaptyv 结局且全部 complete-case；Weak+ 为 287 阳性，Medium+ 为 254 阳性。配对性仅属于外部队列内部。</span>
- [x] <span style="color:#1f4e79">已按冻结协议完成分数偏移、相关、逐靶点 AUC／AP、原始 Brier／log loss 和每靶点 2,000 次设计层配对 bootstrap；未做新的 U／I／F、固定映射跨预测器迁移或最佳模型搜索。</span>
- [x] <span style="color:#1f4e79">主执行线 32 项结构／完整性检查和 55 组数值比较全部通过，800 条重叠 bootstrap 区间、22 个交付文件哈希及 11 个受保护输入均已核验。正式结果：`analysis/codex_outputs/RQ6_execution_2026-11-15/run_20261115_151147/RQ6_results.md`；验收：`analysis/audit_outputs/RQ6_acceptance_run_20261115_151147/ACCEPTANCE_RESULT.md`。</span>
- [x] <span style="color:#1f4e79">RQ6 显示预测器间存在分数偏移和有限排序一致性，逐靶点相对表现随预测器改变，且较低原始 Brier 不能写成“已经校准”。整合表述见 `analysis/paper_results_RQ1_RQ6_draft_2026-11-15.md` 第 8 节。</span>

---

## 14. 结论措辞边界

允许：

- “在与 Overath 镜像精确序列不重叠的当前主队列中……”
- “六个靶点并非六个完全独立设计流程……”
- “AAH-I 使用 SPR 主 assay，BLI 为选择性确认子集……”
- “原始 ipTM 的靶点内校准表现为……”
- “未适配 LOTO 与描述性适配分解分别显示……”

禁止：

- “验证了 Overath 完整 pipeline”；
- “证明了 Boltz‑1 → Boltz‑2 校准漂移”；
- “六个独立外部队列”；
- “assay 间结果稳健”；
- “设计随机效应已解决序列或来源相关”；
- “样本内 I／F 改善证明适配后预测更好”；
- “蛋白特征导致校准差异”；
- “所有 minibinder 的通用概率或阈值”。

---

## 15. 版本冻结规则

以下变更必须另起新版本并重建全部结果：

- 六靶点范围；
- AAH-I 主 assay；
- 精确去重键；
- 主结局多数规则；
- strength 聚合规则；
- 严格结局规则；
- ipTM 配对规则。

以下可以作为版本化分析注释追加，但不得改变已冻结标签：

- campaign／source 字典；
- 序列簇和 identity 审计结果；
- de novo 来源注释；
- ipTM 是否参与设计筛选的注释；
- 统计模型配置和随机种子。

**后续所有 RQ1–RQ4 从 `analysis/v3_outputs/pair_master_v3.csv` 读取核心数据，不再各自解析原始 CSV 或重新定义结局。**

---

## <span style="color:#1f4e79">16. 2026-11-09 讨论决议：结局审计与探索性扩展</span>

> <span style="color:#1f4e79">本节为新增分析注释，不改变第 1 节已经冻结的主队列、主 assay、主结局、严格结局及原始 Boltz‑2 ipTM 主暴露。</span>

### <span style="color:#1f4e79">16.1 主 assay 与结局规则的最终处理</span>

- <span style="color:#1f4e79">主 assay 继续按覆盖完整性确定，不按 strength、阳性率、敏感性或哪个 assay 使 ipTM 表现更好来选择；不得形成用预测指标反向选择真结局的循环论证。</span>
- <span style="color:#1f4e79">当前数据不能公平比较 BLI 与 SPR 的技术可信度。AAH-I 的 18 条 BLI 是选择性确认子集，不能用于估计技术敏感性、特异性或 assay 稳健性。</span>
- <span style="color:#1f4e79">`splitgfp-assay` 对应 `expression-yield`，不是结合结局，不能与 BLI／SPR binding 标签互换。</span>
- <span style="color:#1f4e79">重复标签聚合与 strength 是不同维度：前者处理同一 assay 内多条记录，后者定义实验成功的强度门槛。</span>

| 结局规则 | 可用配对 | 阳性数 | 处理决定 |
|---|---:|---:|---|
| 多数 binding，含 Weak | 2,084 | 748 | 冻结主结局 |
| 全部记录阳性 | 2,084 | 747 | 标签审计，不单独建模 |
| 任一记录阳性 | 2,084 | 748 | 标签审计，不单独建模 |
| 仅保留标签一致配对 | 2,083 | 747 | 标签审计，不单独建模 |
| 严格 Medium／Strong | 2,084 | 678 | 必做严格结局敏感性分析 |

<span style="color:#1f4e79">精确去重后仅 1／2,084 个配对存在 binding 标签不一致；多数与任一阳完全相同，全阳仅改变 1 个配对。因此 S3 不再分别平行拟合三套模型，只报告上表并保留审计字段。严格结局仍用于检验结局强度定义的实际影响。</span>

### <span style="color:#1f4e79">16.2 探索性 RQ5：少量界面指标的增量价值及其与 RQ3 的承接</span>

<span style="color:#1f4e79">**状态：Pinc 技术准备、结果性分析和主执行线验收均已完成（2026-11-15）。** RQ5 已检验 Pinc 单指标及 ipTM＋Pinc 组合能否改善 Boltz‑2 ipTM 的折外区分、概率误差和跨靶点迁移；ipSAE_min、LIS 和少量已有结构特征仍仅作后续探索。所有结果均为探索性，不回填为 RQ1–RQ4 的预设主要假设。正式结果与验收索引见第 13.E 节，整合表述见 `analysis/paper_results_RQ1_RQ6_draft_2026-11-15.md` 第 7 节。</span>

<span style="color:#1f4e79">**与 RQ3 的关系预先固定：** 无论 RQ3 的 U／I／F 结果是否理想，均不修改原始 ipTM 主分析。若零标签迁移或适配不足，RQ5 回答“加入 Pinc 是否能在相同信息层级下产生折外增量”；若 RQ3 表现良好，RQ5 仍只检验增量价值。不得看到 RQ3 结果后临时把组合模型改成主模型。</span>

#### <span style="color:#1f4e79">Pinc 数据基础与共同样本</span>

- <span style="color:#1f4e79">扩展表为 `analysis/pinc_batch/run_20261110_181500/pinc_pair_extension.csv`，按 `pair_id` 连接，不改写冻结主表。2,084 行中 Pinc 成功 2,074，`missing_input=10`，全部为缺少 CIF／PAE 的 PD-L1；其余下载、QC 和 Pinc 失败均为 0。</span>
- <span style="color:#1f4e79">10 个缺失配对中有 8 个主结局阳性、4 个严格结局阳性；不得插补。Pinc 基础分析和 ipTM 对照均使用相同的 2,074 个 complete-case 配对，不能直接与原始 2,084 对结果作非同样本优劣比较。</span>
- <span style="color:#1f4e79">独立验收覆盖全部 2,074 个文件／结果哈希、108 个深度 QC 独立复跑和 Nipah 681×681 强制复跑，均一致。验收记录：`analysis/audit_outputs/C03_pinc_batch_acceptance_2026-11-11/acceptance_audit.json`。冻结主表 SHA-256 仍为 `581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91`。</span>
- <span style="color:#1f4e79">批量统一使用官方 legacy `Pinc.R`，commit `57afb08e29fe97379db5b155da5fb78c4c8715d1`，脚本 SHA-256 `caf054ea54198a28013173422cdd057caf7e9488132842cff06494a1c7a94687`。仅使用已闭环的对称化 Pinc，不解释 Pinc1／Pinc2 的链方向。</span>
- <span style="color:#1f4e79">若其他指标造成额外缺失，每一个模型比较都另行取共同 complete-case 并报告各靶点及阳性覆盖；不得让不同模型使用不同样本后直接比较指标。</span>

#### <span style="color:#1f4e79">先补基础诊断，再做组合与迁移</span>

1. <span style="color:#1f4e79">在 2,074 个配对上，按 RQ1／RQ2 的同一口径计算 Pinc 的逐靶点分布、AP／ROC-AUC、原始 Brier、log loss、校准参数和可靠性曲线；ipTM 在同一批配对重算作可比对照。原始 2,084 对的 RQ1／RQ2 结果及 RQ3 主分析设定保持不变。</span>
2. <span style="color:#1f4e79">Pinc 不替代原始 ipTM 主暴露，不直接套用外部 0.45 阈值。单指标 Pinc 的跨靶点映射同样遵守 RQ3 的 LOTO、靶点等权主分析、配对加权敏感性、序列簇隔离和五靶点／严格结局对照。</span>
3. <span style="color:#1f4e79">ipTM＋Pinc 组合复用 RQ3 的 U／I／F：U 只在源靶点拟合组合映射并原样迁移；I 只用允许的基率信息更新截距；F 在新靶点内部折外评价。组合模型与同折、同样本的 ipTM-only 模型比较，不与不同划分的结果比较。</span>
4. <span style="color:#1f4e79">组合模型以少量特征的惩罚 logistic 为默认。所有特征选择、变换、正则化、超参数调节和校准均只在训练折内完成；需要调参时使用训练折内部嵌套验证，测试折不参与任何选择。</span>
5. <span style="color:#1f4e79">同时报告折外 AP／ROC-AUC、Brier、log loss、校准截距／斜率和曲线，以及逐靶点方向和覆盖量；ML 输出不自动代表良好校准。</span>
6. <span style="color:#1f4e79">随机森林等灵活 ML 仅作为可选探索基准，必须嵌套训练并另行校准，不替代可解释的少量特征主扩展模型，也不进入 RQ3 主分析。</span>

### <span style="color:#1f4e79">16.3 探索性 RQ6：外部队列内的同配对跨预测器描述——已完成并验收</span>

<span style="color:#1f4e79">**状态：正式计算、独立参考核验和主执行线验收均已完成（2026-11-15）。** RQ6 降级回答：在 Anthropic 相同外部设计、相同靶点构建体和相同湿实验结局上，10 种完整预测器管线的 ipTM 数值偏移、设计排序、靶点内区分和原始概率误差是否一致。RQ5 关注同一 Boltz-2 体系内的指标组合，RQ6 只描述更换预测器后的结论一致性，不承担固定映射原样迁移。</span>

1. <span style="color:#1f4e79">正式队列为 RQ3-E 的 8 靶点／720 个共同 `uuid`；10 个预测器在每个设计上均有现成 ipTM。Weak+ 为 287 阳性，Medium+ 为 254 阳性。</span>
2. <span style="color:#1f4e79">已报告配对分数偏移、相关、AUC／AP、原始 Brier／log loss 及靶点等权摘要；未做新的 U／I／F、逐预测器校准模型搜索或最佳预测器选择。</span>
3. <span style="color:#1f4e79">每靶点 2,000 次设计层配对 bootstrap 均完成且有 2,000 次有效重采样；800 条重叠区间与对齐后的独立参考逐值一致，最大端点差为 3.43×10⁻⁷。</span>
4. <span style="color:#1f4e79">具体选取、指标、区间和停止规则仍见 `analysis/RQ6_minimal_descriptive_protocol_2026-11-15.md`。正式结果与验收分别见 `analysis/codex_outputs/RQ6_execution_2026-11-15/run_20261115_151147/RQ6_results.md`、`analysis/audit_outputs/RQ6_acceptance_run_20261115_151147/ACCEPTANCE_RESULT.md`。</span>
5. <span style="color:#1f4e79">核心结论：预测器间既有分数偏移，也有有限排序一致性；Weak+／Medium+ 的靶点等权平均 AUC 范围分别为 0.6429–0.7182／0.6341–0.7239。没有单一预测器在所有靶点一致最佳，且较低原始 Brier 不等于已校准。整合结果见 `analysis/paper_results_RQ1_RQ6_draft_2026-11-15.md` 第 8 节。</span>

### <span style="color:#1f4e79">16.4 Anthropic 外部数据：实验终点、CRO 与预测器口径</span>

<span style="color:#1f4e79">Anthropic 已用于 RQ3-E 外部验证，并提供 RQ6 所需的队列内同配对十预测器数据。它与原 Proteinbase 主队列没有精确序列重叠，因此不能替代主队列内的跨预测器验证；但预测器在其内部共享同一 `uuid` 和结局，配对性成立。不同预测器的 MSA、recycles、采样和实现设置共同构成管线效应，不能只归因于模型权重。上游计算筛选、单一设计生成流程及基率覆盖限制必须保留。</span>

### <span style="color:#1f4e79">16.5 更新后的执行顺序</span>

1. <span style="color:#1f4e79">**RQ1／RQ2：已完成并验收。** 保留冻结结果、必要对照和结论边界。</span>
2. <span style="color:#1f4e79">**RQ3／RQ3‑E：已完成并验收。** 保留零标签 U 与有标签 I／F 的证据等级，不把外部分数生成差异解释为纯靶点效应。</span>
3. <span style="color:#1f4e79">**RQ4：已完成并验收。** 保留固定预算、原始 operating points 与计入成本的共享预算 pilot 结果。</span>
4. <span style="color:#1f4e79">**RQ5：已完成并通过主执行线验收。** 2,074 对基础诊断及相同划分下的单指标／组合折外增量分析均已完成；结论限定为探索性、指标与权重依赖的增量，不宣称通用迁移已解决。</span>
5. <span style="color:#1f4e79">**C04／RQ6 可用性审计：已完成。** 主队列第二预测器不可用；Overath 为非共同样本；Pinc 不是独立预测器。Codex 数据核查通过，但其“只能非配对说明”的建议未被采纳，分歧已留痕。</span>
6. <span style="color:#1f4e79">**最小 RQ6：已完成并通过主执行线验收。** Anthropic 内部 720 个共同 `uuid`、10 个预测器、两种结局和配对 bootstrap 均已核验；结果只作队列内描述性补充，不改写 RQ1–RQ5，也不宣称固定映射跨预测器迁移已经解决。</span>
7. <span style="color:#1f4e79">**截至 2026-11-15 的结果整合：已完成。** 当时主稿为 `analysis/paper_results_RQ1_RQ6_draft_2026-11-15.md`。此后已进入正常英文论文 v2 阶段；最新接续状态见第 17 节，本节不再作为当前稿件状态。</span>

---

## <span style="color:#1f4e79">17. 主线程接续状态（2026-11-17）</span>

> <span style="color:#1f4e79">本节同步本轮新增工作，覆盖第 13.E 节“灵活 ML 仅作后续探索”和第 16.5 节“下一步转化为短篇论文”的旧状态；不覆盖任何冻结数据、旧模型结果或历史验收记录。</span>

### <span style="color:#1f4e79">17.1 当前论文与稿件状态</span>

- <span style="color:#1f4e79">**论文定位：**按导师建议扩展为正常研究论文，不再以短篇报告为当前目标。现用英文主稿为 `manuscript/PEDS_Research_Article_English_draft_v2_2026-11-17.md`；合并审阅 PDF 为 `manuscript/review/PEDS_Main_and_Supplement_Review_Draft_v2_2026-11-17.pdf`。</span>
- <span style="color:#1f4e79">**RQ4E 已入稿：**冻结 outer-test 结果上的事后截距诊断显示，I 相对 U 的靶点等权 Brier 在六靶点下降 45.1%，独立重拟合的五靶点分析下降 42.9%，分别为 6／6、5／5 靶点改善。这只量化 intercept-addressable error，不证明基率失配是最大因果机制。结果：`analysis/rq4e_outputs/run_20261117_162929/RQ4E_RESULTS.md`。</span>
- <span style="color:#1f4e79">**RQ4-C 已入稿：**共享预算模拟保留 pilot 成本和 pilot hits；有利情形的收益不推广为普遍节省，较小预算或较大 pilot 可削弱或逆转收益。RQ4E 与 RQ4-C 均保持事后探索性证据等级。</span>
- <span style="color:#1f4e79">**RQ5-RF 尚未入稿：**本节下述随机森林结果没有写入英文 v2、正式补充材料或 v2 PDF。后续整合须另存新稿件版本，不覆盖 v2。</span>
- <span style="color:#1f4e79">**仍是科学初稿而非最终投稿包：**外部数据集正式引用、Proteinbase 精确版本、公开仓库／DOI、作者声明和期刊排版仍需最终处理。</span>

### <span style="color:#1f4e79">17.2 RQ5-RF：协议、执行与核对均已完成</span>

- <span style="color:#1f4e79">草案：`analysis/RQ5RF_model_complexity_protocol_draft_2026-11-17.md`；唯一一轮 Codex 审查：`analysis/codex_outputs/RQ5RF_protocol_review_2026-11-17/CODEX_REVIEW.md`。</span>
- <span style="color:#1f4e79">冻结协议：`analysis/RQ5RF_model_complexity_frozen_protocol_2026-11-17.md`，SHA-256 `d9db1c3a09d1264580ca56b797922d55fce0e563c04ecb9dab76ac4744a25e28`。审查提出的 1 项 CRITICAL、5 项 MAJOR 和 2 项 MINOR 已一次性纳入；不再开第二轮方法审查。</span>
- <span style="color:#1f4e79">正式运行：`analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_20261117_220855/`，状态 `complete`。34 个 outer source tasks、160 个 inner tasks 和 3,492 次 500-tree RF 拟合均完成。</span>
- <span style="color:#1f4e79">16／16 单元测试、26／26 QC 和独立验证通过；未解决失败为 0。60 个不可评价项均为事后校准截距／斜率诊断的近分离或边界问题，已保留且未纳入这些诊断参数的汇总；主要预测指标完整。</span>
- <span style="color:#1f4e79">主执行线已核对状态、质控、独立复算、失败台账、主要结果表和四张图；这是交付文件与关键结果核对，不宣称重新执行全部拟合的完整第三方复现验收。</span>

### <span style="color:#1f4e79">17.3 RQ5-RF 的结果与措辞边界</span>

主分析 `primary_six × target_equal` 的 U 层靶点等权宏平均：

| 模型 | Brier | Log loss | BSS-transport | AUC | AP |
|---|---:|---:|---:|---:|---:|
| pooled logistic ipTM | 0.282607 | 0.759039 | −0.000133 | 0.651489 | 0.600158 |
| pooled logistic combo | 0.322258 | 0.856314 | −0.139389 | 0.617534 | 0.572235 |
| RF ipTM | 0.262443 | 0.725385 | +0.068965 | 0.568358 | 0.553606 |
| RF combo | 0.289445 | 0.786984 | −0.023982 | 0.593172 | 0.557379 |

- <span style="color:#1f4e79">**Pinc 的 RF 家族内增量不稳定且宏平均方向不利：**在六个配置中，RF-combo 相对 RF-ipTM 的宏平均 Brier、log loss 均升高，BSS-transport 均降低。这是六配置宏平均方向，不得写成每个靶点都变差。</span>
- <span style="color:#1f4e79">**模型复杂度没有稳定救出组合：**RF-combo 相对 pooled-logistic-combo 只在部分配置改善；排除 AAH-I 后方向反转。冻结 10.2 筛查未通过，条件 1–3 失败、条件 4 通过。</span>
- <span style="color:#1f4e79">**保留概率误差与区分度分离：**RF-ipTM 相对 pooled-logistic-ipTM 在六配置中均降低宏平均 Brier／log loss、提高 BSS-transport，但 AUC／AP 均下降；BSS-target-oracle 仍全为负，主分析为 −0.715228。因此不能写成“随机森林完全无效”，也不能宣称 RF 已校准或解决迁移。</span>
- <span style="color:#1f4e79">**对照差异必须明示：**本次 pooled logistic 不含源靶点截距，并按每个实际 outer／inner 拟合行集重新计算源权重；它不同于原 target-intercept RQ5。RQ5-RF 只回答模型复杂度敏感性，不直接改写 RQ3 或已验收 RQ5。</span>
- <span style="color:#1f4e79">本分析不识别纯交互，不能分离 Pinc 主效应、ipTM／Pinc 非线性和交互；没有外部 Pinc 队列，不作未来靶点总体推断。</span>

### <span style="color:#1f4e79">17.4 唯一下一步与停止规则</span>

1. <span style="color:#1f4e79">**下一步：**把 RQ5-RF 作为 RQ5 的事后探索性模型复杂度补充写入英文 Methods、Results、Discussion 限制和 Supplement，保留“RF-ipTM 概率误差改善但排序未改善”的细微结果；随后统一表图编号并生成新审阅 PDF。</span>
2. <span style="color:#1f4e79">**当前不再做：**不追加第二轮协议审查，不重跑已完成的 3,492 次拟合，不继续尝试 XGBoost、神经网络、SVM 或其他模型寻找有利结果，也不把 RQ5-RF 升级为新主分析。</span>
3. <span style="color:#1f4e79">**研究阶段：**主体分析至此停止扩张；当前主线是整合、核查措辞和完成投稿材料，而不是继续结果驱动地增加 RQ。</span>
