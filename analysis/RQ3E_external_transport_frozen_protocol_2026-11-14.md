# RQ3-E：Boltz-2 ipTM 概率校准的外部迁移验证——冻结协议

**冻结日期：2026-11-14**  
**协议状态：已在外部模型性能计算前冻结；尚未运行 RQ3-E 源模型拟合、外部预测或性能评价。**  
**研究定位：原六靶点 RQ3 的独立外部迁移分析，不重写原 RQ3，不扩展为 RQ5／RQ6。**

---

## 0. 冻结声明与本轮边界

本协议固定 RQ3-E 的队列、分析单位、结局、外部 ipTM 分数、源模型、运输公式、参照基率、评价指标、分层和敏感性分析。外部标签已在数据审计阶段用于样本量与阳性数核查，但**尚未用于外部模型拟合、模型选择、阈值选择或性能评价**。

本文件是分析协议，不是结果报告。下列状态必须严格区分：

- [x] 外部数据来源、样本量、标签覆盖、模型覆盖和序列独立性审计完成。
- [x] 外部队列、主结局、主 score、U 公式和评价规则冻结。
- [ ] 在原六靶点上拟合一次部署用 U。
- [ ] 生成外部 720 行 U 预测。
- [ ] 计算外部 Brier、log loss、BSS、校准和可靠性结果。
- [ ] 执行任何 I／F 描述性诊断。

本轮只做协议冻结，**不运行上述未勾选项目**。不得修改 v1、v2、v3、原 RQ3 结果、原 LOTO 折划分、Codex 原始交接文件或冻结主表 `analysis/v3_outputs/pair_master_v3.csv`。

---

## 1. 研究问题与可回答范围

### 1.1 主要问题

> 在原 Proteinbase 六靶点主队列上训练的零标签运输模型 U，能否在不使用外部标签调参的前提下，迁移到 Anthropic 公开队列的 8 个达标靶点，并给出具有可解释概率校准的结合概率？

这里评价的是**已发布 Boltz-2 ipTM 分数及原 U 条件映射的外部运输**，不是重新发现“ipTM 能预测结合”或“校准具有靶点依赖性”。

### 1.2 预定义分层

1. **外部整体：8 个靶点、720 个设计—靶点配对。**这是主外部队列。
2. **相对原六靶点的新增靶点：4 个。**BHRF1、IL-7Ra、RBX1、TrkA，共 360 行。
3. **同一生物靶点的新 campaign：4 个。**Nipah-G、PD-L1、TREM2、VEGF-A，共 360 行。

“同靶点新 campaign”不自动代表相同构建体、表位或实验条件，也不增加独立生物靶点数。两层结果分别报告；不得只展示表现有利的一层。

### 1.3 不可外推的范围

- 结论只适用于经过上游计算筛选并进入实验的候选设计。上游 ipSAE 等筛选造成 range restriction，不能外推到随机生成序列。
- 确认性外部队列在主结局 Adaptyv Weak+ 下的靶点基率支持范围仅为 **21.1%–72.2%**；整合结局为 **21.1%–80.0%**。它没有覆盖原 RQ3 中 Nipah 的 **9.8%（101／1,028）**和 SARS-CoV-2 RBD 的 **84.3%（226／268）**两端。因此，确认性结果只能解释为该合格靶点基率范围内的外部运输表现，不能据此声称对极低或极高基率靶点同样有效。第 3.4 节的低事件压力测试只补充低基率端的描述性证据，不扩大确认性队列，也不消除高基率端未覆盖的限制。
- “外部”指相对当前冻结队列的来源和序列独立，不表示实验时间更晚、Boltz-2 训练完全未见、结构或表位独立，也不表示该数据从未被其他评分研究使用。
- 新旧队列的 ipTM 发布口径不同，因此不能把结果解释为完全同一 score-generation protocol 下的纯靶点效应。

---

## 2. 与原 RQ3 的关系

原 RQ3 的主口径保持为：

`primary_six + sequence_cluster90 + target_equal + linear`

原 LOTO 的 U／I／F 结果、折外预测和五靶点敏感性分析均保持不变。RQ3-E 新增的是一次**全六靶点源域内部署拟合＋外部队列零标签评价**，不是：

- 把外部靶点加入原 LOTO 后重选模型；
- 把原六靶点分析事后改写成十靶点或十四靶点分析；
- 对五个 LOTO 模型简单平均后当作部署模型；
- 用外部标签更新 U 后仍称其为零标签外部验证；
- 启动 Pinc／RQ5、跨预测器 RQ6、随机森林或其他新模型搜索。

---

## 3. 队列构建与冻结纳入流程

### 3.1 原始数据与筛选流程

外部队列按以下顺序形成：

1. Anthropic 原始发布包含 **1,440 个唯一设计、16 个靶点**，均有 binder 序列及设计靶点下的 Boltz-2 结果。
2. Mature GDF-8 的 **120 个设计**因抗原聚集及非特异传感器结合而没有可解释湿实验结局，排除；剩余 **1,320 个具有整合结局的设计**。
3. 在具有整合结局及 Boltz-2 ipTM 的设计中，按既定靶点质量门保留：**每靶点 ≥80 个可用配对且 ≥15 个整合阳性**。
4. 8 个靶点达标，每靶点 90 行，共 **720 行**；这 720 行也同时具有 Adaptyv 主结局和主 ipTM，无额外缺失。
5. 质量门不是按 ipTM 高低筛选，也不是按外部模型表现筛选。所有 8 个达标靶点均冻结保留。

**资格门槛等价性核查：**若把质量门改为“每靶点 ≥80 个具有 Adaptyv Weak+ 结局及 Boltz-2 ipTM 的配对，且 ≥15 个 Adaptyv Weak+ 阳性”，选出的仍是完全相同的 8 个靶点；未达标靶点的 Adaptyv Weak+ 阳性数均不超过 14。因此，本队列不存在“按整合结局入选、再按主结局改变靶点集合”的实际选择差异。该等价性只用于透明说明纳入口径，**不消除 ≥15 阳性门所造成的低基率覆盖限制**。

### 3.2 冻结靶点、结局人数与基率

下表在外部性能计算前固定。百分比仅用于描述队列，不用于平衡样本或选择靶点。

| 靶点 | 相对原六靶点的分层 | 原队列对应 | 整合 n／阳性 | 整合基率 | Adaptyv Weak+ n／阳性 | 主结局基率 | Adaptyv Medium+ 阳性 |
|---|---|---|---:|---:|---:|---:|---:|
| BHRF1 | 新增靶点 | 无 | 90／23 | 25.6% | 90／20 | 22.2% | 18 |
| IL-7Ra | 新增靶点 | 无；不属于原六靶点 | 90／49 | 54.4% | 90／47 | 52.2% | 36 |
| RBX1 | 新增靶点 | 无 | 90／28 | 31.1% | 90／23 | 25.6% | 20 |
| TrkA | 新增靶点 | 无 | 90／20 | 22.2% | 90／19 | 21.1% | 19 |
| Nipah-G | 同靶点新 campaign | Nipah | 90／19 | 21.1% | 90／19 | 21.1% | 17 |
| PD-L1 | 同靶点新 campaign | PD-L1 | 90／39 | 43.3% | 90／39 | 43.3% | 35 |
| TREM2 | 同靶点新 campaign | TREM2 | 90／72 | 80.0% | 90／65 | 72.2% | 62 |
| VEGF-A | 同靶点新 campaign | VEGF110 | 90／54 | 60.0% | 90／55 | 61.1% | 47 |
| **合计** | **8 靶点** | — | **720／304** | **42.2%** | **720／287** | **39.9%** | **254** |

新增靶点子集为 360 行、120 个整合阳性、109 个 Adaptyv Weak+ 阳性；同靶点新 campaign 子集为 360 行、184 个整合阳性、178 个 Adaptyv Weak+ 阳性。

### 3.3 分析单位与重复测量

- 一行一个 `uuid`，对应一个唯一设计—靶点配对。
- 同一设计在 Adaptyv 与 Twist 的测量是两个结局视角，**不得复制成两个独立统计样本**。
- 主队列不按两家 CRO 的共同测量交集缩减；每个结局使用其自身有效分母。
- 未测量、未交付结果或不可解释的实验不当作阴性；逐结局报告可用 n 与阳性数。
- 不重新对各敏感性结局套用 80／15 门槛，不因其阳性较少而删除靶点。例如 TrkA 的 Twist QC 阳性较少，仍须保留并标明估计不稳定。

### 3.4 非达标低事件靶点：描述性压力测试

为直接检查确认性质量门遗漏的低基率端，在外部性能计算前增加以下**非确认性、零标签描述性压力测试**。这些靶点不进入 8 靶点／720 行主队列，不与主队列合并计算宏平均，也不用于选择或修改模型。

| 靶点 | Adaptyv Weak+ n／阳性 | 主结局基率 | 预定解释限制 |
|---|---:|---:|---|
| TNFα | 150／12 | 8.0% | 低事件，逐箱与校准估计不稳定 |
| Cas9 | 82／10 | 12.2% | 低事件；RNP／apo 与测定形式复杂性须单列 |
| EGFR | 90／10 | 11.1% | 低事件 |
| BBF-14 | 79／1 | 1.3% | 单阳性；不作校准斜率的实质解释 |
| 15-PGDH | 30／1 | 3.3% | 极小样本、单阳性；不作校准斜率的实质解释 |
| MBP | 90／0 | 0.0% | 零阳性；不拟合校准截距／斜率，BSS-target-oracle 不可评价 |

上表使用与主分析一致的 **Adaptyv Weak+** 人数，而不是整合结局人数；因此 Cas9 是 82／10，BBF-14 是 79／1。Latent GDF-8 虽有样本量和阳性数不足（55／14），但主结局基率约为 25.5%，不属于本项低基率压力测试，仍不进入确认性队列；Mature GDF-8 因没有可解释结局而继续排除。

压力测试冻结以下规则：

1. 使用与主分析完全相同的全六靶点 target-equal 线性 U、主 ipTM score、变换和 Adaptyv Weak+ 结局；不得在压力测试靶点上重训、更新截距／斜率或调整阈值。
2. 逐靶点报告 n、阳性数、观察基率及 Wilson 95% 区间、平均／中位预测概率、Brier、log loss、BSS-transport，以及原固定分箱的可靠性表和各箱 Wilson 95% 区间。
3. 校准截距／斜率只在两类结局均存在且原数值 QC 通过时列出，并标注为低事件描述性诊断。MBP 不拟合校准诊断；BBF-14 和 15-PGDH 即使数值上可拟合，也不作强结论。
4. BSS-target-oracle 仅在参照分母非零时报告；MBP 的零阳性 oracle 参照使该指标不可评价，不得改成 0 或删除说明。
5. 压力测试不设置“通过”门槛，不汇总成总体成功声明。其结果只能揭示低基率端的失效风险或提供有限支持，不能证明极低基率校准可靠，也不能补救确认性队列的基率覆盖限制。

主确认性队列每靶点也仅 90 行、19–65 个主结局阳性，校准斜率标准误可能较大。无论主分析还是压力测试，均不得把单个靶点的斜率点估计或其置信区间包含 1 作为强结论；主分析关注逐靶点效应量、方向一致性和整体可靠性。

---

## 4. 序列独立性与源训练行处理

已完成的审计结果为：

- 1,440 个外部设计与冻结 Proteinbase 4,280 条序列的精确重叠为 **0**。
- 与本地 Overath HF 3,669 条唯一 binder 序列的精确重叠为 **0**。
- Parasail 1.3.4 全对全 Smith–Waterman 局部比较中，以“identity ≥80% 且 query／reference 双向覆盖均 ≥80%”为近重复边，对 Proteinbase 和 Overath 均为 **0**；identity ≥90% 时也为 0。

这些结果支持外部—源队列的序列非重叠，但不证明生成模型训练独立、结构趋异或表位独立。

RQ3-E 的全六靶点部署拟合没有外部留出靶点参与训练，因此**不套用原 LOTO 的 train-only 划分或留出靶点排除表**。源模型使用原 `primary_six` 的全部 **2,084 行**；外部审计不用于重写原主表、原序列簇或原折间排除记录。

---

## 5. 结局定义

### 5.1 主结局：Adaptyv Weak+

在 `adaptyv_valid=True` 的设计中：

- `adaptyv_binding_strength` 为 `weak`、`medium` 或 `strong`：阳性，`y=1`；对应派生列 `adaptyv_main_positive=True`。
- 其他已测试非结合结果：阴性，`y=0`。
- **已测试但不表达者按原流程计为 tested non-binder，即阴性。**
- 未测试、无交付结局或不可解释者排除，不插补为阴性。

主分析固定为 8 靶点、720 行、287 阳性。

### 5.2 预定义敏感性结局

| 结局 | 有效样本规则 | 阳性规则 | 解释边界 |
|---|---|---|---|
| Adaptyv Medium+ | `adaptyv_valid=True` | `medium`／`strong`；派生列 `adaptyv_strict_positive=True` | 更严格的强度口径；未测试仍排除 |
| 整合 `binder_final` | `integrated_valid=True` | `integrated_positive=True` | 固定规则结合两家证据并经 Claude 曲线复核；**不是第三种独立实验** |
| Twist QC | `twist_valid=True`，即该结局的已测量有效集 | `twist_hit_any_human_qc_pass=True` | 只把至少一个人源／靶抗原 QC-pass kinetic fit 计为阳性 |
| Twist permissive | `twist_valid=True` | `twist_hit_any_human=True` | 至少一个人源／靶抗原 kinetic fit；是宽松检测旗标，不是亲和力阈值 |

缺失实验不得在任何结局中转成阴性。所有敏感性结局均按预先规定报告，不根据结果选择“最合理”的结局。

---

## 6. 外部 Boltz-2 ipTM 主分数

### 6.1 主规则

外部主 score 使用派生列 **`boltz2_iptm_ipsae_best_seed`**，其生成规则固定为：

1. 只取设计靶点，`target_form` 为空的设计靶点记录；不使用 counter-screen 或替代 assay-form 记录。
2. 使用设计复合物自身化学计量：`1to1`、`1to2`、`1to3` 或 `rnp`，不跨化学计量取值。
3. 使用 Boltz **2.2.1**、`boltz2_conf.ckpt` 的 5 个 seed（标签 0–4）；无模板、靶点非配对 MSA、binder 单序列、每 seed 一个 diffusion sample、3 recycles。
4. 在 5 个 seed 中选择 **`ipsae_min` 最大**的 seed，再读取该 seed 的 **`iptm_pae`**。
5. 若 `ipsae_min` 并列，取**数值最小的 seed 标签**。
6. 这不是 max-ipTM；不得改用 `boltz2_iptm_seed_max` 作为主分数。

8 个靶点的 720 行均具有该主分数和完整 5-seed 记录，无需为本轮重新跑结构预测。

### 6.2 预定义 score 敏感性分析

每次只替换主 score，其他设置不变：

1. 五 seed ipTM 均值：`boltz2_iptm_seed_mean`。
2. 五 seed ipTM 中位数：`boltz2_iptm_seed_median`。
3. 固定 seed 0 的 `iptm_pae`：执行时从原始 `anthropic_cofold_predictions.parquet` 按同一设计靶点和化学计量规则提取；不得用“最优 seed”代替。

max-ipTM 不进入确认性分析，因为它引入额外的择优方向，且不是本协议的主发布口径。

### 6.3 与原队列的不可消除差异

原冻结 Proteinbase 主表只有单一 `boltz2_iptm`，缺少完整 seed 元数据，无法证明与外部“最佳 ipSAE seed 对应 ipTM”的生成口径完全一致。因此：

- 本研究可以称为**已发布分数的外部迁移验证**。
- 不能称为完全同协议、纯靶点迁移。
- score-generation protocol 是运输条件的一部分，也是局限；若未来要分离该因素，需新旧队列统一重跑，而不能只重跑外部。

---

## 7. 部署用 U：源域内一次性拟合

### 7.1 源训练集

- 输入：`analysis/v3_outputs/pair_master_v3.csv`。
- 分析集：`primary_six`。
- 结局：`primary_outcome`。
- 行数：2,084；靶点数：6。
- score：原表 `boltz2_iptm`。
- 主源权重：按原 RQ3 prep 的冻结规则从六个源靶点行数派生，而不是读取主表中不存在的权重列。对靶点 \(j\) 的每个源配对，原始权重为 \(1/n_j\)，使每个靶点的原始权重和均为 1；再将全部 2,084 行原始权重除以其行均值，使归一化权重的总体行均值为 1。该规则与原 LOTO 的 `weight_target_equal_raw`／`weight_target_equal` 完全一致，不使用任何外部信息。
- 任何外部行及其标签——包括 720 行确认性主队列和 521 行低事件压力测试——均不得进入拟合、权重构造、参数初始化、收敛补救或模型选择。

### 7.2 变换与模型形式

对源和外部 ipTM 使用同一变换：

\[
x=\operatorname{logit}\{\operatorname{clip}(s,10^{-6},1-10^{-6})\}.
\]

源域内对六个靶点分别设置截距，并共享一个斜率：

\[
\eta_{ij}=\alpha_j+\beta x_{ij},\qquad
\Pr(y_{ij}=1)=\operatorname{expit}(\eta_{ij}),
\]

其中 \(j=1,\ldots,6\)。设计矩阵由六个靶点截距指示列和一列 \(x\) 构成，**无额外全局截距**。

使用原 `fit_source_linear` 的未惩罚加权最大似然实现：目标函数为加权负对数似然之和除以权重和，L-BFGS-B 与既定确定性起点；参数边界 \([-40,40]\) 仅作数值保护。不得静默加入 ridge、改变权重或删除源靶点。

### 7.3 外部运输公式

主分析的运输截距为六个源靶点截距的**未加权算术平均**：

\[
\bar\alpha_U=\frac{1}{6}\sum_{j=1}^{6}\hat\alpha_j.
\]

共享斜率原样运输。对外部配对 \(i\)：

\[
\eta_{U,i}=\bar\alpha_U+\hat\beta x_i,
\qquad
\hat p_{U,i}=\operatorname{expit}(\eta_{U,i}).
\]

这与原 LOTO 的 target-equal U 运输定义一致；区别仅在于本次源训练使用全部六靶点，而不是五个源靶点。不得用外部基率更新 \(\bar\alpha_U\)，不得重新估计 \(\hat\beta\)，也不得对五个 LOTO 模型取平均替代该公式。

### 7.4 配对加权源模型：辅助分析

辅助源权重使用原表 `weight_pair`，模型形式不变。其运输截距为：

\[
\bar\alpha_{U,pair}=\frac{\sum_j n_j\hat\alpha_j}{\sum_j n_j},
\]

其中 \(n_j\) 为本次全六靶点源训练行数。共享斜率同样原样运输。该分析用于检验源权重选择敏感性，不与主分析竞争择优。

### 7.5 数值 QC 与停止规则

出现下列任一情况，部署 U 记为不可评价并停止外部性能解释：

- 优化不收敛、非有限参数或预测；
- Hessian 奇异／不可识别或严重病态；
- 参数触及数值边界；
- 最终梯度不满足原实现的 QC；
- 输入哈希、行数、靶点数或字段不符合本协议。

不得通过外部数据、改变模型形式或临时正则化“修复”后继续称为冻结外部验证。

---

## 8. 基率参照与评价指标

### 8.1 两个常数参照

**运输参照（可部署、只用源标签）：**主 target-equal 分析使用六个源靶点主结局患病／阳性率的等权平均：

\[
\pi_{transport}=\frac{1}{6}\sum_{j=1}^{6}\pi_j.
\]

本次全六靶点部署拟合没有 LOTO 留出排除，因此 \(\pi_j\) 用各源靶点在本次源训练集中的阳性率。所有外部行的参照预测均为该常数。pair-weighted 辅助分析的运输参照为源训练集总阳性数除以总配对数。

**靶点 oracle 参照（不可部署）：**对外部靶点 \(t\)：

\[
\pi_{oracle,t}=\frac{1}{n_t}\sum_{i\in t}y_i.
\]

该参照只用于诊断“若已知外部基率，常数模型能达到什么水平”，绝不能描述为零标签预测或可部署校准器。

### 8.2 逐靶点指标

每个外部靶点均报告：

- 可用配对数 \(n\)、阳性数、实际基率；
- Brier score：\(\operatorname{mean}(p-y)^2\)；
- log loss，概率仅在计算 log loss 时裁剪至 \([10^{-6},1-10^{-6}]\)；
- BSS-transport；
- BSS-target-oracle；
- 校准截距、校准斜率及各自标准误；
- 收敛／不可评价状态；
- 固定分箱可靠性表与 Wilson 95% 区间。

BSS 对任一参照预测 \(r_i\) 定义为：

\[
BSS=1-\frac{\sum_i(p_i-y_i)^2}{\sum_i(r_i-y_i)^2}.
\]

分母为 0 或非有限值时记为不可评价，不静默改成 0。BSS>0 仅表示优于**所指明的参照**，不等于获得通用校准。

### 8.3 校准诊断

在 U 外部预测全部固定后，按靶点拟合未惩罚、未加权逻辑回归：

\[
\operatorname{logit}\{\Pr(y=1)\}=a+b\operatorname{logit}\{\operatorname{clip}(p_U,10^{-6},1-10^{-6})\}.
\]

这是评价阶段的校准诊断，不是重新训练 U。若结局单类、预测非有限或模型不可识别，则该校准项记为不可评价；不得因此删除该靶点的 Brier／log loss 结果。校准斜率置信区间包含 1 不能单独作为“校准成功”的证据。

### 8.4 可靠性图

沿用原 RQ3 固定分箱边界：

`[0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0]`

采用原 `np.digitize(p, edges[1:-1], right=False)` 规则。逐靶点报告每箱 n、阳性数、平均预测、观察阳性率及 Wilson 95% 区间；空箱也保留。不得查看结果后更改分箱。

---

## 9. 汇总层级与主／辅分析

### 9.1 主汇总

- 对 8 个靶点先分别计算指标，再取**未加权靶点宏平均**。
- 同时报告靶点数、各指标最小值、最大值及 BSS 正／负方向计数。
- 宏平均 Brier 或 BSS 不得由先合并全部配对再计算来替代。
- 若某项校准诊断不可评价，明确报告该项可评价靶点数；不得把整行靶点从其他指标中静默删除。

### 9.2 预定义分层汇总

分别对以下两组做靶点等权汇总，并保留逐靶点原始结果：

1. 4 个新增靶点；
2. 4 个同靶点新 campaign。

8 靶点整体仍是主外部队列；两组分层用于区分“新靶点运输”和“同靶点跨 campaign 运输”，不作事后优胜子集选择。

### 9.3 配对加权辅助汇总

将全部可用外部配对合并后按每行等权计算辅助指标，与靶点等权结果分开报告：

- Brier 和 log loss 使用合并行；
- BSS-transport 使用固定的源运输常数；
- BSS-target-oracle 保留每个靶点自身 oracle 基率的逐行参照，不能改成全外部队列总基率；
- 合并校准结果可能受靶点间基率差异驱动，只作辅助描述。

本次不声称 8 个靶点足以估计未来靶点总体的不确定性，不新增以靶点为抽样单位的总体置信区间。

### 9.4 主分析与敏感性分析组合

**主分析只有一套：**

`8靶点外部队列 + Adaptyv Weak+ + 最佳ipSAE seed对应ipTM + 全六靶点target-equal线性U + 靶点等权评价`

敏感性分析采用**相对主分析一次只改一个维度**：

1. 结局：Medium+、整合、Twist QC、Twist permissive；
2. score：五 seed 均值、中位数、固定 seed 0；
3. 源权重：pair-weighted U 及其对应运输参照；
4. 汇总：4 新增／4 同靶点分层、配对加权辅助汇总。

不做全因子大规模结果搜索，不新增 RCS、Pinc、随机森林、其他机器学习或跨预测器比较。所有分析均由本协议预先指定，不按结果增删。

---

## 10. I／F 只允许作为标签依赖诊断

RQ3-E 的零标签外部结论只来自 U。若在 U 结果保存并冻结后继续做 I／F，仅限主结局、主 score、target-equal U，并明确标为描述性诊断：

- **Oracle I：**在每个外部靶点全部标签上求截距修正 \(a\)，使 \(\operatorname{mean}\{\operatorname{expit}(a+\eta_U)\}\) 等于该靶点实际基率；预测为 \(\operatorname{expit}(a+\eta_U)\)。
- **Descriptive F：**在同一外部靶点全部行上拟合未惩罚 \(\eta_F=a+\gamma\eta_U\)，再评价同一批行。

二者可报告 Brier、log loss、BSS-target-oracle 和校准诊断，但必须标注为 **oracle／in-sample apparent adaptation**。它们：

- 不是零标签迁移；
- 不是折外泛化；
- 不能用于证明部署后性能改善；
- 不得成为本协议的主结论。

本轮不新增外部靶点内的 predictive I／F 或嵌套交叉验证。如未来需要，必须另立协议并处理序列簇、campaign 分块及每靶点仅 90 行的精度限制。

---

## 11. Campaign、generator 与 assay 分析边界

### 11.1 Campaign／generator

逐靶点、逐 campaign 报告 n、阳性数、观察基率、平均预测和 Brier；同时描述 design model、generator、single-target／multi-target 流程及序列设计方法。campaign 是预定义分块与异质性描述变量，不把每靶点 90 行视为同质随机样本。

**设计流程的准确边界：**外部 720 行的高层 `design_model` 只有 Mythos Preview（480 行）和 Opus 4.8（240 行），但底层 `generator` 并非单一来源，而是 9 种：PXDesign 178、RFdiffusion3 138、Genie3 91、FreeBindCraft（BindCraft）85、BoltzGen 84、Proteina-Complexa 69、RFdiffusion 64、FoldCraft 10、BoltzDesign1 1。因此不能声称外部队列“未检验 RFdiffusion／BindCraft”或“只有单一生成器”。本研究评价的是一个不同且与原队列部分重叠的 generator 混合分布；由于 generator 与靶点、campaign、design model 相互关联，不能识别独立的 generator 因果效应，也不能外推至未代表的设计流程。

若某 campaign 单类或样本过少，校准斜率记为不可评价，不合并到有利 campaign。campaign 结果用于解释分布差异，不用于事后修改主模型。

### 11.2 两家 CRO／assay

Adaptyv 与 Twist 的对比只通过预定义结局敏感性、标签一致性和分层描述完成。不得把同一设计的两家结果当作独立重复样本。

Twist 使用 Expi293 表达的 Fc 融合二价配体及多循环 SPR；Adaptyv 使用细胞游离表达的单价配体，并包含 SPR／BLI。部分靶点还存在多聚抗原、不同抗原形式或表位遮挡。因此二者差异不能被解释为纯仪器因果效应，也不能单独识别“靶点本身”与“测定体系”的贡献。

同靶点新 campaign 与原靶点内结果的比较同样是描述性的；构建体、表位、campaign 和 score-generation protocol 均是候选运输条件。

---

## 12. 查看外部结果后的不可变规则

一旦查看任何外部性能结果，以下行为均不允许继续称为本冻结协议的确认性分析：

- 删除或替换表现不利的靶点、campaign 或结局；
- 更换 seed 规则、ipTM 变换、源权重、运输截距或斜率；
- 调整可靠性分箱、参照基率或指标以得到更有利结论；
- 在外部数据上重训后仍称零标签验证；
- 只报告有利的新增靶点或同靶点子集；
- 因 U 表现不佳而临时加入 Pinc、非线性或其他 ML，并仍归入本次外部确认。

若确需修改，须写入新的、带日期的修订协议，并将当前外部队列降级为**开发／探索集**；任何新的确认性声明需要新的独立 holdout。原协议和原结果不得覆盖。

本研究不设置事后“必须做成功”的目标。外部迁移失败、仅部分靶点可运输或明显基率漂移均是有效结果。

---

## 13. 输入哈希门

执行前必须逐一核对 SHA-256。任一不匹配即停止，不得静默重下载、重生成或继续分析。

### 13.1 原队列与复用实现

| 文件 | 预期内容 | SHA-256 |
|---|---|---|
| `analysis/v3_outputs/pair_master_v3.csv` | 2,084 行冻结主表 | `581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91` |
| `analysis/codex_tasks/RQ3_EXECUTION_CODEX_HANDOFF_2026-11-13/rq3_execution_config.json` | 原 RQ3 模型、参照和指标定义 | `a416ecd357eabb4c4d95701136c412cf23651b730bf39feab7524fbc07eeef51` |
| `analysis/codex_outputs/RQ3_execution_2026-11-13/run_20261113_175022/rq3_common.py` | 复用模型与指标函数 | `276f46ad6d068ecead50a9236d1dcf6953428e415eef79da2424ad850c471f3e` |

### 13.2 外部原始输入

均位于 `data/external_target_audit_2026-11-13/`。

| 文件 | SHA-256 |
|---|---|
| `anthropic_design_summary.parquet` | `3e796869b8dc0c90d7ab0e12daf25308c0cf272d04c6a5e27281f903ecd18744` |
| `anthropic_wetlab_summary.parquet` | `f80c3cfce67cb7da53576a020d476c5b74efed44e61c6fceaec9ee7de5722f88` |
| `anthropic_cofold_predictions.parquet` | `ba48010610daa9b9d184aa1ddb419beebdf728ccf745f5a8f3b48e4041e1a87d` |

### 13.3 已审计派生输入

均位于 `data/external_target_audit_2026-11-13/`。

| 文件 | 预期内容 | SHA-256 |
|---|---|---|
| `anthropic_audited_pair_table.csv` | 1,440 行审计配对表 | `1f2a4ceab9de06cb75f4aae9a889fc5164e144232101417653944d03ee33d0c6` |
| `anthropic_target_qualification.csv` | 16 靶点门槛表 | `bb3ac9e2a0f3bc253d7b500dcc4909f1aa2291026e56d2960f028f8911e50579` |
| `anthropic_audit_summary.json` | 原始输入哈希与审计汇总 | `df377fc5b755e41a5c2585b14d87f5bf73325b0c8b2d8f94c17381b59c066b81` |
| `anthropic_sequence_independence_summary.json` | 全对全序列独立性汇总 | `8549f59b535357dbe0003a9aa1fa7ffcdadee423bed6de3d9dcf18d0beaf1316` |

执行时另存完整运行配置快照、代码副本和输入／输出哈希清单。不得修改上述已冻结文件来满足哈希门。

---

## 14. 未来执行输出与验收标准

### 14.1 计划输出

未来经确认后，在新目录 `analysis/rq3e_outputs/` 下建立带时间戳的运行目录，不覆盖任何既有结果。至少应生成：

- `run_config.snapshot.json`：完整协议、输入路径、哈希和软件版本快照；
- `rq3e_external_pairs.csv`：720 行主外部队列及各结局有效性旗标；
- `rq3e_source_model_parameters.csv`：六个源截距、共享斜率、运输截距、SE 与 QC；
- `rq3e_u_pair_predictions.csv`：主分析及预定义辅助分析的逐行固定预测与参照；
- `rq3e_target_metrics.csv`：逐靶点指标；
- `rq3e_group_metrics.csv`：8 靶点、4 新增、4 同靶点及配对加权汇总；
- `rq3e_reliability.csv`：固定分箱可靠性表；
- `rq3e_campaign_metrics.csv`：campaign 描述与异质性，包括两种 design model 与九种 generator 的构成；
- `rq3e_sensitivity_metrics.csv`：结局、score 和源权重敏感性；
- `rq3e_low_event_stress_test_pairs.csv`：六个低事件靶点的 521 行描述性压力测试配对；
- `rq3e_low_event_stress_test_metrics.csv`、`rq3e_low_event_stress_test_reliability.csv`：低事件逐靶点指标、固定分箱和不可评价原因；
- `RQ3E_U_external_target_metrics.png`、`RQ3E_U_external_reliability.png`、`RQ3E_U_low_event_stress_test.png`；
- `qc_checks.json`、失败项清单、完整日志和 `final_hash_manifest.csv`。

这些文件当前均不代表已生成或已验收。

### 14.2 最低验收条件

1. 所有第 13 节哈希通过，冻结主表运行前后哈希不变。
2. 主外部队列恰为 720 个唯一 `uuid`、8 靶点、每靶点 90 行；主结局合计 287 阳性。
3. 新增／同靶点分层各 360 行；没有按性能删除靶点。
4. 主 ipTM 与五 seed 规则一致；固定 seed 0 敏感性来自原始逐 seed 表而非最佳 seed 列。
5. U 拟合只读取原六靶点 2,084 行；外部标签未参与源模型、参数初始化或调参。
6. 源模型通过数值 QC，720 行主预测均有限且位于 0–1。
7. 逐靶点、宏平均和配对加权结果按本协议分开计算；两个 BSS 参照逐行保存且没有混用。
8. 所有空箱、缺失结局、不可评价项和失败项显式保留，不作静默完整案例删除。
9. 敏感性分析只改变本协议指定的一个维度；没有新增结果驱动的模型搜索。
10. 结果报告明确区分 U、oracle I 和 descriptive F，且不包含超出数据支持范围的声明。
11. 低事件压力测试恰为 TNFα、Cas9、EGFR、BBF-14、15-PGDH、MBP 六靶点，合计 521 个唯一 `uuid`、34 个主结局阳性；逐靶点人数与第 3.4 节一致，不与 720 行主队列混合，不参与 U 拟合或调参。
12. 结果文件明确写出确认性主队列的 21.1%–72.2% 基率范围、原 RQ3 两端未覆盖、压力测试的非确认性属性，以及两种 design model／九种 generator 的真实构成。

---

## 15. 结果表述边界

无论结果方向如何，论文表述均遵守以下边界：

**可以说：**

> 在来源与既定序列近重复阈值下独立于原训练队列的公开数据上，对已发布 Boltz-2 ipTM 分数及预冻结 U 映射进行外部概率校准与条件映射运输评价。确认性结论限于 Adaptyv Weak+ 基率约 21%–72% 的 8 个达标靶点；低事件压力测试另行描述，不作为广泛基率覆盖的确认性证据。

**不能说：**

- 首次证明 ipTM 能预测实验成功或首次发现靶点依赖；
- 外部队列在实验时间上必然晚于 Overath；
- 数据从未用于任何评分或校准研究；
- 已证明对随机生成设计、极低／极高基率靶点或所有新 campaign 通用；
- 因低事件压力测试结果有利，就声称已确认极低基率校准可靠；
- 外部队列只有单一底层 generator，或没有包含 RFdiffusion／BindCraft 家族；
- 已把靶点、assay、campaign、generator 和 score-generation protocol 的效应完全分离；
- oracle 或同批数据适配后的表现等于零标签或折外泛化性能。

---

## 16. 审计与复现依据

本协议依赖且不重写以下已验收材料：

- `analysis/external_target_audit_2026-11-13.md`
- `analysis/audit_anthropic_external_2026-11-13.py`
- `analysis/audit_anthropic_sequence_similarity_2026-11-13.py`
- `data/external_target_audit_2026-11-13/anthropic_target_qualification.csv`
- `data/external_target_audit_2026-11-13/anthropic_audited_pair_table.csv`
- `data/external_target_audit_2026-11-13/anthropic_campaign_counts.csv`
- `data/external_target_audit_2026-11-13/anthropic_model_coverage.csv`
- `data/external_target_audit_2026-11-13/anthropic_sequence_independence_summary.json`
- 原 RQ3 配置与 `rq3_common.py`，路径和哈希见第 13 节。

**冻结结论：开跑前补充内容已并入原协议，且仍早于任何外部性能查看。RQ3-E 可以进入“源域内拟合一次部署 U，并在外部 720 行确认性主队列及 521 行低事件描述性压力测试队列上分别评价”的执行阶段；本文件本身不构成任何外部性能结果。**
