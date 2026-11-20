# RQ5-RF：模型复杂度敏感性分析冻结协议

日期：2026-11-17  
版本：`RQ5RF_v1_frozen`  
工作区：`<REPRODUCTION_ROOT>`  
状态：**已在唯一一轮 Codex 方法审查后一次性修订并冻结；不得根据外测结果修改。**

唯一正式输出根目录：`analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_YYYYMMDD_HHMMSS/`。

## 1. 目的与证据边界

已验收 RQ5 在相同 2,074 个 complete-case 配对、相同外测行上比较了 ipTM、Pinc 及其线性组合。组合的部分概率误差改善不稳定，且对源权重敏感。已有 RQ3 的限制性立方样条只检验了 ipTM 的平滑非线性，没有检验随机森林所能表达的更一般非线性关系。

本分析回答：

> 在完全相同的 complete-case 样本、外层验证、当前行集源权重和本地适配规则下，随机森林是否相对同特征的 pooled ridge logistic 产生稳定的折外概率预测增量；并检验在 RF 家族内加入 Pinc 是否相对 RF-ipTM 产生稳定的联合预测增量。

估计对象是**给定特征集下的模型家族条件预测差异**和**给定模型家族下加入 Pinc 的联合预测增量**。本设计不识别纯交互，也不可分离归因 Pinc 主效应、ipTM 非线性、Pinc 非线性及二者交互。

本分析是 **post hoc、exploratory model-complexity sensitivity**，不是算法竞赛、验证性比较或通用校准方案开发，不自动改写 RQ1–RQ4 或已验收 RQ5。只新增随机森林一个算法家族；无论结果如何，不再依次尝试 XGBoost、神经网络、SVM 或其他模型寻找有利结果。没有外部 Pinc 队列，因此任何改善仍缺少独立外部验证。

## 2. 不可变输入与旧结果

所有输入、行数、键、哈希和 folds 以已验收配置为准：

- RQ5 配置：`analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/run_config.snapshot.json`。
- 主表：`analysis/v3_outputs/pair_master_v3.csv`，2,084 行，SHA-256 `581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91`。
- Pinc：`analysis/pinc_batch/run_20261110_181500/pinc_pair_extension.csv`；仅 `pinc_status=success` 的对称化 `pinc`，2,074 成功，10 个 PD-L1 缺失不插补。
- 序列簇、LOTO、target-CV5 和 cross-target exclusions：使用上述配置中的冻结文件和哈希。
- 既有 RQ5：`analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720`；只读复用其预测、参数和实现，不覆盖或重写。
- 最终 RQ3 修正版实现：仅按既有 RQ5 的引用方式复用公式与适配规则。

预期共同样本：

| 分析集 | 行数 | 阳性数 | 结局 |
|---|---:|---:|---|
| `primary_six` | 2,074 | 740 | `primary_outcome` |
| `primary_five` | 1,941 | 721 | `primary_outcome` |
| `strict_six` | 2,074 | 674 | `strict_outcome` |

任何不符立即停止。`boltz2_iptm` 与 `pinc` 必须有限且位于 `[0,1]`；越界输入在变换前即失败，不得裁剪后静默接受。

不得修改 v1/v2/v3 流程、稿件、冻结主表、Pinc 扩展表、序列簇、folds 或旧结果。不得联网、安装或升级依赖。

## 3. 分析集、独立重拟合和外层验证

- 分析集：`primary_six`、`primary_five`、`strict_six`。
- 主分析：`primary_six × target_equal × sequence_cluster90`。
- 必做敏感性：`{primary_six, primary_five, strict_six} × {target_equal, pair}` 六个配置全部完成。
- 三个分析集各自从对应 complete-case 行、结局和冻结 folds **独立重拟合全部四个主分析模型**。`primary_five` 不得由 `primary_six` 的预测子集得到。
- 外层复用冻结 LOTO assignments 和 cross-target-cluster exclusions；仅过滤 10 个 Pinc 缺失行，不重分折。
- 主聚类规则为 `sequence_cluster90`；本分析不增加其他聚类规则。
- 结局、样本、外测标签、folds 或模型网格不得因新结果改变。

## 4. 每个实际拟合行集重新计算源权重

对每个 outer 或 inner 的实际拟合行集 \(S\)，先完成 Pinc complete-case 过滤、held-out target 移除和跨靶点序列簇排除，再从 \(S\) 重新计算源权重；**不得直接沿用冻结 fold 文件中的 `weight_target_equal` 或 `weight_pair` 数值**。

令 \(T_S\) 为 \(S\) 中仍存在的源靶点集合，\(n_t\) 为靶点 \(t\) 在 \(S\) 中的当前行数：

- `target_equal`：\(w_i=1/(|T_S|n_{t(i)})\)；
- `pair`：\(w_i=1/|S|\)。

两种权重均满足 \(\sum_i w_i=1\)。`target_equal` 中每个当前源靶点总权重严格为 \(1/|T_S|\)，数值误差不得超过 `1e-12`。

若不能形成预期靶点集合、任一预期靶点被排空、权重非有限或非正，则该拟合明确失败；不得静默删除一个靶点后继续。每次拟合保存行集键哈希、\(|T_S|\)、各靶点 \(n_t\)、各靶点权重总和、权重向量 SHA-256 和归一化误差。四个主分析模型在相同任务上使用完全相同的重算权重。

验证和外测指标不加权；靶点等权宏平均也不使用训练权重。

RF 的 `fit()` 只接收上述权重一次作为 `sample_weight`。不得按权重复制／删除行，不得另行过采样或欠采样，`class_weight=None`。配置必须冻结 scikit-learn 版本；已审查环境为 **1.9.0**，版本不符为运行阻断条件，不得跨版本拼接结果。

此“当前行集重算权重”是本新分析对概念规则的明确实现，只用于四个 pooled 模型。旧 target-intercept RQ5 锚点保留原验收预测，不重算、不纳入第 10.2 节门槛，也不宣称与新模型只有函数形式差异。

## 5. 同信息 2×2 模型对照

### 5.1 模型暴露

| 模型暴露 | 特征 | 函数形式 |
|---|---|---|
| `pooled_logistic_ipTM` | `x_ipTM` | pooled ridge logistic |
| `pooled_logistic_combo` | `x_ipTM`, `x_Pinc` | pooled ridge logistic |
| `rf_ipTM` | `x_ipTM` | random forest |
| `rf_combo` | `x_ipTM`, `x_Pinc` | random forest |

`x_ipTM=logit(clip(ipTM,1e-6,1-1e-6))`，`x_Pinc=logit(clip(Pinc,1e-6,1-1e-6))`。两个模型家族使用相同输入特征。禁止加入 target ID、源靶点阳性率、assay、campaign、序列、ipSAE、LIS 或其他协变量。

四个模型均不接收源靶点截距指示符。pooled 设计避免把“是否有源靶点截距”与“linear/RF”同时改变；代价是其结果不能与原 target-intercept RQ5 作纯函数形式分解。

既有 `ipTM_only`、`Pinc_only`、`ipTM_plus_Pinc` target-intercept U/I/F 预测只作次要描述锚点，不能作为模型家族效应的唯一对照。

### 5.2 预设对比及含义

U 层主要对比：

1. `pooled_logistic_combo − pooled_logistic_ipTM`：logistic 家族内加入 Pinc 的联合预测增量；
2. `rf_ipTM − pooled_logistic_ipTM`：单变量 ipTM 下的模型家族条件预测差异；
3. `rf_combo − pooled_logistic_combo`：双变量特征集下的模型家族条件预测差异；
4. `rf_combo − rf_ipTM`：RF 家族内加入 Pinc 的联合预测增量；
5. `rf_combo − pooled_logistic_ipTM`：完整 RF＋Pinc 与单变量 pooled logistic 的条件预测差异。

Brier／log loss 差值为左减右，负值利于左侧；AUC／AP／BSS 差值正值利于左侧。U、I、F 各层内分别比较，不混淆信息层级。

`rf_combo − rf_ipTM` 不得命名为“纯交互效应”。若结果可由 ipTM 非线性、Pinc 非线性或交互中的任一路径产生，必须同时保留这些解释。

## 6. Pooled ridge logistic 规则

- 一个未惩罚全局截距＋相应共同斜率；不得加入源靶点截距。
- 每个 outer／inner 训练集内，仅用当前训练行及第 4 节重算权重计算加权均值和总体 SD，标准化 x；验证／外测行沿用该变换。SD=0 记失败，不得用外测信息修复。
- 目标为加权平均 NLL＋`0.5*lambda*sum(beta^2)`；截距不惩罚。
- lambda 网格固定 `[0.001,0.01,0.1,1,10,100]`。
- 每个 outer 源训练集内部用 leave-one-source-target-out 选择 lambda；每个 inner fold 重新应用冻结跨靶点序列簇排除、第 4 节权重、标准化和拟合。
- 选择各 inner 留出靶点 log loss 的等权宏平均最小者；1e-12 内并列取较大 lambda。不得使用 outer 留出靶点标签。
- 最终用全部当前 outer 源训练行重拟合，直接预测 outer 留出靶点。
- 保存每折 lambda、全部 inner metrics、变换、权重审计字段和参数。失败不得静默改成无惩罚或其他模型。

## 7. 随机森林规则

使用已安装 scikit-learn 1.9.0 的 `RandomForestClassifier`，不得安装新包。

### 7.1 有效网格

基础网格：`max_depth ∈ {2,4,None}`、`min_samples_leaf ∈ {10,30}`、`max_features ∈ {1,2}`。

对维度为 \(d\) 的暴露，使用 `max_features_eff=min(max_features,d)`，并按键 `(max_depth,min_samples_leaf,max_features_eff)` 去重：

- `rf_ipTM`：6 个有效候选；
- `rf_combo`：12 个有效候选。

单特征暴露的 `max_features=2` 不得计作独立候选。

### 7.2 共享确定性 seed

同一 `analysis_set × source_weighting × outer fold × inner fold × 拟合阶段` 下的所有超参数候选、以及可配对的两个 RF 暴露，使用同一个固定 `random_state`。

seed 由 SHA-256 对规范字符串

`RQ5RF|v1|analysis_set|weighting|outer_fold|inner_fold|inner_or_final`

取前 4 字节并转为 `[0,2^32-1]` 的整数。字符串不包含超参数、标签值或任何结果值。最终 outer 重拟合使用 `inner_or_final=final` 的任务 seed，不复用 inner candidate seed。

相同输入、参数和 seed 重跑必须得到逐元素相同概率。`n_jobs` 只控制资源；`n_jobs>1` 仅在当前环境通过重复拟合逐元素一致性单元测试后允许。

### 7.3 冻结参数与选参

显式冻结：

- `n_estimators=500`、`criterion="gini"`、`bootstrap=True`、`max_samples=None`；
- `class_weight=None`、`min_samples_split=2`、`min_weight_fraction_leaf=0.0`；
- `max_leaf_nodes=None`、`min_impurity_decrease=0.0`、`oob_score=False`；
- `warm_start=False`、`ccp_alpha=0.0`、`verbose=0`；
- `max_depth`、`min_samples_leaf`、`max_features=max_features_eff` 按候选设置；
- `n_jobs` 和 seed 写入配置及每个候选记录；其他参数不得用于结果驱动调整。

选择规则：

- 与 logistic 共用当前 outer 源训练集内部的 leave-one-source-target-out 划分、序列簇排除和评价行。
- 每组参数在每个 inner 训练集重新拟合，使用第 4 节当前行集权重作为 `sample_weight`。
- 以各 inner 留出靶点 log loss 的等权宏平均选择最小者；不得用 AUC、Brier 或 outer 标签选参。
- 1e-12 内并列时依次选择较浅深度（`None` 视为最深）、较大 `min_samples_leaf`、较小 `max_features_eff`。
- 选定后在全部当前 outer 源训练行上以 final seed 重拟合，再预测 untouched outer 靶点。
- 保存全部有效候选的 inner metrics、状态、参数和 seed，不得只保存最佳组。

`predict_proba` 前硬性要求 `classes_` 恰为 `[0,1]`，按 `classes_==1` 取阳性列，不得使用固定列号。训练集只有一类、classes 不符或预测含非有限值均明确失败，不得补 0.5。

U 层不额外做 Platt、isotonic 或其他校准。RF 原始输出只是候选概率，不预设已校准。

## 8. U／I／F 信息层级与适配

对四个 pooled 模型暴露分别保存 outer LOTO `p_U`，并令 `eta_U=logit(clip(p_U,1e-6,1-1e-6))`。

本地适配复用已验收 RQ5 的 target-CV5、概率截断、F lambda 网格和失败规则，并冻结以下实现细节：

- 所有 source、outer-test 和 adaptation 数据在进入切分器前按 `pair_id` 的 UTF-8 字符串升序稳定排序；不得依赖 merge、CSV 或 groupby 的偶然顺序。
- 对每个 `analysis_set × held_out_target × outer target-CV5 fold`，只构造一次三折 `StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=20261113)` 和 inner assignments；同一 outer 训练行集上的四个暴露、两种 source weighting 及全部 F lambda 必须复用相同 inner assignments。保存 assignments SHA-256。
- U：不使用 destination-target 标签。
- I：`eta_I=a+eta_U`，只用适配训练折更新截距。
- F：`eta_F=a+gamma*eta_U=a+eta_U+delta*eta_U`，只更新总 logit 的截距和一个校准斜率；不重新训练 RF，不重新估计 ipTM／Pinc 的相对作用。
- 每个适配外测行不得参与对应适配训练或调参；source weighting 只影响固定 U base，适配训练对权重仍为 1。
- Pinc 过滤后重新检查每个 target-CV5 outer train/test 非空、类别可拟合、训练测试簇无交集；不满足只能记失败，不得重新分折。
- 每个 outer adaptation test 行只出现一次；逐靶点指标由该靶点五个 outer test fold 拼成，禁止重复或遗漏。
- U、I、F 及四个暴露在相同 outer-test rows 上比较。
- `eta_U`、`eta_I`、`eta_F`、`p_U`、`p_I`、`p_F` 必须全部有限；任一失败保留为不可评价项。

不得因 RF 输出形态改变适配定义，不增加新的校准算法。

## 9. 指标、完整样本汇总和图表

逐靶点、逐模型、逐层输出：n／阳性数、Brier、log loss、BSS-transport、BSS-target-oracle、BSS-adaptation-train（仅适配层）、ROC-AUC、AP、校准截距／斜率及状态、固定 bins 可靠性。

- 靶点等权宏平均为主，pair-pooled 汇总另列；AAH-I 单独可见，五靶点不可隐藏。
- 所有 BSS reference 与对应源训练行和当前权重规则一致。target-oracle 只作诊断，不是可部署比较器。
- 输出全部预设 pairwise differences、绝对差值、逐靶点方向、最小／最大差和改善靶点数。
- Brier 使用模型原始概率；log loss 使用 `1e-6` 数值裁剪。
- 不新增显著性检验或“最优模型”排行榜；不把六靶点方向计数解释为未来靶点总体推断。

完整样本硬门：

1. 每个 `analysis_set × weighting × model exposure × layer` 的靶点等权宏平均，只有在预期全部靶点均可评价时才成立；任一靶点不可评价时，该宏平均及依赖它的第 10.2 节条件均不可评价。禁止 `dropna()` 后对剩余靶点取均值。
2. pair-pooled 对比只在左右 `pair_id` 集合逐元素相等、无重复且预测均有限时计算；任何不符记为失败，不得 inner join 后静默缩小样本。
3. 同层跨模型比较也必须使用完全相同的外测 pair_id 集合。

至少生成：U 层四模型逐靶点 Brier／BSS-transport 图、主要 pairwise difference 图、I/F 同层比较图、primary_six target_equal 可靠性图。图必须实际打开检查，空 bins 和不可评价项保留。

## 10. 冻结解释规则

### 10.1 不支持改变当前结论的情形

以下任一情况不得宣称 RF 解决迁移或校准：

- 仅随机划分、靶点内拟合、训练集表现或 feature importance 改善；
- 仅 AUC／AP 改善，而 Brier／log loss／BSS-transport 未改善；
- 平均改善只出现在一个结局、一个权重或 AAH-I 主导的设定；
- 仅 I/F 使用本地标签后改善，而 U 层未改善；
- RF 相对一个对照改善，但相对另一个同信息对照没有一致方向。

### 10.2 探索性方向一致性筛查

数值方向先按以下中性带定义：

- Brier／log loss：差值 `< -1e-12` 为左侧改善，`> 1e-12` 为反转，其余为中性；
- BSS：差值 `> 1e-12` 为左侧改善，`< -1e-12` 为反转，其余为中性。

只有同时满足以下条件，才可写为：**“出现方向一致的探索性增量，值得在独立 complete-case 数据中重复核验。”**

1. `rf_combo − pooled_logistic_combo` 与 `rf_combo − rf_ipTM` 在 `primary_six × target_equal` 的宏平均 Brier、log loss 均严格改善，BSS-transport 均严格改善；
2. 上述两个对比在主六靶点分别至少 5／6 靶点 Brier 严格改善，在五靶点分别至少 4／5 靶点 Brier 严格改善；
3. 在 `{primary_six,primary_five,strict_six} × {target_equal,pair}` 六个配置中，上述两个对比的宏平均 Brier、log loss 和 BSS-transport 均无反转；中性不算反转；
4. 完整逐靶点结果、外测预测、所有预期靶点和全部候选调参记录均通过 QC，无静默删行或删靶点。

该筛查不是统计显著性、实际效应重要性或外部验证优先级门槛。必须同时报告绝对效应量；不得写为“效果显著”“具有实际重要性”“RF/Pinc 更优”。BSS-transport 与 Brier 同时报告，但不得计作两项独立证据。

若未全部满足，按实际结果描述为方向不一致、局部增量或无明确增量；不得移动门槛，也不得把极小但同号的差异包装成实质价值。

### 10.3 与原结论的关系

- RF 无稳定增量：支持“已评估的线性与所选 RF 网格下均未发现稳定 Pinc 增量”，不得推广为“任何复杂模型都无效”。
- RF 仅改善同特征模型家族对比、但 `rf_combo − rf_ipTM` 不稳定：只能描述模型家族条件差异，不能归因 Pinc 增量。
- 达到 10.2：需要收窄现有 RQ5“未发现稳定增量”的措辞，并明确缺少外部 Pinc 验证；不能据此改写 RQ3 的 ipTM 零标签失败结果。
- 无论结果如何，不建立未来靶点总体推断或因果机制。

## 11. 实现、计算恢复、QC 和交付

预计 34 个 final outer source tasks、160 个 inner source-target tasks；去重后共 3,492 个 500-tree RF fits。计算量不得通过减少树数、网格、fold、分析集或权重来压缩。

允许不改变规则的优化：共享只读行集／权重／split、去重有效网格、单次 fit 内固定并行且不做外层进程 oversubscription、及时释放森林、按 outer fold checkpoint 和分批恢复。正式全量前记录代表性 fit 的耗时与峰值内存；若资源不足，保持规则并分批，不得改模型。

执行方必须完成：

1. 保存本冻结协议、执行交接、配置快照、环境、输入／脚本哈希、命令和日志。
2. 单元测试至少覆盖：当前行集权重及靶点权重和、logit、概率列映射、确定性 seed、重复拟合一致性、相同 outer rows、LOTO／inner 隔离、标准化只来自训练行、无 target ID 泄漏、ridge penalty、I 平均约束、F 的 delta=0 重现 U、适配行序与 split 复用。
3. 每个 inner candidate 和 final fit 保存：行集 `pair_id` 集合哈希、各源靶点 `n_t` 与权重总和、权重 SHA-256、有效参数、seed、训练阳性数、`classes_`、`n_in_bag`、状态和验证 log loss。
4. 独立验证不得只调用生产汇总函数：从 inner metrics 重算所有选参，从预测重算所有主要指标和 pairwise 汇总；至少对一个单特征和一个双特征 fold 重算权重、相同 seed 重拟合并逐元素核对预测；检查行数、哈希、fold 隔离、完整样本硬门和图像可读性。
5. 所有失败进入 `failed_item_ledger.csv`；未解决失败不得隐藏。
6. 交付完整 pair predictions、fold／target／macro／pair metrics、pairwise comparisons、全部 inner metrics 和参数、可靠性、图、单元测试、独立验证、`qc_checks.json`、`RUN_STATUS.json`、`HANDOFF_RESULT.md`、`handoff_result.json` 和最终 hash manifest。
7. 结果报告同时列出主分析、五靶点、严格结局和两种权重；报告是否满足 10.2，但不得把该筛查包装成验证性成功。

## 12. 唯一一轮审查的处置记录

审查文件：`analysis/codex_outputs/RQ5RF_protocol_review_2026-11-17/CODEX_REVIEW.md`。

- C-01：已改为每个实际拟合行集重算权重并设硬断言。
- M-01：已统一降格为条件预测增量，删除纯交互归因。
- M-02：已去重单特征网格、共享任务 seed 并冻结 RF 全参数与版本。
- M-03：已冻结 pair_id 排序、inner split 复用、有限性和完整外测行规则。
- M-04：已将 10.2 定义为方向一致性筛查，加入中性带并明确 Brier/BSS 非独立证据。
- M-05：已要求三个分析集独立重拟合及完整样本硬门。
- m-01、m-02：已加入输入范围、概率／裁剪、classes 映射和审计字段要求。

不再进行第二轮方法审查。后续只允许修正不改变科学规则的代码错误；任何规则变化都必须停止运行并重新立项，不得静默纳入本冻结版本。
