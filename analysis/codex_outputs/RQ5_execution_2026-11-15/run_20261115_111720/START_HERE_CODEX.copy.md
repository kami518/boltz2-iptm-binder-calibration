# RQ5 执行交接：Pinc 基础诊断与组合迁移

日期：2026-11-15。工作区：`<REPRODUCTION_ROOT>`。

你是本轮的计算执行方。请直接检查文件、编写脚本、运行、完成自检并交付；不要只给建议。所有新文件只能写入：
- `analysis/codex_outputs/RQ5_execution_2026-11-15/run_YYYYMMDD_HHMMSS/`（唯一正式运行目录）；
- 当前交接目录中的执行日志／配置快照（如确有需要）。

不得修改 v1／v2／v3、冻结主表、序列簇、RQ1–RQ4 结果、Pinc 扩展表、既有 fold 或旧 Codex 运行。不得联网、下载、安装或升级依赖。不得启动 RQ6，不加入随机森林、其他 ML、ipSAE／LIS 或额外特征搜索。

## 1. 已验收输入

1. 主表：`analysis/v3_outputs/pair_master_v3.csv`，2,084 行；SHA-256 `581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91`。
2. Pinc：`analysis/pinc_batch/run_20261110_181500/pinc_pair_extension.csv`，按 `pair_id` 一对一连接；只使用 `pinc_status=success` 的对称化 `pinc` 列。2,074 成功、10 个 `missing_input`，缺失全在 PD-L1；不得插补，不解释 Pinc1／Pinc2。
3. 序列簇：`analysis/rq12_outputs/sequence_clusters_rq12.csv`。
4. RQ3 冻结配置：`analysis/codex_tasks/RQ3_EXECUTION_CODEX_HANDOFF_2026-11-13/rq3_execution_config.json`。
5. 冻结 folds／exclusions：`analysis/codex_outputs/RQ3_prep_2026-11-13/run_20261113_105301/` 中配置指定的 LOTO、target CV5 和 cross-target exclusion 文件。
6. 可复用实现：最终修正版 `analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/` 的 `rq3_common.py`、stage 0–4。复用公式和 folds，但不得直接沿用 2,084 行预测作为本轮模型比较。

开机先核验输入哈希、行数、唯一键和 Pinc 验收状态，输出 `final_input_hashes.csv`。2,074 行共同 complete-case 中 PD-L1 为 267 行、190 个主结局阳性；缺失的 10 行含 8 个主结局阳性、4 个严格阳性。若不一致立即停止。

## 2. 研究定位

RQ5 全为探索性增量分析，不改写原始 ipTM 的 RQ1–RQ4 主结果。回答：
1. 在相同 2,074 行上，Pinc 与 ipTM 的逐靶点区分和原始概率误差如何？
2. Pinc 单指标及 ipTM＋Pinc 是否在相同信息层级、相同 folds、相同外测行上产生折外增量？

不根据结果更换结局、特征、lambda、权重或 folds。

## 3. 冻结分析集与通用口径

- `primary_six`：共同 complete-case，主结局 `primary_outcome`；预期 2,074 行、740 阳性。
- `primary_five`：再排除 AAH-I；预期 1,941 行。
- `strict_six`：共同 complete-case，结局改为 `strict_outcome`；预期 674 阳性（678−4；2026-11-15 预检算术纠错，不改变结局定义）。
- 主权重 target_equal；pair source weighting 为预设敏感性。权重公式完全复用 RQ3。
- 主聚类 sequence_cluster90；复用冻结 fold 和 cross-target exclusions 后，仅过滤掉 10 个 Pinc 缺失行，不重新分折。每个模型比较必须使用相同分析行和相同外测行。
- `x_ipTM=logit(clip(ipTM,1e-6,1-1e-6))`，`x_Pinc=logit(clip(Pinc,1e-6,1-1e-6))`。
- 原始 score 的 Brier／log loss 只用于“把数值误当概率”的诊断；Pinc 原始数值同样不自动是概率。
- 逐靶点及靶点等权宏平均为主，配对加权汇总另列；AAH-I 单独可见，五靶点不可省略。
- 不采用“校准斜率 CI 包含 1”作为成功标准。六靶点不提供未来靶点总体不确定性声明。

## 4. Stage A：同样本 Pinc／ipTM 基础诊断

对 primary_six、primary_five、strict_six，分别评价 ipTM 与 Pinc：

1. 逐靶点 n、阳性数、score 分布、中位数和范围、ROC-AUC、AP。
2. 原始 score 的 Brier、log loss、相对各靶点恒定阳性率的 BSS。
3. 逐靶点未汇聚 logistic 校准截距／斜率及固定 bins 可靠性表；bins 复用 RQ3 `[0,.2,.4,.6,.8,.9,.95,1]`。本阶段不新增层级部分汇聚模型。
4. 每靶点、每结局、每 score 做 2,000 次 sequence_cluster90 整簇 bootstrap，固定种子 20261115。对 ipTM 与 Pinc 使用同一组抽取簇，输出 AUC、AP、Brier 的配对差及 95% percentile 区间；不得用两套独立 bootstrap 拼差。
5. 输出靶点等权与配对加权汇总，明确 AP／Brier 的基率依赖。

Stage A 至少生成：`rq5_target_score_metrics.csv`、`rq5_score_bootstrap_comparisons.csv`、`rq5_score_reliability.csv`、`RQ5_score_discrimination.png`、`RQ5_raw_reliability.png`。

## 5. Stage B：U／I／F 折外迁移与组合增量

三个模型暴露：
- `ipTM_only`：单变量 x_ipTM；作为同样本、同折对照，必须重新拟合。
- `Pinc_only`：单变量 x_Pinc。
- `ipTM_plus_Pinc`：两个特征 x_ipTM、x_Pinc。

### 5.1 单指标 U

完全复用 RQ3 linear U：源靶点截距指示符＋一个共同 score 斜率，未惩罚加权 ML；target_equal 的运输截距为源靶点截距算术平均，pair 敏感性按源靶点行数加权。只用 outer-fold 源训练行及冻结 cross-target exclusions。

### 5.2 组合 U（固定、少量特征的 ridge logistic）

- 设计矩阵为源靶点截距指示符＋x_ipTM、x_Pinc 两个共同斜率；截距不惩罚，只惩罚两个斜率。
- 每个 outer LOTO  fold 内，两个 x 特征的均值和 SD 只由当前源训练行计算；用当前 source-weighting 权重，SD=0 则该 fold 记失败，不得用外测数据。用相同中心／尺度变换留出靶点。
- 目标为加权平均 NLL＋`0.5*lambda*(beta_ipTM²+beta_Pinc²)`；lambda grid 固定 `[0.001,0.01,0.1,1,10,100]`。
- lambda 只在当前 outer 源训练靶点内部用 leave-one-source-target-out 选择。每个 inner fold 重新估计中心、SD、源靶点截距和斜率；inner train／validation 间同序列簇行必须排除。按各 inner 留出靶点 log loss 的等权宏平均选最小；1e-12 内并列取较大 lambda。不得使用 outer 留出靶点标签。
- 最终组合 U 用全部当前源训练行重拟合；运输截距按当前 weighting 规则汇总源靶点截距，斜率迁移不变。
- 保存每折 lambda、inner metrics、特征变换和系数。若组合失败，不得静默改成单指标或无惩罚模型。

### 5.3 I／F 本地适配

对每个模型暴露分别得到 outer LOTO `eta_U` 后，复用冻结 target-CV5 外折及 RQ3 最终修正版适配：
- U：不使用当前靶点标签。
- I：`eta_I=a+eta_U`，只用适配训练折更新截距。
- F：`eta_F=a+gamma*eta_U=a+eta_U+delta*eta_U`；只更新组合总 logit 的截距和一个校准斜率，**不重新估计 ipTM／Pinc 的相对特征权重**。
- F penalty、lambda grid、3 折 StratifiedGroupKFold、inner seeds、log-loss 选择和 failure rules 完全复用 RQ3；只过滤 Pinc 缺失行，不重新分折。
- U、I、F 以及三个暴露在相同 outer-test rows 上比较。外测标签不得参与拟合、lambda 选择或补救。

### 5.4 评价与比较

逐靶点、逐层输出 n／阳性、Brier、log loss、BSS-transport、BSS-target-oracle、BSS-adaptation-train、校准截距／斜率、可靠性和状态；另输出折外 ROC-AUC／AP。

预设比较：
1. Pinc_only − ipTM_only；
2. combo − ipTM_only；
3. combo − Pinc_only；
4. 各暴露 I−U、F−U、F−I。

报告靶点等权宏平均、逐靶点方向、最小／最大差、配对汇总和五靶点／严格结局／pair weighting 敏感性。Brier／log loss 差值负值表示前者更好，AUC／AP／BSS 差值正值表示前者更好。不得只展示有利靶点或有利指标。

Stage B 至少生成完整 pair predictions、fold metrics、target metrics、macro／pair summaries、pairwise comparisons、lambda／parameters、reliability、`RQ5_U_transport.png`、`RQ5_adaptation_increment.png`。

## 6. QC、独立自检与交付

1. 先保存 `run_config.snapshot.json` 和本交接文件副本，再计算结果；记录所有输入／脚本哈希、软件版本、命令和日志。
2. 单元测试至少覆盖：logit 变换、target weights、LOTO 泄漏、inner 标准化只用训练行、combo penalty、delta=0 重现 U、I 平均约束、F 梯度／目标函数、相同 outer-test rows。
3. 所有失败写 `failed_item_ledger.csv`，不得静默改模型。空 reliability bins 保留。
4. 独立 `independent_verification.py` 不得只调用生产汇总函数：至少从冻结输入重算 Stage A 点估计和抽样配对 bootstrap 代表项，重算主要 U 预测／指标及所有 pairwise 汇总，检查 fold 隔离、行数、哈希、相同外测行和图像可读性。
5. 输出 `qc_checks.json`、`RUN_STATUS.json`、`HANDOFF_RESULT.md`、`handoff_result.json`、完整脚本、命令、日志、最终 hash manifest。图必须实际打开检查，不仅确认文件存在。
6. 结果报告只称“探索性增量”；不得声称 Pinc 或组合已经解决通用迁移，也不得用 RQ5 改写原 RQ3 主结论。

若全部执行预计过长，优先完整完成 Stage A 并交付，再执行 Stage B；不得为了缩短时间减少 bootstrap、fold、敏感性或改变规则。完成后在最终消息列出运行目录、主要数字、失败项和自检状态。
