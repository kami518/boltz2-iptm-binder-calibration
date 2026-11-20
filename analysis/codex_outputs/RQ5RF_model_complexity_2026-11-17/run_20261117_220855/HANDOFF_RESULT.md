# RQ5-RF 执行交付

- 运行目录：`<REPRODUCTION_ROOT>\analysis\codex_outputs\RQ5RF_model_complexity_2026-11-17\run_20261117_220855`
- 状态：`complete`
- 证据边界：post hoc exploratory model-complexity sensitivity；不是验证性比较、通用校准或未来靶点推断。
- scikit-learn：`1.9.0`；RF：`3492` 次 500-tree 拟合。

## 主分析 U 层

`primary_six × target_equal`，靶点等权宏平均，未使用 destination-target 标签：

| model | Brier | log loss | BSS-transport | ROC-AUC | AP |
|---|---:|---:|---:|---:|---:|
| pooled_logistic_combo | 0.322258 | 0.856314 | -0.139389 | 0.617534 | 0.572235 |
| pooled_logistic_ipTM | 0.282607 | 0.759039 | -0.000133 | 0.651489 | 0.600158 |
| rf_combo | 0.289445 | 0.786984 | -0.023982 | 0.593172 | 0.557379 |
| rf_ipTM | 0.262443 | 0.725385 | 0.068965 | 0.568358 | 0.553606 |

## 六配置结果

六个配置均完成并保存在 `rq5rf_six_configuration_summary.csv`。两种预设核心对比的 U 层宏平均差值如下（Brier/log loss 为左减右，负值利于左侧；BSS 为正值利于左侧）：

- `primary_six × target_equal`：RF-combo − logistic-combo：Brier -0.0328131，log loss -0.0693298，BSS 0.115407；RF-combo − RF-ipTM：Brier 0.0270021，log loss 0.061599，BSS -0.0929467。
- `primary_six × pair`：RF-combo − logistic-combo：Brier -0.0316154，log loss -0.0540149，BSS 0.104483；RF-combo − RF-ipTM：Brier 0.0286833，log loss 0.0731139，BSS -0.0800347。
- `primary_five × target_equal`：RF-combo − logistic-combo：Brier 0.0466245，log loss 0.172516，BSS -0.129417；RF-combo − RF-ipTM：Brier 0.0622966，log loss 0.201819，BSS -0.173207。
- `primary_five × pair`：RF-combo − logistic-combo：Brier 0.00491417，log loss 0.0518272，BSS -0.0128654；RF-combo − RF-ipTM：Brier 0.0275508，log loss 0.0856496，BSS -0.0666497。
- `strict_six × target_equal`：RF-combo − logistic-combo：Brier -0.0359785，log loss -0.0692363，BSS 0.149603；RF-combo − RF-ipTM：Brier 0.0175213，log loss 0.0423095，BSS -0.0751858。
- `strict_six × pair`：RF-combo − logistic-combo：Brier -0.0092657，log loss 0.00894916，BSS 0.0486844；RF-combo − RF-ipTM：Brier 0.0454303，log loss 0.120459，BSS -0.151498。

## 10.2 方向一致性筛查

结论：**未满足冻结的 10.2 全部条件；只能按实际结果报告方向不一致、局部增量或无明确增量。**

| condition | passed | description |
|---:|:---:|---|
| 1 | FAIL | primary_six target_equal: both comparisons strictly improve Brier, log loss, BSS-transport |
| 2 | FAIL | primary six at least 5/6 and five-target at least 4/5 Brier-improved targets for both comparisons |
| 3 | FAIL | all six configurations show no reversal in the three required macro metrics |
| 4 | PASS | all predictions, expected targets, candidate tuning records, and QC gates pass |

该筛查不是显著性检验或效应重要性门槛；即使满足，也只能表述为需要在独立 complete-case 数据中重复核验。绝对效应与逐靶点方向见 pairwise 文件。

## 验证与失败项

- 独立验证：`PASS`，包含选参重算、指标重算、单/双特征 final RF 同 seed 复跑和逐元素预测比对。
- 未解决失败项：`0`；失败台帐：`<REPRODUCTION_ROOT>\analysis\codex_outputs\RQ5RF_model_complexity_2026-11-17\run_20261117_220855\failed_item_ledger.csv`。
- 已记录但不可评价的诊断项：`60`（保留在 target 表，不纳入预测指标宏平均）。
- 单元测试：`16` 项，successful=`True`。
- AAH-I、五靶点、严格结局与两种权重均保留在 target/macro 表中，未删除不可评价项或不可见靶点。

## 解释限制

- 不识别纯交互；RF-combo − RF-ipTM 可来自 ipTM 非线性、Pinc 非线性或二者交互。
- 不声称 RF 概率已校准，不声称 RF/Pinc 普遍更优，不声称解决跨靶点迁移。
- 缺少独立外部 Pinc 队列；本结果不建立未来靶点总体推断。
