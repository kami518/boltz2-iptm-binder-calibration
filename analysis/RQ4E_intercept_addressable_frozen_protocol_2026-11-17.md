# RQ4E：截距可修正误差的事后诊断（冻结协议）

日期：2026-11-17  
状态：执行前冻结；仅派生既有冻结结果，不重新拟合模型。

## 1. 目的与定位

RQ4E 服务于 RQ4 的 pilot 逻辑，量化在既有留一靶点外测中：

1. 本地截距更新（I）相对零标签迁移（U）能减少多少预测误差；
2. 在截距更新之后进一步更新斜率（F）还有多少增量收益；
3. 截距更新后仍剩余多少误差。

这是对冻结 RQ3 结果的事后诊断，不是因果贡献分解。它不能证明基率、靶点性质、测定方法或任何其他机制是迁移误差的“最大原因”，也不能单独证明 pilot 实验值得推荐。

## 2. 冻结输入

唯一目标级取数来源：

`analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/rq3_target_adaptation_metrics.csv`

跨表一致性核对来源：

- `rq3_loto_metrics.csv`：独立 U 层结果；
- `rq3_adaptation_comparisons.csv`：冻结的逐对比宏平均与改善／变差靶点数；
- `rq3_adaptation_macro_summary.csv`：冻结的 U／I／F 宏平均；
- `qc_checks.json` 与 `final_hash_gate_result.json`：原分析 QC 与输入哈希门状态。

所有输入只读；运行前后记录并核对 SHA-256。

## 3. 分析范围

- 主分析：`primary_six + primary_outcome + sequence_cluster90 + target_equal`。
- 敏感性分析：`primary_five` 的完整重训结果，使用相同 outcome、cluster rule 和 source weighting；不能通过从六靶点结果删除 AAH-I 获得。
- 模型层：`U_transport_zero_label`（U）、`predictive_I`（I）、`predictive_F`（F）。
- 每个分析集内，U／I／F 必须对应相同靶点及相同外测行；以 `n`、`positives` 和 `prevalence` 完全一致为核对条件。
- 不混入 RQ3E 外部队列，不加入 Pinc，不生成新外折，不拟合新模型，不开展 pilot 模拟。

## 4. 指标与计算

Brier score 为主，log loss 为辅。逐靶点计算：

- `delta_brier_I_minus_U = brier_I - brier_U`；负值表示改善；
- `delta_brier_F_minus_I = brier_F - brier_I`；负值表示改善；
- `delta_log_loss_I_minus_U = log_loss_I - log_loss_U`；
- `delta_log_loss_F_minus_I = log_loss_F - log_loss_I`；
- `relative_brier_reduction_I_vs_U = (brier_U - brier_I) / brier_U`，仅命名为“相对 Brier 误差下降”，不是“基率解释比例”；
- `relative_brier_reduction_F_vs_I = (brier_I - brier_F) / brier_I`，仅表示斜率更新的相对增量。

同时保留 U／I／F 的 Brier、log loss、`bss_transport`、`bss_target_oracle`、校准截距与校准斜率。I 层的 Brier 与 BSS 只能描述为“截距更新后的剩余误差／相对参考模型的表现”，不得据此作机制归因。

汇总采用靶点等权宏平均，逐靶点值与宏平均同时报告。宏平均的相对误差下降由宏平均 Brier 计算，同时保留逐靶点相对下降，二者不得混称。

## 5. 推断边界

- 不新增置信区间、bootstrap、显著性检验或未来靶点总体推断。
- 不把改善／变差靶点计数当作随机符号检验。
- 不根据六靶点与五靶点结果选择更有利的口径。
- RQ4E 可支持的最强结论是：在当前队列和已有适配设计中，截距更新是否具有一致的可修正误差成分，以及斜率更新是否提供稳定增量。
- 是否将 pilot 纳入主文及其措辞强度，必须在查看 RQ4E 结果后另行决定。

## 6. 输出与验收

输出到新的 `analysis/rq4e_outputs/run_<实际时间>/`：

- `rq4e_target_diagnostic.csv`；
- `rq4e_macro_summary.csv`；
- `RQ4E_RESULTS.md`；
- `rq4e_qc_results.json`；
- `input_hashes.csv`；
- `run_config.snapshot.json`。

必须通过以下 QC：

1. 六靶点 18 行（6×3 层），五靶点 15 行；靶点集合唯一且完整；
2. 六靶点每层合计 2,084／748，五靶点每层合计 1,951／729；五靶点不含 AAH-I；
3. 每个靶点 U／I／F 的 `n`、`positives`、`prevalence` 完全一致；
4. 所有纳入结果可评估，核心指标数值有限；
5. 所有派生差值算术一致；
6. 派生宏平均、改善／变差计数与 `rq3_adaptation_comparisons.csv` 在绝对容差 1e-12 内一致；
7. U／I／F 的 Brier 与 log loss 宏平均与 `rq3_adaptation_macro_summary.csv` 在绝对容差 1e-12 内一致；
8. 主六靶点 U 与 `rq3_loto_metrics.csv` 的核心指标一致；样本量、结局数、Brier、log loss 与 BSS 使用绝对容差 1e-12。首次执行发现 VEGF110 的重复数值优化仅使校准截距／斜率产生约 1e-8–1e-9 差异，因此仅这两个诊断参数允许 1e-7 容差，并保留失败运行作为审计记录；
9. 原 RQ3 QC 为全部通过，原输入哈希门为通过；
10. 所有来源文件运行前后 SHA-256 不变。
