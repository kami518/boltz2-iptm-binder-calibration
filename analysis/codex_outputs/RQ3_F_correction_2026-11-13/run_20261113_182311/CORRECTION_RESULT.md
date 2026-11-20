# RQ3 F 层实现修正结果

运行目录：`<REPRODUCTION_ROOT>\analysis\codex_outputs\RQ3_F_correction_2026-11-13\run_20261113_182311`

状态：**Codex 执行完成，仍需主执行线独立验收。**

本结果不替代主执行线的独立 Newton/IRLS 复核，也不自行宣称正式科学验收通过。

## 1. 根因

旧运行 `fit_offset_slope()` 使用设计矩阵 `[1, eta_U]`，实际优化的是 `eta_F = a + theta_1 * eta_U`，却把 `theta_1` 保存为 `delta`，并按 `gamma = 1 + theta_1` 预测。因此拟合目标与保存/预测使用了两条不同的线性预测子，二者恰好相差一个 `eta_U`。旧单测只检查 lambda 变化时的参数趋势，没有从保存参数重构训练目标，因而漏检。

冻结定义未改变：`eta_F = a + eta_U + delta * eta_U = a + gamma * eta_U`，`gamma = 1 + delta`；目标为未裁剪 `mean_NLL(eta_F, y) + 0.5*lambda*delta^2`，截距不惩罚。

## 2. 修改位置

- `rq3_common.py`：通用 `fit_weighted_logistic()` 增加显式 `offset`；目标、梯度、Hessian 与收敛后概率均使用 `X @ theta + offset`。F 专用 `fit_offset_slope()` 固定 `offset=eta_U`、设计 `[1, eta_U]`、惩罚索引 1，保存 `delta=theta[1]` 与 `gamma=1+delta`。
- `stage_2_oracle_descriptive.py`：参数表新增描述性 F 解析梯度最大值。
- `stage_3_predictive_I_F.py`：扩展合成单测，覆盖保存目标重构、恒等关系、解析梯度、中心有限差分、lambda=0 独立拟合与大 lambda 极限；原 lambda 网格、`1e-12` 并列取更大 lambda 规则和概率裁剪均未改变。
- `rq3_f_regression.py`：新增全部 88 个描述性 F 与 440 个最终外折 F 的逐记录目标—预测一致性检查，并生成新旧 U/I/F 数值对照。
- `stage_4_finalize.py`：把修正回归检查纳入最终 QC，并生成本报告。

## 3. 修正后回归测试

逐记录检查覆盖 528 个 F 拟合（88 个描述性 F、440 个最终外折 F）。容差在测试 JSON 与 QC JSON 中明确记录，优化器停止条件没有替代目标一致性验证。

| check | status | maximum | tolerance |
| --- | --- | --- | --- |
| target-to-prediction objective reconstruction | passed | 2.220446e-16 | 1.0e-10 |
| gamma = 1 + delta for all F records | passed | 1.776357e-15 | 1.0e-14 |
| eta_F = a + gamma*eta_U | passed | 3.552714e-15 | 1.0e-12 |
| p_F = expit(eta_F) | passed | 4.440892e-16 | 1.0e-12 |
| saved analytical gradients | passed | 1.649816e-07 | 1.0e-05 |
| lambda=0 independent gamma fit | passed | 2.438230e-08 | 1.0e-07 |
| central finite-difference gradient | passed | 3.268789e-11 | 1.0e-06 |
| large-lambda delta compression | passed | 9.994923e-10 | 1.0e-05 |

全部 88 个描述性 F 的 lambda=0 直接 `y ~ intercept + gamma*eta_U` 独立拟合均成功；有限大 lambda 单测仅要求 `delta` 被压缩，不把“严格等于 I”作为要求。

## 4. 旧结果失效范围

旧运行的 F 层拟合参数、目标、lambda 候选与选择、F 预测、F 指标、比较及含 F 的图表/报告均不能用于验收。不失效的范围是旧 U、oracle I 与预测性 I；新运行仍从头重算它们并做数值对照。

| old F artifact | count |
| --- | --- |
| descriptive F fits | 88 |
| descriptive F pair predictions | 32280 |
| predictive F final outer fits | 440 |
| predictive F pair predictions | 32280 |
| lambda candidates | 2640 |

旧 `RQ3_adaptation.png` 和旧 HANDOFF 中包含 F 的面板或结论同时失效；U 层图不因本次根因而失效，但仍以新运行文件为交付。

## 5. 新旧 U／I 数值一致性

| layer | new_rows | old_rows | joined_rows | unmatched_rows | max_absolute_eta_difference | max_absolute_probability_difference | material_change_rows |
| --- | --- | --- | --- | --- | --- | --- | --- |
| U | 64560 | 64560 | 64560 | 0 | 0.0000 | 0.0000 | 0 |
| oracle_I | 32280 | 32280 | 32280 | 0 | 0.0000 | 0.0000 | 0 |
| predictive_I | 32280 | 32280 | 32280 | 0 | 0.0000 | 0.0000 | 0 |
| descriptive_F | 32280 | 32280 | 32280 | 0 | 4.4686 | 0.6496 | 32280 |
| predictive_F | 32280 | 32280 | 32280 | 0 | 3.7666 | 0.6430 | 32280 |

`material_change_rows` 阈值：绝对 eta 差 `>1e-10` 或绝对概率差 `>1e-12`。U、oracle I、预测性 I 若存在实质变化会同时写入 `rq3_old_new_material_changes.csv` 并触发最终 QC 失败；本表直接报告实测最大差。

## 6. 修正后关键结果

- 主 U 口径：Brier `0.308627`，BSS-transport `-0.080917`，BSS-target-oracle `-1.120489`。
- 预测性 I 相对 U：Brier 宏平均差 `-0.139279`，6/6 个靶点改善；BSS-adaptation-train `0.015318`。
- 预测性 F 相对 I：Brier 宏平均差 `-0.007625`，4/6 个靶点改善、2/6 个变差；BSS-adaptation-train `0.080209`。
- 2,640 个 lambda 候选全部重新计算；440 个最终外折各选一个 lambda。选择频数：`{"0.001": 145, "0.01": 106, "0.1": 25, "1.0": 6, "10.0": 2, "100.0": 156}`。
- 上述为六个留出靶点的描述性执行结果，不构成对未来靶点的总体不确定性声明。

## 7. QC 与失败概况

- 冻结配置 SHA-256：`a416ecd357eabb4c4d95701136c412cf23651b730bf39feab7524fbc07eeef51`。
- 冻结方案 SHA-256：`36dc1503e86960b5f05d0eb273d3b68091eb2072975fb51e5556ee821e73cfd8`。
- 最终 QC：40/40 通过。
- 失败台账：28 条，其中校准诊断边界 `28` 条，其他模型/回归失败 `0` 条。
- 所有校准边界缺失均按原规则记账，没有替换方法、删除折、放宽裁剪或使用外测标签修补。

**Codex 执行完成，仍需主执行线独立验收。**
