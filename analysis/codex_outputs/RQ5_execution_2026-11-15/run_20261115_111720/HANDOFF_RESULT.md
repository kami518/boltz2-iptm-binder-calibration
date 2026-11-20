# RQ5 探索性增量结果

运行目录：`<REPRODUCTION_ROOT>\analysis\codex_outputs\RQ5_execution_2026-11-15\run_20261115_111720`

状态：**完成**。所有结果仅定位为探索性增量分析；不声称 Pinc 或组合模型已经解决通用迁移，也不改写 RQ3 主结论。

## 输入与分析集

- 主表 2,084 行，冻结哈希通过；Pinc 成功 2,074 行，10 个缺失全部在 PD-L1，无插补。
- `primary_six`：2,074 行，主结局 740 阳性。
- `primary_five`：排除 AAH-I 后 1,941 行；AAH-I 在图中始终单独可见。
- `strict_six`：2,074 行，严格结局 674 阳性。
- 主权重 `target_equal`，主聚类 `sequence_cluster90`；pair weighting 为预设敏感性。

## Stage A：同一样本原始分数诊断

以下是 `primary_six` 六靶点等权宏平均。原始 score 的 Brier/log loss 只表示“把数值直接当作概率”的诊断误差；它们不说明对应 score 已校准。

| score | metric | mean | min | max |
| --- | --- | --- | --- | --- |
| Pinc | roc_auc | 0.6366 | 0.5403 | 0.8500 |
| Pinc | ap | 0.5533 | 0.1315 | 0.9186 |
| Pinc | brier | 0.3834 | 0.2425 | 0.5861 |
| Pinc | log_loss | 1.0750 | 0.6980 | 1.5889 |
| Pinc | bss_target_constant | -1.7332 | -4.4077 | -0.0380 |
| ipTM | roc_auc | 0.6770 | 0.5766 | 0.9501 |
| ipTM | ap | 0.6132 | 0.1245 | 0.9168 |
| ipTM | brier | 0.4154 | 0.2400 | 0.6299 |
| ipTM | log_loss | 1.2224 | 0.6650 | 1.8276 |
| ipTM | bss_target_constant | -2.0569 | -6.1097 | -0.0404 |

Pinc 的宏平均 ROC-AUC/AP 低于 ipTM，但其原始数值的 Brier/log loss 更低。逐靶点方向并不一致：例如 Pinc 的 Brier 在 AAH-I、nipah、PD-L1 和 TREM2 更低，在 human-VEGF110 与 SARS-CoV-2 RBD 不低于 ipTM。代表性 bootstrap 使用同一组 sequence-cluster90 整簇抽取。

## Stage B：U 运输

`primary_six`、`target_equal` 的 LOTO 宏平均：

| model_exposure | metric | mean | min | max |
| --- | --- | --- | --- | --- |
| Pinc_only | roc_auc | 0.6366 | 0.5403 | 0.8500 |
| Pinc_only | ap | 0.5533 | 0.1315 | 0.9186 |
| Pinc_only | brier | 0.3144 | 0.2251 | 0.5029 |
| Pinc_only | log_loss | 0.8695 | 0.6428 | 1.4555 |
| Pinc_only | bss_transport | -0.1038 | -0.4989 | 0.2260 |
| Pinc_only | bss_target_oracle | -1.1113 | -2.8056 | -0.0968 |
| ipTM_only | roc_auc | 0.6770 | 0.5766 | 0.9501 |
| ipTM_only | ap | 0.6132 | 0.1245 | 0.9168 |
| ipTM_only | brier | 0.3139 | 0.2132 | 0.5094 |
| ipTM_only | log_loss | 0.9204 | 0.6178 | 1.7749 |
| ipTM_only | bss_transport | -0.0991 | -0.5183 | 0.2670 |
| ipTM_only | bss_target_oracle | -1.1561 | -2.8548 | -0.0298 |
| ipTM_plus_Pinc | roc_auc | 0.6662 | 0.5659 | 0.9331 |
| ipTM_plus_Pinc | ap | 0.5903 | 0.1271 | 0.9173 |
| ipTM_plus_Pinc | brier | 0.2949 | 0.2441 | 0.3933 |
| ipTM_plus_Pinc | log_loss | 0.7859 | 0.6812 | 0.9944 |
| ipTM_plus_Pinc | bss_transport | -0.0381 | -0.1722 | 0.0338 |
| ipTM_plus_Pinc | bss_target_oracle | -0.9988 | -2.2142 | -0.0031 |

预设 U 对比（左减右；Brier/log loss 负值更好，AUC/AP/BSS 正值更好）：

| comparison | metric | macro_difference | min_target_difference | max_target_difference | targets_left_better | target_count |
| --- | --- | --- | --- | --- | --- | --- |
| Pinc_only_minus_ipTM_only | brier | 0.0005 | -0.0532 | 0.0279 | 2.0000 | 6.0000 |
| Pinc_only_minus_ipTM_only | log_loss | -0.0509 | -0.3194 | 0.0569 | 2.0000 | 6.0000 |
| Pinc_only_minus_ipTM_only | roc_auc | -0.0404 | -0.1002 | 0.0170 | 2.0000 | 6.0000 |
| Pinc_only_minus_ipTM_only | ap | -0.0598 | -0.1832 | 0.0070 | 2.0000 | 6.0000 |
| Pinc_only_minus_ipTM_only | bss_transport | -0.0047 | -0.1019 | 0.1869 | 2.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_ipTM_only | brier | -0.0189 | -0.1161 | 0.0963 | 4.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_ipTM_only | log_loss | -0.1344 | -0.7805 | 0.1966 | 4.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_ipTM_only | roc_auc | -0.0108 | -0.0282 | 0.0107 | 2.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_ipTM_only | ap | -0.0229 | -0.0559 | 0.0025 | 2.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_ipTM_only | bss_transport | 0.0611 | -0.3312 | 0.3460 | 4.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_Pinc_only | brier | -0.0194 | -0.1096 | 0.0844 | 4.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_Pinc_only | log_loss | -0.0835 | -0.4611 | 0.1715 | 4.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_Pinc_only | roc_auc | 0.0297 | -0.0063 | 0.0831 | 4.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_Pinc_only | ap | 0.0369 | -0.0044 | 0.1273 | 4.0000 | 6.0000 |
| ipTM_plus_Pinc_minus_Pinc_only | bss_transport | 0.0658 | -0.2902 | 0.3267 | 4.0000 | 6.0000 |

组合相对 ipTM 的宏平均概率误差更低，BSS-transport 更高，但 ROC-AUC 和 AP 的宏平均下降；逐靶点方向也混合。因此这里支持的是“探索性概率误差与运输基线改善”，不是通用判别力改善。

## Stage B：I/F 适配增量

所有比较均在相同 target-CV5 外测行上完成，U 不更新当前靶点标签；I 只更新截距；F 更新截距和一个校准斜率。

| exposure | comparison | metric | macro_difference | min_target_difference | max_target_difference | targets_left_better |
| --- | --- | --- | --- | --- | --- | --- |
| ipTM_only | predictive_I_minus_U_transport_zero_label | brier | -0.1453 | -0.3620 | -0.0115 | 6.0000 |
| ipTM_only | predictive_I_minus_U_transport_zero_label | log_loss | -0.4080 | -1.2829 | -0.0251 | 6.0000 |
| ipTM_only | predictive_I_minus_U_transport_zero_label | roc_auc | -0.0009 | -0.0068 | 0.0028 | 2.0000 |
| ipTM_only | predictive_I_minus_U_transport_zero_label | ap | 0.0063 | -0.0101 | 0.0502 | 2.0000 |
| ipTM_only | predictive_F_minus_U_transport_zero_label | brier | -0.1526 | -0.3782 | -0.0088 | 6.0000 |
| ipTM_only | predictive_F_minus_U_transport_zero_label | log_loss | -0.4382 | -1.3479 | -0.0138 | 6.0000 |
| ipTM_only | predictive_F_minus_U_transport_zero_label | roc_auc | -0.0176 | -0.0530 | 0.0004 | 1.0000 |
| ipTM_only | predictive_F_minus_U_transport_zero_label | ap | -0.0218 | -0.0962 | 0.0034 | 1.0000 |
| ipTM_only | predictive_F_minus_predictive_I | brier | -0.0073 | -0.0422 | 0.0135 | 3.0000 |
| ipTM_only | predictive_F_minus_predictive_I | log_loss | -0.0302 | -0.1516 | 0.0320 | 2.0000 |
| ipTM_only | predictive_F_minus_predictive_I | roc_auc | -0.0167 | -0.0462 | 0.0004 | 1.0000 |
| ipTM_only | predictive_F_minus_predictive_I | ap | -0.0281 | -0.0860 | 0.0036 | 1.0000 |
| Pinc_only | predictive_I_minus_U_transport_zero_label | brier | -0.1442 | -0.3696 | -0.0221 | 6.0000 |
| Pinc_only | predictive_I_minus_U_transport_zero_label | log_loss | -0.3571 | -1.0233 | -0.0467 | 6.0000 |
| Pinc_only | predictive_I_minus_U_transport_zero_label | roc_auc | -0.0062 | -0.0192 | 0.0003 | 1.0000 |
| Pinc_only | predictive_I_minus_U_transport_zero_label | ap | -0.0014 | -0.0155 | 0.0103 | 3.0000 |
| Pinc_only | predictive_F_minus_U_transport_zero_label | brier | -0.1461 | -0.3711 | -0.0220 | 6.0000 |
| Pinc_only | predictive_F_minus_U_transport_zero_label | log_loss | -0.3676 | -1.0275 | -0.0465 | 6.0000 |
| Pinc_only | predictive_F_minus_U_transport_zero_label | roc_auc | -0.0150 | -0.0478 | -0.0024 | 0.0000 |
| Pinc_only | predictive_F_minus_U_transport_zero_label | ap | -0.0156 | -0.0319 | -0.0001 | 0.0000 |
| Pinc_only | predictive_F_minus_predictive_I | brier | -0.0019 | -0.0140 | 0.0031 | 2.0000 |
| Pinc_only | predictive_F_minus_predictive_I | log_loss | -0.0105 | -0.0681 | 0.0073 | 2.0000 |
| Pinc_only | predictive_F_minus_predictive_I | roc_auc | -0.0087 | -0.0481 | 0.0162 | 1.0000 |
| Pinc_only | predictive_F_minus_predictive_I | ap | -0.0142 | -0.0306 | 0.0003 | 1.0000 |
| ipTM_plus_Pinc | predictive_I_minus_U_transport_zero_label | brier | -0.1229 | -0.2631 | -0.0008 | 6.0000 |
| ipTM_plus_Pinc | predictive_I_minus_U_transport_zero_label | log_loss | -0.2666 | -0.5680 | -0.0015 | 6.0000 |
| ipTM_plus_Pinc | predictive_I_minus_U_transport_zero_label | roc_auc | -0.1051 | -0.3601 | -0.0051 | 0.0000 |
| ipTM_plus_Pinc | predictive_I_minus_U_transport_zero_label | ap | -0.0775 | -0.2991 | 0.0012 | 1.0000 |
| ipTM_plus_Pinc | predictive_F_minus_U_transport_zero_label | brier | -0.1229 | -0.2631 | -0.0008 | 6.0000 |
| ipTM_plus_Pinc | predictive_F_minus_U_transport_zero_label | log_loss | -0.2670 | -0.5696 | -0.0016 | 6.0000 |
| ipTM_plus_Pinc | predictive_F_minus_U_transport_zero_label | roc_auc | -0.1098 | -0.3523 | -0.0371 | 0.0000 |
| ipTM_plus_Pinc | predictive_F_minus_U_transport_zero_label | ap | -0.0809 | -0.2894 | -0.0065 | 0.0000 |
| ipTM_plus_Pinc | predictive_F_minus_predictive_I | brier | -0.0000 | -0.0001 | 0.0001 | 5.0000 |
| ipTM_plus_Pinc | predictive_F_minus_predictive_I | log_loss | -0.0003 | -0.0015 | -0.0000 | 6.0000 |
| ipTM_plus_Pinc | predictive_F_minus_predictive_I | roc_auc | -0.0047 | -0.0320 | 0.0078 | 2.0000 |
| ipTM_plus_Pinc | predictive_F_minus_predictive_I | ap | -0.0034 | -0.0172 | 0.0097 | 1.0000 |

主分析中，I/F 相对 U 的 Brier 和 log loss 在六靶点宏平均上改善；AUC/AP 方向往往略有下降。`bss_adaptation_train` 使用适配训练折的靶点内基线，不是零标签指标。F 相对 I 的额外增量较小且逐靶点混合。

## 敏感性

`primary_five`、`strict_six` 和 pair weighting 均已完成并保留在 `rq5_pairwise_comparisons.csv`、`rq5_sensitivity_summary.csv`。完整敏感性结果不裁剪为有利方向；以下另列 Brier/log loss 的宏平均差值：

| analysis_set | source_weighting | comparison_family | comparison | left_model | metric | mean_difference | positive_difference_targets | negative_difference_targets | target_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary_five | pair | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | brier | -0.0095 | 2.0000 | 3.0000 | 5.0000 |
| primary_five | pair | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | log_loss | -0.0408 | 2.0000 | 3.0000 | 5.0000 |
| primary_five | pair | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | brier | 0.0200 | 3.0000 | 2.0000 | 5.0000 |
| primary_five | pair | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | log_loss | 0.0254 | 3.0000 | 2.0000 | 5.0000 |
| primary_five | pair | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | brier | 0.0105 | 2.0000 | 3.0000 | 5.0000 |
| primary_five | pair | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | log_loss | -0.0154 | 2.0000 | 3.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | brier | -0.2015 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | log_loss | -0.4804 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | brier | -0.2120 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | log_loss | -0.5232 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.2219 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.5065 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | brier | 0.0010 | 5.0000 | 0.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | log_loss | 0.0036 | 5.0000 | 0.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | brier | 0.0020 | 5.0000 | 0.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | log_loss | 0.0052 | 5.0000 | 0.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | brier | 0.0000 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | log_loss | 0.0003 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | brier | -0.2025 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | log_loss | -0.4839 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | brier | -0.2140 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | log_loss | -0.5284 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.2219 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.5068 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | brier | -0.0094 | 3.0000 | 2.0000 | 5.0000 |
| primary_five | target_equal | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | log_loss | -0.0478 | 3.0000 | 2.0000 | 5.0000 |
| primary_five | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | brier | -0.0033 | 2.0000 | 3.0000 | 5.0000 |
| primary_five | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | log_loss | -0.0141 | 2.0000 | 3.0000 | 5.0000 |
| primary_five | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | brier | -0.0126 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | log_loss | -0.0619 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | brier | -0.1324 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | log_loss | -0.2969 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | brier | -0.1431 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | log_loss | -0.3477 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.1286 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.2816 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | brier | 0.0010 | 5.0000 | 0.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | log_loss | 0.0034 | 5.0000 | 0.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | brier | 0.0007 | 4.0000 | 1.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | log_loss | 0.0012 | 4.0000 | 1.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | brier | 0.0001 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | log_loss | 0.0007 | 1.0000 | 4.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | brier | -0.1333 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | log_loss | -0.3003 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | brier | -0.1438 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | log_loss | -0.3488 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.1287 | 0.0000 | 5.0000 | 5.0000 |
| primary_five | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.2822 | 0.0000 | 5.0000 | 5.0000 |
| primary_six | pair | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | brier | -0.0043 | 3.0000 | 3.0000 | 6.0000 |
| primary_six | pair | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | log_loss | -0.0337 | 3.0000 | 3.0000 | 6.0000 |
| primary_six | pair | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | brier | 0.0127 | 3.0000 | 3.0000 | 6.0000 |
| primary_six | pair | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | log_loss | 0.0030 | 3.0000 | 3.0000 | 6.0000 |
| primary_six | pair | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | brier | 0.0084 | 3.0000 | 3.0000 | 6.0000 |
| primary_six | pair | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | log_loss | -0.0307 | 3.0000 | 3.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | brier | -0.1686 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | log_loss | -0.4266 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | brier | -0.1802 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | log_loss | -0.4811 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.1793 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.4187 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | brier | -0.0013 | 5.0000 | 1.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | log_loss | -0.0081 | 5.0000 | 1.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | brier | -0.0062 | 3.0000 | 3.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | log_loss | -0.0236 | 4.0000 | 2.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | brier | -0.0010 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | log_loss | -0.0034 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | brier | -0.1673 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | log_loss | -0.4185 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | brier | -0.1740 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | log_loss | -0.4575 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.1783 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.4153 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | brier | 0.0005 | 4.0000 | 2.0000 | 6.0000 |
| primary_six | target_equal | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | log_loss | -0.0509 | 4.0000 | 2.0000 | 6.0000 |
| primary_six | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | brier | -0.0194 | 2.0000 | 4.0000 | 6.0000 |
| primary_six | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | log_loss | -0.0835 | 2.0000 | 4.0000 | 6.0000 |
| primary_six | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | brier | -0.0189 | 2.0000 | 4.0000 | 6.0000 |
| primary_six | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | log_loss | -0.1344 | 2.0000 | 4.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | brier | -0.1461 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | log_loss | -0.3676 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | brier | -0.1526 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | log_loss | -0.4382 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.1229 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.2670 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | brier | -0.0019 | 4.0000 | 2.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | log_loss | -0.0105 | 4.0000 | 2.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | brier | -0.0073 | 3.0000 | 3.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | log_loss | -0.0302 | 4.0000 | 2.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | brier | -0.0000 | 1.0000 | 5.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | log_loss | -0.0003 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | brier | -0.1442 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | log_loss | -0.3571 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | brier | -0.1453 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | log_loss | -0.4080 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.1229 | 0.0000 | 6.0000 | 6.0000 |
| primary_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.2666 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | pair | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | brier | -0.0030 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | log_loss | -0.0264 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | brier | 0.0080 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | log_loss | -0.0043 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | brier | 0.0050 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | log_loss | -0.0307 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | brier | -0.1334 | 1.0000 | 5.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | log_loss | -0.3518 | 1.0000 | 5.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | brier | -0.1448 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | log_loss | -0.4009 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.1390 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.3359 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | brier | -0.0022 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | log_loss | -0.0110 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | brier | -0.0090 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | log_loss | -0.0303 | 4.0000 | 2.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | brier | -0.0012 | 1.0000 | 5.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | log_loss | -0.0044 | 1.0000 | 5.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | brier | -0.1312 | 1.0000 | 5.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | log_loss | -0.3409 | 1.0000 | 5.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | brier | -0.1358 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | log_loss | -0.3706 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.1378 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | pair | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.3314 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | brier | 0.0017 | 4.0000 | 2.0000 | 6.0000 |
| strict_six | target_equal | preset_U_exposure | Pinc_only_minus_ipTM_only | Pinc_only | log_loss | -0.0463 | 4.0000 | 2.0000 | 6.0000 |
| strict_six | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | brier | -0.0182 | 2.0000 | 4.0000 | 6.0000 |
| strict_six | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_Pinc_only | ipTM_plus_Pinc | log_loss | -0.0781 | 2.0000 | 4.0000 | 6.0000 |
| strict_six | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | brier | -0.0165 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | target_equal | preset_U_exposure | ipTM_plus_Pinc_minus_ipTM_only | ipTM_plus_Pinc | log_loss | -0.1244 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | brier | -0.1156 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | Pinc_only | log_loss | -0.3017 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | brier | -0.1230 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_only | log_loss | -0.3728 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.0935 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.2066 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | brier | -0.0029 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | Pinc_only | log_loss | -0.0130 | 3.0000 | 3.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | brier | -0.0123 | 4.0000 | 2.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_only | log_loss | -0.0439 | 4.0000 | 2.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | brier | -0.0003 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_F_minus_predictive_I | ipTM_plus_Pinc | log_loss | -0.0010 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | brier | -0.1127 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | Pinc_only | log_loss | -0.2887 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | brier | -0.1107 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_only | log_loss | -0.3289 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | brier | -0.0931 | 0.0000 | 6.0000 | 6.0000 |
| strict_six | target_equal | preset_adaptation_increment | predictive_I_minus_U_transport_zero_label | ipTM_plus_Pinc | log_loss | -0.2055 | 0.0000 | 6.0000 | 6.0000 |

## QC 与失败项

- 输入哈希、行数、唯一键、Pinc 验收状态：通过。
- 15/15 单元测试：通过。
- 独立重算 18/18：通过，包括全部 Stage A 点估计、代表性 2,000 次配对 bootstrap、主 U 预测/指标，以及全部 3,672 条靶点差值、648 条 macro 和 648 条 pair-pooled pairwise 汇总。
- U 与适配预测各 36,534 行；相同外测行检查全部通过。
- 四张图均已实际打开检查。
- 失败台账 12 行，均为已解决的开发期数值/中间状态失败；最终 unresolved 失败为 0。

## 解释边界

六靶点结果不提供未来靶点总体不确定性声明。原始 ipTM/Pinc 数值不自动是概率，校准截距/斜率只作描述。当前结果不改变原始 ipTM 的 RQ1–RQ4 结论。
