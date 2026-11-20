# RQ3 F 层修正后的正式执行交接结果

运行目录：`<REPRODUCTION_ROOT>\analysis\codex_outputs\RQ3_F_correction_2026-11-13\run_20261113_182311`

总体状态：**Codex 执行完成，仍需主执行线独立验收。**

本运行针对已确认的 F 层实现错误从冻结输入完整重跑 stage 0–4；根因、修改、回归测试和旧结果失效范围见 `CORRECTION_RESULT.md`。

本运行严格复用已验收的 44 个 LOTO 折与 220 个 CV5 折；没有重新生成外折，没有修改任何冻结输入或旧结果，没有安装或升级依赖，也没有启动 RQ4/RQ5/RQ6。所有新文件均位于本运行目录。

## 1. 主分析：U 层零标签迁移

唯一主口径：`primary_six + sequence_cluster90 + target_equal + linear`。U 层是零标签主问题；oracle I、描述性 F、预测性 I/F 均不替代该结论。

| held_out_target | n | positives | brier | log_loss | bss_transport | bss_target_oracle | calibration_intercept | calibration_slope |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aah-i | 133 | 19 | 0.2655 | 0.7244 | 0.0327 | -1.1685 | -7.8767 | 16.0453 |
| human-vegf110 | 268 | 158 | 0.3226 | 0.8567 | -0.2248 | -0.3330 | 1.2168 | 0.9391 |
| nipah-glycoprotein-g | 1028 | 101 | 0.3001 | 0.7948 | -0.0521 | -2.3878 | -2.4565 | 0.8542 |
| pd-l1 | 277 | 198 | 0.2263 | 0.6525 | 0.2230 | -0.1103 | 0.8377 | 0.2418 |
| sars-cov-2-rbd | 268 | 226 | 0.4899 | 1.4699 | -0.4626 | -2.7070 | 2.3340 | 0.3946 |
| trem2 | 110 | 46 | 0.2473 | 0.6876 | -0.0017 | -0.0163 | -1.0364 | 3.3008 |

主分析宏平均与逐靶点范围：

| metric | mean | min | max | n_positive_sign | n_negative_sign |
| --- | --- | --- | --- | --- | --- |
| brier | 0.3086 | 0.2263 | 0.4899 | 6 | 0 |
| log_loss | 0.8643 | 0.6525 | 1.4699 | 6 | 0 |
| bss_transport | -0.0809 | -0.4626 | 0.2230 | 2 | 4 |
| bss_target_oracle | -1.1205 | -2.7070 | -0.0163 | 0 | 6 |
| calibration_slope | 3.6293 | 0.2418 | 16.0453 | 6 | 0 |

配对汇总（与靶点等权宏平均分开报告）：Brier=0.312622，log loss=0.860474，BSS-transport=-0.089528，BSS-target-oracle=-1.239891。

U 层初步判定：主分析靶点等权宏平均 `BSS-transport=-0.080917`，平均未改善，按方案 6.6 初步判为**当前零标签映射未解决迁移**；逐靶点方向为 2 个正、4 个负，显示明显靶点依赖。AAH-I 的 19 个阳性保留在主分析中并单独可见。

## 2. 预设敏感性

U 层五项预设敏感性均保留，不因结果删减：

| dimension | analysis_set | cluster_rule | source_weighting | model_form | brier | bss_target_oracle | bss_transport | log_loss |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| main | primary_six | sequence_cluster90 | target_equal | linear | 0.3086 | -1.1205 | -0.0809 | 0.8643 |
| pair_source_weighting | primary_six | sequence_cluster90 | pair | linear | 0.3395 | -1.3191 | -0.0397 | 0.9342 |
| primary_five_excluding_aah_i | primary_five | sequence_cluster90 | target_equal | linear | 0.3140 | -1.1930 | -0.0921 | 0.8411 |
| restricted_cubic_spline | primary_six | sequence_cluster90 | target_equal | restricted_cubic_spline | 0.2909 | -0.9203 | -0.0240 | 0.8032 |
| sequence_cluster80 | primary_six | sequence_cluster80 | target_equal | linear | 0.3086 | -1.1200 | -0.0804 | 0.8638 |
| strict_outcome | strict_six | sequence_cluster90 | target_equal | linear | 0.2929 | -0.7550 | -0.0670 | 0.8504 |

预测性适配的对应敏感性（预测性 F 不扩展到 RCS，因为 RCS 只用于 U 层预设形态敏感性）：

| dimension | analysis_set | cluster_rule | source_weighting | model_form | model_layer_label | brier | bss_target_oracle | log_loss |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| main | primary_six | sequence_cluster90 | target_equal | linear | U_transport_zero_label | 0.3086 | -1.1205 | 0.8643 |
| main | primary_six | sequence_cluster90 | target_equal | linear | predictive_F | 0.1617 | 0.0798 | 0.4835 |
| main | primary_six | sequence_cluster90 | target_equal | linear | predictive_I | 0.1693 | 0.0148 | 0.5110 |
| pair_source_weighting | primary_six | sequence_cluster90 | pair | linear | U_transport_zero_label | 0.3395 | -1.3191 | 0.9342 |
| pair_source_weighting | primary_six | sequence_cluster90 | pair | linear | predictive_F | 0.1622 | 0.0766 | 0.4848 |
| pair_source_weighting | primary_six | sequence_cluster90 | pair | linear | predictive_I | 0.1680 | 0.0241 | 0.5065 |
| primary_five_excluding_aah_i | primary_five | sequence_cluster90 | target_equal | linear | U_transport_zero_label | 0.3140 | -1.1930 | 0.8411 |
| primary_five_excluding_aah_i | primary_five | sequence_cluster90 | target_equal | linear | predictive_F | 0.1807 | 0.0060 | 0.5378 |
| primary_five_excluding_aah_i | primary_five | sequence_cluster90 | target_equal | linear | predictive_I | 0.1793 | 0.0125 | 0.5338 |
| sequence_cluster80 | primary_six | sequence_cluster80 | target_equal | linear | U_transport_zero_label | 0.3086 | -1.1200 | 0.8638 |
| sequence_cluster80 | primary_six | sequence_cluster80 | target_equal | linear | predictive_F | 0.1599 | 0.0871 | 0.4790 |
| sequence_cluster80 | primary_six | sequence_cluster80 | target_equal | linear | predictive_I | 0.1699 | 0.0108 | 0.5128 |
| strict_outcome | strict_six | sequence_cluster90 | target_equal | linear | U_transport_zero_label | 0.2929 | -0.7550 | 0.8504 |
| strict_outcome | strict_six | sequence_cluster90 | target_equal | linear | predictive_F | 0.1723 | 0.0940 | 0.5073 |
| strict_outcome | strict_six | sequence_cluster90 | target_equal | linear | predictive_I | 0.1831 | 0.0156 | 0.5444 |

## 3. Oracle I 与描述性 F（非部署诊断）

以下两层使用留出靶点标签，分别标为 `oracle_base_rate_scenario` 与 `descriptive_in_sample_adaptation`，不是零标签或折外预测结果。

| held_out_target | model_layer_label | brier | log_loss | bss_target_oracle | calibration_status |
| --- | --- | --- | --- | --- | --- |
| aah-i | descriptive_in_sample_adaptation | 0.0635 | 0.1914 | 0.4813 | evaluable |
| aah-i | oracle_base_rate_scenario | 0.1115 | 0.3674 | 0.0897 | evaluable |
| human-vegf110 | descriptive_in_sample_adaptation | 0.2359 | 0.6640 | 0.0252 | evaluable |
| human-vegf110 | oracle_base_rate_scenario | 0.2359 | 0.6641 | 0.0250 | evaluable |
| nipah-glycoprotein-g | descriptive_in_sample_adaptation | 0.0877 | 0.3149 | 0.0101 | evaluable |
| nipah-glycoprotein-g | oracle_base_rate_scenario | 0.0878 | 0.3151 | 0.0095 | evaluable |
| pd-l1 | descriptive_in_sample_adaptation | 0.2030 | 0.5958 | 0.0040 | evaluable |
| pd-l1 | oracle_base_rate_scenario | 0.2116 | 0.6159 | -0.0380 | evaluable |
| sars-cov-2-rbd | descriptive_in_sample_adaptation | 0.1290 | 0.4192 | 0.0242 | evaluable |
| sars-cov-2-rbd | oracle_base_rate_scenario | 0.1380 | 0.4479 | -0.0439 | evaluable |
| trem2 | descriptive_in_sample_adaptation | 0.2190 | 0.6283 | 0.0997 | evaluable |
| trem2 | oracle_base_rate_scenario | 0.2302 | 0.6525 | 0.0541 | evaluable |

Oracle I 的截距严格求解为 `mean(expit(a+eta_U)) = target prevalence`；描述性 F 使用留出靶点全样本拟合，未被称为预测性改善。

## 4. 预测性 I 与 F（同一 CV5 外测行）

每个预测性结果均来自对应 LOTO 折固定 `eta_U` 后在适配训练折上的更新；U、I、F 使用完全相同的外测配对行。lambda 只由外折训练侧内部的 3 折 `StratifiedGroupKFold` 选择。

| held_out_target | model_layer_label | brier | log_loss | bss_target_oracle | bss_adaptation_train | calibration_status |
| --- | --- | --- | --- | --- | --- | --- |
| aah-i | U_transport_zero_label | 0.2655 | 0.7244 | -1.1685 | -1.1632 | evaluable |
| aah-i | predictive_F | 0.0687 | 0.2150 | 0.4388 | 0.4402 | evaluable |
| aah-i | predictive_I | 0.1117 | 0.3682 | 0.0880 | 0.0902 | evaluable |
| human-vegf110 | U_transport_zero_label | 0.3226 | 0.8567 | -0.3330 | -0.3328 | evaluable |
| human-vegf110 | predictive_F | 0.2360 | 0.6644 | 0.0245 | 0.0247 | evaluable |
| human-vegf110 | predictive_I | 0.2360 | 0.6643 | 0.0247 | 0.0248 | evaluable |
| nipah-glycoprotein-g | U_transport_zero_label | 0.3001 | 0.7948 | -2.3878 | -2.3877 | evaluable |
| nipah-glycoprotein-g | predictive_F | 0.0877 | 0.3151 | 0.0098 | 0.0098 | evaluable |
| nipah-glycoprotein-g | predictive_I | 0.0878 | 0.3151 | 0.0095 | 0.0095 | evaluable |
| pd-l1 | U_transport_zero_label | 0.2263 | 0.6525 | -0.1103 | -0.1102 | evaluable |
| pd-l1 | predictive_F | 0.2055 | 0.6012 | -0.0082 | -0.0081 | evaluable |
| pd-l1 | predictive_I | 0.2116 | 0.6160 | -0.0382 | -0.0381 | evaluable |
| sars-cov-2-rbd | U_transport_zero_label | 0.4899 | 1.4699 | -2.7070 | -2.7062 | evaluable |
| sars-cov-2-rbd | predictive_F | 0.1319 | 0.4286 | 0.0018 | 0.0021 | evaluable |
| sars-cov-2-rbd | predictive_I | 0.1385 | 0.4492 | -0.0482 | -0.0480 | evaluable |
| trem2 | U_transport_zero_label | 0.2473 | 0.6876 | -0.0163 | -0.0155 | evaluable |
| trem2 | predictive_F | 0.2404 | 0.6770 | 0.0119 | 0.0127 | evaluable |
| trem2 | predictive_I | 0.2305 | 0.6531 | 0.0527 | 0.0534 | evaluable |

主分析目标级配对比较：

| comparison | metric | macro_mean_difference | n_improved | n_worsened | min_difference | max_difference |
| --- | --- | --- | --- | --- | --- | --- |
| predictive_I_minus_U | brier | -0.1393 | 6 | 0 | -0.3514 | -0.0147 |
| predictive_I_minus_U | log_loss | -0.3533 | 6 | 0 | -1.0207 | -0.0345 |
| predictive_I_minus_U | bss_target_oracle | 1.1352 | 6 | 0 | 0.0690 | 2.6588 |
| predictive_F_minus_U | brier | -0.1469 | 6 | 0 | -0.3580 | -0.0069 |
| predictive_F_minus_U | log_loss | -0.3807 | 6 | 0 | -1.0413 | -0.0106 |
| predictive_F_minus_U | bss_target_oracle | 1.2003 | 6 | 0 | 0.0282 | 2.7088 |
| predictive_F_minus_I | brier | -0.0076 | 4 | 2 | -0.0430 | 0.0099 |
| predictive_F_minus_I | log_loss | -0.0274 | 3 | 3 | -0.1532 | 0.0239 |
| predictive_F_minus_I | bss_target_oracle | 0.0650 | 4 | 2 | -0.0408 | 0.3508 |

预测性结果按六靶点逐项报告：I 相对 U 的 Brier 宏平均差为 -0.139279，6/6 个靶点改善；I 相对训练折基率预测的 `BSS-adaptation-train` 宏平均为 0.015318。F 相对 I 的 Brier 宏平均差为 -0.007625，4/6 改善、2/6 变差，`BSS-adaptation-train=0.080209`。方向计数按靶点等权统计，不能作为未来靶点上的总体不确定性声明。

## 5. 阶段状态与失败台账

| Stage | Status | Evidence |
| --- | --- | --- |
| stage_0_input_and_fold_qc | complete | 19/19 structural checks |
| stage_1_U_linear_and_RCS_with_predictions_metrics_reliability | complete | 176/176 evaluable |
| stage_2_oracle_I_and_descriptive_F | complete | 176/176 evaluable diagnostics |
| stage_3_predictive_I_F_with_nested_inner_cv | complete | 2640/2640 lambda candidates evaluable |
| stage_4_figures_tables_qc_and_handoff | complete | 40/40 final checks passed |

失败台账共 28 条，其中预测性阶段逐外折校准诊断触及 `±40` 数值边界 28 条，其他模型或回归失败 0 条；相应校准指标保留为缺失并标记 `failed`，未替换算法、未删除外折、未使用外测标签补救。预测模型与 lambda 候选可估计情况以 `stage3_summary.json` 和失败台账为准。

| stage | failure_type | rows |
| --- | --- | --- |
| stage_3_predictive_I_F_with_nested_inner_cv | calibration_failure | 28 |

## 6. QC、冻结与交付

- 输入 SHA-256 在开机和结束各核对一次，全部一致。
- 2,084 个 `pair_id`、44 个 LOTO 折、220 个 CV5 折、364 条跨靶点同簇剔除记录均与冻结清单一致。
- U：176 个模型-靶点评估、64,560 条测试配对预测，全部可估计；RCS 公式与尾部线性单元测试通过。
- Oracle/描述性：176 个诊断全覆盖；oracle 平均预测概率与真实基率误差小于 `1e-10`。
- 预测性：440 个目标×权重×外折评估、32,280 条相同外测行 U/I/F 预测、2,640 个候选 lambda、387,360 条内折分配；所有内折训练/测试序列簇隔离通过。
- F 修正回归：状态 `passed`；528 个 F 拟合的目标重构、恒等关系、预测重构与梯度结果详见 `rq3_f_regression_tests.json`。
- 最终 QC：40/40 项通过。
- `failed_item_ledger.csv`、`qc_checks.json`、脚本、命令、阶段日志、图与总表均已保存；`final_hash_manifest.csv` 和 `manifest_sha256.txt` 在本脚本退出、日志关闭后生成。

信息层级边界：不得把 `BSS-target-oracle` 当作可部署零标签参照；不得把 oracle I 或描述性 F 称为预测改善；不得用“校准斜率区间包含 1”作为成功标准。六靶点结果不提供对未来靶点的总体不确定性声明。

**正式 RQ3 结论待主执行线独立验收。**
