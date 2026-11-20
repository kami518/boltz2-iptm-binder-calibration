# RQ3-E 720 行确认性外部迁移结果

**范围：8 靶点／720 行主外部队列及预定义敏感性分析；target-equal U 从上一运行原样复用，未使用外部标签训练或调参，未拟合 I／F。**

主 U：运输截距 **-1.139290**，共享 logit-ipTM 斜率 **0.510542**，源靶点等权基率参照 **0.467821**。

## 主分析：逐靶点

| 靶点 | 分层 | n／阳性 | 基率 | 平均预测 | Brier | Log loss | BSS-transport | BSS-target-oracle | 校准截距 | 校准斜率 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BHRF1 | new_target | 90／20 | 0.2222 | 0.6432 | 0.3613 | 0.9273 | -0.5495 | -1.0902 | -1.0722 | -0.3074 |
| IL-7Ra | new_target | 90／47 | 0.5222 | 0.5192 | 0.2347 | 0.6621 | 0.0705 | 0.0594 | -0.0195 | 1.3266 |
| RBX1 | new_target | 90／23 | 0.2556 | 0.6123 | 0.3098 | 0.8178 | -0.3166 | -0.6285 | -2.0197 | 1.8415 |
| TrkA | new_target | 90／19 | 0.2111 | 0.5238 | 0.2505 | 0.6928 | -0.0777 | -0.5042 | -1.8705 | 2.6018 |
| Nipah-G | same_target_new_campaign | 90／19 | 0.2111 | 0.4984 | 0.2465 | 0.6831 | -0.0606 | -0.4802 | -1.3782 | 1.0375 |
| PD-L1 | same_target_new_campaign | 90／39 | 0.4333 | 0.6085 | 0.2466 | 0.6857 | 0.0007 | -0.0042 | -4.7644 | 8.3242 |
| TREM2 | same_target_new_campaign | 90／65 | 0.7222 | 0.6137 | 0.1854 | 0.5598 | 0.3013 | 0.0758 | -1.4771 | 5.7244 |
| VEGF-A | same_target_new_campaign | 90／55 | 0.6111 | 0.5179 | 0.2370 | 0.6672 | 0.0822 | 0.0029 | 0.3665 | 1.3864 |

## 主分析：预定义汇总

| 汇总 | n 靶点／配对 | 阳性 | Brier | Log loss | BSS-transport | BSS-target-oracle | 校准截距 | 校准斜率 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| target_equal_macro／all_8 | 8／720 | 287 | 0.2590 | 0.7120 | -0.0687 | -0.3211 | -1.5294 | 2.7419 |
| target_equal_macro／new_target_4 | 4／360 | 109 | 0.2891 | 0.7750 | -0.2184 | -0.5409 | -1.2455 | 1.3656 |
| target_equal_macro／same_target_campaign_4 | 4／360 | 178 | 0.2289 | 0.6489 | 0.0809 | -0.1014 | -1.8133 | 4.1181 |
| pair_pooled／all_available_pairs | 8／720 | 287 | 0.2590 | 0.7120 | -0.0591 | -0.2714 | -0.7129 | 1.0046 |

## 敏感性分析：8 靶点等权宏平均

| 配置 | 类型 | n 配对 | Brier | Log loss | BSS-transport | BSS-target-oracle | BSS-T 为正靶点数 |
|---|---|---:|---:|---:|---:|---:|---:|
| primary_target_equal_weak_plus_best_ipsae | primary | 720 | 0.2590 | 0.7120 | -0.0687 | -0.3211 | 4／8 |
| outcome_medium_plus | outcome_sensitivity | 720 | 0.2639 | 0.7219 | -0.1011 | -0.3982 | 2／8 |
| outcome_integrated | outcome_sensitivity | 720 | 0.2532 | 0.6999 | -0.0388 | -0.2660 | 4／8 |
| outcome_twist_qc | outcome_sensitivity | 718 | 0.2620 | 0.7176 | -0.0961 | -0.5933 | 2／8 |
| outcome_twist_permissive | outcome_sensitivity | 718 | 0.2541 | 0.7015 | -0.0473 | -0.4516 | 3／8 |
| score_seed_mean | score_sensitivity | 720 | 0.2438 | 0.6796 | -0.0038 | -0.2355 | 6／8 |
| score_seed_median | score_sensitivity | 720 | 0.2471 | 0.6863 | -0.0177 | -0.2533 | 6／8 |
| score_seed0 | score_sensitivity | 720 | 0.2464 | 0.6840 | -0.0148 | -0.2495 | 5／8 |
| source_pair_weight | source_weight_sensitivity | 720 | 0.2311 | 0.6531 | 0.0266 | -0.1417 | 6／8 |

## 强制解释边界

- 确认性结论仅适用于 Adaptyv Weak+ 基率约 21.1%–72.2% 的 8 个达标靶点；低事件压力测试另行报告，高基率端仍未覆盖。
- 本研究评价已发布 score＋冻结 U 的运输，不是完全同一 score-generation protocol 下的纯靶点验证。
- 外部队列有两种高层 design model、九种底层 generator；不识别独立 generator 因果效应。
- 靶点等权宏平均、4 新增／4 同靶点分层与配对加权汇总不得混用；不同结局的绝对 Brier 不作同基率直接比较。
- BSS>0 只表示优于对应参照；校准斜率置信区间包含 1 不能单独作为校准成功证据。
- 未运行 Pinc、跨预测器、其他 ML 或任何外部标签适配。
