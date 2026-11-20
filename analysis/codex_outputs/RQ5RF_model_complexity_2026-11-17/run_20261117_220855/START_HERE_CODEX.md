# RQ5-RF 全量执行交接

日期：2026-11-17  
工作区：`<REPRODUCTION_ROOT>`  
角色：Codex 是本轮唯一计算执行方，负责实现、全量运行、自检、独立验证和结果交付；不是第二轮方法审查方。

## 1. 唯一规则来源

必须完整阅读并严格执行：

`<REPRODUCTION_ROOT>\analysis\RQ5RF_model_complexity_frozen_protocol_2026-11-17.md`

SHA-256：

`d9db1c3a09d1264580ca56b797922d55fce0e563c04ecb9dab76ac4744a25e28`

草案和一轮审查仅供追踪：

- `analysis/RQ5RF_model_complexity_protocol_draft_2026-11-17.md`
- `analysis/codex_outputs/RQ5RF_protocol_review_2026-11-17/CODEX_REVIEW.md`

不得重新评审、增删模型、修改门槛或提出第二轮审查。冻结协议优先于草案及本交接中的任何概括；如有疑义，按冻结协议原文执行。

## 2. 写入边界

唯一正式运行目录：

`analysis/codex_outputs/RQ5RF_model_complexity_2026-11-17/run_YYYYMMDD_HHMMSS/`

所有脚本、配置、日志、中间 checkpoint、预测、图表和结果只能写入该新目录。当前交接目录只允许保存执行日志或交接副本（如确有需要）。

不得修改：

- 冻结协议及其哈希；
- v1/v2/v3 流程、任何稿件或 PDF；
- 冻结主表、Pinc 扩展表、序列簇和 folds；
- RQ1–RQ4、原 RQ5、RQ3E、RQ4C/E 或任何旧 Codex 运行；
- 文献、数据库和外部数据。

不得联网、下载、安装或升级依赖；不得启动其他 RQ 或额外 ML。

## 3. 必须实际完成的工作

不要只给建议或示例代码。请直接完成以下全部阶段：

### Stage 0：预检与冻结快照

1. 创建唯一 run 目录，复制冻结协议和本交接。
2. 核验全部输入路径、SHA-256、行数、唯一键、complete-case 数量和阳性数。
3. 核验 Python 与 scikit-learn；必须为审查确认的 scikit-learn 1.9.0，否则阻断，不得换环境拼接。
4. 保存 `run_config.snapshot.json`、`final_input_hashes.csv`、环境、命令和日志。
5. 在任何新外测结果生成前完成配置快照。

### Stage 1：共享行集、当前权重、folds 和单元测试

1. 实现并测试第 4 节当前行集权重重算，不直接沿用冻结权重列。
2. 复用冻结 LOTO、cross-target exclusions 和 target-CV5；实现规定排序与 inner split 复用。
3. 建立可恢复的逐 outer fold checkpoint。
4. 先通过单元测试再进入正式拟合；测试须覆盖冻结协议第 11 节全部项目。

### Stage 2：U 层全量拟合

1. 四个 pooled 暴露：`pooled_logistic_ipTM`、`pooled_logistic_combo`、`rf_ipTM`、`rf_combo`。
2. 完成 `{primary_six,primary_five,strict_six} × {target_equal,pair}` 六配置。
3. 严格执行 source-only leave-one-source-target-out 内层选参、有效 RF 网格去重、共享 seed、最终重拟合和外测预测。
4. 保存所有 inner candidates，不只保存最佳候选；保存权重、seed、classes 和行集审计字段。
5. 原 target-intercept RQ5 结果只读复用为次要锚点，不重跑、不纳入 10.2。

### Stage 3：I／F 全量适配

对四个新暴露、六配置执行冻结 target-CV5 适配。固定各自 `eta_U`，复用相同 inner assignments，不重训 RF，不改变 F penalty／lambda／失败规则。保存完整折外预测和适配参数。

### Stage 4：指标、图表和汇总

1. 从完整 pair predictions 生成逐折、逐靶点、靶点等权宏平均、pair-pooled、全部预设 pairwise comparisons 和可靠性。
2. 完整样本硬门不得被 `dropna()` 或 inner join 绕过。
3. 生成冻结协议要求的全部图，并实际打开检查；AAH-I 和五靶点结果均可见。
4. 按冻结中性带和 10.2 规则计算方向一致性筛查结果，同时报告绝对效应量。

### Stage 5：独立验证和最终交付

1. 编写并运行独立验证，不得只调用生产汇总函数。
2. 重算全部选参和主要汇总；至少复跑一个单特征与一个双特征 outer fold，逐元素核对相同 seed 预测。
3. 汇总单元测试、独立验证、失败台账和 QC；所有失败均透明记录。
4. 生成最终 hash manifest、`RUN_STATUS.json`、`HANDOFF_RESULT.md`、`handoff_result.json`。
5. 最终报告列出运行目录、主要数字、六配置结果、10.2 是否满足、失败项和自检状态。

## 4. 计算与恢复

冻结规模约为 3,492 个 500-tree RF fits。不得为缩短时间减少树数、候选、fold、分析集或权重。

允许：共享只读数据与 split、有效网格去重、单次 fit 内固定并行、及时释放森林、逐 outer fold checkpoint、分批恢复。正式全量前记录代表性 fit 的耗时与峰值内存；资源不足时保持规则并分批，不改变科学规则。

不得在候选拟合外层再开启进程池造成 oversubscription。若并行不能通过重复拟合逐元素一致性测试，改为串行。

## 5. 失败和结果驱动修改禁令

- 代码错误可修复并记录，但不得改变冻结科学规则。
- 若发现必须改协议才能继续的问题，停止受影响阶段，在 `RUN_STATUS.json` 和失败台账标记 blocked，并完整交付已完成部分；不得静默改规则。
- 不因中间结果有利／不利提前停止、追加模型、调网格、改权重、改结局或隐藏敏感性。
- 不用 outer 标签调参，不用 target-oracle 作可部署规则，不把 RF 概率称为已校准。

## 6. 交付口径

所有结果仅为 post hoc exploratory model-complexity sensitivity。即使满足 10.2，也只能表述为“方向一致的探索性增量，值得在独立 complete-case 数据中重复核验”；不得声称通用校准、纯交互、RF/Pinc 普遍优越或已解决跨靶点迁移。

开始执行后持续推进至全部完成，除非遇到冻结协议规定的真正阻断。最终消息只汇报真实完成情况和文件路径，不把计划当作结果。
