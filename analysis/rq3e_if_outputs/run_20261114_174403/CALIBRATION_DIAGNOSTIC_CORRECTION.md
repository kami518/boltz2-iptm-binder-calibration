# 校准诊断数值修正记录

最终有效运行：`analysis/rq3e_if_outputs/run_20261114_174403`。

原运行 `run_20261114_173316` 的 OOF 校准诊断使用未带回溯的 Newton 迭代。BHRF1／F 因 oversized step 后出现奇异 Hessian，被记为 `singular_or_separation`。独立复核同时发现：

1. 独立 BFGS 优化可得到有限解：校准截距约 −2.5027、斜率约 −0.9453；最终 Hessian 条件数约 61.7。
2. 对 `y ~ logit(p_F)` 的线性规划检查不支持完全分离。
3. 因而原标记是数值优化失败，不应解释为完全分离。

最终运行只将 OOF 校准诊断改为带解析梯度的 BFGS，并以有限参数、梯度和未归一化 Hessian 条件数验收；没有加入惩罚、伪逆或性能驱动的模型改动。I／F 的模型、折分、lambda 选择、OOF 概率、Brier、log loss、BSS、可靠性记录和两张图均与原运行逐字节一致。校准标准误仍定义为 `inv(X'WX)`。

原运行还暴露出日志哈希顺序问题：完成记录在哈希清单之后写入。最终运行已把完成日志置于哈希清单之前；24 个原始交付文件的大小与 SHA‑256 均通过独立核验。

`independent_verification.py` 和 `independent_verification.json` 是冻结运行完成后新增的审计文件，不属于上述 24 文件清单。冻结主表未修改，SHA‑256 仍为 `581c90e80e3cec3a6a974254c4ee0b151ce74387ea5b5298087860b02fbe5c91`。
