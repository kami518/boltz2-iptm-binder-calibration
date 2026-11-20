# RQ4 冻结协议：执行一致性澄清

日期：2026-11-14。

本澄清不改变任何队列、结局、选择规则、命中数、预算、pilot 参数或模型设定。

原冻结协议第 3 节把空选集合的 precision、recall、enrichment 和单位成本统一写为不可评价。执行脚本将空选时的 **recall 定义为 0**，因为此时 hits=0、候选库阳性数大于 0，召回率具有明确的数学值 0／总阳性数；只有 precision、enrichment 和单位 hit 成本因零入选或零命中而记为 NA。

因此，正式结果采用：

- 空选：`n_selected=0`、`hits=0`、`recall=0`；
- 空选：`precision=NA`、`enrichment=NA`、`cost_per_hit=NA`。

这是一项文字口径纠错，不是看到结果后的分析规则调整。原冻结协议及其哈希保持不变，本文件单独留痕。
