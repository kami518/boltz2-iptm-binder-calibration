"""Build the human-readable RQ1/RQ2 report only from accepted output tables."""
from pathlib import Path
import json
import pandas as pd

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
OUT = ROOT / "analysis" / "rq12_outputs"
manifest = json.loads((OUT / "rq12_run_manifest.json").read_text(encoding="utf-8"))
assert manifest["status"] == "complete"
assert all(manifest["quality_gates"].values())
metrics = pd.read_csv(OUT / "rq12_target_metrics.csv")
summary = pd.read_csv(OUT / "rq12_summary_metrics.csv")
ht = pd.read_csv(OUT / "rq2_hierarchical_targets.csv")
hh = pd.read_csv(OUT / "rq2_hierarchical_hyperparameters.csv").set_index("analysis_set")
reliability = pd.read_csv(OUT / "rq2_reliability_bins.csv")

labels = {
    "aah-i": "AAH-I",
    "human-vegf110": "VEGF110",
    "nipah-glycoprotein-g": "Nipah G",
    "pd-l1": "PD-L1",
    "sars-cov-2-rbd": "SARS-CoV-2 RBD",
    "trem2": "TREM2",
}
set_labels = {
    "primary_six": "主结局／六靶点",
    "primary_five": "主结局／排除 AAH-I",
    "strict_six": "严格结局／六靶点",
    "strict_five": "严格结局／排除 AAH-I",
}
weight_labels = {"pair_weighted": "配对数加权", "target_equal": "靶点等权"}
mi = metrics.set_index(["analysis_set", "target", "metric"])
si = summary.set_index(["analysis_set", "weighting", "metric"])
hi = ht.set_index(["analysis_set", "target"])

def value(target, metric, analysis_set="primary_six"):
    return float(mi.loc[(analysis_set, target, metric), "estimate"])

def ci(target, metric, analysis_set="primary_six"):
    row = mi.loc[(analysis_set, target, metric)]
    return f'{row.estimate:.3f}（{row.ci_low:.3f}–{row.ci_high:.3f}）'

def summary_value(analysis_set, weighting, metric):
    return float(si.loc[(analysis_set, weighting, metric), "estimate"])

def summary_ci(analysis_set, weighting, metric):
    row = si.loc[(analysis_set, weighting, metric)]
    return f'{row.estimate:.3f}（{row.ci_low:.3f}–{row.ci_high:.3f}）'

lines = []
add = lines.append
add("# RQ1／RQ2 结果与验收记录")
add("")
add("日期：2026-11-09。分析对象：冻结主表中的 2,084 个设计—靶点配对、6 个靶点。此报告只覆盖原始 Boltz‑2 ipTM 的靶点内区分（RQ1）与原始校准／靶点异质性（RQ2），不包含 LOTO、筛选阈值、Pinc 或跨模型迁移结论。")
add("")
add("## 一、先看结论")
add("")
add(f'1. **ipTM 有靶点内排序信息，但强度并不统一。** 六靶点等权的逐靶点 ROC-AUC 均值为 {summary_ci("primary_six", "target_equal", "roc_auc")}；排除 AAH-I 后为 {summary_ci("primary_five", "target_equal", "roc_auc")}。PD-L1 的区间包含 0.5，其他五靶点区间高于 0.5，但这不等于都有足够的实际筛选价值。')
add(f'2. **AAH-I 必须单独解释。** 其 AUC 为 {ci("aah-i", "roc_auc")}，远高于其余靶点；19 个 SPR 阳性中 18 个位于选择性 BLI 确认子集。它明显抬高六靶点等权 AUC，因此五靶点结果是必要对照，而不是可选补充。')
add(f'3. **原始 ipTM 不能直接解释为实验结合成功概率。** 主结局下靶点等权 Brier 为 {summary_value("primary_six", "target_equal", "brier"):.3f}、log loss 为 {summary_value("primary_six", "target_equal", "log_loss"):.3f}；相对各靶点恒定阳性率参照的汇总 BSS 为 {summary_value("primary_six", "target_equal", "bss_target_prevalence"):.3f}，配对数加权 BSS 为 {summary_value("primary_six", "pair_weighted", "bss_target_prevalence"):.3f}，均为负。')
add(f'4. **必要对照没有改变主要判断。** 严格结局的六靶点等权 AUC 为 {summary_value("strict_six", "target_equal", "roc_auc"):.3f}，五靶点为 {summary_value("strict_five", "target_equal", "roc_auc"):.3f}；相应 BSS 仍为 {summary_value("strict_six", "target_equal", "bss_target_prevalence"):.3f} 和 {summary_value("strict_five", "target_equal", "bss_target_prevalence"):.3f}。改变 Weak 的处理或排除 AAH-I，都不会把原始 ipTM 变成已校准概率。')
add("5. **层级模型支持靶点异质性，但不支持用一个新映射直接外推。** AAH-I 的极端校准参数明显收缩；排除 AAH-I 后斜率的靶点间标准差大幅下降。当前结果是描述性校准诊断，不是折外适配性能。")
add("")
add("## 二、数据与验收")
add("")
add("| 项目 | 已验收结果 |")
add("|---|---|")
add("| 主结局／严格结局 | 六靶点 N=2,084，阳性 748／678；五靶点 N=1,951，阳性 729／659 |")
add("| 主表 SHA-256 | `" + manifest["input_master_sha256"] + "` |")
add("| 序列簇映射 SHA-256 | `" + manifest["cluster_mapping_sha256"] + "` |")
add("| 序列簇 | parasail 1.3.4 全对局部比对；90% identity、双向覆盖 ≥80% 的主规则；同时保存 80% 严格簇 |")
add("| Bootstrap | 固定种子 20261109；每靶点、每结局均保留 2,000 次有效序列簇重采样；失败抽样补抽，不减少有效次数 |")
add(f'| 层级模型 | 四个分析集均保留随机截距＋随机斜率，无模型回退；均收敛且信息矩阵正定 |')
add(f'| 积分稳定性 | 自适应 Gauss–Hermite 9／11 节点；最大靶点估计差 {manifest["hierarchical_stability"]["max_target_estimate_absolute_difference"]:.2e}，最大超参数差 {max(manifest["hierarchical_stability"]["hyperparameter_absolute_differences"].values()):.2e}，均小于 0.02 |')
add("| 独立数值核验 | 6 项测试通过，包括与独立密集二维网格直接积分比较 |")
add("")
add("普通非自适应 Gauss–Hermite 积分未通过预设稳定性检查，失败清单和模型结果已单独归档；正式结果没有采用那些估计。积分实现的修正不改变主表、结局、先验、bootstrap 或模型结构。")
add("")
add("## 三、RQ1：逐靶点区分能力")
add("")
add("区间均为 90% 序列簇 bootstrap 的 95% 百分位区间。AP 为 average precision，是 PR 曲线的摘要，不应与阳性率或跨靶点 pooled AUC 混淆。")
add("")
add("| 靶点 | 主 assay | N | 主／严格阳性 | 阳性率 | ROC-AUC（95% CI） | AP（95% CI） | 序列簇数 |")
add("|---|---|---:|---:|---:|---:|---:|---:|")
for target in labels:
    row = mi.loc[("primary_six", target, "roc_auc")]
    strict_n = int(mi.loc[("strict_six", target, "roc_auc"), "positives"])
    add(f'| {labels[target]} | {row.assay.upper()} | {int(row.n)} | {int(row.positives)}／{strict_n} | {value(target, "positive_rate"):.3f} | {ci(target, "roc_auc")} | {ci(target, "average_precision")} | {int(row.sequence_clusters)} |')
add("")
add("![逐靶点 ROC-AUC 与 AP](RQ1_target_discrimination.png)")
add("")
add("![逐靶点 ipTM 分布](RQ1_ipTM_distributions.png)")
add("")
add("### 汇总结论")
add("")
add("下表的 AUC／AP 是**逐靶点指标的加权平均**，不是把所有配对混合后计算的 pooled AUC／AP。配对数加权回答当前送测配对的平均表现；靶点等权避免 Nipah 的大样本主导。")
add("")
add("| 分析集 | 权重 | ROC-AUC（95% CI） | AP（95% CI） |")
add("|---|---|---:|---:|")
for analysis_set in ["primary_six", "primary_five"]:
    for weighting in ["pair_weighted", "target_equal"]:
        add(f'| {set_labels[analysis_set]} | {weight_labels[weighting]} | {summary_ci(analysis_set, weighting, "roc_auc")} | {summary_ci(analysis_set, weighting, "average_precision")} |')
add("")
add("除 AAH-I 外，AUC 点估计约为 0.56–0.69：存在排序信号，但并非强区分。AP 必须结合阳性率读，例如 SARS-CoV-2 RBD 的高 AP 部分来自 0.843 的高阳性率基线，Nipah 的 AP 只有 0.125、基线为 0.098。")
add("")
add("## 四、RQ2：原始校准")
add("")
add("校准诊断模型为 `logit(P(Y=1)) = α + β × logit(clip(ipTM))`，理想参照为 α=0、β=1。α、β 是联合拟合参数，不能把 α 单独解释为固定斜率后的 calibration-in-the-large。下表同时给出未汇聚和层级部分汇聚点估计；Brier、log loss 和可靠性图始终直接使用原始 ipTM，而非拟合后的概率。")
add("")
add("| 靶点 | Brier | BSS | log loss | α 未汇聚 | β 未汇聚 | α 部分汇聚 | β 部分汇聚 |")
add("|---|---:|---:|---:|---:|---:|---:|---:|")
for target in labels:
    row = hi.loc[("primary_six", target)]
    add(f'| {labels[target]} | {value(target, "brier"):.3f} | {value(target, "bss_target_prevalence"):.3f} | {value(target, "log_loss"):.3f} | {row.alpha_unpooled:.3f} | {row.beta_unpooled:.3f} | {row.alpha_partial:.3f} | {row.beta_partial:.3f} |')
add("")
add("![原始 ipTM 可靠性曲线](RQ2_raw_reliability.png)")
add("")
add("Nipah 是最直观的反例：阳性／阴性 ipTM 中位数分别为 0.910／0.895，虽然 AUC 为 0.603，但观测阳性率只有 0.098，原始 Brier 达 0.630。AAH-I、TREM2 也明显过自信；VEGF110 的原始 Brier／log loss 相对较好，但其 BSS 点估计仍为 −0.040，不能据此宣布理想校准。可靠性图中稀疏箱只作低精度展示。")
add("")
add("### 权重与必要敏感性")
add("")
add("汇总 BSS 先按指定权重合并逐靶点 Brier 和各靶点恒定阳性率参照的 Brier，再计算 `1 − BS／BS_reference`；不是六个逐靶点 BSS 的简单平均。")
add("")
add("| 分析集 | 权重 | Brier | log loss | BSS | ECE（次要） |")
add("|---|---|---:|---:|---:|---:|")
for analysis_set in ["primary_six", "primary_five", "strict_six", "strict_five"]:
    for weighting in ["pair_weighted", "target_equal"]:
        add(f'| {set_labels[analysis_set]} | {weight_labels[weighting]} | {summary_value(analysis_set, weighting, "brier"):.3f} | {summary_value(analysis_set, weighting, "log_loss"):.3f} | {summary_value(analysis_set, weighting, "bss_target_prevalence"):.3f} | {summary_value(analysis_set, weighting, "ece"):.3f} |')
add("")
add("严格结局下阳性率降低，因此 Brier／log loss 的绝对值不能与主结局作同基率下的性能比较；关键对照是区分指标相近、BSS 仍为负，以及逐靶点可靠性关系没有被修复。ECE 受固定分箱和稀疏箱影响，只作辅助描述。")
add("")
add("## 五、层级收缩与异质性")
add("")
add("![校准参数的部分汇聚](RQ2_partial_pooling.png)")
add("")
add("| 分析集 | 总体 α 均值 | 总体 β 均值 | τ(α) | τ(β) | 模型 |")
add("|---|---:|---:|---:|---:|---|")
for analysis_set in ["primary_six", "primary_five", "strict_six", "strict_five"]:
    row = hh.loc[analysis_set]
    add(f'| {set_labels[analysis_set]} | {row.mu_intercept:.3f} | {row.mu_slope:.3f} | {row.tau_intercept:.3f} | {row.tau_slope:.3f} | 随机截距＋随机斜率 |')
add("")
aah = hi.loc[("primary_six", "aah-i")]
trem = hi.loc[("primary_five", "trem2")]
add(f'AAH-I 的 α／β 从未汇聚的 {aah.alpha_unpooled:.3f}／{aah.beta_unpooled:.3f} 收缩至 {aah.alpha_partial:.3f}／{aah.beta_partial:.3f}。排除 AAH-I 后 τ(β) 从 {hh.loc["primary_six", "tau_slope"]:.3f} 降为 {hh.loc["primary_five", "tau_slope"]:.3f}，说明六靶点中明显的斜率异质性很大程度受 AAH-I 影响。此时 TREM2 的 α／β 收缩至 {trem.alpha_partial:.3f}／{trem.beta_partial:.3f}。')
add("")
add("这不等于证明五靶点共享同一个真实校准斜率：靶点数少、来源块更少，τ 的不确定性仍大。靶点后验 SD 以估计后的超参数为条件，超参数不确定性另行列于 CSV，不能将其当作完整的新靶点预测区间。")
add("")
add("## 六、可以写入论文与不能写入论文的表述")
add("")
add("**可以：** 在与 Overath 镜像精确序列不重叠的当前主队列中，Boltz‑2 ipTM 在部分靶点具有中等靶点内排序信息；原始数值不能直接解释为实验结合成功概率，区分和校准均呈明显上下文依赖。AAH-I 的特殊阳性结构显著影响六靶点汇总。")
add("")
add("**暂时不可以：** 已验证通用再校准映射；已证明跨靶点零标签可迁移；已验证 Pinc 或组合指标更优；已证明 Boltz‑1→Boltz‑2 漂移；已证明 assay 间稳健；已覆盖所有 minibinder。上述分别需要 RQ3、RQ5／RQ6、统一模型对照及额外数据审计。")
add("")
add("还需保留的限制：六靶点并非六个完全独立设计流程；当前序列簇 bootstrap 不消除 campaign／source 混杂；AAH-I 的选择性确认可能放大表观区分；区间不能超越有限独立靶点和来源块；当前没有折外预测评价。")
add("")
add("## 七、可复现文件")
add("")
add("| 文件 | 内容 |")
add("|---|---|")
for name, description in [
    ("rq12_run_manifest.json", "最终质量门、输入哈希、bootstrap 审计、积分稳定性和软件版本"),
    ("rq12_analysis_config.json", "冻结分析集、种子、bootstrap、分箱和层级模型配置"),
    ("rq1_target_metrics.csv", "逐靶点 ROC-AUC、AP、阳性率及 bootstrap 区间"),
    ("rq2_target_metrics.csv", "逐靶点 Brier、BSS、log loss、ECE、校准参数及区间"),
    ("rq12_summary_metrics.csv", "四种分析集、两种权重的汇总指标"),
    ("rq2_reliability_bins.csv", "固定分箱计数、平均 ipTM、阳性率和 Wilson 区间"),
    ("rq2_hierarchical_targets.csv", "未汇聚与部分汇聚的靶点校准参数"),
    ("rq2_hierarchical_hyperparameters.csv", "层级模型超参数、标准误及收敛诊断"),
    ("sequence_clusters_rq12.csv", "90%／80% 序列簇映射"),
    ("sequence_cluster_summary.json", "全对比对及簇规模审计"),
]:
    add(f'| [{name}]({name}) | {description} |')
add("")
add("主执行脚本：`analysis/run_rq12_analysis_2026-11-09.py`；数值核心：`analysis/rq12_core.py`、`analysis/rq12_adaptive.py`；测试：`analysis/test_rq12_core.py`。配置文件的运行副本也嵌入最终 manifest。")
add("")

(OUT / "RQ1_RQ2_results_2026-11-09.md").write_text("\n".join(lines), encoding="utf-8")
print(f"Wrote {OUT / 'RQ1_RQ2_results_2026-11-09.md'} with {len(lines)} lines")
