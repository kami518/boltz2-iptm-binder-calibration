# Evaluating Boltz-2 ipTM as an experimental-success probability across binder-design targets

**Article type:** Original Article  
**Draft status:** Version 2 structural revision following advisor feedback; not yet formatted for submission. The title, Introduction, terminology, and Discussion have been revised, and post hoc intercept-error and shared-budget analyses have been integrated. Author metadata, the verified external-dataset citation, repository links, and final journal-specific sizing remain unresolved. Composite Figures 1–7 retain the frozen sources listed in `manuscript/figures/figure_source_manifest_2026-11-15.md`.  
**Authors:** [To be supplied]  
**Affiliations:** [To be supplied]  
**Corresponding author:** [To be supplied]

## Abstract


Computational binder-design scores are used to rank candidates, estimate experimental success, and allocate assays, although these uses require different evidence. We evaluated Boltz-2 interface predicted TM-score (ipTM) in 2,084 experimentally assayed design–target pairs spanning six targets and 748 primary binding-positive outcomes. Analyses separated within-target discrimination, probability error, application to held-out targets, local adaptation, and experimental selection. A frozen mapping was also evaluated without external-label refitting in 720 pairs across eight targets. The target-equal mean area under the receiver operating characteristic curve was 0.675, decreasing to 0.620 after exclusion of the influential AAH-I target. Raw ipTM was not directly interpretable as a binding probability (target-equal Brier score 0.416; Brier skill score versus target-specific constant-prevalence references −1.417). Mappings learned from other targets had a transport-reference Brier skill score of −0.0809, positive in only two of six targets. A post hoc diagnostic found that local intercept adaptation reduced target-equal Brier by 45.1% (42.9% without AAH-I), whereas additional slope adaptation was inconsistent. The frozen external mapping also showed context-dependent performance under a different score-generation protocol. Common raw-score thresholds produced markedly different selected-set sizes. An exploratory retrospective shared-budget simulation suggested budget- and cost-dependent pilot value, not prospective assay savings. In these contexts, ipTM was a target-local ranking aid rather than a portable success probability or universal selection threshold. Absolute probabilities and shared-budget decisions require representative local measurements and cost-aware evaluation.

**Keywords:** binder design; Boltz-2; ipTM; probability calibration; external validation; experimental allocation

## Introduction


De novo protein design now combines target selection, backbone generation, sequence design, structural prediction, and experimental testing in a single development chain. Structure-guided methods can generate small binding proteins against selected target surfaces, while ProteinMPNN provides sequence designs for generated or existing backbones (Cao et al., 2022; Dauparas et al., 2022). RFdiffusion extended generative backbone design across several protein-design tasks, including binders, and BindCraft integrated binder generation with AlphaFold2-based design evaluation (Watson et al., 2023; Pacesa et al., 2025). These advances have changed both the scale and the diversity of candidate pools (Yang et al., 2026). They have not removed the experimental bottleneck: a researcher must still decide which designs to express and assay, how many to test, and how to distribute a finite experimental budget among targets. Confidence scores from predicted complexes are attractive at this decision point because they are available before the wet-lab work begins.

AlphaFold2 and AlphaFold3 established increasingly general frameworks for predicting protein and biomolecular structures, and Boltz-2 provides predicted complex structures, confidence measures, and affinity-related outputs (Jumper et al., 2021; Abramson et al., 2024; Passaro et al., 2025). In binder design, however, structural confidence can be asked to serve a further purpose: predicting whether a candidate will succeed in an experimental binding assay. Overath et al. (2025) analysed 3,766 experimentally characterised binders, compared structural, energetic, and confidence-derived features, and reported target-dependent success prediction and filtering behaviour. Badonyi and Toth-Petroczy (2026) compared Pinc with ipTM-type measures against experimental outcomes. These studies show that computational features can contain useful information about experimental success. They also make clear that the present study should not claim to discover target dependence or to introduce a better ranking score.

A distinct question concerns the numerical interpretation of a score when the experimental setting changes. Ranking asks whether binders receive higher scores than non-binders; calibration asks whether predicted probabilities match observed frequencies (Steyerberg et al., 2010; Van Calster et al., 2019). We use “application to a new target” to mean fitting a probability mapping on other targets and applying it without experimental labels from the destination target. “Local adaptation” then means updating that mapping after some destination-target measurements become available. These definitions distinguish two situations that can otherwise be conflated under the general term “generalization.” Sagadinov et al. (2026) addressed related calibration questions for DockQ-defined interface quality, whereas Lee et al. (2026) studied target-held-out top-K shortlisting. Neither structural interface quality nor shortlist ordering is interchangeable with an absolute probability of wet-lab binding. Moreover, a locally updated probability model has practical value only if its benefit exceeds the cost of obtaining the local measurements.

Here, we evaluate Boltz-2 ipTM as a candidate experimental-success probability in two retrospectively assembled, experimentally assayed collections. We first measure within-target ranking and raw probability error. We then hold out entire targets to test whether mappings learned elsewhere transfer without local labels, and separately quantify how much local intercept and slope updates reduce error. A frozen mapping is evaluated in a distinct external cohort. Finally, we ask how raw scores and locally updated probabilities affect candidate selection and experimentally confirmed hits when pilot assays consume the same finite budget as subsequent testing. Exploratory analyses examine whether Pinc or alternative predictor pipelines resolve these limitations. The aim is not to tune a universal mapping or propose a superior shortlisting algorithm, but to determine which probability and experimental-allocation claims are supported by each layer of evidence.

## Methods

### Study design and provenance

We conducted a retrospective, multi-target evaluation comprising a frozen Proteinbase main cohort, target-blocked internal–external validation, a separately frozen external evaluation without destination-target labels, and exploratory extensions. The main source was Proteinbase dataset version `v2026-11-01` (Hugging Face dataset `adaptyvbio/proteinbase`, revision `2d2825070da5ba77e0c7e1a0942146e8307e9282`); the local CSV was saved on 5 November 2026 and accessed for analysis on 7 November 2026 (Adaptyv Bio, 2026). The external source was the Claude protein binder design data release v1.0 (revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`), which contained experimentally assayed designs, cofolding predictions, and outcome records (The protein binder design data release authors, 2026). Analyses followed frozen protocols defining cohort construction, outcomes, folds, model forms, references, and sensitivity analyses before the corresponding performance results were examined. No new wet-lab experiments were performed. Models for explicitly labelled post hoc extensions were fitted under separately frozen protocols; the original analysis outputs were not overwritten. Frozen analysis tables, input/output hashes, software versions, and analysis-to-display mappings are designated for Supplementary Table S15.

### Main-cohort construction and analysis unit

The broader main-source audit contained 4,280 designs and 4,943 design–target pairs, of which 2,893 had a released Boltz-2 ipTM value. The frozen primary analysis retained 2,084 eligible pairs across six targets: AAH-I, VEGF110, Nipah G, PD-L1, SARS-CoV-2 receptor-binding domain (RBD), and TREM2 (Table 1). One row represented one design–target pair with one prespecified primary assay. Every design occurred once in this primary table, so no repeated cross-target structure was available for a design-level random intercept. Sequence near-duplicates and campaign/source dependence were instead addressed through clustering, blocked resampling, and sensitivity descriptions. All design classes were eligible; the cohort was not restricted to confirmed de novo minibinders. Miniprotein-only and identifiable de novo-source analyses were treated as scope sensitivity analyses rather than used to relabel the full cohort.

### Primary assays and outcome definitions

Primary assays were selected by measurement coverage, not by which assay yielded stronger prediction. BLI was the primary assay for five targets. AAH-I used SPR because all 133 designs had SPR outcomes, whereas only 18 had selectively acquired BLI confirmations; all 18 were BLI-positive and also SPR-positive. This selective confirmation structure precludes an assay-robustness claim and motivates separate reporting of AAH-I.

Exact records were deduplicated within design × target × primary assay, without assuming that distinct timestamps represented independent experimental replicates. The primary outcome was majority binding after deduplication, including Weak calls. The strict outcome was the subset of primary positives whose highest observed strength was Medium or Strong. The main cohort contained 748 primary and 678 strict positives. Alternative label rules were audited to expose definition sensitivity, not fitted as a family of parallel models from which favorable results could be selected.

### External cohort and transport conditions

The external release contained 1,440 unique designs across 16 targets. We excluded 120 Mature GDF-8 designs because antigen aggregation and nonspecific sensor binding made their wet-lab outcomes uninterpretable, leaving 1,320 designs with integrated outcomes. The frozen eligibility rule required at least 80 usable pairs and at least 15 integrated positives per target. Eight targets qualified, each contributing 90 pairs, yielding 720 pairs with complete primary scores and Adaptyv outcomes. The primary external outcome, Adaptyv Weak+, comprised 287 positives; Medium+ comprised 254. Four targets were absent from the main cohort, and four were new campaigns for represented biological targets (Table 1). The latter stratification does not assume identical target constructs, epitopes, or assay systems. Using Adaptyv Weak+ rather than integrated positives for the 15-event gate selected the same eight targets; this equivalence documents the inclusion rule but does not remove event-based selection effects.

The confirmatory cohort covered primary-outcome prevalences of 21.1%–72.2%. Six low-event targets (521 pairs) were therefore retained outside that cohort as a descriptive stress test without destination-target labels, without refitting or updating the mapping. MBP had zero positives, making its target-oracle Brier reference and calibration slope non-evaluable; BBF-14 and 15-PGDH each had one positive and were not used for substantive slope inference. The high-prevalence extreme of the main cohort remained uncovered externally.

External predictions used Boltz 2.2.1 with five seeds, no templates, an unpaired target MSA, a single binder sequence, one diffusion sample per seed, and three recycles. For each design, the seed with maximum ipSAE_min was selected and its iptm_pae value used; ties were resolved by the smallest numerical seed label. This was not a maximum-ipTM rule. The main cohort contained one released Boltz-2 ipTM value per pair without complete seed metadata. External transport thus evaluated the released score plus frozen mapping under jointly changed target/campaign and score-generation conditions, not a pure target effect under an identical protocol.

The external designs were not generated by one pipeline. Two high-level design-model categories were Mythos Preview (480 designs) and Opus 4.8 (240). Nine underlying generators were represented: PXDesign (178), RFdiffusion3 (138), Genie3 (91), FreeBindCraft/BindCraft (85), BoltzGen (84), Proteina-Complexa (69), RFdiffusion (64), FoldCraft (10), and BoltzDesign1 (1). Generator, target, campaign, and design model were associated, so their independent causal effects were not identifiable. Adaptyv and Twist outcomes were treated as alternative outcome views, not independent replicates. Differences in expression, valency, antigen format, and assay setup also prevented a pure instrument-effect interpretation.

### Relationship to previous data and sequence dependence

Exact-sequence matching was distinguished from near-sequence similarity. Against the accessible Overath Hugging Face mirror, the full main-source audit matched 345 of 4,280 designs (8.1%) and 342 of 4,943 pairs (6.9%); the frozen Boltz-2 main analysis had zero exact-sequence matches. The mirror contained 3,669 records rather than the 3,766 reported by Overath et al. (2025), leaving a 97-record provenance discrepancy. We therefore did not describe the observed percentages as lower bounds for the complete original dataset. Maximum sequence-identity distributions were retained in addition to thresholded overlap indicators. The measuredAt field was not used to claim temporal independence or a prospective post-Overath cohort.

For the external release, exact overlap was zero with both the 4,280 main-source sequences and the accessible Overath mirror. All-pairs local Smith–Waterman comparisons using Parasail 1.3.4 found no links at either 80% or 90% identity when query and reference coverage were both required to be at least 80%. These audits support sequence non-overlap under the stated rules, not independence of generative-model training data, structural dissimilarity, or epitope independence.

### Discrimination and raw probability error

Within-target discrimination was quantified using the area under the receiver operating characteristic curve (AUC) and average precision (AP). AP was interpreted relative to observed target prevalence. Raw ipTM was deliberately evaluated as a candidate experimental-success probability to test that stronger interpretation; this stress test does not imply that the model developers defined ipTM as a wet-lab success probability (Passaro et al., 2025).

For outcome \(y_i\) and candidate probability \(\hat p_i\), Brier score was \(n^{-1}\sum_i(\hat p_i-y_i)^2\); logarithmic loss was also calculated, with lower values preferable for both. Brier skill score was defined as \(1-BS_{model}/BS_{reference}\), with the reference named for every analysis. Target-constant prevalence references were used for raw within-target diagnostics. Transport references used source-target constant prevalence, target-oracle references used held-out-target prevalence for diagnosis only, and adaptation-train references used prevalence in the corresponding adaptation-training folds. Target-oracle references are not deployable comparators. Target-equal macro summaries were primary, while pair-weighted summaries assessed sensitivity to the much larger Nipah G contribution.

### Calibration heterogeneity

For score \(s\), calibration diagnostics used \(x=\operatorname{logit}\{\operatorname{clip}(s)\}\) with frozen clipping constants. Target-specific logistic models estimated the joint intercept and slope on this scale, and reliability plots used fixed intervals of [0,0.2), [0.2,0.4), [0.4,0.6), [0.6,0.8), [0.8,0.9), [0.9,0.95), and [0.95,1]. A hierarchical logistic model with target-specific random intercepts and slopes, \(\operatorname{logit}\{P(y_i=1)\}=\alpha_{t[i]}+\beta_{t[i]}x_i\), described between-target heterogeneity through partial pooling (Van Calster et al., 2016). It was not evaluated as a mapping for an unseen target. All four primary/strict six-/five-target models converged with positive-definite information matrices. Adaptive Gauss–Hermite integration was checked using 9 and 11 nodes; the maximum target-estimate difference was \(2.71\times10^{-6}\). Target-protein covariates were not added because six target-level contexts and confounding by assay and source could not identify their effects.

### Target-held-out validation and local adaptation

Leave-one-target-out (LOTO) analysis asked whether a mapping learned from other targets could predict experimental outcomes for a completely held-out target without using its labels, while respecting target and sequence structure (Roberts et al., 2017). For each held-out target, only source targets were used for fitting, transformation, and tuning, and same-cluster training rows were removed where required. The source model included target-specific intercepts and a common slope on clipped logit-ipTM. Each source target received equal total weight in the primary analysis, with equal within-target pair weights; pair weighting was a sensitivity analysis. The unadapted mapping (U) used the common source slope and the weighted mean source intercept without labels from the held-out target.

Intercept adaptation (I) retained the U linear predictor as an offset and estimated only a local intercept, thereby addressing baseline-level mismatch without changing the score relationship. Full adaptation (F) estimated an intercept and one multiplicative slope on the U linear predictor. Predictive I/F used five-fold, 90%-identity sequence-cluster-separated out-of-fold evaluation within each target. Intercepts were unpenalized, slope adaptation was shrunk toward no change from U, and regularization was selected within adaptation-training folds using nested grouped validation. This prevented tuning on the rows used for evaluation (Varma and Simon, 2006). Descriptive and oracle adaptations were kept separate from predictive claims. A three-knot restricted cubic spline on clipped logit-ipTM was the sole prespecified nonlinear sensitivity analysis, with knots estimated only from source training data; its results did not trigger further model searching.

A post hoc diagnostic quantified how much of the observed probability error was addressable by local intercept updating. Using the already evaluated outer-test rows, we calculated Brier(I)−Brier(U), Brier(F)−Brier(I), corresponding log-loss differences, and the relative Brier reduction (Brier(U)−Brier(I))/Brier(U). Brier differences were primary, and target-level results were averaged with equal target weights. The five-target sensitivity analysis used the separately refitted primary-five results rather than deleting AAH-I from six-target fits. No model was refitted for this diagnostic, and no new confidence intervals or significance tests were introduced. Relative error reduction was not interpreted as a causal proportion explained by baseline prevalence.

### Frozen external evaluation and post hoc adaptation

For external U, the deployment mapping was fitted once using all six main-cohort targets with target-equal source weighting and frozen before external performance was examined. It was applied unchanged to the released external score. The six LOTO models were not averaged, external prevalence was not used to update the intercept, and no score rule was selected after seeing external performance. The primary analysis comprised the eight eligible targets, with prespecified strata for new targets and for new campaigns involving represented biological targets, together with outcome and score/seed sensitivity analyses and source-weight sensitivity analyses.

External I/F was added after the external U results had been examined and is therefore post hoc and exploratory. It used the frozen U logit as an offset, within-target five-fold sequence-cluster-separated out-of-fold predictions, and three-fold grouped tuning inside adaptation-training folds. These label-dependent results were evaluated on the same outer rows as U but were not combined with the confirmatory evidence obtained without destination-target labels.

### Experimental-selection analyses

Decision analyses were designed around explicit experimental rules rather than inferred from AUC alone. Fixed proportional budgets selected the top 5%, 10%, and 20% by raw ipTM within each target. Precision, recall, enrichment over same-target random selection, and absolute hit gains were reported. Raw operating points at ipTM 0.80, 0.90, and 0.95 described selected-set size and observed precision; they were not interpreted as success probabilities or recommended cutoffs. Empty selections retained zero recall and non-evaluable precision, enrichment, and cost per hit.

A post hoc exploratory retrospective shared-budget simulation examined whether representative local pilots could improve allocation across targets. Total budgets \(B=120\) or 300 assays and pilot sizes \(m=10\), 20, or 40 per target were evaluated where the total pilot did not exceed the budget. Strategies were raw-ipTM allocation, equal target budgets with ipTM ranking, random allocation, pilot-I, and pilot-F. Pilot assays counted against the total budget, pilot-identified hits counted toward total hits, and already measured pairs were not selected again. Pilots were sampled by sequence cluster, with related rows removed from subsequent selection. Small-pilot models used prespecified Jeffreys smoothing and fixed ridge regularization and were distinct from the larger-fold I/F models in the transport analysis. Each held-out-target U served only as a local adaptation starting point; different LOTO mappings were not presented as one deployed mapping. The primary six-target analysis used 2,000 simulation replicates, and five-target and strict-outcome analyses used 1,000. Simulation quantiles describe variation under the observed pools and sampling mechanism, not confidence intervals for future targets.

### Exploratory Pinc analyses

Pinc analyses were prespecified exploratory incremental evaluations. All comparators used the same 2,074 complete-case pairs, frozen folds, held-out rows, and weighting schemes; ten PD-L1 pairs lacking Pinc were not imputed. Raw-score discrimination and probability error were compared on identical rows. For LOTO U, ipTM-only, Pinc-only, and combined models were refitted on those rows. The combined model used logit-transformed scores with ridge-penalized common slopes, and feature scaling and tuning remained inside source training data. For I/F, only the intercept or intercept plus one calibration slope on the combined linear predictor was updated; relative feature weights were not silently re-estimated.

The symmetric Pinc score was calculated using the official implementation. Batch processing used official legacy R for throughput after small-sample agreement with the compiled official C implementation; this accepted toolchain deviation, implementation version, and script hashes were recorded in the provenance files. Chain-specific Pinc1/Pinc2 interpretation was not attempted because released PAE conversion metadata were insufficient to close chain orientation.

### Random-forest model-complexity sensitivity

The random-forest analysis was added after the original Pinc results had been examined and was therefore post hoc and exploratory. A two-by-two comparison crossed pooled ridge logistic regression and random forests with ipTM alone or ipTM plus Pinc. All four models were refitted independently for the primary-six, primary-five, and strict-six analysis sets under target-equal and source-pair weighting, using identical complete-case rows and frozen outer folds. Source weights were recalculated within each current training row set after complete-case filtering and sequence-cluster exclusions. Unlike the original target-intercept RQ5 models, the pooled logistic comparators used one global intercept, allowing matched model-family comparisons without simultaneously changing the source-target intercept structure. The earlier target-intercept models were therefore retained as descriptive anchors rather than pure model-family comparators.

Random forests used scikit-learn 1.9.0 with 500 trees and a prespecified grid over maximum depth {2, 4, unrestricted}, minimum samples per leaf {10, 30}, and maximum features {1, 2}, with duplicate one-feature candidates removed. Hyperparameters were selected using leave-one-source-target-out logarithmic loss within each outer source-training set, with sequence-cluster exclusions and deterministic seeds. U used no destination-target labels. I and F applied the same frozen local intercept or intercept-plus-slope updates to each base model on identical outer-test rows. A prespecified directional-consistency screen required agreement across metrics, targets, analysis sets, and weighting schemes; it was a descriptive screen, not a significance test. No additional algorithm families were tried, and no independent external Pinc cohort was available. Complete parameters, candidate-level tuning, evaluation metrics, and provenance are assigned to Supplementary Table S16.

### Paired cross-predictor description

Ten complete predictor pipelines were compared on the same 720 external designs and wet-lab outcomes: Boltz-2, AFM3, Chai-1, EF2-fast, EF2-full, PTXv2, OF3, ODDE, RF3, and af3of3. For each predictor, the released run/model with maximum ipSAE_min supplied its corresponding iptm_pae under the frozen rule. The analysis described score shifts, within-target Spearman agreement, AUC/AP, and raw Brier/logarithmic loss under Weak+ and Medium+ outcomes. Predictors were treated as complete pipelines because implementations, MSA construction, recycles, and sampling could differ. af3of3 used AlphaFold 3 architecture code with OpenFold3 weights and was not an AlphaFold 3-weights result. No U/I/F mapping was fitted across predictors, no superiority test was performed, and this analysis did not constitute cross-model validation in the Proteinbase main cohort.

### Uncertainty, sensitivity analyses, and reproducibility

Unless otherwise stated, uncertainty used at least 2,000 valid within-target sequence-cluster bootstrap replicates. Paired external-predictor analyses used design-level paired bootstrap within target. These intervals quantify conditional precision in the observed targets, not sampling uncertainty over a population of future targets. Core sensitivity analyses comprised strict outcomes, exclusion and complete refitting without AAH-I, source pair weighting, an 80%-identity cluster rule, absolute budgets, miniprotein-only scope, external outcome/score definitions, and campaign/source descriptions. Counts of improving targets were descriptive, not significance tests or estimated success probabilities for future targets. The reproducibility package will provide code, frozen derived tables where licensing permits, seeds, versions, and input/output hashes before submission.

## Results

### The two cohorts represent different experimental settings

The frozen main cohort comprised 2,084 unique design–target pairs across six targets, with 748 primary and 678 strict binding-positive outcomes (Table 1; Fig. 1). Target prevalence ranged from 9.8% for Nipah G to 84.3% for SARS-CoV-2 RBD. AAH-I differed structurally from the other settings: SPR provided the primary outcome for all 133 designs, whereas 18 selectively acquired BLI confirmations were all positive and were also SPR-positive. AAH-I therefore required a prespecified sensitivity analysis rather than interpretation as an assay comparison.

The external confirmatory cohort contained 720 pairs across eight targets, including 287 Adaptyv Weak+ and 254 Medium+ positives. Four targets were new relative to the main cohort and four were new campaigns for represented biological targets. Its five-seed, best-ipSAE-seed scoring rule differed from the single released score available for the main cohort. The external designs also spanned two high-level design models and nine underlying generators. These features provide a changed deployment context, but they do not isolate target effects from campaign, generator, or protocol effects. A further 521 pairs from six low-event targets were evaluated separately as a descriptive stress test.

### ipTM contains target-dependent ranking information

Within targets, ipTM contained information about experimental binding, but the strength of discrimination varied substantially (Table 2; Fig. 2). The target-equal mean AUC was 0.675 (95% sequence-cluster bootstrap interval 0.642–0.705), compared with 0.629 (0.593–0.667) when targets were weighted by pair count. AAH-I had an AUC of 0.950 (0.902–0.982), whereas the other five targets ranged from 0.565 to 0.689; the PD-L1 interval included 0.5. Excluding AAH-I reduced the target-equal mean AUC to 0.620 (0.582–0.655). Strict-outcome analyses gave corresponding six- and five-target means of 0.677 and 0.622, showing that the principal ranking pattern was not created by inclusion of Weak calls.

Average precision required the same target-specific interpretation. SARS-CoV-2 RBD had an AP of 0.917, but its baseline prevalence was already 0.843. Nipah G had an AP of 0.125 against a prevalence of 0.098. These values do not make the score uniformly uninformative or uniformly useful: they show why a pooled performance number cannot describe the ranking task without its target context. The influence of AAH-I further demonstrates that even a target-equal average can be affected strongly by one unusual setting.

### Ranking signal does not make raw ipTM a binding probability

Treating raw ipTM as an experimental-success probability produced substantial error (Table 2; Fig. 3). The six-target, target-equal Brier score was 0.416, logarithmic loss was 1.222, and BSS versus target-specific constant-prevalence references was −1.417. Pair weighting increased the corresponding values to 0.474, 1.386, and −2.397, respectively. Excluding AAH-I gave a target-equal Brier score of 0.374 and BSS of −1.053; the strict six-target outcome gave 0.443 and −1.386. Neither sensitivity analysis supported direct probability interpretation.

Nipah G illustrates why ranking and calibration must be separated. Its AUC was 0.603, yet only 101 of 1,028 designs were primary-positive. Median ipTM was 0.910 among positives and 0.895 among negatives, both far above the observed prevalence of 0.098. The resulting raw Brier score was 0.630 and BSS was −6.110. Thus, a modest ability to order candidates coexisted with severe overstatement of absolute binding frequency.

The hierarchical model likewise indicated target-dependent calibration. Between-target standard deviations were 3.021 for intercepts and 0.824 for slopes on the clipped logit-ipTM scale. Excluding AAH-I reduced the slope standard deviation to 0.167, indicating that this target influenced estimated slope heterogeneity; it did not establish a common true slope among the other five targets. These are failures of a proposed probability interpretation of a structural-confidence score, not evidence that ipTM lacks value for its original structural or ranking purposes.

### Probability mappings are setting-dependent, and local measurements reduce error

#### Applying mappings to a completely held-out target

When each target was excluded in turn, mappings fitted to the other targets did not provide reliable probability estimates for the excluded target (Table 3; Fig. 4). The target-equal Brier score was 0.3086, logarithmic loss was 0.8643, and BSS against the source-prevalence reference was −0.0809, with improvement in only two of six targets. BSS against the held-out-target prevalence reference was −1.1205 and negative in every target. The latter reference is diagnostic rather than deployable because it uses the held-out target's observed prevalence. The five-target analysis excluding AAH-I had BSS-transport −0.0921, positive in only one of five targets. Strict-outcome and 80%-identity cluster sensitivity analyses gave six-target values of −0.0670 and −0.0804, respectively.

Restricted cubic splines improved the six-target BSS-transport from −0.0809 to −0.0240, but the mean remained negative. Allowing a nonlinear score relationship therefore did not, by itself, resolve application to a new target.

#### Local measurements improved probability estimates

Updating the overall predicted level with local measurements reduced out-of-fold Brier in all six targets, with a mean intercept-minus-unadapted difference of −0.1393 (Table 3). Mean performance against the corresponding adaptation-training prevalence references remained modest: BSS was 0.0153 and was negative for PD-L1 and SARS-CoV-2 RBD. Thus, correcting the overall probability level did not guarantee reliable design-level probabilities in every setting.

The post hoc intercept-error diagnostic quantified the size of this intercept-addressable reduction on the same outer-test results. Intercept updating reduced target-equal macro Brier from 0.3086 to 0.1693, a 45.1% relative reduction. In the separately refitted five-target sensitivity analysis, Brier decreased from 0.3140 to 0.1793, a 42.9% reduction, with improvement in five of five targets. Target-level relative reductions in the six-target analysis ranged from 6.5% to 71.7%, showing that a consistent direction did not imply a uniform benefit. These values quantify error reduced by an intercept update; they do not estimate the causal contribution of baseline prevalence or establish it as the largest mechanism of failure.

Additional slope updating had smaller and inconsistent effects. The mean F-minus-I Brier difference was −0.0076, with four of six targets improving and two worsening. For logarithmic loss, three targets improved and three worsened. After AAH-I was excluded and all mappings were refitted, the mean F-minus-I Brier difference changed to +0.0014, with two of five targets improving. The contrast between consistent intercept gains and unstable incremental slope gains argues against treating local recalibration as a single, uniformly effective procedure.

#### Evaluation in the separate external cohort

The frozen external U mapping also showed context-dependent performance (Table 3; Fig. 4). Across eight targets, target-equal Brier was 0.2590, logarithmic loss was 0.7120, BSS-transport was −0.0687, and BSS-target-oracle was −0.3211. The four new targets had mean BSS-transport −0.2184. The four new campaigns for represented biological targets had mean BSS-transport +0.0809, yet their BSS-target-oracle remained −0.1014. Improvement over a source-prevalence rule in some campaigns therefore did not establish useful probability prediction once local prevalence was taken into account.

The stratum averages conceal markedly different target-level behavior. Among the new targets, only IL-7Ra had positive values for both BSS-transport and BSS-target-oracle; its mean prediction of 0.519 was close to its observed prevalence of 0.522. BHRF1, RBX1, and TrkA instead showed probability overstatement, with mean predictions of 0.643, 0.612, and 0.524 against prevalences of 0.222, 0.256, and 0.211, respectively. Among new campaigns for represented biological targets, TREM2 improved over both references, whereas PD-L1 was close to the transport reference and slightly worse than its target-oracle reference. These contrasts are more informative about deployment than a single pooled external score because they show where the frozen mapping was directionally helpful and where it was not.

Prespecified sensitivity analyses did not yield a uniformly calibrated alternative. BSS-transport was negative under Medium+, integrated, Twist QC, and Twist permissive outcomes, ranging from −0.1011 to −0.0388. Using the five-seed mean, median, or fixed seed 0 reduced Brier to 0.2438–0.2471 and improved BSS-transport to −0.0177 to −0.0038, but BSS-target-oracle remained negative. Source-pair weighting produced BSS-transport +0.0266 while BSS-target-oracle remained −0.1417. Score summarization and source weighting can therefore alter error magnitude without establishing a stable probability mapping.

The confirmatory external cohort covered prevalences from 21.1% to 72.2%, not the main cohort's 9.8% and 84.3% extremes. In the separate low-event stress test, mean predicted probabilities exceeded observed prevalence for all six targets. MBP's zero positives made target-oracle BSS and calibration-slope estimation non-evaluable, and the one-event targets did not support slope inference. These analyses extend the descriptive evidence toward low prevalence without converting it into confirmatory coverage.

Post hoc external I/F analyses provided a separate, label-dependent check. Target-equal Brier values were 0.2590 for U, 0.1918 for I, and 0.1861 for F. Intercept adaptation reduced Brier in seven of eight targets, while F improved on I only in PD-L1 and TREM2. No new target had lower Brier under F than under I. Because these analyses were added after external U was examined and used local labels in out-of-fold adaptation, they do not alter the confirmatory result obtained without destination-target labels. Absolute Brier values were also not compared directly across main and external cohorts with different outcomes and prevalence distributions.

### Ranking and calibration have different consequences for experimental selection

Within-target ranking produced enrichment at some targets and budgets but not others (Fig. 5). AAH-I, VEGF110, SARS-CoV-2 RBD, and TREM2 showed enrichment at multiple proportional budgets. Nipah G top-5% selection had enrichment of 0.783 relative to random selection, while PD-L1 enrichment point estimates were no greater than 1 at all three budgets. SARS-CoV-2 RBD illustrates the difference between precision and enrichment: selected-set precision was high, but baseline prevalence was already 0.843, limiting possible improvement over random selection.

Common raw-score thresholds did not define common experimental operating points. An ipTM threshold of 0.80 selected 786 Nipah G designs and all 110 TREM2 designs. A threshold of 0.95 selected no VEGF110 designs. Empty, complete, and very large selected sets are incompatible with interpreting 0.80, 0.90, or 0.95 as portable experimental-success rules, even though such numerical cutoffs may be convenient within a particular campaign.

We next asked whether the intercept-addressable error identified above translated into more experimentally confirmed hits after charging pilot measurements to the same budget (Table 4; Fig. 5). This post hoc exploratory retrospective simulation evaluated net experimental yield, rather than assuming that better probability estimates automatically save assays. With \(B=300\) and \(m=10\) per target, pilot-I yielded 225.08 mean total hits versus 207.00 for raw-ipTM allocation, an I-minus-raw difference of +18.08. The 95% simulation interval for the paired difference was −1 to +29, and 97.2% of simulation differences were positive. The direction persisted in five-target and strict-outcome analyses, with mean differences of +19.87 and +22.66, respectively, although the five-target interval still crossed zero.

That favorable scenario was not general. With \(B=120\) and \(m=10\), the mean I-minus-raw difference was −3.39. At \(B=300\), increasing the pilot to \(m=20\) or 40 changed the difference to −0.81 or −38.64. Pilot-F did not add a stable benefit over pilot-I. Pilot measurements can therefore improve allocation among targets in some observed settings by correcting budget allocation, but their assays consume the same limited budget. These retrospective simulations neither prescribe a universal pilot size nor estimate prospective assay savings.

### Pinc does not provide a stable incremental remedy in complete-case analyses

The complete-case Pinc cohort comprised 2,074 pairs and 740 primary positives; all ten missing Pinc values were in PD-L1. On identical rows, target-equal AUC was 0.6366 for Pinc and 0.6770 for ipTM, and AP was 0.5533 and 0.6132, respectively (Fig. 6). Pinc had slightly higher AUC only in Nipah G and SARS-CoV-2 RBD; its AUC was lower in AAH-I, VEGF110, PD-L1, and TREM2, with paired sequence-cluster bootstrap intervals for those four contrasts not crossing zero. Raw Pinc had lower Brier than raw ipTM (0.3834 versus 0.4154), but both scores had negative target-constant BSS in every target. Lower raw error therefore did not amount to achieved calibration.

Under target-equal U weighting, the combined ipTM+Pinc model had Brier 0.2949 versus 0.3139 for ipTM-only and BSS-transport −0.0381 versus −0.0991. It did not improve AUC or AP, however, and mean BSS-transport remained negative. The Brier improvement also varied by context: the combined-model difference relative to ipTM was −0.1161 in SARS-CoV-2 RBD but +0.0963 in PD-L1. Under source-pair weighting, the Brier/BSS advantage reversed across analysis sets. Intercept adaptation remained the dominant consistent gain, with little stable increment from additional slope adaptation. These analyses do not establish that Pinc is generally inferior; they show no stable incremental solution under the evaluated rows, folds, and weighting schemes.

### Random forests alter probability error but do not consistently outperform pooled logistic regression on held-out targets

To examine whether the unstable Pinc increment could be explained by the linear functional form, we performed a post hoc two-by-two sensitivity analysis comparing pooled ridge logistic regression and random forests using ipTM alone or ipTM plus Pinc. All models used the same complete-case pairs, outer folds, source weights, and adaptation rules. The pooled logistic comparators omitted source-target intercepts so that matched comparisons varied model family without simultaneously changing the intercept structure.

Under six-target, target-equal leave-one-target-out evaluation without labels from the held-out target, the ipTM-only random forest lowered Brier from 0.2826 to 0.2624 and log loss from 0.7590 to 0.7254 relative to its logistic comparator, while BSS-transport increased from −0.0001 to 0.0690. Brier was lower in five of six targets, although the difference in SARS-CoV-2 RBD was negligible. These probability-error gains were accompanied by poorer discrimination: AUC decreased from 0.6515 to 0.5684 and AP from 0.6002 to 0.5536. BSS-target-oracle also remained negative (−0.7152). Thus, lower average probability error did not establish improved ranking or transportable calibration.

The two-feature random forest likewise did not yield a stable improvement. Relative to the two-feature logistic model, it lowered Brier and log loss in the primary six-target analysis, but this direction reversed after refitting on the five-target analysis set (Brier difference, +0.0466; log-loss difference, +0.1725; BSS-transport difference, −0.1294). Within the random-forest family, adding Pinc worsened macro-averaged Brier, log loss, and BSS-transport in all six analysis-set-by-weighting configurations. In the primary analysis, Brier increased from 0.2624 for ipTM alone to 0.2894 for ipTM plus Pinc, log loss increased from 0.7254 to 0.7870, and BSS-transport decreased from 0.0690 to −0.0240. In the primary analysis, the full random-forest-plus-Pinc model also failed to outperform the simplest pooled logistic ipTM model. With destination-target labels, the fully adapted two-feature random forest achieved Brier 0.1656 and BSS-target-oracle 0.0452 in the new pooled comparison; this local-label result does not alter the interpretation for prediction without destination-target labels and is not a pure model-family contrast with the previous target-intercept analysis. Overall, these exploratory findings show metric- and target-dependent changes under the evaluated random-forest grid, with no consistent improvement over pooled logistic regression when predicting held-out targets without their experimental labels, and no stable incremental contribution from Pinc (Supplementary Table S16).

### Predictor pipelines differ in both score scale and design ordering

The paired external-cohort description compared ten complete predictor pipelines on the same 720 designs (Fig. 7). Mean within-target Spearman correlations with Boltz-2 ranged from 0.4561 to 0.6146, so pipeline differences were not limited to a constant numerical shift. RF3, OF3, and Chai-1 had target-equal median ipTM shifts of −0.0905, −0.0716, and −0.0513 relative to Boltz-2, respectively, while also changing candidate ordering. Target-equal mean AUC ranges were 0.6429–0.7182 for Weak+ and 0.6341–0.7239 for Medium+. No predictor had the highest point-estimate AUC across all targets; Weak+ maxima involved five predictors and Medium+ maxima six. RF3 had positive Weak+ AUC differences relative to Boltz-2 in seven of eight targets, but this directional pattern did not make it uniformly best across targets or outcomes.

Some pipelines had lower raw Brier scores, but this did not establish calibration: Boltz-2, Chai-1, and OF3 had mean ipTM above observed prevalence in all eight targets under Weak+, and RF3 did so in seven. Because implementations, MSA construction, recycles, and sampling could differ, these observations are not attributable to model weights alone. They describe predictor-related context within the external cohort and do not test cross-predictor calibration or transport.

## Discussion


This study separates several tasks that are often encountered together in a computational binder-design workflow. Boltz-2 ipTM contained moderate but target-dependent information for ordering candidates within a target. Its raw numerical value was nevertheless a poor estimate of experimental-binding probability, and a probability mapping learned from other targets did not transfer reliably without destination-target measurements. Local intercept updating helped more consistently than additional slope updating, while the experimental value of a pilot depended on its size and the remaining assay budget. These findings are compatible: a score can retain useful ordering information while lacking the absolute scale needed for probability-based selection or allocation.

The distinction matters because modern design pipelines connect generative models, sequence design, complex prediction, and experimental testing rather than treating them as isolated tools (Cao et al., 2022; Dauparas et al., 2022; Watson et al., 2023). BindCraft illustrates the practical value of integrating design and prediction, while broader reviews describe rapid growth in the scope of de novo design (Pacesa et al., 2025; Yang et al., 2026). Our results do not challenge those design capabilities or propose a better generator. They concern the narrower decision of how a released confidence value should be interpreted when it controls experimental work.

The findings also complement, rather than duplicate, previous analyses of binder success. Overath et al. (2025) established the value of large experimentally characterised collections and investigated confidence, structural, and energetic predictors across targets. Badonyi and Toth-Petroczy (2026) showed why interface scores themselves should be compared. The present study instead asks whether a numerical score can be read as an absolute experimental-success probability when applied to a target with no local labels, and whether local measurements improve that interpretation. Lee et al. (2026) addressed target-held-out top-K shortlisting; our fixed-budget analyses likewise concern selection, but the probability mappings and pilot-cost analysis answer different questions from optimising the order of a shortlist.

The post hoc intercept diagnostic sharpens the interpretation of local adaptation. Intercept-only updating reduced target-equal Brier by 45.1% in the six-target analysis and 42.9% after complete refitting without AAH-I, with improvement in every analysed target. This establishes a substantial intercept-addressable error component in these data. It does not establish baseline prevalence as the largest causal mechanism, because target, assay, campaign, design source, and score-generation conditions vary together. Calibration research distinguishes mean-level adjustment, calibration slope, and more complex forms of miscalibration (Van Calster et al., 2016; Van Calster et al., 2019). Consistent with that distinction, residual error remained after intercept updating, and further slope updating was unstable. The appropriate conclusion is therefore not that a single local base-rate estimate solves transport, but that representative local measurements can correct an important part of the observed error.

The shared-budget simulation adds the necessary experimental qualification. A small pilot produced more total hits than raw-ipTM allocation in the favourable B=300, m=10 scenario, but larger pilots or a smaller total budget reduced or reversed that advantage. Better probability estimates therefore cannot be equated automatically with assay savings: the measurements used for adaptation consume part of the same budget. This is why we report net experimentally confirmed hits under explicit budgets rather than recommending a universal pilot size. The result is a conditional decision boundary for the observed candidate pools, not evidence of prospective monetary savings.

The secondary extensions reinforce the need to test each proposed remedy. Adding Pinc did not yield a stable increment across complete-case weighting schemes, despite the rationale for comparing interface-confidence measures (Badonyi and Toth-Petroczy, 2026). The post hoc random-forest analysis produced metric- and target-dependent changes, but no consistent improvement over pooled logistic regression when predicting held-out targets without their experimental labels and no stable incremental contribution from Pinc. Alternative predictor pipelines changed both numerical scale and within-target ordering, so changing predictors introduces another context requiring validation rather than an automatic correction. Related work on monomeric design and antibody–antigen specificity likewise shows that the relationship between computational confidence and experimental or structural outcomes depends on the task (Garcia et al., 2026; Smorodina et al., 2026). Those different endpoints do not directly prove our mechanism, and our results do not show that every possible calibration mapping or predictor must fail.

Several limitations bound these conclusions. The study was retrospective, and six main targets are not six independent design or sampling cohorts. Design classes and sources varied, and sequence/campaign dependence remains relevant despite cluster-separated validation. AAH-I had only 19 positives and a selective BLI-confirmation structure; target and assay effects cannot be separated, and its exclusion is a sensitivity analysis rather than a new primary cohort. External evaluation combined four new targets and four new campaigns for represented biological targets with a different score-generation protocol; identical target constructs, epitopes, and assay systems were not assumed. Its two high-level design models and nine generators were structurally associated with targets and campaigns, each target contained only 90 pairs, and external I/F was post hoc. Confirmatory prevalence coverage was 21.1%–72.2%; low-event stress tests cannot support slope inference at zero or one event, and the main cohort's high-prevalence extreme was not reproduced externally. Pinc analyses were complete-case, with missingness concentrated in PD-L1. The random-forest extension was post hoc, limited to the evaluated parameter grid, and lacked independent external Pinc validation; its combined-feature contrasts did not identify pure interaction effects. Predictor analyses concerned complete pipelines in an upstream-filtered external cohort. Shared-budget simulations used observed candidate pools and retrospective sampling, so their intervals are not future-target confidence intervals. This diagnostic was a post hoc descriptive contrast of frozen outer-test metrics and offered no causal decomposition or future-target inference. Finally, the Overath mirror discrepancy and unresolved repository/citation metadata limit provenance statements, and temporal independence has not been established.

In the evaluated contexts, Boltz-2 ipTM is best regarded as a context-dependent, target-local ranking aid, not as a portable experimental-success probability or universal selection threshold. When absolute probabilities or shared-budget allocation matter, representative local measurements should be evaluated together with their experimental cost. This conclusion does not preclude transportable calibration with richer future data; it specifies the evidence that such a claim would require.

## Data and code availability

The main data source was Proteinbase version `v2026-11-01` (Hugging Face dataset `adaptyvbio/proteinbase`, revision `2d2825070da5ba77e0c7e1a0942146e8307e9282`), distributed under ODC-BY. The dataset card requires the attribution statement: “This work used Proteinbase by Adaptyv Bio under ODC-BY license.” The external data source was the Claude protein binder design data release v1.0 (revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`); its data and documentation are distributed under CC BY 4.0 and included scripts under MIT. Third-party structure-prediction outputs, reference sequences, vendor images, and other third-party material remain subject to their own terms. Analysis code, frozen derived tables permitted by these licenses, software versions, random seeds, and input/output hashes will be deposited in a public reproducibility package before submission [REPOSITORY AND DOI TO BE ADDED]. The package will provide acquisition instructions rather than re-host any third-party material whose terms do not permit redistribution. Until the repository link is supplied, the local frozen run manifests and accepted result tables constitute the numerical source of record for this draft.

## Declarations

**Funding:** [To be supplied.]  
**Conflicts of interest:** [To be supplied.]  
**Author contributions:** [To be supplied from actual contributions.]  
**Acknowledgements:** [To be supplied.]  
**Ethics:** No new human or animal experiments were conducted; the study reanalysed publicly released computational and experimental design records.

## References

Adaptyv Bio (2026) *Proteinbase* (version `v2026-11-01`; Hugging Face revision `2d2825070da5ba77e0c7e1a0942146e8307e9282`) [Data set]. Hugging Face. https://huggingface.co/datasets/adaptyvbio/proteinbase (local CSV saved 5 November 2026; accessed for analysis 7 November 2026; licensed under ODC-BY).

Abramson, J. et al. (2024) Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*. https://doi.org/10.1038/s41586-024-07487-w

Badonyi, M. and Toth-Petroczy, A. (2026) A simple probabilistic AlphaFold interaction score. *Protein Science*, **35**(9). https://doi.org/10.1002/pro.70760

Cao, L. et al. (2022) Design of protein-binding proteins from the target structure alone. *Nature*. https://doi.org/10.1038/s41586-022-04654-9

Dauparas, J. et al. (2022) Robust deep learning-based protein sequence design using ProteinMPNN. *Science*. https://doi.org/10.1126/science.add2187

Garcia, M., Dixit, S.M. and Rocklin, G.J. (2026) Evaluating zero-shot prediction of monomeric protein design success by AlphaFold, ESMFold, and ProteinMPNN. *Protein Science*, **35**(2). https://doi.org/10.1002/pro.70453

Jumper, J. et al. (2021) Highly accurate protein structure prediction with AlphaFold. *Nature*. https://doi.org/10.1038/s41586-021-03819-2

Lee, G., Yoo, K., Seo, J., Kim, J., Hur, K. and Choi, E. (2026) Natural-Language-Guided Generator-Agnostic Shortlisting for Protein Binder Design. *arXiv*, arXiv:2608.20755v2. https://arxiv.org/abs/2608.20755 (ICML 2026 Workshop on Generative and Agentic AI for Biology).

Overath, M.D., Rygaard, A.S.H., Jacobsen, C.P., Brasas, V., Morell, O., Sormanni, P. and Jenkins, T.P. (2025) Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders. *bioRxiv*. https://doi.org/10.1101/2025.08.14.670059

Pacesa, M. et al. (2025) One-shot design of functional protein binders with BindCraft. *Nature*. https://doi.org/10.1038/s41586-025-09429-6

Passaro, S., Corso, G., Wohlwend, J., Reveiz, M., Thaler, S., Somnath, V.R., Getz, N., Portnoi, T., Roy, J., Stark, H., Kwabi-Addo, D., Beaini, D., Jaakkola, T. and Barzilay, R. (2025) Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction. *bioRxiv*. https://doi.org/10.1101/2025.06.14.659707

Roberts, D.R., Bahn, V., Ciuti, S., Boyce, M.S., Elith, J., Guillera-Arroita, G., Hauenstein, S., Lahoz-Monfort, J.J., Schröder, B., Thuiller, W., Warton, D.I., Wintle, B.A., Hartig, F. and Dormann, C.F. (2017) Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, **40**(8), 913–929. https://doi.org/10.1111/ecog.02881

Sagadinov, S., Viliuga, V. and Elofsson, A. (2026) Calibrating AlphaFold confidence scores for interface quality. *arXiv*, arXiv:2610.17108v1. https://arxiv.org/abs/2610.17108

Smorodina, E., Ali, M., Kropivšek Brumat, K., Salicari, L., Miklavc, S., Kappassov, A., Fu, C., Sormanni, P., de Marco, A. and Greiff, V. (2026) Structural Plausibility Without Binding Specificity: Limits of AI-Based Antibody-Antigen Structure Prediction Confidence Scores. *bioRxiv*. https://doi.org/10.64898/2026.03.02.709004

Steyerberg, E.W. and Harrell, F.E. (2016) Prediction models need appropriate internal, internal–external, and external validation. *Journal of Clinical Epidemiology*, **69**, 245–247. https://doi.org/10.1016/j.jclinepi.2015.04.005

Steyerberg, E.W., Vickers, A.J., Cook, N.R., Gerds, T., Gonen, M., Obuchowski, N., Pencina, M.J. and Kattan, M.W. (2010) Assessing the Performance of Prediction Models: A Framework for Traditional and Novel Measures. *Epidemiology*, **21**(1), 128–138. https://doi.org/10.1097/EDE.0b013e3181c30fb2

The protein binder design data release authors (2026) *Protein binder design data release: de novo miniprotein binders against 16 targets with two-vendor binding data and co-folded structure predictions* (version v1.0) [Data set]. Hugging Face. https://huggingface.co/datasets/Anthropic/claude-protein-binder-design (revision `9e1b81696da46835e9e9cde9a3da976e0abc92ab`; accessed 13 November 2026; data and documentation licensed under CC BY 4.0).

Van Calster, B., McLernon, D.J., van Smeden, M., Wynants, L. and Steyerberg, E.W. (2019) Calibration: the Achilles heel of predictive analytics. *BMC Medicine*, **17**(1). https://doi.org/10.1186/s12916-019-1466-7

Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina, M.J. and Steyerberg, E.W. (2016) A calibration hierarchy for risk models was defined: from utopia to empirical data. *Journal of Clinical Epidemiology*, **74**, 167–176. https://doi.org/10.1016/j.jclinepi.2015.12.005

Varma, S. and Simon, R. (2006) Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics*, **7**, 91. https://doi.org/10.1186/1471-2105-7-91

Watson, J.L. et al. (2023) De novo design of protein structure and function with RFdiffusion. *Nature*. https://doi.org/10.1038/s41586-023-06415-8

Yang, W. et al. (2026) The past, present and future of de novo protein design. *Nature*. https://doi.org/10.1038/s41586-026-10328-7

## Main tables

### Table 1. Cohorts, assays, outcomes, and transport strata

**Panel A. Frozen Proteinbase main cohort**

| Target | Primary assay | Pairs, n | Primary positives, n (%) | Strict positives, n |
|---|---|---:|---:|---:|
| AAH-I | SPR | 133 | 19 (14.3) | 19 |
| VEGF110 | BLI | 268 | 158 (59.0) | 155 |
| Nipah G | BLI | 1,028 | 101 (9.8) | 100 |
| PD-L1 | BLI | 277 | 198 (71.5) | 159 |
| SARS-CoV-2 RBD | BLI | 268 | 226 (84.3) | 204 |
| TREM2 | BLI | 110 | 46 (41.8) | 41 |
| **Total: six targets** | — | **2,084** | **748 (35.9)** | **678** |

**Panel B. Frozen external confirmatory cohort**

| Target | Transport stratum | Pairs, n | Adaptyv Weak+ positives, n (%) | Adaptyv Medium+ positives, n |
|---|---|---:|---:|---:|
| BHRF1 | New target | 90 | 20 (22.2) | 18 |
| IL-7Ra | New target | 90 | 47 (52.2) | 36 |
| RBX1 | New target | 90 | 23 (25.6) | 20 |
| TrkA | New target | 90 | 19 (21.1) | 19 |
| Nipah G | Same biological target, new campaign | 90 | 19 (21.1) | 17 |
| PD-L1 | Same biological target, new campaign | 90 | 39 (43.3) | 35 |
| TREM2 | Same biological target, new campaign | 90 | 65 (72.2) | 62 |
| VEGF-A | Same biological target, new campaign | 90 | 55 (61.1) | 47 |
| **Total: eight targets** | Four new targets; four new campaigns | **720** | **287 (39.9)** | **254** |

Primary main-cohort outcomes include Weak calls; strict outcomes require Medium or Strong. Primary assays were chosen by measurement coverage. For AAH-I, all 133 designs had SPR outcomes; 18 selective BLI confirmations were all BLI-positive and SPR-positive. For Panel B, “Same biological target, new campaign” denotes a campaign involving a biological target represented in the main cohort; it does not assume identical target constructs, epitopes, or assay systems, including between external VEGF-A and main-cohort VEGF110. The external cohort contains Mythos Preview (480 designs) and Opus 4.8 (240), with nine underlying generators: PXDesign (178), RFdiffusion3 (138), Genie3 (91), FreeBindCraft/BindCraft (85), BoltzGen (84), Proteina-Complexa (69), RFdiffusion (64), FoldCraft (10), and BoltzDesign1 (1). Generator, target, campaign, and design model are associated and are not independently identifiable effects. The separate low-event stress test contains 521 pairs and is not included in Panel B.

### Table 2. Within-target discrimination and raw probability error under the primary main-cohort outcome

| Analysis target/set | AUC (95% interval) | AP (95% interval) | Raw Brier | BSS vs target-constant prevalence |
|---|---:|---:|---:|---:|
| AAH-I | 0.950 (0.902–0.982) | 0.601 (0.447–0.835) | 0.628 | −4.126 |
| VEGF110 | 0.602 (0.529–0.667) | 0.686 (0.608–0.766) | 0.252 | −0.040 |
| Nipah G | 0.603 (0.541–0.671) | 0.125 (0.092–0.172) | 0.630 | −6.110 |
| PD-L1 | 0.565 (0.480–0.645) | 0.732 (0.661–0.812) | 0.258 | −0.267 |
| SARS-CoV-2 RBD | 0.642 (0.564–0.717) | 0.917 (0.877–0.947) | 0.240 | −0.816 |
| TREM2 | 0.689 (0.573–0.792) | 0.618 (0.470–0.769) | 0.488 | −1.004 |
| **Six targets, target-equal** | **0.675 (0.642–0.705)** | — | **0.416** | **−1.417** |
| **Six targets, pair-weighted** | **0.629 (0.593–0.667)** | — | **0.474** | **−2.397** |
| **Five targets excluding AAH-I, target-equal** | **0.620 (0.582–0.655)** | — | **0.374** | **−1.053** |

AUC, area under the receiver operating characteristic curve; AP, average precision; BSS, Brier skill score. Intervals are 95% within-target sequence-cluster bootstrap intervals and quantify conditional precision in the observed targets, not uncertainty over future targets. Six-target target-equal logarithmic loss was 1.222, pair-weighted logarithmic loss was 1.386, and five-target logarithmic loss was 1.117. Strict-outcome target-equal AUCs were 0.677 for six targets and 0.622 for five targets; corresponding raw Brier/BSS values were 0.443/−1.386 and 0.406/−1.048. AP should be interpreted against the target prevalences in Table 1.

### Table 3. Target-held-out probability prediction and local adaptation: key target-equal summaries

**Panel A. Main-cohort six-target LOTO evaluation and predictive out-of-fold adaptation**

| Layer or contrast | Brier or ΔBrier | Log loss or Δlog loss | BSS-transport | Additional reference/result |
|---|---:|---:|---:|---|
| Unadapted mapping (U) | 0.3086 | 0.8643 | −0.0809 | BSS-target-oracle −1.1205; positive BSS-transport in 2/6 targets |
| Intercept update minus unadapted (I − U) | −0.1393 | −0.3533 | — | Brier improved in 6/6; relative macro-Brier reduction 45.1%; I BSS-adaptation-train +0.0153 |
| F − I | −0.0076 | −0.0274 | — | Brier improved in 4/6; log loss improved in 3/6 |
| Five-target I − U after complete refitting | −0.1347 | −0.3073 | — | Brier improved in 5/5; relative macro-Brier reduction 42.9% |
| Five-target F − I after complete refitting | +0.0014 | +0.0040 | — | Brier and log loss improved in 2/5 targets |

**Panel B. Frozen external U and post hoc exploratory external I/F**

| Evaluation/layer | Brier | Log loss | BSS-transport | Additional reference/result |
|---|---:|---:|---:|---|
| External U: all eight targets | 0.2590 | 0.7120 | −0.0687 | BSS-target-oracle −0.3211 |
| External U: four new targets | — | — | −0.2184 | BSS-target-oracle −0.5409 |
| External U: four new campaigns for represented biological targets | — | — | +0.0809 | BSS-target-oracle −0.1014 |
| External I: all eight targets | 0.1918 | 0.5679 | +0.2170 | BSS-adaptation-train +0.0561; Brier improved versus U in 7/8 |
| External F: all eight targets | 0.1861 | 0.5562 | +0.2391 | BSS-adaptation-train +0.0807; Brier improved versus I in 2/8 |

U is the unadapted mapping fitted without destination-target labels, I uses local measurements to update only the overall predicted level, and F additionally updates the slope. The relative error reductions are post hoc descriptive contrasts of frozen outer-test metrics, not causal proportions explained by prevalence. Differences are calculated as the named left layer minus the right layer; negative Brier/log-loss differences indicate improvement. BSS-transport uses the source-prevalence transport reference, BSS-target-oracle uses held-out-target prevalence for diagnosis only, and BSS-adaptation-train uses the corresponding adaptation-training prevalence. Main I/F and external I/F use local labels and do not constitute evidence for prediction without destination-target labels. External I/F was added after external U results were examined and is post hoc exploratory. Absolute Brier values should not be compared naively across cohorts with different outcomes and prevalence distributions. Complete per-target results and sensitivity analyses are assigned to Supplementary Tables S5, S8, and S10. The intercept-error derivation is archived in `analysis/rq4e_outputs/run_20261117_162929/`.

### Table 4. Post hoc exploratory retrospective shared-budget pilot outcomes: representative scenarios

| Total budget B | Pilot m per target | Raw ipTM hits | Equal-budget + ipTM hits | Random-allocation hits | Pilot-I hits | Pilot-F hits | Pilot I − raw hits (95% simulation interval) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 120 | 10 | 88.00 | 76.00 | 42.89 | 84.61 | 82.44 | −3.39 (−18 to +4) |
| 300 | 10 | 207.00 | 173.00 | 107.86 | 225.08 | 222.65 | +18.08 (−1 to +29) |
| 300 | 20 | 207.00 | 173.00 | 107.64 | 206.19 | 205.46 | −0.81 (−10 to +8) |
| 300 | 40 | 207.00 | 173.00 | 107.23 | 168.36 | 167.97 | −38.64 (−50 to −28) |

Values are mean total experimentally confirmed hits in the six-target main cohort under 2,000 retrospective simulation replicates. Pilot assays count against the total budget and pilot-identified hits count toward total hits. The simulation interval is the central 95% range of paired differences under the observed candidate pools and sampling mechanism, not a confidence interval for future targets. These analyses were added after the transport results and are post hoc exploratory; they do not establish prospective assay savings or a universally recommended pilot size. The complete B/m grid and sensitivity analyses are assigned to Supplementary Table S12.

## Figure files and legends

> Composite Figures 1–7 have been assembled at 300 dpi from the frozen panels. The figures below are review artwork; final journal-specific sizing and any publisher-requested format conversion remain submission-stage tasks. Exact source paths and SHA-256 hashes are recorded in `manuscript/figures/figure_source_manifest_2026-11-15.md`. Original internal panel labels are retained for source traceability.

### Figure 1. Study design and cohort flow

![Figure 1. Study design and cohort flow](figures/Figure_1_study_design.png)

The main branch follows the Proteinbase audit from 4,280 designs and 4,943 design–target pairs to 2,893 pairs with Boltz-2 ipTM and the frozen 2,084-pair/six-target primary cohort. The external branch follows 1,440 designs across 16 targets to 1,320 designs with interpretable integrated outcomes and the 720-pair/eight-target confirmatory cohort. A separate branch identifies the 521-pair/six-target low-event stress test. Labels distinguish within-target ranking, raw calibration, target-blocked U/I/F, external confirmatory U without destination-target labels, post hoc external I/F, and decision analyses. Four external targets are new and four are new campaigns for represented biological targets; external target/campaign and score-generation conditions change together.

### Figure 2. Target-dependent ipTM distributions and discrimination

![Figure 2. Target-dependent ipTM distributions and discrimination](figures/Figure_2_ipTM_discrimination.png)

(A) ipTM distributions by experimental outcome and target. (B) Within-target AUC (left) and AP (right) with 95% sequence-cluster bootstrap intervals; diamonds in the AP plot mark observed target prevalence. The main outcome includes Weak calls; n=2,084 pairs and 748 positives. AAH-I is shown explicitly because of its influence on the target-equal summary.  
*Frozen panel sources:* `analysis/rq12_outputs/RQ1_ipTM_distributions.png`; `analysis/rq12_outputs/RQ1_target_discrimination.png`.  
*Numerical source:* `analysis/rq12_outputs/rq12_target_metrics.csv`; `analysis/rq12_outputs/rq12_summary_metrics.csv`.

### Figure 3. Raw calibration and hierarchical heterogeneity

![Figure 3. Raw calibration and hierarchical heterogeneity](figures/Figure_3_raw_calibration.png)

(A) Target-specific raw-ipTM reliability curves using the frozen common bins. (B) Target-specific intercept and slope estimates before and after hierarchical partial pooling, with between-target heterogeneity summaries. The hierarchy is descriptive and is not a validated mapping for future targets. Raw-ipTM analysis treats the score as a candidate experimental-success probability solely to test that interpretation; n=2,084 pairs and 748 primary positives.  
*Frozen panel sources:* `analysis/rq12_outputs/RQ2_raw_reliability.png`; `analysis/rq12_outputs/RQ2_partial_pooling.png`.  
*Numerical source:* accepted RQ2 target and hierarchical result tables in `analysis/rq12_outputs/` and `analysis/rq12_outputs/rq12_run_manifest.json`.

### Figure 4. Prediction without destination-target labels and local adaptation

![Figure 4. Prediction without destination-target labels and local adaptation](figures/Figure_4_transport_adaptation.png)

(A) Main-cohort LOTO U performance by held-out target. (B) Main-cohort predictive out-of-fold I/F adaptation on identical outer rows. (C) Frozen external U performance by target and transport stratum. (D) External U reliability curves. U uses no destination-target labels; main I/F uses local labels and out-of-fold adaptation. External U is confirmatory for the frozen released-score-plus-mapping evaluation but combines target/campaign and score-generation changes. Post hoc external I/F is intentionally reserved for Supplementary Figures 3–4 and is not presented as confirmatory evidence here.  
*Frozen panel sources:* `analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/RQ3_U_LOTO.png`; `analysis/codex_outputs/RQ3_F_correction_2026-11-13/run_20261113_182311/RQ3_adaptation.png`; `analysis/rq3e_outputs/run_20261114_164950/RQ3E_U_external_target_metrics.png`; `analysis/rq3e_outputs/run_20261114_164950/RQ3E_U_external_reliability.png`.

### Figure 5. Experimental-selection consequences

![Figure 5. Experimental-selection consequences](figures/Figure_5_experimental_selection.png)

(A) Precision and enrichment for top 5%, 10%, and 20% within-target ipTM selections, relative to same-target random selection. (B) Selected-set sizes and observed precision at raw ipTM operating points 0.80, 0.90, and 0.95; these are not probability-derived or recommended universal cutoffs. (C) Mean total hits and paired strategy differences in the **post hoc exploratory retrospective shared-budget simulation**, including pilot costs and pilot hits. Simulation intervals describe the observed pools and sampling mechanism, not prospective deployment uncertainty.  
*Frozen panel sources:* `analysis/rq4_outputs/run_20261114_231925/RQ4A_fixed_budget.png`; `analysis/rq4_outputs/run_20261114_231925/RQ4B_raw_thresholds.png`; `analysis/rq4_outputs/run_20261114_231925/RQ4C_shared_budget.png`.  
*Numerical sources:* accepted RQ4 result CSV files in the same run directory, including `rq4c_summary.csv` and `rq4c_comparisons.csv`.

### Figure 6. Exploratory Pinc increment on common complete cases

![Figure 6. Exploratory Pinc increment on common complete cases](figures/Figure_6_Pinc_increment.png)

(A) Same-row target-specific AUC/AP for ipTM and Pinc. (B) LOTO U transport for ipTM-only, Pinc-only, and ipTM+Pinc models refitted on identical rows and folds. (C) Local-adaptation increments on identical outer-test rows. All analyses use n=2,074 complete-case pairs and 740 primary positives; ten PD-L1 pairs without Pinc are not imputed. Improvements in selected probability-error summaries do not imply restored universal transport or stable incremental discrimination.  
*Frozen panel sources:* `analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/RQ5_score_discrimination.png`; `analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/RQ5_U_transport.png`; `analysis/codex_outputs/RQ5_execution_2026-11-15/run_20261115_111720/RQ5_adaptation_increment.png`.

### Figure 7. Exploratory paired description across predictor pipelines

![Figure 7. Exploratory paired description across predictor pipelines](figures/Figure_7_predictor_pipelines.png)

(A) Target-specific ipTM shifts and within-target Spearman correlations with Boltz-2. (B) Target-specific AUC and raw Brier under Adaptyv Weak+ and Medium+ outcomes. Ten complete pipelines are compared on the same 720 external designs and outcomes; this is a paired description, not a superiority test or cross-predictor U/I/F transport analysis. Pipeline differences may involve implementation, MSA construction, recycles, and sampling, not weights alone. af3of3 uses AlphaFold 3 architecture code with OpenFold3 weights. The cohort contains two high-level design models and nine generators listed in Table 1. Lower raw Brier does not establish calibration.  
*Accepted panel sources:* `analysis/audit_outputs/RQ6_acceptance_run_20261115_151147/RQ6_figure17_score_shift.png`; `analysis/audit_outputs/RQ6_acceptance_run_20261115_151147/RQ6_figure18_discrimination_error.png`.  
*Evidence boundary:* no main-cohort cross-model validation and no claim that context dependence is unique to Boltz-2.

## Supplementary materials

The complete scientific supplement draft consists of the [supplementary guide](supplement/Supplementary_Materials_PEDS_2026-11-16.md), the [Supplementary Tables S1–S16 workbook](supplement/Supplementary_Tables_S1_S16.xlsx), and the [source and SHA-256 manifest](supplement/supplement_source_manifest_2026-11-16.md). The workbook contains Cover and Sources sheets plus 16 supplementary-table sheets with 75 Excel tables; S16 contains six logical panels. Frozen configurations and source precision are retained; source-level non-evaluable cells remain blank and are explained by status fields and panel notes. Journal-specific file conversion and public-release review remain submission-stage tasks.

Supplementary Tables S1–S16 contain, respectively: operational outcome rules and alternative-label audit; assay composition and AAH-I selective confirmation; complete main-cohort discrimination/calibration results; hierarchical model diagnostics; full LOTO U/I/F and sensitivity results; Overath exact/near-sequence audits; campaign/source and generator composition; full external U and low-event results; external outcome/score/source-weight sensitivity; post hoc external I/F; fixed-budget and raw-operating-point results; complete shared-budget simulations; Pinc complete-case analyses and implementation provenance; full predictor × target results; reproducibility metadata with analysis-to-display mappings; and random-forest model-complexity sensitivity on common complete cases.

### Supplementary figure files

- [Supplementary Figure S1](supplement/figures/Supplementary_Figure_S1.png): primary/strict and six-/five-target discrimination and raw probability error.
- [Supplementary Figure S2](supplement/figures/Supplementary_Figure_S2.png): linear versus restricted-cubic-spline prediction without destination-target labels.
- [Supplementary Figure S3](supplement/figures/Supplementary_Figure_S3.png): post hoc exploratory external U/I/F out-of-fold metrics.
- [Supplementary Figure S4](supplement/figures/Supplementary_Figure_S4.png): post hoc exploratory external U/I/F reliability.
- [Supplementary Figure S5](supplement/figures/Supplementary_Figure_S5.png): final-run nonconfirmatory low-event stress test.
- [Supplementary Figure S6](supplement/figures/Supplementary_Figure_S6.png): external outcome, score-rule, seed, and source-weight sensitivity.
- [Supplementary Figure S7](supplement/figures/Supplementary_Figure_S7.png): shared-budget pilot sensitivity across analysis sets.
- [Supplementary Figure S8](supplement/figures/Supplementary_Figure_S8.png): raw Pinc and ipTM reliability on common complete cases.
- [Supplementary Figure S9](supplement/figures/Supplementary_Figure_S9.png): Pinc increment by analysis set and source weighting.
- [Supplementary Figure S10](supplement/figures/Supplementary_Figure_S10.png): target-equal summaries and observed target ranges for ten complete predictor pipelines.
