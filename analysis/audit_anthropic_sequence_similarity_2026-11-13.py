"""Coverage-constrained nearest-neighbour audit for external sequence independence.

All query-reference pairs are aligned (no heuristic prefilter). A cluster edge requires
both directional alignment coverages >=0.80 plus identity >= threshold; this avoids
calling a short exact motif a near-duplicate. Does not fit outcome models.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json
import re
import sys
import time
import pandas as pd
import duckdb

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
DATA = ROOT / "data" / "external_target_audit_2026-11-13"
sys.path.insert(0, str(ROOT / "analysis" / "rq12_lib"))
import parasail
CIGAR_RE = re.compile(r"(\d+)([=XIDM])")
REFS = {}

def values(cigar, qlen, tlen):
    if isinstance(cigar, bytes):
        cigar = cigar.decode("ascii")
    matches = qc = tc = cols = 0
    for number, op in CIGAR_RE.findall(str(cigar)):
        n = int(number)
        cols += n
        if op == "=":
            matches += n; qc += n; tc += n
        elif op == "X":
            qc += n; tc += n
        elif op == "I":
            qc += n
        elif op == "D":
            tc += n
        else:
            raise ValueError(f"Unexpected CIGAR operation {op}")
    if not cols or not qc or not tc:
        return 0., 0., 0., 0, 0
    return matches / cols, qc / qlen, tc / tlen, cols, matches

def initialise(pb_path, ov_path):
    pb = pd.read_csv(pb_path, usecols=["sequence"])
    pb["sequence"] = pb.sequence.astype(str).str.strip().str.upper()
    con = duckdb.connect()
    # HF SeqA is binder (3,669 unique); SeqB is target (19 unique).
    ov = con.execute(f"select SeqA from '{ov_path}'").df()
    ov["sequence"] = ov.SeqA.astype(str).str.strip().str.upper()
    assert ov.sequence.nunique() == len(ov) == 3669
    for name, frame in [("proteinbase", pb), ("overath", ov)]:
        seqs = frame.sequence.drop_duplicates().tolist()
        REFS[name] = seqs

def process(item):
    qi, uuid, full_name, target, query = item
    row = {"query_index":qi, "uuid":uuid, "full_name":full_name, "target":target,
           "query_length":len(query)}
    for refname, targets in REFS.items():
        best_key = (-1., -1., -1.)
        best = None
        edge80 = edge90 = 0
        for ti, ref in enumerate(targets):
            r = parasail.sw_trace_striped_16(query, ref, 11, 1, parasail.blosum62)
            identity, qcov, tcov, cols, matches = values(r.cigar.decode, len(query), len(ref))
            cov80 = qcov >= .8 and tcov >= .8
            edge80 += int(cov80 and identity >= .8)
            edge90 += int(cov80 and identity >= .9)
            # Rank eligible high-coverage alignments before any short local match.
            key = (identity if cov80 else -1., min(qcov, tcov), identity)
            if key > best_key:
                best_key = key
                best = {"reference_row":ti, "identity":identity, "query_coverage":qcov,
                        "target_coverage":tcov, "alignment_columns":cols, "matches":matches,
                        "reference_length":len(ref), "both_coverage_ge_80":cov80}
        for k, v in best.items():
            row[f"{refname}_{k}"] = v
        row[f"{refname}_n_edges_identity80_coverage80"] = edge80
        row[f"{refname}_n_edges_identity90_coverage80"] = edge90
    return row

if __name__ == "__main__":
    start = time.time()
    x = pd.read_csv(DATA / "anthropic_audited_pair_table.csv", usecols=["uuid", "full_name", "target", "sequence"])
    items = [(i, r.uuid, r.full_name, r.target, r.sequence) for i, r in enumerate(x.itertuples(index=False))]
    args = (str(ROOT / "data" / "proteinbase_latest.csv"), str(ROOT / "data" / "overath_hf.parquet"))
    initialise(*args)
    print(f"Aligning {len(items):,} queries against " + ", ".join(f"{k}: {len(v):,}" for k,v in REFS.items()), flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=6, initializer=initialise, initargs=args) as pool:
        for n, row in enumerate(pool.map(process, items, chunksize=1), start=1):
            rows.append(row)
            if n % 120 == 0 or n == len(items):
                print(f"Completed {n}/{len(items)}; {time.time()-start:.1f}s", flush=True)
    out = pd.DataFrame(rows).sort_values("query_index")
    out.to_csv(DATA / "anthropic_sequence_nearest_neighbour.csv", index=False)
    summary = {"status":"complete", "n_queries":len(out), "software":f"parasail {parasail.__version__}",
        "method":"Smith-Waterman local alignment; BLOSUM62; gap open 11; extend 1; all query-reference pairs, no prefilter",
        "identity":"exact matches / local-alignment columns including gaps",
        "cluster_edge_rule":"identity threshold AND query coverage >=0.80 AND reference coverage >=0.80",
        "interpretation_limit":"Sequence non-overlap does not establish independent target constructs, assay platforms, or design pipelines.",
        "comparisons":{}}
    for refname in REFS:
        elig = out[out[f"{refname}_both_coverage_ge_80"]]
        ids = elig[f"{refname}_identity"]
        summary["comparisons"][refname] = {
            "reference_unique_sequences":len(REFS[refname]),
            "pairwise_alignments":len(out)*len(REFS[refname]),
            "n_with_any_both_coverage_ge_80":len(elig),
            "n_edges_identity80_coverage80":int((out[f"{refname}_n_edges_identity80_coverage80"]>0).sum()),
            "n_edges_identity90_coverage80":int((out[f"{refname}_n_edges_identity90_coverage80"]>0).sum()),
            "eligible_max_identity_quantiles":ids.quantile([0,.05,.25,.5,.75,.95,1]).to_dict() if len(elig) else {},
            "unconstrained_short_local_identity_is_not_a_duplicate_measure":True}
    summary["elapsed_seconds"] = time.time()-start
    (DATA / "anthropic_sequence_independence_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
