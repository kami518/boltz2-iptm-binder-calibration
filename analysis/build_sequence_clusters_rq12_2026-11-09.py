"""Exhaustive all-pairs local-alignment clustering for the frozen RQ1/RQ2 cohort.

No length or k-mer prefilter is used. Passing alignments are linked into connected
components at 90% and 80% identity, each requiring >=80% query and target coverage.
Completed per-query blocks are checkpointed so an interrupted run can resume.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / ".release-root").is_file())
LIB = ROOT / "analysis" / "rq12_lib"
sys.path.insert(0, str(LIB))

import numpy as np
import pandas as pd
import parasail

MASTER = ROOT / "analysis" / "v3_outputs" / "pair_master_v3.csv"
WORK = ROOT / "analysis" / "rq12_intermediate" / "alignment_blocks"
OUT = ROOT / "analysis" / "rq12_outputs"
CIGAR_RE = re.compile(r"(\d+)([=XMID])")


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.size[a] < self.size[b]:
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]


def alignment_values(result, qlen: int, tlen: int) -> tuple[float, float, float]:
    matches = mismatches = query_consumed = target_consumed = columns = 0
    cigar = result.cigar.decode
    if isinstance(cigar, bytes):
        cigar = cigar.decode("ascii")
    for count_text, operation in CIGAR_RE.findall(cigar):
        count = int(count_text)
        columns += count
        if operation == "=":
            matches += 1 * count
            query_consumed += count
            target_consumed += count
        elif operation in ("X", "M"):
            if operation == "M":
                matches += count
            else:
                mismatches += count
            query_consumed += count
            target_consumed += count
        elif operation == "I":
            query_consumed += count
        elif operation == "D":
            target_consumed += count
    if columns == 0 or query_consumed == 0 or target_consumed == 0:
        return 0.0, 0.0, 0.0
    return matches / columns, query_consumed / qlen, target_consumed / tlen


def process_block(i: int, sequences: list[str], overwrite: bool = False) -> dict:
    WORK.mkdir(parents=True, exist_ok=True)
    final_path = WORK / f"block_{i:04d}.csv"
    metadata_path = WORK / f"block_{i:04d}.json"
    if final_path.exists() and metadata_path.exists() and not overwrite:
        return json.loads(metadata_path.read_text(encoding="utf-8"))
    started = time.time()
    temp_path = final_path.with_suffix(".csv.tmp")
    n80 = n90 = 0
    with temp_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["query_index", "target_index", "identity", "query_coverage", "target_coverage"])
        query = sequences[i]
        for j in range(i + 1, len(sequences)):
            target = sequences[j]
            result = parasail.sw_trace_striped_16(query, target, 11, 1, parasail.blosum62)
            identity, qcov, tcov = alignment_values(result, len(query), len(target))
            if identity >= 0.80 and qcov >= 0.80 and tcov >= 0.80:
                writer.writerow([i, j, f"{identity:.8f}", f"{qcov:.8f}", f"{tcov:.8f}"])
                n80 += 1
                if identity >= 0.90:
                    n90 += 1
    os.replace(temp_path, final_path)
    metadata = {
        "query_index": i,
        "alignments": len(sequences) - i - 1,
        "edges_80": n80,
        "edges_90": n90,
        "seconds": time.time() - started,
        "status": "complete",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def combine(df: pd.DataFrame) -> None:
    n = len(df)
    expected = [WORK / f"block_{i:04d}.json" for i in range(n - 1)]
    missing = [p.name for p in expected if not p.exists()]
    if missing:
        raise RuntimeError(f"Cannot combine: {len(missing)} blocks missing; first={missing[:5]}")
    uf80, uf90 = UnionFind(n), UnionFind(n)
    edges80 = edges90 = 0
    alignments = 0
    for i in range(n - 1):
        metadata = json.loads((WORK / f"block_{i:04d}.json").read_text(encoding="utf-8"))
        alignments += metadata["alignments"]
        block_path = WORK / f"block_{i:04d}.csv"
        with block_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                a, b = int(row["query_index"]), int(row["target_index"])
                uf80.union(a, b)
                edges80 += 1
                if float(row["identity"]) >= 0.90:
                    uf90.union(a, b)
                    edges90 += 1
    mapping = df[["pair_id", "design_id", "target", "sequence"]].copy()
    for threshold, uf in ((90, uf90), (80, uf80)):
        roots = [uf.find(i) for i in range(n)]
        ordered_roots = {}
        for root in roots:
            if root not in ordered_roots:
                ordered_roots[root] = len(ordered_roots) + 1
        mapping[f"sequence_cluster{threshold}"] = [f"pb{threshold}_{ordered_roots[r]:04d}" for r in roots]
    OUT.mkdir(parents=True, exist_ok=True)
    mapping_path = OUT / "sequence_clusters_rq12.csv"
    mapping.to_csv(mapping_path, index=False)
    summary = {
        "status": "complete",
        "method": "Exhaustive all-pairs Smith-Waterman local alignment; connected components",
        "software": f"parasail {getattr(parasail, '__version__', 'unknown')}",
        "matrix": "BLOSUM62",
        "gap_open": 11,
        "gap_extend": 1,
        "identity_definition": "exact matches / all local-alignment columns, including gaps",
        "coverage_definition": "aligned residues consumed / full sequence length; both directions >=0.80",
        "prefilter": "none",
        "pairs": n,
        "pairwise_alignments": alignments,
        "expected_pairwise_alignments": n * (n - 1) // 2,
        "edges_80": edges80,
        "edges_90": edges90,
        "mapping_file": str(mapping_path),
    }
    for threshold in (90, 80):
        column = f"sequence_cluster{threshold}"
        counts = mapping.groupby(column).size()
        by_target = mapping.groupby("target")[column].nunique().to_dict()
        summary[f"clusters_{threshold}"] = {
            "count": int(counts.size),
            "singletons": int((counts == 1).sum()),
            "max_size": int(counts.max()),
            "median_size": float(counts.median()),
            "clusters_by_target": by_target,
        }
    (OUT / "sequence_cluster_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--combine", action="store_true")
    parser.add_argument("--benchmark-blocks", type=int, default=0)
    args = parser.parse_args()
    df = pd.read_csv(MASTER)
    if len(df) != 2084 or df.pair_id.nunique() != 2084 or df.sequence.nunique() != 2084:
        raise RuntimeError("Frozen-cohort gate failed: expected 2084 unique pair IDs and sequences")
    sequences = df.sequence.astype(str).tolist()
    if args.combine:
        combine(df)
        return
    stop = min(len(sequences) - 1, args.benchmark_blocks) if args.benchmark_blocks else len(sequences) - 1
    started = time.time()
    completed = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(process_block, i, sequences): i for i in range(stop)}
        for future in as_completed(futures):
            metadata = future.result()
            completed += 1
            if completed % 25 == 0 or completed == stop:
                elapsed = time.time() - started
                print(f"completed_blocks={completed}/{stop} elapsed={elapsed:.1f}s last={metadata}", flush=True)
    if args.benchmark_blocks:
        print(f"Benchmark complete in {time.time() - started:.2f}s")
    else:
        combine(df)


if __name__ == "__main__":
    main()
