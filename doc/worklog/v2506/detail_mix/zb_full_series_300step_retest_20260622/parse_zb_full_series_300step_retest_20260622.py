#!/usr/bin/env python3
import csv
import re
import statistics
from pathlib import Path

ROOT = Path("/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb")
OUTDIR = ROOT / "doc/worklog/v2506/detail_mix/zb_full_series_300step_retest_20260622"
MANIFEST = OUTDIR / "run_manifest.tsv"
RESULTS = OUTDIR / "results.csv"
AGG = OUTDIR / "aggregate.csv"

FIELDS = [
    "case", "mode", "np", "omp_threads", "start", "end", "exit",
    "log", "controlDict", "expected_iterations", "log_exists",
    "bad_markers", "end_main", "real", "user", "sys", "iterations",
    "collisions", "candidates", "particles", "stuck", "total_energy",
    "move_collide", "move", "build_occ", "collision_phase", "post",
    "full_evolve", "openmp_enabled", "openmp_threads_reported",
    "migration_calls", "migration_wall", "dlb_checks", "dlb_rebalances",
    "dlb_triggered", "dlb_wall_max", "dlb_check_max", "parmetis_max",
    "dlb_migration_max", "particles_rank_min", "particles_rank_max",
    "particles_rank_ratio", "rank_wall_steps", "rank_wall_min",
    "rank_wall_max", "rank_wall_ratio", "correctness_ok",
]

FLOAT_PATTERNS = {
    "real": r"^real\s+([0-9.eE+-]+)",
    "user": r"^user\s+([0-9.eE+-]+)",
    "sys": r"^sys\s+([0-9.eE+-]+)",
    "total_energy": r"Total energy\s+=\s+([0-9.eE+-]+)",
    "move_collide": r"move\+collide wall \[s\]\s+=\s+([0-9.eE+-]+)",
    "move": r"move only \[s\]\s+=\s+([0-9.eE+-]+)",
    "build_occ": r"buildCellOccupancy \[s\]\s+=\s+([0-9.eE+-]+)",
    "collision_phase": r"collision phase \[s\]\s+=\s+([0-9.eE+-]+)",
    "post": r"post fields/output \[s\]\s+=\s+([0-9.eE+-]+)",
    "full_evolve": r"full evolve wall \[s\]\s+=\s+([0-9.eE+-]+)",
    "migration_wall": r"migration wall time \[s\]\s+=\s+([0-9.eE+-]+)",
    "dlb_wall_max": r"Phase C auto DLB wall max \[s\]\s*=\s*([0-9.eE+-]+)",
    "dlb_check_max": r"Phase C auto DLB check max \[s\]\s*=\s*([0-9.eE+-]+)",
    "parmetis_max": r"Phase C ParMETIS max \[s\]\s*=\s*([0-9.eE+-]+)",
    "dlb_migration_max": r"Phase C migration max \[s\]\s*=\s*([0-9.eE+-]+)",
    "particles_rank_ratio": r"particles per rank\s+=\s+min\s+[0-9.eE+-]+\s+max\s+[0-9.eE+-]+\s+max/min\s+([0-9.eE+-]+)",
    "rank_wall_ratio": r"rank wall time \(evolve,\s+[0-9]+\s+steps\):\s+min\s+[0-9.eE+-]+\s+max\s+[0-9.eE+-]+\s+max/min\s+([0-9.eE+-]+)",
}

INT_PATTERNS = {
    "iterations": r"Total Iterations\s+=\s+([0-9]+)",
    "collisions": r"Collisions\s+=\s+([0-9]+)",
    "candidates": r"Collision candidates\s+=\s+([0-9]+)",
    "particles": r"Number of DSMC particles\s+=\s+([0-9]+)",
    "stuck": r"Number of stuck particles\s+=\s+([0-9]+)",
    "openmp_enabled": r"OpenMP enabled\s+=\s+([0-9]+)",
    "openmp_threads_reported": r"OpenMP max threads\s+=\s+([0-9]+)",
    "migration_calls": r"migration calls\s+=\s+([0-9]+)",
    "dlb_checks": r"Phase C auto DLB checks\s+=\s+([0-9]+)",
    "dlb_rebalances": r"Phase C auto DLB rebalances\s+=\s+([0-9]+)",
    "dlb_triggered": r"Phase C auto DLB triggered checks\s+=\s+([0-9]+)",
}

def last_match(pattern, text, cast):
    matches = re.findall(pattern, text, re.MULTILINE)
    if not matches:
        return ""
    return cast(matches[-1])

def parse_rank_wall(text):
    m = re.findall(
        r"rank wall time \(evolve,\s+([0-9]+)\s+steps\):\s+min\s+([0-9.eE+-]+)\s+max\s+([0-9.eE+-]+)\s+max/min\s+([0-9.eE+-]+)",
        text,
    )
    if not m:
        return "", "", "", ""
    steps, min_v, max_v, ratio = m[-1]
    return int(steps), float(min_v), float(max_v), float(ratio)

def parse_particles_rank(text):
    m = re.findall(
        r"particles per rank\s+=\s+min\s+([0-9.eE+-]+)\s+max\s+([0-9.eE+-]+)\s+max/min\s+([0-9.eE+-]+)",
        text,
    )
    if not m:
        return "", "", ""
    min_v, max_v, ratio = m[-1]
    return float(min_v), float(max_v), float(ratio)

def parse_row(row):
    log = Path(row["log"])
    data = {
        "case": row["case"],
        "mode": row["mode"],
        "np": int(row["np"]),
        "omp_threads": int(row["omp_threads"]),
        "start": row["start"],
        "end": row["end"],
        "exit": int(row["exit"]),
        "log": row["log"],
        "controlDict": row["controlDict"],
        "expected_iterations": 300,
        "log_exists": log.exists(),
    }

    if not log.exists():
        for field in FIELDS:
            data.setdefault(field, "")
        data["correctness_ok"] = False
        return data

    text = log.read_text(errors="ignore")
    markers = []
    for marker in ("FOAM FATAL", "Segmentation", "BAD TERMINATION", "Killed", "Floating point exception"):
        if marker in text:
            markers.append(marker)
    data["bad_markers"] = ",".join(markers)
    data["end_main"] = "End main" in text

    for key, pat in FLOAT_PATTERNS.items():
        data[key] = last_match(pat, text, float)
    for key, pat in INT_PATTERNS.items():
        data[key] = last_match(pat, text, int)

    pr_min, pr_max, pr_ratio = parse_particles_rank(text)
    data["particles_rank_min"] = pr_min
    data["particles_rank_max"] = pr_max
    data["particles_rank_ratio"] = pr_ratio

    rw_steps, rw_min, rw_max, rw_ratio = parse_rank_wall(text)
    data["rank_wall_steps"] = rw_steps
    data["rank_wall_min"] = rw_min
    data["rank_wall_max"] = rw_max
    data["rank_wall_ratio"] = rw_ratio

    data["correctness_ok"] = (
        data["exit"] == 0
        and data["end_main"]
        and not markers
        and data.get("iterations") == 300
        and data.get("stuck") == 0
    )

    for field in FIELDS:
        data.setdefault(field, "")
    return data

def main():
    with MANIFEST.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    parsed = [parse_row(row) for row in rows]

    with RESULTS.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(parsed)

    metric_fields = [
        "real", "full_evolve", "move_collide", "move", "build_occ",
        "collision_phase", "post", "particles", "collisions", "candidates",
        "total_energy", "migration_wall", "particles_rank_ratio",
        "rank_wall_ratio", "dlb_rebalances", "dlb_triggered",
    ]
    agg_fields = ["case", "mode", "runs", "ok_runs"]
    for m in metric_fields:
        agg_fields += [f"{m}_mean", f"{m}_stdev", f"{m}_min", f"{m}_max"]

    grouped = {}
    for row in parsed:
        grouped.setdefault((row["case"], row["mode"]), []).append(row)

    agg_rows = []
    for (case, mode), group in grouped.items():
        ok_group = [r for r in group if r["correctness_ok"]]
        out = {
            "case": case,
            "mode": mode,
            "runs": len(group),
            "ok_runs": len(ok_group),
        }
        for m in metric_fields:
            vals = [r[m] for r in ok_group if r.get(m) not in ("", None)]
            vals = [float(v) for v in vals]
            if vals:
                out[f"{m}_mean"] = statistics.mean(vals)
                out[f"{m}_stdev"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
                out[f"{m}_min"] = min(vals)
                out[f"{m}_max"] = max(vals)
            else:
                out[f"{m}_mean"] = ""
                out[f"{m}_stdev"] = ""
                out[f"{m}_min"] = ""
                out[f"{m}_max"] = ""
        agg_rows.append(out)

    with AGG.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=agg_fields)
        writer.writeheader()
        writer.writerows(agg_rows)

if __name__ == "__main__":
    main()
