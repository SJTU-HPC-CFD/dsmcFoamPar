#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import re
import statistics
import sys
from pathlib import Path


ROOT = Path("/home/superxcx/code/OpenFoam/OF-1706/hyStrath_dlb")
OUTDIR = ROOT / "doc/worklog/v2506/detail_mix/repeat_nowrite_compute_20260610_rerun"
MANIFEST = OUTDIR / "run_manifest.tsv"
RESULTS_CSV = OUTDIR / "results.csv"
AGG_CSV = OUTDIR / "aggregate.csv"
REPORT_MD = OUTDIR / "performance_correctness_summary_20260610.md"

EXPECTED_ITERS = {"ourmesh": 500, "zb": 300}
MODE_ORDER = ["OMP8", "MPI2xOMP4", "MPI4xOMP2", "MPI8", "MPI8origin"]
CASE_ORDER = ["ourmesh", "zb"]

BAD_MARKERS = [
    "FOAM FATAL",
    "Segmentation",
    "SIGSEGV",
    "MPI_ABORT",
    "BAD TERMINATION",
    "Killed",
]


def last_float(text: str, pattern: str) -> float | None:
    vals = re.findall(pattern, text, flags=re.MULTILINE)
    if not vals:
        return None
    value = vals[-1]
    if isinstance(value, tuple):
        value = value[-1]
    try:
        return float(value)
    except ValueError:
        return None


def last_int(text: str, pattern: str) -> int | None:
    value = last_float(text, pattern)
    if value is None:
        return None
    return int(value)


def marker_hits(text: str) -> list[str]:
    hits = [m for m in BAD_MARKERS if m in text]
    if re.search(r"(^|[^A-Za-z])nan([^A-Za-z]|$)", text):
        hits.append("nan")
    if re.search(r"(^|[^A-Za-z])NaN([^A-Za-z]|$)", text):
        hits.append("NaN")
    return hits


def parse_log(row: dict[str, str]) -> dict[str, object]:
    log_path = Path(row["log"])
    text = log_path.read_text(errors="replace") if log_path.exists() else ""
    hits = marker_hits(text)

    particles_rank = re.findall(
        r"particles per rank\s+=\s+min\s+([0-9.+\-Ee]+)\s+max\s+([0-9.+\-Ee]+)\s+max/min\s+([0-9.+\-Ee]+)",
        text,
    )
    rank_wall = re.findall(
        r"rank wall time \(evolve,\s+([0-9]+)\s+steps\):\s+min\s+([0-9.+\-Ee]+)\s+max\s+([0-9.+\-Ee]+)\s+max/min\s+([0-9.+\-Ee]+)",
        text,
    )

    parsed: dict[str, object] = dict(row)
    parsed.update(
        {
            "expected_iterations": EXPECTED_ITERS.get(row["case"], None),
            "log_exists": log_path.exists(),
            "bad_markers": ",".join(hits),
            "end_main": "End main" in text,
            "real": last_float(text, r"^real\s+([0-9.+\-Ee]+)$"),
            "user": last_float(text, r"^user\s+([0-9.+\-Ee]+)$"),
            "sys": last_float(text, r"^sys\s+([0-9.+\-Ee]+)$"),
            "iterations": last_int(text, r"Total Iterations\s+=\s+([0-9]+)"),
            "collisions": last_int(text, r"Collisions\s+=\s+([0-9]+)"),
            "candidates": last_int(text, r"Collision candidates\s+=\s+([0-9]+)"),
            "particles": last_int(text, r"Number of DSMC particles\s+=\s+([0-9]+)"),
            "stuck": last_int(text, r"Number of stuck particles\s+=\s+([0-9]+)"),
            "total_energy": last_float(text, r"Total energy\s+=\s+([0-9.+\-Ee]+)"),
            "move_collide": last_float(text, r"move\+collide wall \[s\]\s+=\s+([0-9.+\-Ee]+)"),
            "move": last_float(text, r"move only \[s\]\s+=\s+([0-9.+\-Ee]+)"),
            "build_occ": last_float(text, r"buildCellOccupancy \[s\]\s+=\s+([0-9.+\-Ee]+)"),
            "collision_phase": last_float(text, r"collision phase \[s\]\s+=\s+([0-9.+\-Ee]+)"),
            "post": last_float(text, r"post fields/output \[s\]\s+=\s+([0-9.+\-Ee]+)"),
            "full_evolve": last_float(text, r"full evolve wall \[s\]\s+=\s+([0-9.+\-Ee]+)"),
            "openmp_enabled": last_int(text, r"OpenMP enabled\s+=\s+([0-9]+)"),
            "openmp_threads": last_int(text, r"OpenMP max threads\s+=\s+([0-9]+)"),
            "migration_calls": last_int(text, r"migration calls\s+=\s+([0-9]+)"),
            "migration_wall": last_float(text, r"migration wall time \[s\]\s+=\s+([0-9.+\-Ee]+)"),
            "dlb_checks": last_int(text, r"Phase C auto DLB checks\s+=\s+([0-9]+)"),
            "dlb_rebalances": last_int(text, r"Phase C auto DLB rebalances\s+=\s+([0-9]+)"),
            "dlb_triggered": last_int(text, r"Phase C auto DLB triggered checks\s+=\s+([0-9]+)"),
            "dlb_wall_max": last_float(text, r"Phase C auto DLB wall max \[s\]=\s*([0-9.+\-Ee]+)"),
            "dlb_check_max": last_float(text, r"Phase C auto DLB check max \[s\]=\s*([0-9.+\-Ee]+)"),
            "parmetis_max": last_float(text, r"Phase C ParMETIS max \[s\]\s+=\s+([0-9.+\-Ee]+)"),
            "dlb_migration_max": last_float(text, r"Phase C migration max \[s\]\s+=\s+([0-9.+\-Ee]+)"),
        }
    )

    if particles_rank:
        pr_min, pr_max, pr_ratio = particles_rank[-1]
        parsed["particles_rank_min"] = float(pr_min)
        parsed["particles_rank_max"] = float(pr_max)
        parsed["particles_rank_ratio"] = float(pr_ratio)
    else:
        parsed["particles_rank_min"] = None
        parsed["particles_rank_max"] = None
        parsed["particles_rank_ratio"] = None

    if rank_wall:
        steps, rw_min, rw_max, rw_ratio = rank_wall[-1]
        parsed["rank_wall_steps"] = int(steps)
        parsed["rank_wall_min"] = float(rw_min)
        parsed["rank_wall_max"] = float(rw_max)
        parsed["rank_wall_ratio"] = float(rw_ratio)
    else:
        parsed["rank_wall_steps"] = None
        parsed["rank_wall_min"] = None
        parsed["rank_wall_max"] = None
        parsed["rank_wall_ratio"] = None

    expected = EXPECTED_ITERS.get(row["case"])
    exit_ok = row.get("exit") == "0"
    iter_ok = parsed["iterations"] == expected
    end_ok = bool(parsed["end_main"])
    marker_ok = not hits
    particles_ok = isinstance(parsed["particles"], int) and parsed["particles"] > 0
    stuck_ok = parsed["stuck"] == 0
    collisions_ok = isinstance(parsed["collisions"], int) and parsed["collisions"] >= 0
    candidates_ok = isinstance(parsed["candidates"], int) and parsed["candidates"] >= 0
    energy_ok = isinstance(parsed["total_energy"], float) and math.isfinite(parsed["total_energy"])
    parsed["correctness_ok"] = all(
        [
            exit_ok,
            iter_ok,
            end_ok,
            marker_ok,
            particles_ok,
            stuck_ok,
            collisions_ok,
            candidates_ok,
            energy_ok,
        ]
    )
    return parsed


def fnum(value: object, digits: int = 2) -> str:
    if value is None or value == "":
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if math.isnan(value):
            return "n/a"
        return f"{value:.{digits}f}"
    return str(value)


def mean(values: list[float]) -> float | None:
    vals = [v for v in values if v is not None]
    return statistics.mean(vals) if vals else None


def stdev(values: list[float]) -> float | None:
    vals = [v for v in values if v is not None]
    return statistics.stdev(vals) if len(vals) > 1 else 0.0 if len(vals) == 1 else None


def aggregate(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    aggs: list[dict[str, object]] = []
    for case in CASE_ORDER:
        for mode in MODE_ORDER:
            group = [r for r in rows if r["case"] == case and r["mode"] == mode]
            if not group:
                continue
            agg: dict[str, object] = {
                "case": case,
                "mode": mode,
                "runs": len(group),
                "ok_runs": sum(1 for r in group if r["correctness_ok"]),
            }
            for key in [
                "real",
                "full_evolve",
                "move_collide",
                "move",
                "build_occ",
                "collision_phase",
                "post",
                "particles",
                "collisions",
                "candidates",
                "total_energy",
                "migration_wall",
                "particles_rank_ratio",
                "rank_wall_ratio",
                "dlb_rebalances",
                "dlb_triggered",
            ]:
                vals = [r[key] for r in group if isinstance(r.get(key), (int, float))]
                agg[f"{key}_mean"] = mean(vals)
                agg[f"{key}_stdev"] = stdev(vals)
                agg[f"{key}_min"] = min(vals) if vals else None
                agg[f"{key}_max"] = max(vals) if vals else None
            aggs.append(agg)
    return aggs


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def agg_lookup(aggs: list[dict[str, object]], case: str, mode: str) -> dict[str, object] | None:
    for agg in aggs:
        if agg["case"] == case and agg["mode"] == mode:
            return agg
    return None


def pct_delta(value: object, base: object) -> str:
    if not isinstance(value, (int, float)) or not isinstance(base, (int, float)) or float(base) == 0.0:
        return "n/a"
    delta = float(value) - float(base)
    pct = delta / float(base) * 100.0
    return f"{delta:.2f} s ({pct:+.2f}%)"


def render_report(rows: list[dict[str, object]], aggs: list[dict[str, object]]) -> str:
    lines: list[str] = []
    lines.append("# No-write compute performance repeat summary - 2026-06-10")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append("This report summarizes three repeat runs for each requested 8-core mode:")
    lines.append("")
    lines.append("- `ourmesh`: 500 steps (`endTime 5.e-05`, `deltaT 1.e-07`)")
    lines.append("- `zb-cylinder-react`: 300 steps (`endTime 1.9920146682e-05`, `deltaT 6.640048894e-08`)")
    lines.append("- output writing disabled for these step ranges by keeping `writeControl runTime` and setting `writeInterval 1.e-3`, which is larger than both case end times")
    lines.append("- `OMP8`, `MPI2xOMP4`, `MPI4xOMP2`, and `MPI8` use the replicated-mesh/raw-MPI path and do not use OpenFOAM `-parallel`")
    lines.append("- `MPI8origin` is the standard OpenFOAM decomposed baseline, run as `mpirun -np 8 dsmcFoam+ -parallel`")
    lines.append("")
    lines.append("Logs and control snapshots are under:")
    lines.append("")
    lines.append("```text")
    lines.append(str(OUTDIR.relative_to(ROOT)))
    lines.append("```")
    lines.append("")

    bad_rows = [r for r in rows if not r["correctness_ok"]]
    lines.append("## Correctness Status")
    lines.append("")
    if bad_rows:
        lines.append(f"{len(bad_rows)} runs failed correctness gates:")
        lines.append("")
        lines.append("| case | mode | rep | exit | iterations | end_main | markers | log |")
        lines.append("|---|---|---:|---:|---:|---|---|---|")
        for r in bad_rows:
            lines.append(
                f"| {r['case']} | {r['mode']} | {r['rep']} | {r['exit']} | "
                f"{r['iterations']} | {r['end_main']} | {r['bad_markers']} | `{Path(str(r['log'])).name}` |"
            )
    else:
        lines.append("All parsed runs passed the correctness gates: exit 0, expected iteration count, `End main`, no fatal/segmentation/MPI abort/NaN markers, nonzero final particles, zero stuck particles, and finite total energy.")
    lines.append("")

    for case in CASE_ORDER:
        expected = EXPECTED_ITERS[case]
        case_rows = [r for r in rows if r["case"] == case]
        lines.append(f"## {case} Performance")
        lines.append("")
        lines.append(f"Expected iterations: `{expected}`.")
        lines.append("")
        lines.append("| mode | ok/runs | real mean | real stdev | real min-max | full evolve mean | move mean | buildOcc mean | collision mean | post mean |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        best_mode = None
        best_real = None
        for mode in MODE_ORDER:
            agg = agg_lookup(aggs, case, mode)
            if not agg:
                continue
            real_mean = agg["real_mean"]
            if isinstance(real_mean, float) and (best_real is None or real_mean < best_real):
                best_real = real_mean
                best_mode = mode
            lines.append(
                f"| {mode} | {agg['ok_runs']}/{agg['runs']} | "
                f"{fnum(agg['real_mean'])} | {fnum(agg['real_stdev'])} | "
                f"{fnum(agg['real_min'])}-{fnum(agg['real_max'])} | "
                f"{fnum(agg['full_evolve_mean'])} | {fnum(agg['move_mean'])} | "
                f"{fnum(agg['build_occ_mean'])} | {fnum(agg['collision_phase_mean'])} | "
                f"{fnum(agg['post_mean'])} |"
            )
        lines.append("")
        if best_mode and best_real is not None:
            lines.append(f"Best mean external `real`: `{best_mode}` at `{best_real:.2f} s`.")
            omp8 = agg_lookup(aggs, case, "OMP8")
            if omp8 and isinstance(omp8.get("real_mean"), float):
                base = float(omp8["real_mean"])
                lines.append("")
                lines.append("| mode | mean real delta vs OMP8 | percent |")
                lines.append("|---|---:|---:|")
                for mode in MODE_ORDER:
                    agg = agg_lookup(aggs, case, mode)
                    if not agg or not isinstance(agg.get("real_mean"), float):
                        continue
                    delta = float(agg["real_mean"]) - base
                    pct = delta / base * 100.0 if base else 0.0
                    lines.append(f"| {mode} | {delta:.2f} s | {pct:.2f}% |")
        lines.append("")

        lines.append(f"### {case} Correctness Metrics")
        lines.append("")
        lines.append("| mode | particles mean | particles min-max | collisions mean | candidates mean | total energy mean | energy min-max | stuck max |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
        for mode in MODE_ORDER:
            agg = agg_lookup(aggs, case, mode)
            if not agg:
                continue
            lines.append(
                f"| {mode} | {fnum(agg['particles_mean'], 0)} | "
                f"{fnum(agg['particles_min'], 0)}-{fnum(agg['particles_max'], 0)} | "
                f"{fnum(agg['collisions_mean'], 0)} | {fnum(agg['candidates_mean'], 0)} | "
                f"{fnum(agg['total_energy_mean'], 10)} | "
                f"{fnum(agg['total_energy_min'], 10)}-{fnum(agg['total_energy_max'], 10)} | "
                f"{max((r['stuck'] for r in case_rows if r['mode'] == mode and isinstance(r.get('stuck'), int)), default='n/a')} |"
            )
        lines.append("")

        lines.append(f"### {case} Load-Balance Diagnostics")
        lines.append("")
        lines.append("| mode | migration wall mean | particles max/min mean | rank wall max/min mean | DLB rebalances mean | triggered checks mean |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        for mode in MODE_ORDER:
            agg = agg_lookup(aggs, case, mode)
            if not agg:
                continue
            lines.append(
                f"| {mode} | {fnum(agg['migration_wall_mean'])} | "
                f"{fnum(agg['particles_rank_ratio_mean'], 4)} | "
                f"{fnum(agg['rank_wall_ratio_mean'], 4)} | "
                f"{fnum(agg['dlb_rebalances_mean'], 2)} | {fnum(agg['dlb_triggered_mean'], 2)} |"
            )
        lines.append("")

    origin_cases = [case for case in CASE_ORDER if agg_lookup(aggs, case, "MPI8origin")]
    if origin_cases:
        lines.append("## MPI8origin Baseline Comparison")
        lines.append("")
        lines.append("`MPI8origin` is the original decomposed OpenFOAM `-parallel` baseline.  It is compared with both the single-process OpenMP baseline and the raw-MPI replicated-mesh `MPI8` mode for every case where the baseline was run.")
        lines.append("")
        lines.append("| case | comparison | real mean | delta | full evolve mean | delta |")
        lines.append("|---|---|---:|---:|---:|---:|")
        for case in origin_cases:
            origin = agg_lookup(aggs, case, "MPI8origin")
            replicated = agg_lookup(aggs, case, "MPI8")
            omp8 = agg_lookup(aggs, case, "OMP8")
            if not origin:
                continue
            if omp8:
                lines.append(
                    f"| {case} | MPI8origin vs OMP8 | {fnum(origin.get('real_mean'))} vs {fnum(omp8.get('real_mean'))} | "
                    f"{pct_delta(origin.get('real_mean'), omp8.get('real_mean'))} | "
                    f"{fnum(origin.get('full_evolve_mean'))} vs {fnum(omp8.get('full_evolve_mean'))} | "
                    f"{pct_delta(origin.get('full_evolve_mean'), omp8.get('full_evolve_mean'))} |"
                )
            if replicated:
                lines.append(
                    f"| {case} | MPI8origin vs replicated-mesh MPI8 | {fnum(origin.get('real_mean'))} vs {fnum(replicated.get('real_mean'))} | "
                    f"{pct_delta(origin.get('real_mean'), replicated.get('real_mean'))} | "
                    f"{fnum(origin.get('full_evolve_mean'))} vs {fnum(replicated.get('full_evolve_mean'))} | "
                    f"{pct_delta(origin.get('full_evolve_mean'), replicated.get('full_evolve_mean'))} |"
                )
        lines.append("")

    lines.append("## Per-Run Logs")
    lines.append("")
    lines.append("| case | mode | rep | exit | real | iterations | correctness | log |")
    lines.append("|---|---|---:|---:|---:|---:|---|---|")
    for case in CASE_ORDER:
        for mode in MODE_ORDER:
            for r in [x for x in rows if x["case"] == case and x["mode"] == mode]:
                lines.append(
                    f"| {r['case']} | {r['mode']} | {r['rep']} | {r['exit']} | "
                    f"{fnum(r['real'])} | {r['iterations']} | {r['correctness_ok']} | "
                    f"`{Path(str(r['log'])).name}` |"
                )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    if not MANIFEST.exists():
        print(f"missing manifest: {MANIFEST}", file=sys.stderr)
        return 2
    manifest_rows = read_manifest()
    rows = [parse_log(row) for row in manifest_rows]
    aggs = aggregate(rows)
    write_csv(RESULTS_CSV, rows)
    write_csv(AGG_CSV, aggs)
    REPORT_MD.write_text(render_report(rows, aggs))
    print(REPORT_MD)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
