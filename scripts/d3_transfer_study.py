#!/usr/bin/env python
"""
D3 parameter transfer study.

Evaluates how well D3-BJ parameters optimized on GMTKN55 subsets
transfer to NENCI-2021 and MPCONF196 benchmarks. Compares against
published parameters, no-D3 baselines, and direct optimization on
the target benchmarks.

Usage:
    python scripts/d3_transfer_study.py \
        --gmtkn-energies dft_data/gmtkn_b3lyp_ma_tzvpp/dft_energies.csv \
        --nenci-xyz nenci2021/xyzfiles \
        --nenci-energies dft_data/nenci_b3lyp_ma_tzvpp/dft_energies.csv \
        --mpconf-xyz mpconf196/xyzfiles \
        --mpconf-energies dft_data/mpconf_b3lyp_ma_tzvpp/dft_energies.csv \
        --output results.csv
"""

import argparse
import csv
import sys
import time

from dftd3.benchmark import (
    load_gmtkn55,
    load_mpconf196,
    load_nenci,
    read_energies_csv,
)
from dftd3.optimize import (
    OptimizeConfig,
    _published_params_array,
    evaluate_no_d3,
    evaluate_params,
    fit_d3_params,
    precompute_d3,
)

# ── NCI-relevant GMTKN55 subsets ──────────────────────────────────────────

NCI_SUBSETS = [
    "RG18", "ADIM6", "S22", "S66", "HEAVY28", "CARBHB12",
    "PNICO23", "HAL59", "AHB21", "CHB6", "IL16", "IDISP",
    "ICONF", "ACONF", "Amino20x4", "PCONF21", "MCONF",
    "SCONF", "UPU23", "BUT14DIOL",
]

# Thematic groupings
THEMATIC_GROUPS = {
    "dispersion": ["S22", "S66", "ADIM6", "RG18", "HEAVY28"],
    "hbond": ["AHB21", "CHB6", "CARBHB12"],
    "mixed-nci": ["S22", "S66", "HEAVY28", "AHB21", "IDISP", "HAL59"],
    "conformer": ["ICONF", "ACONF", "PCONF21", "MCONF", "SCONF", "BUT14DIOL"],
    "ionic": ["IL16", "PNICO23", "HAL59"],
    "all-nci": list(NCI_SUBSETS),
}


def _fmt_subsets(subsets):
    """Format subset list for display."""
    if len(subsets) <= 4:
        return "+".join(subsets)
    return "+".join(subsets[:3]) + f"+...({len(subsets)})"


def _load_and_filter(dataset, dft_energies):
    """Filter dataset to reactions with available DFT energies."""
    missing = set(dataset.species.keys()) - set(dft_energies.keys())
    if missing:
        n_before = len(dataset.reactions)
        dataset.reactions = [
            rxn for rxn in dataset.reactions
            if not (set(rxn.stoichiometry.keys()) & missing)
        ]
        n_dropped = n_before - len(dataset.reactions)
        if n_dropped:
            print(f"   Warning: skipped {n_dropped}/{n_before} reactions "
                  f"({len(missing)} missing species)")
    return dataset


def run_baselines(nenci_ds, nenci_energies, nenci_ims,
                  mpconf_ds, mpconf_energies, mpconf_ims,
                  functional, damp, workers):
    """Run Phase 1: baselines on NENCI and MPCONF."""
    rows = []
    pub_params = _published_params_array(functional, damp)

    # Published parameters
    if pub_params is not None:
        nenci_pub = evaluate_params(nenci_ds, nenci_energies, damp, pub_params,
                                    intermediates=nenci_ims)
        mpconf_pub = evaluate_params(mpconf_ds, mpconf_energies, damp, pub_params,
                                     intermediates=mpconf_ims)
        rows.append({
            "label": "published",
            "training": "-",
            "n_train": 0,
            "train_wmad": None,
            "nenci_wmad": nenci_pub["wmad"],
            "nenci_mad": nenci_pub["mad"],
            "mpconf_wmad": mpconf_pub["wmad"],
            "mpconf_mad": mpconf_pub["mad"],
            "s8": pub_params[0],
            "a1": pub_params[1] if damp == "bj" else None,
            "a2": pub_params[2] if damp == "bj" else None,
            "elapsed": 0,
        })
        print(f"   Published:    NENCI WMAD={nenci_pub['wmad']:.3f}  "
              f"MPCONF WMAD={mpconf_pub['wmad']:.3f}")

    # No D3
    nenci_nod3 = evaluate_no_d3(nenci_ds, nenci_energies)
    mpconf_nod3 = evaluate_no_d3(mpconf_ds, mpconf_energies)
    rows.append({
        "label": "no_d3",
        "training": "-",
        "n_train": 0,
        "train_wmad": None,
        "nenci_wmad": nenci_nod3["wmad"],
        "nenci_mad": nenci_nod3["mad"],
        "mpconf_wmad": mpconf_nod3["wmad"],
        "mpconf_mad": mpconf_nod3["mad"],
        "s8": None, "a1": None, "a2": None,
        "elapsed": 0,
    })
    print(f"   No D3:        NENCI WMAD={nenci_nod3['wmad']:.3f}  "
          f"MPCONF WMAD={mpconf_nod3['wmad']:.3f}")

    # Fit directly on NENCI
    print("   Fitting on NENCI...")
    config = OptimizeConfig(damp=damp, test_fraction=0.0, workers=workers)
    t0 = time.time()
    nenci_result = fit_d3_params(nenci_ds, nenci_energies, config)
    elapsed = time.time() - t0
    nenci_opt_params = _result_to_array(nenci_result, damp)
    mpconf_from_nenci = evaluate_params(mpconf_ds, mpconf_energies, damp,
                                        nenci_opt_params, intermediates=mpconf_ims)
    rows.append({
        "label": "direct_nenci",
        "training": "NENCI",
        "n_train": nenci_result.n_train,
        "train_wmad": nenci_result.weighted_mad,
        "nenci_wmad": nenci_result.weighted_mad,
        "nenci_mad": nenci_result.mad,
        "mpconf_wmad": mpconf_from_nenci["wmad"],
        "mpconf_mad": mpconf_from_nenci["mad"],
        "s8": nenci_result.params["s8"],
        "a1": nenci_result.params.get("a1"),
        "a2": nenci_result.params.get("a2"),
        "elapsed": elapsed,
    })
    print(f"   Direct NENCI: NENCI WMAD={nenci_result.weighted_mad:.3f}  "
          f"MPCONF WMAD={mpconf_from_nenci['wmad']:.3f}  ({elapsed:.0f}s)")

    # Fit directly on MPCONF
    print("   Fitting on MPCONF...")
    t0 = time.time()
    mpconf_result = fit_d3_params(mpconf_ds, mpconf_energies, config)
    elapsed = time.time() - t0
    mpconf_opt_params = _result_to_array(mpconf_result, damp)
    nenci_from_mpconf = evaluate_params(nenci_ds, nenci_energies, damp,
                                        mpconf_opt_params, intermediates=nenci_ims)
    rows.append({
        "label": "direct_mpconf",
        "training": "MPCONF",
        "n_train": mpconf_result.n_train,
        "train_wmad": mpconf_result.weighted_mad,
        "nenci_wmad": nenci_from_mpconf["wmad"],
        "nenci_mad": nenci_from_mpconf["mad"],
        "mpconf_wmad": mpconf_result.weighted_mad,
        "mpconf_mad": mpconf_result.mad,
        "s8": mpconf_result.params["s8"],
        "a1": mpconf_result.params.get("a1"),
        "a2": mpconf_result.params.get("a2"),
        "elapsed": elapsed,
    })
    print(f"   Direct MPCONF: NENCI WMAD={nenci_from_mpconf['wmad']:.3f}  "
          f"MPCONF WMAD={mpconf_result.weighted_mad:.3f}  ({elapsed:.0f}s)")

    return rows


def _result_to_array(result, damp):
    """Extract parameter array from FitResult."""
    if damp == "bj":
        return [result.params["s8"], result.params["a1"], result.params["a2"]]
    else:
        return [result.params["rs6"], result.params["s8"]]


def run_transfer(combinations, gmtkn_energies,
                 nenci_ds, nenci_energies, nenci_ims,
                 mpconf_ds, mpconf_energies, mpconf_ims,
                 damp, workers):
    """Run Phase 2: fit on GMTKN55 subsets, evaluate on NENCI + MPCONF."""
    rows = []
    config = OptimizeConfig(damp=damp, test_fraction=0.0, workers=workers)

    for label, subsets in combinations:
        print(f"\n   [{label}] Fitting on {_fmt_subsets(subsets)}...")
        try:
            gmtkn_ds = load_gmtkn55(subsets=subsets)
        except Exception as e:
            print(f"   Error loading {subsets}: {e}")
            continue

        # Filter to available energies
        gmtkn_ds = _load_and_filter(gmtkn_ds, gmtkn_energies)
        if len(gmtkn_ds.reactions) == 0:
            print(f"   Skipping {label}: no reactions after filtering")
            continue

        t0 = time.time()
        result = fit_d3_params(gmtkn_ds, gmtkn_energies, config)
        elapsed = time.time() - t0

        params = _result_to_array(result, damp)
        nenci_stats = evaluate_params(nenci_ds, nenci_energies, damp, params,
                                      intermediates=nenci_ims)
        mpconf_stats = evaluate_params(mpconf_ds, mpconf_energies, damp, params,
                                       intermediates=mpconf_ims)

        rows.append({
            "label": label,
            "training": "+".join(subsets),
            "n_train": result.n_train,
            "train_wmad": result.weighted_mad,
            "nenci_wmad": nenci_stats["wmad"],
            "nenci_mad": nenci_stats["mad"],
            "mpconf_wmad": mpconf_stats["wmad"],
            "mpconf_mad": mpconf_stats["mad"],
            "s8": result.params["s8"],
            "a1": result.params.get("a1"),
            "a2": result.params.get("a2"),
            "elapsed": elapsed,
        })
        print(f"   [{label}] Train WMAD={result.weighted_mad:.3f}  "
              f"NENCI={nenci_stats['wmad']:.3f}  MPCONF={mpconf_stats['wmad']:.3f}  "
              f"({elapsed:.0f}s)")

    return rows


def build_incremental(single_rows, metric="nenci_wmad"):
    """Build incremental combinations ranked by single-subset performance."""
    # Sort individual subsets by their NENCI WMAD
    ranked = sorted(single_rows, key=lambda r: r.get(metric, float("inf")))
    combos = []
    accumulated = []
    for row in ranked[:10]:  # top 10
        subset = row["training"]
        accumulated.append(subset)
        if len(accumulated) >= 2:
            label = f"top{len(accumulated)}_{metric.split('_')[0]}"
            combos.append((label, list(accumulated)))
    return combos


def print_summary(all_rows, damp):
    """Print formatted summary table."""
    print("\n")
    print("   " + "=" * 100)
    print("   Transfer Study Summary")
    print("   " + "=" * 100)
    print()

    if damp == "bj":
        header = (f"   {'Label':<25s} {'N':>5s} {'Train':>7s} {'NENCI':>7s} {'MPCONF':>7s}"
                  f"  {'s8':>7s} {'a1':>7s} {'a2':>7s} {'Time':>5s}")
    else:
        header = (f"   {'Label':<25s} {'N':>5s} {'Train':>7s} {'NENCI':>7s} {'MPCONF':>7s}"
                  f"  {'rs6':>7s} {'s8':>7s} {'Time':>5s}")
    print(header)
    print("   " + "-" * (len(header) - 3))

    for row in all_rows:
        n = row["n_train"] if row["n_train"] else "-"
        train = f"{row['train_wmad']:.3f}" if row["train_wmad"] is not None else "-"
        nenci = f"{row['nenci_wmad']:.3f}"
        mpconf = f"{row['mpconf_wmad']:.3f}"
        elapsed = f"{row['elapsed']:.0f}s" if row["elapsed"] else "-"

        if damp == "bj":
            s8 = f"{row['s8']:.4f}" if row["s8"] is not None else "-"
            a1 = f"{row['a1']:.4f}" if row["a1"] is not None else "-"
            a2 = f"{row['a2']:.4f}" if row["a2"] is not None else "-"
            print(f"   {row['label']:<25s} {str(n):>5s} {train:>7s} {nenci:>7s} {mpconf:>7s}"
                  f"  {s8:>7s} {a1:>7s} {a2:>7s} {elapsed:>5s}")
        else:
            rs6 = f"{row['s8']:.4f}" if row["s8"] is not None else "-"
            s8 = f"{row.get('a1', '-')}" if row.get("a1") is not None else "-"
            print(f"   {row['label']:<25s} {str(n):>5s} {train:>7s} {nenci:>7s} {mpconf:>7s}"
                  f"  {rs6:>7s} {s8:>7s} {elapsed:>5s}")

    print("   " + "=" * (len(header) - 3))


def save_csv(all_rows, output_path, damp):
    """Save results to CSV."""
    fieldnames = ["label", "training", "n_train", "train_wmad",
                  "nenci_wmad", "nenci_mad", "mpconf_wmad", "mpconf_mad",
                  "s8", "a1", "a2", "elapsed"]
    if damp == "zero":
        fieldnames = [f for f in fieldnames if f not in ("a1", "a2")]
        fieldnames.insert(fieldnames.index("s8") + 1, "rs6")

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)
    print(f"\n   Saved results to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="D3 parameter transfer study: GMTKN55 → NENCI/MPCONF")
    parser.add_argument("--gmtkn-energies", required=True,
                        help="GMTKN55 DFT energies CSV")
    parser.add_argument("--nenci-xyz", required=True,
                        help="NENCI-2021 XYZ directory")
    parser.add_argument("--nenci-energies", required=True,
                        help="NENCI DFT energies CSV")
    parser.add_argument("--nenci-subsets", default="S101,IonPi",
                        help="NENCI subsets (default: S101,IonPi)")
    parser.add_argument("--mpconf-xyz", required=True,
                        help="MPCONF196 XYZ directory")
    parser.add_argument("--mpconf-energies", required=True,
                        help="MPCONF DFT energies CSV")
    parser.add_argument("--damp", default="bj", choices=("bj", "zero"))
    parser.add_argument("--func", default="B3LYP",
                        help="Functional for published params (default: B3LYP)")
    parser.add_argument("--workers", type=int, default=-1,
                        help="DE parallel workers (default: -1 = all cores)")
    parser.add_argument("--output", default="transfer_results.csv",
                        help="Output CSV path (default: transfer_results.csv)")
    parser.add_argument("--skip-phase1", action="store_true",
                        help="Skip Phase 1 baselines")
    parser.add_argument("--skip-phase2", action="store_true",
                        help="Skip Phase 2 GMTKN55 transfer")
    args = parser.parse_args()

    print("\n   " + "=" * 72)
    print("   D3 Parameter Transfer Study")
    print("   " + "=" * 72)

    # ── Load validation datasets ──────────────────────────────────────

    nenci_subsets = [s.strip() for s in args.nenci_subsets.split(",")]
    print(f"\n   Loading NENCI ({','.join(nenci_subsets)})...")
    nenci_ds = load_nenci(args.nenci_xyz, subsets=nenci_subsets)
    nenci_energies, _, _ = read_energies_csv(args.nenci_energies)
    nenci_ds = _load_and_filter(nenci_ds, nenci_energies)
    print(f"   NENCI: {len(nenci_ds.reactions)} reactions, {len(nenci_ds.species)} species")

    print("\n   Loading MPCONF196...")
    mpconf_ds = load_mpconf196(args.mpconf_xyz)
    mpconf_energies, _, _ = read_energies_csv(args.mpconf_energies)
    mpconf_ds = _load_and_filter(mpconf_ds, mpconf_energies)
    print(f"   MPCONF: {len(mpconf_ds.reactions)} reactions, {len(mpconf_ds.species)} species")

    # ── Pre-compute intermediates ─────────────────────────────────────

    print("\n   Pre-computing D3 intermediates for validation sets...")
    t0 = time.time()
    nenci_ims = {name: precompute_d3(sp) for name, sp in nenci_ds.species.items()}
    mpconf_ims = {name: precompute_d3(sp) for name, sp in mpconf_ds.species.items()}
    print(f"   Done ({time.time() - t0:.1f}s)")

    all_rows = []

    # ── Phase 1: Baselines ────────────────────────────────────────────

    if not args.skip_phase1:
        print("\n   " + "-" * 72)
        print("   Phase 1: Baselines")
        print("   " + "-" * 72)
        baseline_rows = run_baselines(
            nenci_ds, nenci_energies, nenci_ims,
            mpconf_ds, mpconf_energies, mpconf_ims,
            args.func, args.damp, args.workers)
        all_rows.extend(baseline_rows)

    # ── Phase 2: GMTKN55 Transfer ─────────────────────────────────────

    if not args.skip_phase2:
        print("\n   " + "-" * 72)
        print("   Phase 2: GMTKN55 Subset Transfer")
        print("   " + "-" * 72)

        # Load GMTKN energies once
        gmtkn_energies, _, _ = read_energies_csv(args.gmtkn_energies)

        # 2a: Individual subsets
        print("\n   --- Individual Subsets ---")
        individual_combos = [(s, [s]) for s in NCI_SUBSETS]
        single_rows = run_transfer(
            individual_combos, gmtkn_energies,
            nenci_ds, nenci_energies, nenci_ims,
            mpconf_ds, mpconf_energies, mpconf_ims,
            args.damp, args.workers)
        all_rows.extend(single_rows)

        # 2b: Thematic groups
        print("\n   --- Thematic Groups ---")
        thematic_combos = list(THEMATIC_GROUPS.items())
        thematic_rows = run_transfer(
            thematic_combos, gmtkn_energies,
            nenci_ds, nenci_energies, nenci_ims,
            mpconf_ds, mpconf_energies, mpconf_ims,
            args.damp, args.workers)
        all_rows.extend(thematic_rows)

        # 2c: Incremental — ranked by NENCI WMAD
        print("\n   --- Incremental (ranked by NENCI WMAD) ---")
        incr_combos = build_incremental(single_rows, "nenci_wmad")
        if incr_combos:
            incr_rows = run_transfer(
                incr_combos, gmtkn_energies,
                nenci_ds, nenci_energies, nenci_ims,
                mpconf_ds, mpconf_energies, mpconf_ims,
                args.damp, args.workers)
            all_rows.extend(incr_rows)

    # ── Summary ───────────────────────────────────────────────────────

    print_summary(all_rows, args.damp)
    save_csv(all_rows, args.output, args.damp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
