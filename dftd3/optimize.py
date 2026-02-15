"""
optimize.py - D3 parameter optimization against benchmark datasets.

Fits D3 damping parameters (s8, a1, a2 for BJ; rs6, s8 for zero) by
minimizing weighted mean absolute deviation against benchmark reference
energies. Uses scipy.optimize.differential_evolution for robust global
optimization.
"""

import logging
import os
import time
from argparse import ArgumentParser
from dataclasses import dataclass, field

import numpy as np

from dftd3.benchmark import (
    generate_orca_inputs,
    load_dataset,
    parse_orca_outputs,
    read_energies_csv,
    write_energies_csv,
)
from dftd3.dftd3 import (
    ALPHA6,
    ALPHA8,
    AUTOANG,
    AUTOKCAL,
    _getc6_all_pairs,
    _ncoord_vectorized,
    _r,
    _r2r4,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# D3 precomputation
# ---------------------------------------------------------------------------

@dataclass
class D3Intermediates:
    """Precomputed D3 intermediates for a single species.

    These quantities depend only on the molecular geometry and C6
    coefficients, not on the damping parameters. Computing them once
    and reusing across optimizer iterations is the key to fast fitting.
    """

    c6_pairs: np.ndarray
    c8_pairs: np.ndarray
    dist_au_6: np.ndarray
    dist_au_8: np.ndarray
    rr_zero: np.ndarray   # _r[A,B] / dist_au — for zero-damping
    rr_bj: np.ndarray     # sqrt(C8/C6) — for BJ-damping


def precompute_d3(species):
    """Precompute parameter-independent D3 quantities for a species.

    Parameters
    ----------
    species : Species
        Molecular species with atomnos and positions.

    Returns
    -------
    D3Intermediates
    """
    atomnos = species.atomnos
    coords = np.array(species.positions)
    elem_idx = (atomnos - 1).astype(int)

    cn, dist_ang = _ncoord_vectorized(coords, elem_idx)
    c6_matrix = _getc6_all_pairs(elem_idx, cn)

    r2r4_atoms = _r2r4[elem_idx]
    c8_matrix = 3.0 * c6_matrix * np.outer(r2r4_atoms, r2r4_atoms)

    j_idx, k_idx = np.triu_indices(len(atomnos), k=1)
    c6_pairs = c6_matrix[j_idx, k_idx]
    c8_pairs = c8_matrix[j_idx, k_idx]
    dist_au = dist_ang[j_idx, k_idx] / AUTOANG

    dist_au_6 = dist_au ** 6
    dist_au_8 = dist_au ** 8

    rr_zero = _r[elem_idx[j_idx], elem_idx[k_idx]] / dist_au

    safe_c6 = np.where(c6_pairs > 0, c6_pairs, 1.0)
    rr_bj = np.sqrt(c8_pairs / safe_c6)

    return D3Intermediates(
        c6_pairs=c6_pairs,
        c8_pairs=c8_pairs,
        dist_au_6=dist_au_6,
        dist_au_8=dist_au_8,
        rr_zero=rr_zero,
        rr_bj=rr_bj,
    )


# ---------------------------------------------------------------------------
# Fast D3 energy functions (parameter-dependent only)
# ---------------------------------------------------------------------------

def d3_energy_bj(im, s6, s8, a1, a2):
    """Compute D3(BJ) energy in kcal/mol from precomputed intermediates."""
    tmp = a1 * im.rr_bj + a2
    damp6 = tmp ** 6
    damp8 = tmp ** 8
    e6 = -s6 * np.sum(im.c6_pairs / (im.dist_au_6 + damp6)) * AUTOKCAL
    e8 = -s8 * np.sum(im.c8_pairs / (im.dist_au_8 + damp8)) * AUTOKCAL
    return e6 + e8


def d3_energy_zero(im, s6, rs6, s8):
    """Compute D3(zero) energy in kcal/mol from precomputed intermediates."""
    damp6 = 1.0 / (1.0 + 6.0 * (rs6 * im.rr_zero) ** ALPHA6)
    damp8 = 1.0 / (1.0 + 6.0 * im.rr_zero ** ALPHA8)  # rs8=1.0 always
    e6 = -s6 * np.sum(im.c6_pairs * damp6 / im.dist_au_6) * AUTOKCAL
    e8 = -s8 * np.sum(im.c8_pairs * damp8 / im.dist_au_8) * AUTOKCAL
    return e6 + e8


# ---------------------------------------------------------------------------
# Optimizer configuration and results
# ---------------------------------------------------------------------------

@dataclass
class OptimizeConfig:
    """Configuration for the D3 parameter optimizer."""

    damp: str = "bj"
    test_fraction: float = 0.2
    random_seed: int = 42
    n_folds: int = 0
    maxiter: int = 1000
    tol: float = 1e-6
    popsize: int = 15
    workers: int = 1
    s8_bounds: tuple = (0.0, 4.0)
    a1_bounds: tuple = (0.0, 1.0)
    a2_bounds: tuple = (2.0, 9.0)
    rs6_bounds: tuple = (0.5, 2.0)


@dataclass
class FitResult:
    """Result of D3 parameter optimization."""

    damp: str
    params: dict
    weighted_mad: float
    mad: float
    rmsd: float
    max_error: float
    test_weighted_mad: float = None
    test_mad: float = None
    test_rmsd: float = None
    test_max_error: float = None
    n_train: int = 0
    n_test: int = 0
    elapsed_seconds: float = 0.0
    per_reaction_errors: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Objective functions
# ---------------------------------------------------------------------------

def _reaction_errors(reactions, intermediates_dict, dft_energies, damp, params):
    """Compute per-reaction errors for given parameters.

    Returns list of (reaction, error, weight) tuples.
    """
    s6 = 1.0
    results = []

    for rxn in reactions:
        dft_rxn = 0.0
        d3_rxn = 0.0
        for sp_name, coeff in rxn.stoichiometry.items():
            dft_rxn += coeff * dft_energies[sp_name] * AUTOKCAL
            im = intermediates_dict[sp_name]
            if damp == "bj":
                s8, a1, a2 = params
                d3_rxn += coeff * d3_energy_bj(im, s6, s8, a1, a2)
            else:
                rs6, s8 = params
                d3_rxn += coeff * d3_energy_zero(im, s6, rs6, s8)

        error = dft_rxn + d3_rxn - rxn.reference_energy
        results.append((rxn, error, rxn.weight))

    return results


def _weighted_mad(reactions, intermediates_dict, dft_energies, damp, params):
    """Compute weighted mean absolute deviation."""
    errors = _reaction_errors(reactions, intermediates_dict, dft_energies, damp, params)
    weight_sum = sum(w for _, _, w in errors)
    if weight_sum == 0:
        return float("inf")
    return sum(w * abs(e) for _, e, w in errors) / weight_sum


def _compute_statistics(reactions, intermediates_dict, dft_energies, damp, params):
    """Compute full statistics (WMAD, MAD, RMSD, max error)."""
    errors = _reaction_errors(reactions, intermediates_dict, dft_energies, damp, params)
    if not errors:
        return 0.0, 0.0, 0.0, 0.0

    abs_errors = [abs(e) for _, e, _ in errors]
    weights = [w for _, _, w in errors]
    weight_sum = sum(weights)

    wmad = sum(w * ae for w, ae in zip(weights, abs_errors)) / weight_sum if weight_sum > 0 else 0.0
    mad = np.mean(abs_errors)
    rmsd = np.sqrt(np.mean([e ** 2 for _, e, _ in errors]))
    max_err = max(abs_errors)

    return wmad, mad, rmsd, max_err


# ---------------------------------------------------------------------------
# Train/test splitting
# ---------------------------------------------------------------------------

def train_test_split(reactions, test_fraction, seed):
    """Stratified train/test split by subset name.

    Parameters
    ----------
    reactions : list of Reaction
    test_fraction : float
        Fraction for test set (0.0 to 0.5).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    (train, test) : tuple of lists
    """
    if test_fraction <= 0:
        return list(reactions), []

    rng = np.random.default_rng(seed)

    by_subset = {}
    for rxn in reactions:
        by_subset.setdefault(rxn.subset, []).append(rxn)

    train, test = [], []
    for subset_name in sorted(by_subset):
        subset_rxns = by_subset[subset_name]
        n_test = round(len(subset_rxns) * test_fraction)
        if n_test == 0:
            # Subset too small for requested fraction — keep all in train
            train.extend(subset_rxns)
            continue
        if n_test >= len(subset_rxns):
            n_test = max(1, len(subset_rxns) // 2)
        indices = rng.permutation(len(subset_rxns))
        for i, idx in enumerate(indices):
            if i < n_test:
                test.append(subset_rxns[idx])
            else:
                train.append(subset_rxns[idx])

    return train, test


def k_fold_split(reactions, n_folds, seed):
    """Generate k-fold cross-validation splits, stratified by subset.

    Returns
    -------
    list of (train, test) tuples
    """
    rng = np.random.default_rng(seed)

    by_subset = {}
    for rxn in reactions:
        by_subset.setdefault(rxn.subset, []).append(rxn)

    # Assign each reaction to a fold
    fold_assignments = {}
    for subset_name in sorted(by_subset):
        subset_rxns = by_subset[subset_name]
        indices = rng.permutation(len(subset_rxns))
        for i, idx in enumerate(indices):
            fold_assignments[id(subset_rxns[idx])] = i % n_folds

    folds = []
    for fold_idx in range(n_folds):
        test = [r for r in reactions if fold_assignments[id(r)] == fold_idx]
        train = [r for r in reactions if fold_assignments[id(r)] != fold_idx]
        folds.append((train, test))

    return folds


# ---------------------------------------------------------------------------
# Main fitting function
# ---------------------------------------------------------------------------

def fit_d3_params(dataset, dft_energies, config=None):
    """Optimize D3 parameters against a benchmark dataset.

    Parameters
    ----------
    dataset : Dataset
        Benchmark dataset with reactions and species.
    dft_energies : dict of str to float
        DFT energies in Hartree, keyed by species name.
    config : OptimizeConfig, optional
        Optimizer configuration.

    Returns
    -------
    FitResult
    """
    from scipy.optimize import differential_evolution

    if config is None:
        config = OptimizeConfig()

    t0 = time.perf_counter()

    # Precompute D3 intermediates for all species
    intermediates = {}
    for sp_name, sp in dataset.species.items():
        intermediates[sp_name] = precompute_d3(sp)

    # Train/test split
    train_rxns, test_rxns = train_test_split(
        dataset.reactions, config.test_fraction, config.random_seed)

    # Set up bounds
    if config.damp == "bj":
        bounds = [config.s8_bounds, config.a1_bounds, config.a2_bounds]
    else:
        bounds = [config.rs6_bounds, config.s8_bounds]

    # Objective: weighted MAD on training set
    def objective(params):
        return _weighted_mad(train_rxns, intermediates, dft_energies, config.damp, params)

    # Run optimizer
    result = differential_evolution(
        objective,
        bounds=bounds,
        maxiter=config.maxiter,
        tol=config.tol,
        popsize=config.popsize,
        workers=config.workers,
        seed=config.random_seed,
    )

    elapsed = time.perf_counter() - t0

    # Extract parameters
    if config.damp == "bj":
        s8, a1, a2 = result.x
        params_dict = {"s6": 1.0, "s8": round(s8, 4), "a1": round(a1, 4), "a2": round(a2, 4)}
    else:
        rs6, s8 = result.x
        params_dict = {"s6": 1.0, "rs6": round(rs6, 4), "s8": round(s8, 4)}

    # Compute statistics
    train_wmad, train_mad, train_rmsd, train_max = _compute_statistics(
        train_rxns, intermediates, dft_energies, config.damp, result.x)

    test_wmad = test_mad = test_rmsd = test_max = None
    if test_rxns:
        test_wmad, test_mad, test_rmsd, test_max = _compute_statistics(
            test_rxns, intermediates, dft_energies, config.damp, result.x)

    # Per-reaction errors on full dataset
    all_errors = _reaction_errors(
        dataset.reactions, intermediates, dft_energies, config.damp, result.x)
    per_rxn = [
        {"subset": rxn.subset, "index": rxn.index, "error": err, "weight": w}
        for rxn, err, w in all_errors
    ]

    return FitResult(
        damp=config.damp,
        params=params_dict,
        weighted_mad=train_wmad,
        mad=train_mad,
        rmsd=train_rmsd,
        max_error=train_max,
        test_weighted_mad=test_wmad,
        test_mad=test_mad,
        test_rmsd=test_rmsd,
        test_max_error=test_max,
        n_train=len(train_rxns),
        n_test=len(test_rxns),
        elapsed_seconds=elapsed,
        per_reaction_errors=per_rxn,
    )


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def _published_params_array(functional, damp):
    """Return published parameter array for the functional, or None."""
    if not functional:
        return None
    from .pars import bj_parms, zero_parms, resolve_functional
    canonical = resolve_functional(functional)
    if canonical is None:
        return None
    if damp == "bj" and canonical in bj_parms:
        p = bj_parms[canonical]
        return [p[2], p[1], p[3]]  # s8, a1, a2
    elif damp == "zero" and canonical in zero_parms:
        p = zero_parms[canonical]
        return [p[1], p[2]]  # rs6, s8
    return None


def compute_published_stats(dataset, dft_energies, functional, damp):
    """Compute statistics for published parameters, if they exist.

    Returns a tuple (wmad, mad, rmsd, max_err, n) or None.
    """
    pub_params = _published_params_array(functional, damp)
    if pub_params is None:
        return None
    intermediates = {}
    for sp_name, sp in dataset.species.items():
        intermediates[sp_name] = precompute_d3(sp)
    wmad, mad, rmsd, max_err = _compute_statistics(
        dataset.reactions, intermediates, dft_energies, damp, pub_params)
    return (wmad, mad, rmsd, max_err, len(dataset.reactions))


def print_results(result, verbose=False, functional=None, basis=None,
                  published_stats=None):
    """Print formatted optimization results."""
    print()
    print("   "+"=" * 72)
    print("   D3 Parameter Optimization Results")
    print("   "+"=" * 72)
    print()

    if functional:
        print(f"   Functional:        {functional}")
    if basis:
        print(f"   Basis set:         {basis}")
    damp_label = "Becke-Johnson (BJ)" if result.damp == "bj" else "Zero"
    print(f"   Damping scheme:    {damp_label}")
    print(f"   Train / Test:      {result.n_train} / {result.n_test}")
    print(f"   Elapsed time:      {result.elapsed_seconds:.1f} s")

    print()
    print("   "+"-" * 72)
    print("   Optimized Parameters")
    print("   "+"-" * 72)

    # Look up published parameters for comparison
    published = None
    if functional:
        from .pars import bj_parms, zero_parms, resolve_functional
        canonical = resolve_functional(functional)
        if canonical is not None and result.damp == "bj" and canonical in bj_parms:
            p = bj_parms[canonical]
            published = {"s6": p[0], "a1": p[1], "s8": p[2], "a2": p[3]}
        elif canonical is not None and result.damp == "zero" and canonical in zero_parms:
            p = zero_parms[canonical]
            published = {"s6": p[0], "rs6": p[1], "s8": p[2]}

    if published:
        print(f"   {'':4s}   {'Optimized':>10s}  {'Published':>10s}")
        for key, val in result.params.items():
            fixed = "  (fixed)" if key == "s6" else ""
            pub_val = published.get(key)
            print(f"   {key:4s}   {val:10.4f}  {pub_val:10.4f}{fixed}")
    else:
        for key, val in result.params.items():
            fixed = "  (fixed)" if key == "s6" else ""
            print(f"   {key:4s} = {val:.4f}{fixed}")
        if functional:
            print(f"\n   No published {damp_label} parameters found for {functional}")

    # CLI command
    cli_parts = " ".join(f"--{k} {v}" for k, v in result.params.items())
    print(f"\n   pydftd3 CLI:  {cli_parts}")

    print()
    stat_label = "Statistics (kcal/mol)"
    stat_header = f"   {stat_label:<24s} {'N':>5s} {'WMAD':>8s} {'MAD':>8s} {'RMSD':>8s} {'Max Err':>8s}"
    print("   " + "-" * 72)
    print(stat_header)
    print("   " + "-" * 72)
    print(f"   {'Train:':<24s} {result.n_train:5d} {result.weighted_mad:8.3f} {result.mad:8.3f} "
          f"{result.rmsd:8.3f} {result.max_error:8.3f}")
    if result.test_weighted_mad is not None:
        print(f"   {'Test:':<24s} {result.n_test:5d} {result.test_weighted_mad:8.3f} {result.test_mad:8.3f} "
              f"{result.test_rmsd:8.3f} {result.test_max_error:8.3f}")
    if published_stats is not None:
        pw, pm, pr, px, pn = published_stats
        print(f"   {'Published:':<24s} {pn:5d} {pw:8.3f} {pm:8.3f} {pr:8.3f} {px:8.3f}")

    # Per-subset breakdown
    if verbose:
        print()
        print("   " + "-" * 72)
        print(f"   {'Per-Subset Breakdown':<24s} {'N':>5s} {'WMAD':>8s} {'MAD':>8s} {'RMSD':>8s} {'Max Err':>8s}")
        print("   " + "-" * 72)

        by_subset = {}
        for entry in result.per_reaction_errors:
            by_subset.setdefault(entry["subset"], []).append(entry)

        for subset_name in sorted(by_subset):
            entries = by_subset[subset_name]
            n = len(entries)
            abs_errors = [abs(e["error"]) for e in entries]
            weights = [e["weight"] for e in entries]
            w_sum = sum(weights)
            wmad = sum(w * ae for w, ae in zip(weights, abs_errors)) / w_sum if w_sum > 0 else 0.0
            mad = np.mean(abs_errors)
            rmsd = np.sqrt(np.mean([e["error"] ** 2 for e in entries]))
            max_err = max(abs_errors)
            print(f"   {subset_name:<24s} {n:5d} {wmad:8.3f} {mad:8.3f} {rmsd:8.3f} {max_err:8.3f}")

    print("   "+"=" * 72)


# ---------------------------------------------------------------------------
# CLI subcommands
# ---------------------------------------------------------------------------

def optimize_main(argv):
    """Entry point for 'pydftd3 optimize' subcommands."""
    parser = ArgumentParser(
        prog="pydftd3 optimize",
        description="D3 parameter optimization tools.",
    )
    subparsers = parser.add_subparsers(dest="command")

    # --- prep subcommand ---
    prep = subparsers.add_parser("prep", help="Generate ORCA input files from benchmark dataset.")
    prep.add_argument("dataset", help="Path to DietGMTKN55 YAML file")
    prep.add_argument("-o", "--output-dir", required=True, help="Directory for ORCA input files")
    prep.add_argument("--func", required=True, help="DFT functional for ORCA")
    prep.add_argument("--basis", default="def2-TZVP", help="Basis set (default: def2-TZVP)")
    prep.add_argument("--nprocs", type=int, default=8, help="ORCA parallel processes (default: 8)")
    prep.add_argument("--maxcore", type=int, default=4000, help="Memory per core in MB (default: 4000)")
    prep.add_argument("--extra", default="", help="Additional ORCA keywords")

    # --- fit subcommand ---
    fit = subparsers.add_parser("fit", help="Optimize D3 parameters against benchmark data.")
    fit.add_argument("dataset", help="Path to DietGMTKN55 YAML file")
    fit.add_argument("--energies", required=True,
                     help="CSV file with DFT energies, or directory of ORCA .out files")
    fit.add_argument("--damp", default="bj", choices=("zero", "bj"),
                     help="Damping scheme (default: bj)")
    fit.add_argument("--test-frac", type=float, default=0.2,
                     help="Fraction for test set (default: 0.2, 0=no split)")
    fit.add_argument("--folds", type=int, default=0,
                     help="Cross-validation folds (0=single split, default: 0)")
    fit.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    fit.add_argument("--maxiter", type=int, default=1000, help="Max optimizer iterations (default: 1000)")
    fit.add_argument("--workers", type=int, default=1, help="Parallel workers (-1=all cores)")
    fit.add_argument("--skip-missing", action="store_true",
                     help="Skip reactions with missing DFT energies instead of aborting")
    fit.add_argument("--save-csv", default=None, help="Save per-reaction errors to CSV")
    fit.add_argument("-v", dest="verbose", action="store_true", default=False,
                     help="Print per-subset breakdown")

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "prep":
        return _cmd_prep(args)
    elif args.command == "fit":
        return _cmd_fit(args)
    return 1


def _cmd_prep(args):
    """Handle 'pydftd3 optimize prep' command."""
    dataset = load_dataset(args.dataset)
    print(f"\n   Loaded {len(dataset.reactions)} reactions, {len(dataset.species)} unique species")

    paths = generate_orca_inputs(
        dataset, args.output_dir, args.func,
        basis=args.basis, nprocs=args.nprocs,
        maxcore=args.maxcore, extra_keywords=args.extra,
    )
    print(f"Generated {len(paths)} ORCA input files in {args.output_dir}/")
    return 0


def _cmd_fit(args):
    """Handle 'pydftd3 optimize fit' command."""
    dataset = load_dataset(args.dataset)
    print(f"\n   Loaded {len(dataset.reactions)} reactions, {len(dataset.species)} unique species")

    # Load DFT energies
    detected_functional = None
    detected_basis = None
    if args.energies.endswith(".csv"):
        dft_energies, detected_functional, detected_basis = read_energies_csv(args.energies)
    else:
        dft_energies, detected_functional, detected_basis = parse_orca_outputs(args.energies)
        csv_path = os.path.join(args.energies, "dft_energies.csv")
        write_energies_csv(dft_energies, csv_path,
                           functional=detected_functional, basis=detected_basis)
        print(f"\n   Saved DFT energies to {csv_path}")

    # Validate completeness
    missing = set(dataset.species.keys()) - set(dft_energies.keys())
    if missing:
        if not args.skip_missing:
            print(f"\n   Error: Missing DFT energies for {len(missing)} species:")
            for name in sorted(missing)[:10]:
                print(f"  - {name}")
            if len(missing) > 10:
                print(f"  ... and {len(missing) - 10} more")
            print("\n   Use --skip-missing to drop reactions that reference these species.")
            return 1

        # Filter out reactions that use missing species
        n_before = len(dataset.reactions)
        dataset.reactions = [
            rxn for rxn in dataset.reactions
            if not (set(rxn.stoichiometry.keys()) & missing)
        ]
        n_dropped = n_before - len(dataset.reactions)
        print(f"   Skipping {len(missing)} missing species, dropped {n_dropped}/{n_before} reactions")

    print(f"   Loaded DFT energies for {len(dft_energies)} species, {len(dataset.reactions)} reactions")

    config = OptimizeConfig(
        damp=args.damp,
        test_fraction=args.test_frac,
        random_seed=args.seed,
        n_folds=args.folds,
        maxiter=args.maxiter,
        workers=args.workers,
    )

    pub_stats = compute_published_stats(
        dataset, dft_energies, detected_functional, args.damp)

    if args.folds > 0:
        _cmd_fit_cv(dataset, dft_energies, config, args,
                    functional=detected_functional, basis=detected_basis,
                    published_stats=pub_stats)
    else:
        result = fit_d3_params(dataset, dft_energies, config)
        print_results(result, verbose=args.verbose,
                      functional=detected_functional, basis=detected_basis,
                      published_stats=pub_stats)

        if args.save_csv:
            _save_errors_csv(result, args.save_csv)

    return 0


def _fit_one_fold(fold_idx, train_rxns, test_rxns, intermediates, dft_energies, config):
    """Fit a single CV fold. Designed to be called in parallel."""
    from scipy.optimize import differential_evolution

    if config.damp == "bj":
        bounds = [config.s8_bounds, config.a1_bounds, config.a2_bounds]
    else:
        bounds = [config.rs6_bounds, config.s8_bounds]

    def objective(params):
        return _weighted_mad(train_rxns, intermediates, dft_energies, config.damp, params)

    result = differential_evolution(
        objective, bounds=bounds,
        maxiter=config.maxiter, tol=config.tol,
        popsize=config.popsize, seed=config.random_seed + fold_idx,
    )

    test_wmad, test_mad, test_rmsd, test_max = _compute_statistics(
        test_rxns, intermediates, dft_energies, config.damp, result.x)
    return (result.x, test_wmad, test_mad, test_rmsd)


def _cmd_fit_cv(dataset, dft_energies, config, args, functional=None, basis=None,
                published_stats=None):
    """Run k-fold cross-validation."""
    from concurrent.futures import ProcessPoolExecutor, as_completed

    folds = k_fold_split(dataset.reactions, config.n_folds, config.random_seed)

    # Precompute intermediates once
    intermediates = {}
    for sp_name, sp in dataset.species.items():
        intermediates[sp_name] = precompute_d3(sp)

    n_workers = args.workers if args.workers > 0 else None  # None = all cores
    if n_workers == 1 or config.n_folds == 1:
        # Sequential fallback
        fold_results = []
        for fold_idx, (train_rxns, test_rxns) in enumerate(folds):
            r = _fit_one_fold(fold_idx, train_rxns, test_rxns,
                              intermediates, dft_energies, config)
            fold_results.append(r)
            print(f"   Fold {fold_idx + 1}/{config.n_folds} complete (WMAD: {r[1]:.3f})")
    else:
        print(f"   Running {config.n_folds} folds in parallel ({n_workers or 'all'} workers)...")
        fold_results = [None] * config.n_folds
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {}
            for fold_idx, (train_rxns, test_rxns) in enumerate(folds):
                fut = executor.submit(
                    _fit_one_fold, fold_idx, train_rxns, test_rxns,
                    intermediates, dft_energies, config)
                futures[fut] = fold_idx
            for fut in as_completed(futures):
                idx = futures[fut]
                fold_results[idx] = fut.result()
                print(f"   Fold {idx + 1}/{config.n_folds} complete (WMAD: {fold_results[idx][1]:.3f})")

    # Print CV results
    print()
    print("   "+"=" * 72)
    print(f"   Cross-Validation Results ({config.n_folds} folds)")
    print("    "+"=" * 72)
    print(f"\n   {'Fold':>4s}  {'WMAD':>8s}  {'MAD':>8s}  {'RMSD':>8s}")

    wmads = [r[1] for r in fold_results]
    mads = [r[2] for r in fold_results]
    rmsds = [r[3] for r in fold_results]

    for i, (_, wmad, mad, rmsd) in enumerate(fold_results):
        print(f"   {i + 1:4d}  {wmad:8.3f}  {mad:8.3f}  {rmsd:8.3f}")

    print(f"   {'Mean':>4s}  {np.mean(wmads):8.3f}  {np.mean(mads):8.3f}  {np.mean(rmsds):8.3f}")
    print(f"   {'Std':>4s}  {np.std(wmads):8.3f}  {np.std(mads):8.3f}  {np.std(rmsds):8.3f}")

    # Final fit on all data
    print("\n   Fitting on all data for final parameters...")
    config_all = OptimizeConfig(
        damp=config.damp,
        test_fraction=0.0,
        random_seed=config.random_seed,
        maxiter=config.maxiter,
        tol=config.tol,
        popsize=config.popsize,
    )
    final_result = fit_d3_params(dataset, dft_energies, config_all)
    print_results(final_result, verbose=args.verbose,
                  functional=functional, basis=basis,
                  published_stats=published_stats)


def _save_errors_csv(result, csv_path):
    """Save per-reaction errors to CSV."""
    import csv
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["subset", "index", "error_kcal", "weight"])
        for entry in result.per_reaction_errors:
            writer.writerow([
                entry["subset"], entry["index"],
                f"{entry['error']:.6f}", f"{entry['weight']:.4f}",
            ])
    print(f"Saved per-reaction errors to {csv_path}")
