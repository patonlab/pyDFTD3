"""Benchmark D3(BJ) dispersion energy and timing as a function of cutoff radius.

Runs pyDFTD3 on large molecules (maitotoxin, 3I40 protein, 1YA5 protein) with
varying cutoffs to show convergence and performance tradeoffs.

Usage:
    python examples/benchmark_cutoff.py
"""

import os
import sys
import time

# Add parent directory to path so we can import dftd3
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from dftd3.dftd3 import CalcD3, read_file

EXAMPLES_DIR = os.path.dirname(__file__)

CUTOFFS = [10, 15, 25, None]

BENCHMARKS = [
    ("maitotoxin.xyz", "B3LYP", 285),
    ("3I40.pdb", "B3LYP", None),
    ("1YA5.pdb", "B3LYP", None),
]


def run_benchmark(filename, functional, cutoffs):
    filepath = os.path.join(EXAMPLES_DIR, filename)
    data = read_file(filepath)
    if data is None:
        print(f"  Could not parse {filepath}, skipping")
        return None

    natom = len(data.atomnos)
    results = []

    for cutoff in cutoffs:
        label = f"{cutoff} A" if cutoff is not None else "None"
        t0 = time.perf_counter()
        result = CalcD3(data, functional, damp="bj", cutoff=cutoff)
        elapsed = time.perf_counter() - t0
        total_kcal = result.attractive_r6_vdw + result.attractive_r8_vdw
        results.append((label, total_kcal, elapsed))

    return natom, results


def main():
    print("# D3(BJ)/B3LYP Cutoff Radius Benchmark\n")

    for filename, functional, expected_natom in BENCHMARKS:
        print(f"## {filename}\n")
        out = run_benchmark(filename, functional, CUTOFFS)
        if out is None:
            continue
        natom, results = out
        print(f"Atoms: {natom}\n")

        # Reference = no cutoff (last entry)
        ref_energy = results[-1][1]

        print(f"| {'Cutoff':>8} | {'Edisp (kcal/mol)':>18} | {'Error (kcal/mol)':>18} | {'Time (s)':>10} |")
        print(f"|{'-' * 10}|{'-' * 20}|{'-' * 20}|{'-' * 12}|")
        for label, energy, elapsed in results:
            error = energy - ref_energy
            print(f"| {label:>8} | {energy:>18.3f} | {error:>+18.3f} | {elapsed:>10.3f} |")
        print()

    print("Reference energy = no cutoff (all pairwise interactions included).")


if __name__ == "__main__":
    main()
