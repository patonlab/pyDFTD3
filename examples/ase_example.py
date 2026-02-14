"""
Example: Using pyDFTD3 with an ASE Atoms object.

This script shows how to compute the D3(BJ) dispersion correction for a
water dimer built with ASE. Since pyDFTD3 expects cclib-style parsed data,
we create a lightweight adapter that wraps ASE coordinates into the
interface that CalcD3 expects (atomnos and atomcoords attributes).

Requirements:
    pip install ase
"""

import numpy as np

from ase.build import molecule


class ASEAdapter:
    """Wrap an ASE Atoms object to look like cclib parsed data."""

    def __init__(self, atoms):
        self.atomnos = np.array(atoms.get_atomic_numbers())
        # cclib stores coordinates as (n_steps, n_atoms, 3); use a single step
        self.atomcoords = np.array([atoms.get_positions()])


def main():
    from dftd3.dftd3 import CalcD3, AUTOKCAL

    # --- Build a water dimer from ASE ---
    water = molecule("H2O")
    # Create a second water molecule shifted 3 Angstrom along x
    water2 = molecule("H2O")
    water2.translate([3.0, 0.0, 0.0])
    dimer = water + water2

    print("Water dimer geometry (ASE):")
    for sym, pos in zip(dimer.get_chemical_symbols(), dimer.get_positions()):
        print(f"  {sym:2s}  {pos[0]:10.5f} {pos[1]:10.5f} {pos[2]:10.5f}")

    # --- Compute D3(BJ) correction for B3LYP ---
    data = ASEAdapter(dimer)
    result = CalcD3(data, "B3LYP", damp="bj")

    total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
    total_kcal = result.attractive_r6_vdw + result.attractive_r8_vdw

    print(f"\nD3(BJ)/B3LYP dispersion correction:")
    print(f"  D3(R6)  = {result.attractive_r6_vdw / AUTOKCAL:12.8f} au")
    print(f"  D3(R8)  = {result.attractive_r8_vdw / AUTOKCAL:12.8f} au")
    print(f"  Total   = {total_au:12.8f} au  ({total_kcal:.4f} kcal/mol)")


if __name__ == "__main__":
    main()
