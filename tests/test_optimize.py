"""Tests for D3 parameter optimization (benchmark + optimize modules)."""

import os

import numpy as np
import pytest

from dftd3.benchmark import (
    Dataset,
    Reaction,
    Species,
    _ELEMENT_TO_Z,
    generate_orca_inputs,
    load_dataset,
    load_mpconf196,
    load_nenci,
    read_energies_csv,
    resolve_dataset,
    write_energies_csv,
)
from dftd3.dftd3 import AUTOKCAL, CalcD3
from dftd3.optimize import (
    OptimizeConfig,
    d3_energy_bj,
    d3_energy_zero,
    fit_d3_params,
    precompute_d3,
    train_test_split,
)

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

WATER_SPECIES = Species(
    name="water",
    charge=0,
    uhf=0,
    elements=["O", "H", "H"],
    positions=[[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]],
)

METHANE_SPECIES = Species(
    name="methane",
    charge=0,
    uhf=0,
    elements=["C", "H", "H", "H", "H"],
    positions=[
        [0.0, 0.0, 0.0],
        [0.6276, 0.6276, 0.6276],
        [-0.6276, -0.6276, 0.6276],
        [-0.6276, 0.6276, -0.6276],
        [0.6276, -0.6276, -0.6276],
    ],
)

MINI_YAML = """\
TestSubsetA:
  1:
    Energy: 3.5
    Weight: 1.0
    Species:
      water_a:
        Count: -1
        Charge: 0
        UHF: 0
        Number: 3
        Elements: [O, H, H]
        Positions:
          - [0.0, 0.0, 0.1173]
          - [0.0, 0.7572, -0.4692]
          - [0.0, -0.7572, -0.4692]
      water_b:
        Count: 1
        Charge: 0
        UHF: 0
        Number: 3
        Elements: [O, H, H]
        Positions:
          - [0.0, 0.0, 0.12]
          - [0.0, 0.76, -0.47]
          - [0.0, -0.76, -0.47]
  2:
    Energy: -1.2
    Weight: 2.0
    Species:
      methane_a:
        Count: -1
        Charge: 0
        UHF: 0
        Number: 5
        Elements: [C, H, H, H, H]
        Positions:
          - [0.0, 0.0, 0.0]
          - [0.6276, 0.6276, 0.6276]
          - [-0.6276, -0.6276, 0.6276]
          - [-0.6276, 0.6276, -0.6276]
          - [0.6276, -0.6276, -0.6276]
      methane_b:
        Count: 1
        Charge: 0
        UHF: 0
        Number: 5
        Elements: [C, H, H, H, H]
        Positions:
          - [0.0, 0.0, 0.0]
          - [0.63, 0.63, 0.63]
          - [-0.63, -0.63, 0.63]
          - [-0.63, 0.63, -0.63]
          - [0.63, -0.63, -0.63]
TestSubsetB:
  1:
    Energy: 0.8
    Weight: 1.5
    Species:
      water_c:
        Count: -1
        Charge: 0
        UHF: 0
        Number: 3
        Elements: [O, H, H]
        Positions:
          - [0.0, 0.0, 0.115]
          - [0.0, 0.755, -0.465]
          - [0.0, -0.755, -0.465]
      water_d:
        Count: 1
        Charge: 0
        UHF: 0
        Number: 3
        Elements: [O, H, H]
        Positions:
          - [0.0, 0.0, 0.118]
          - [0.0, 0.758, -0.468]
          - [0.0, -0.758, -0.468]
"""


@pytest.fixture
def mini_yaml(tmp_path):
    yaml_file = tmp_path / "test_dataset.yaml"
    yaml_file.write_text(MINI_YAML)
    return str(yaml_file)


# ---------------------------------------------------------------------------
# Species dataclass
# ---------------------------------------------------------------------------

class TestSpecies:
    def test_atomnos(self):
        assert list(WATER_SPECIES.atomnos) == [8, 1, 1]

    def test_atomcoords_shape(self):
        assert WATER_SPECIES.atomcoords.shape == (1, 3, 3)

    def test_multiplicity_singlet(self):
        assert WATER_SPECIES.multiplicity == 1

    def test_multiplicity_doublet(self):
        sp = Species("rad", 0, 1, ["O", "H"], [[0, 0, 0], [1, 0, 0]])
        assert sp.multiplicity == 2

    def test_natom(self):
        assert WATER_SPECIES.natom == 3
        assert METHANE_SPECIES.natom == 5


# ---------------------------------------------------------------------------
# Dataset loading
# ---------------------------------------------------------------------------

class TestDatasetLoading:
    def test_load_reaction_count(self, mini_yaml):
        ds = load_dataset(mini_yaml)
        assert len(ds.reactions) == 3

    def test_load_species_count(self, mini_yaml):
        ds = load_dataset(mini_yaml)
        # 2 per reaction × 3 reactions, but names are subset-qualified so all unique
        assert len(ds.species) == 6

    def test_reaction_stoichiometry(self, mini_yaml):
        ds = load_dataset(mini_yaml)
        rxn = [r for r in ds.reactions if r.subset == "TestSubsetA" and r.index == 1][0]
        assert rxn.stoichiometry["TestSubsetA-water_a"] == -1
        assert rxn.stoichiometry["TestSubsetA-water_b"] == 1

    def test_reference_energy(self, mini_yaml):
        ds = load_dataset(mini_yaml)
        rxn = [r for r in ds.reactions if r.subset == "TestSubsetA" and r.index == 1][0]
        assert rxn.reference_energy == 3.5

    def test_weight(self, mini_yaml):
        ds = load_dataset(mini_yaml)
        rxn = [r for r in ds.reactions if r.subset == "TestSubsetA" and r.index == 2][0]
        assert rxn.weight == 2.0

    def test_species_atomnos(self, mini_yaml):
        ds = load_dataset(mini_yaml)
        sp = ds.species["TestSubsetA-water_a"]
        assert list(sp.atomnos) == [8, 1, 1]

    def test_species_positions(self, mini_yaml):
        ds = load_dataset(mini_yaml)
        sp = ds.species["TestSubsetA-water_a"]
        assert sp.positions[0] == [0.0, 0.0, 0.1173]


# ---------------------------------------------------------------------------
# D3 precomputation
# ---------------------------------------------------------------------------

class TestPrecomputation:
    def test_bj_matches_calcd3(self):
        """D3(BJ) via precomputed intermediates must match CalcD3."""
        im = precompute_d3(WATER_SPECIES)
        s6, a1, s8, a2 = 1.0, 0.3981, 1.9889, 4.4211
        fast = d3_energy_bj(im, s6, s8, a1, a2)

        result = CalcD3(WATER_SPECIES, "B3LYP", damp="bj")
        calcd3 = result.attractive_r6_vdw + result.attractive_r8_vdw

        assert fast == pytest.approx(calcd3, abs=1e-8)

    def test_zero_matches_calcd3(self):
        """D3(zero) via precomputed intermediates must match CalcD3."""
        im = precompute_d3(WATER_SPECIES)
        s6, rs6, s8 = 1.0, 1.261, 1.703
        fast = d3_energy_zero(im, s6, rs6, s8)

        result = CalcD3(WATER_SPECIES, "B3LYP", damp="zero")
        calcd3 = result.attractive_r6_vdw + result.attractive_r8_vdw

        assert fast == pytest.approx(calcd3, abs=1e-8)

    def test_bj_methane_matches_calcd3(self):
        """Methane: D3(BJ) via intermediates matches CalcD3."""
        im = precompute_d3(METHANE_SPECIES)
        s6, a1, s8, a2 = 1.0, 0.3981, 1.9889, 4.4211
        fast = d3_energy_bj(im, s6, s8, a1, a2)

        result = CalcD3(METHANE_SPECIES, "B3LYP", damp="bj")
        calcd3 = result.attractive_r6_vdw + result.attractive_r8_vdw

        assert fast == pytest.approx(calcd3, abs=1e-8)

    def test_zero_methane_matches_calcd3(self):
        """Methane: D3(zero) via intermediates matches CalcD3."""
        im = precompute_d3(METHANE_SPECIES)
        s6, rs6, s8 = 1.0, 1.261, 1.703
        fast = d3_energy_zero(im, s6, rs6, s8)

        result = CalcD3(METHANE_SPECIES, "B3LYP", damp="zero")
        calcd3 = result.attractive_r6_vdw + result.attractive_r8_vdw

        assert fast == pytest.approx(calcd3, abs=1e-8)

    def test_pbe_bj_matches_calcd3(self):
        """PBE/BJ via intermediates matches CalcD3."""
        im = precompute_d3(WATER_SPECIES)
        # PBE BJ params: s6=1.0, a1=0.4289, s8=0.7875, a2=4.4407
        from dftd3.pars import bj_parms
        s6, a1, s8, a2 = bj_parms["PBE"]
        fast = d3_energy_bj(im, s6, s8, a1, a2)

        result = CalcD3(WATER_SPECIES, "PBE", damp="bj")
        calcd3 = result.attractive_r6_vdw + result.attractive_r8_vdw

        assert fast == pytest.approx(calcd3, abs=1e-8)


# ---------------------------------------------------------------------------
# Train/test split
# ---------------------------------------------------------------------------

class TestTrainTestSplit:
    def test_sizes(self):
        rxns = [Reaction("A", i, 1.0, 1.0, {"x": 1}) for i in range(20)]
        train, test = train_test_split(rxns, 0.2, seed=42)
        assert len(train) + len(test) == 20
        assert len(test) >= 1

    def test_reproducibility(self):
        rxns = [Reaction("A", i, 1.0, 1.0, {"x": 1}) for i in range(20)]
        t1, _ = train_test_split(rxns, 0.2, seed=42)
        t2, _ = train_test_split(rxns, 0.2, seed=42)
        assert [r.index for r in t1] == [r.index for r in t2]

    def test_different_seeds(self):
        rxns = [Reaction("A", i, 1.0, 1.0, {"x": 1}) for i in range(20)]
        t1, _ = train_test_split(rxns, 0.2, seed=42)
        t2, _ = train_test_split(rxns, 0.2, seed=99)
        assert [r.index for r in t1] != [r.index for r in t2]

    def test_stratification(self):
        rxns = (
            [Reaction("A", i, 1.0, 1.0, {"x": 1}) for i in range(10)]
            + [Reaction("B", i, 1.0, 1.0, {"x": 1}) for i in range(10)]
        )
        train, test = train_test_split(rxns, 0.2, seed=42)
        test_subsets = {r.subset for r in test}
        assert "A" in test_subsets and "B" in test_subsets

    def test_zero_fraction_no_test(self):
        rxns = [Reaction("A", i, 1.0, 1.0, {"x": 1}) for i in range(10)]
        train, test = train_test_split(rxns, 0.0, seed=42)
        assert len(train) == 10
        assert len(test) == 0


# ---------------------------------------------------------------------------
# CSV I/O
# ---------------------------------------------------------------------------

class TestCSVIO:
    def test_roundtrip(self, tmp_path):
        energies = {"mol_A": -76.123456789012, "mol_B": -152.789012345678}
        csv_path = str(tmp_path / "energies.csv")
        write_energies_csv(energies, csv_path)
        loaded, _, _ = read_energies_csv(csv_path)
        assert loaded["mol_A"] == pytest.approx(-76.123456789012, abs=1e-10)
        assert loaded["mol_B"] == pytest.approx(-152.789012345678, abs=1e-10)

    def test_csv_has_header(self, tmp_path):
        energies = {"mol_A": -76.0}
        csv_path = str(tmp_path / "energies.csv")
        write_energies_csv(energies, csv_path)
        with open(csv_path) as f:
            header = f.readline().strip()
        assert header == "species_name,energy_hartree"


# ---------------------------------------------------------------------------
# ORCA input generation
# ---------------------------------------------------------------------------

class TestOrcaInputGeneration:
    def test_generates_files(self, tmp_path):
        ds = Dataset("test", [], {"water": WATER_SPECIES})
        paths = generate_orca_inputs(ds, str(tmp_path), "PBE", basis="def2-SVP")
        assert len(paths) == 1
        assert os.path.exists(paths[0])

    def test_correct_format(self, tmp_path):
        ds = Dataset("test", [], {"water": WATER_SPECIES})
        paths = generate_orca_inputs(ds, str(tmp_path), "PBE", basis="def2-SVP")
        content = open(paths[0]).read()
        assert "! PBE def2-SVP" in content
        assert "* xyz 0 1" in content
        assert "O" in content

    def test_no_d3_keywords(self, tmp_path):
        ds = Dataset("test", [], {"water": WATER_SPECIES})
        paths = generate_orca_inputs(ds, str(tmp_path), "B3LYP")
        content = open(paths[0]).read().upper()
        assert "D3ZERO" not in content
        assert "D3BJ" not in content
        assert "D3" not in content

    def test_extra_keywords(self, tmp_path):
        ds = Dataset("test", [], {"water": WATER_SPECIES})
        paths = generate_orca_inputs(ds, str(tmp_path), "PBE", extra_keywords="TightSCF RIJCOSX")
        content = open(paths[0]).read()
        assert "TightSCF RIJCOSX" in content

    def test_charge_and_multiplicity(self, tmp_path):
        sp = Species("radical", 1, 1, ["O", "H"], [[0, 0, 0], [1, 0, 0]])
        ds = Dataset("test", [], {"radical": sp})
        paths = generate_orca_inputs(ds, str(tmp_path), "PBE")
        content = open(paths[0]).read()
        assert "* xyz 1 2" in content  # charge=1, mult=uhf+1=2


# ---------------------------------------------------------------------------
# End-to-end parameter recovery (synthetic data)
# ---------------------------------------------------------------------------

class TestParameterRecovery:
    def test_recover_bj_params(self, mini_yaml):
        """Optimizer should recover known BJ params from synthetic DFT energies."""
        ds = load_dataset(mini_yaml)

        # Known BJ parameters to recover
        known_s8, known_a1, known_a2 = 1.9889, 0.3981, 4.4211
        s6 = 1.0

        # Precompute intermediates
        intermediates = {}
        for sp_name, sp in ds.species.items():
            intermediates[sp_name] = precompute_d3(sp)

        # Construct synthetic DFT energies: E_DFT = (E_ref - D3_rxn) / stoich
        # We set arbitrary base energies, then ensure reaction energies are consistent
        base_energies = {}
        for sp_name in ds.species:
            # Arbitrary base energy in Hartree
            base_energies[sp_name] = -76.0 + hash(sp_name) % 1000 * 0.001

        # For each reaction, compute the D3 reaction contribution
        # E_ref = E_DFT_rxn + E_D3_rxn  =>  E_DFT_rxn = E_ref - E_D3_rxn
        # We need to adjust base_energies so that the DFT reaction energies are correct
        dft_energies = dict(base_energies)

        # The optimizer minimizes |E_DFT_rxn + E_D3_rxn(params) - E_ref|
        # If we set E_DFT for species such that E_DFT_rxn = E_ref - E_D3_rxn(known_params),
        # then the optimizer should recover known_params.

        # For simplicity with multiple reactions sharing species, we use a least-squares
        # approach: pick one "free" species per reaction and solve.
        # Actually, simpler: just set DFT energies so that every reaction error is 0
        # for the known parameters. We can do this by adjusting one species per reaction.

        # Even simpler: generate fake DFT energies by computing D3 for each species
        # and setting E_DFT_species = E_base - D3(species) so that:
        # E_DFT_rxn + D3_rxn = sum(coeff * E_base) for all reactions
        # Then set E_ref = sum(coeff * E_base).
        # This makes E_ref - E_DFT_rxn - D3_rxn = 0 exactly.

        for sp_name, sp in ds.species.items():
            im = intermediates[sp_name]
            d3_sp = d3_energy_bj(im, s6, known_s8, known_a1, known_a2)
            # DFT energy = base - D3/AUTOKCAL (convert D3 from kcal to Hartree)
            dft_energies[sp_name] = base_energies[sp_name] - d3_sp / AUTOKCAL

        # Update reference energies to match
        for rxn in ds.reactions:
            rxn.reference_energy = sum(
                coeff * base_energies[sp_name] * AUTOKCAL
                for sp_name, coeff in rxn.stoichiometry.items()
            )

        # Now fit — the optimizer should recover the known parameters
        config = OptimizeConfig(
            damp="bj",
            test_fraction=0.0,
            random_seed=42,
            maxiter=500,
        )
        result = fit_d3_params(ds, dft_energies, config)

        # Check recovered parameters
        assert result.params["s8"] == pytest.approx(known_s8, abs=0.05)
        assert result.params["a1"] == pytest.approx(known_a1, abs=0.05)
        assert result.params["a2"] == pytest.approx(known_a2, abs=0.05)
        assert result.weighted_mad < 0.01  # Should be near-zero


# ---------------------------------------------------------------------------
# NENCI-2021 loading
# ---------------------------------------------------------------------------

# Helper to write a synthetic NENCI XYZ file
def _write_nenci_xyz(path, elements_A, positions_A, elements_B, positions_B,
                     ref_energy, charge_A=0, mult_A=1, charge_B=0, mult_B=1):
    """Write one synthetic NENCI-format XYZ file."""
    natoms = len(elements_A) + len(elements_B)
    dimer_charge = charge_A + charge_B
    dimer_mult = 1  # assume singlet dimer
    comment = (
        f"{dimer_charge} {dimer_mult} "
        f"{charge_A} {mult_A} {charge_B} {mult_B} "
        f"{len(elements_A)} {len(elements_B)} "
        f"{ref_energy}"
    )
    lines = [str(natoms), comment]
    for el, pos in zip(elements_A, positions_A):
        lines.append(f"{el}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")
    for el, pos in zip(elements_B, positions_B):
        lines.append(f"{el}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")
    lines.append("")
    with open(path, "w") as f:
        f.write("\n".join(lines))


@pytest.fixture
def nenci_dir(tmp_path):
    """Create a temp directory with 3 synthetic NENCI XYZ files."""
    xyz_dir = tmp_path / "nenci_xyz"
    xyz_dir.mkdir()

    # Water dimer (config 1)
    _write_nenci_xyz(
        str(xyz_dir / "water_water_01.xyz"),
        elements_A=["O", "H", "H"],
        positions_A=[[0.0, 0.0, 0.0], [0.0, 0.76, -0.47], [0.0, -0.76, -0.47]],
        elements_B=["O", "H", "H"],
        positions_B=[[3.0, 0.0, 0.0], [3.0, 0.76, -0.47], [3.0, -0.76, -0.47]],
        ref_energy=-4.97,
    )

    # Water dimer (config 2 — stretched)
    _write_nenci_xyz(
        str(xyz_dir / "water_water_02.xyz"),
        elements_A=["O", "H", "H"],
        positions_A=[[0.0, 0.0, 0.0], [0.0, 0.76, -0.47], [0.0, -0.76, -0.47]],
        elements_B=["O", "H", "H"],
        positions_B=[[5.0, 0.0, 0.0], [5.0, 0.76, -0.47], [5.0, -0.76, -0.47]],
        ref_energy=-1.23,
    )

    # Methane-water
    _write_nenci_xyz(
        str(xyz_dir / "methane_water_01.xyz"),
        elements_A=["C", "H", "H", "H", "H"],
        positions_A=[
            [0.0, 0.0, 0.0], [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63],
        ],
        elements_B=["O", "H", "H"],
        positions_B=[[4.0, 0.0, 0.0], [4.0, 0.76, -0.47], [4.0, -0.76, -0.47]],
        ref_energy=-0.65,
    )

    return str(xyz_dir)


class TestLoadNENCI:
    """Tests for loading NENCI-2021 XYZ files."""

    def test_load_reaction_count(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        assert len(ds.reactions) == 3

    def test_species_count(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        # 3 dimers × 3 species each = 9
        assert len(ds.species) == 9

    def test_species_naming(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        assert "NENCI-water_water_01" in ds.species
        assert "NENCI-water_water_01_monA" in ds.species
        assert "NENCI-water_water_01_monB" in ds.species

    def test_monomer_atom_counts(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        # Water dimer: monA = 3 atoms (O,H,H), monB = 3 atoms (O,H,H)
        assert ds.species["NENCI-water_water_01_monA"].natom == 3
        assert ds.species["NENCI-water_water_01_monB"].natom == 3
        # Methane-water: monA = 5 atoms (C,H,H,H,H), monB = 3 atoms (O,H,H)
        assert ds.species["NENCI-methane_water_01_monA"].natom == 5
        assert ds.species["NENCI-methane_water_01_monB"].natom == 3

    def test_dimer_atom_count(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        assert ds.species["NENCI-water_water_01"].natom == 6
        assert ds.species["NENCI-methane_water_01"].natom == 8

    def test_stoichiometry(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        # Find the water_water_01 reaction (files sorted alphabetically)
        rxn = [r for r in ds.reactions
               if "NENCI-water_water_01" in r.stoichiometry][0]
        # dimer = +1, monA = -1, monB = -1
        assert rxn.stoichiometry["NENCI-water_water_01"] == 1
        assert rxn.stoichiometry["NENCI-water_water_01_monA"] == -1
        assert rxn.stoichiometry["NENCI-water_water_01_monB"] == -1

    def test_reference_energy(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        refs = {list(r.stoichiometry.keys())[0]: r.reference_energy for r in ds.reactions}
        assert refs["NENCI-water_water_01"] == pytest.approx(-4.97)
        assert refs["NENCI-water_water_02"] == pytest.approx(-1.23)

    def test_subset_label(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        assert all(r.subset == "NENCI" for r in ds.reactions)

    def test_dataset_name(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        assert ds.name == "NENCI-2021"

    def test_element_symbols_valid(self, nenci_dir):
        ds = load_nenci(nenci_dir)
        for sp in ds.species.values():
            for el in sp.elements:
                assert el in _ELEMENT_TO_Z

    def test_d3_computable(self, nenci_dir):
        """D3 intermediates can be computed for all NENCI species."""
        ds = load_nenci(nenci_dir)
        for sp in ds.species.values():
            im = precompute_d3(sp)
            energy = d3_energy_bj(im, 1.0, 1.9889, 0.3981, 4.4211)
            assert np.isfinite(energy)

    def test_missing_dir_raises(self):
        with pytest.raises(FileNotFoundError):
            load_nenci("/nonexistent/path")

    def test_empty_dir_raises(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        with pytest.raises(ValueError, match="No valid NENCI"):
            load_nenci(str(empty_dir))

    def test_resolve_dataset_nenci(self, nenci_dir):
        ds = resolve_dataset(f"nenci:{nenci_dir}")
        assert ds.name == "NENCI-2021"
        assert len(ds.reactions) == 3


# ---------------------------------------------------------------------------
# NENCI subdirectory subsets
# ---------------------------------------------------------------------------

@pytest.fixture
def nenci_subdir(tmp_path):
    """Create NENCI files organized into subdirectories (S66, IonPi)."""
    root = tmp_path / "nenci_sub"
    root.mkdir()

    s66 = root / "S66"
    s66.mkdir()
    _write_nenci_xyz(
        str(s66 / "water_water_01.xyz"),
        elements_A=["O", "H", "H"],
        positions_A=[[0.0, 0.0, 0.0], [0.0, 0.76, -0.47], [0.0, -0.76, -0.47]],
        elements_B=["O", "H", "H"],
        positions_B=[[3.0, 0.0, 0.0], [3.0, 0.76, -0.47], [3.0, -0.76, -0.47]],
        ref_energy=-4.97,
    )
    _write_nenci_xyz(
        str(s66 / "water_water_02.xyz"),
        elements_A=["O", "H", "H"],
        positions_A=[[0.0, 0.0, 0.0], [0.0, 0.76, -0.47], [0.0, -0.76, -0.47]],
        elements_B=["O", "H", "H"],
        positions_B=[[5.0, 0.0, 0.0], [5.0, 0.76, -0.47], [5.0, -0.76, -0.47]],
        ref_energy=-1.23,
    )

    ionpi = root / "IonPi"
    ionpi.mkdir()
    _write_nenci_xyz(
        str(ionpi / "methane_water_01.xyz"),
        elements_A=["C", "H", "H", "H", "H"],
        positions_A=[
            [0.0, 0.0, 0.0], [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63],
        ],
        elements_B=["O", "H", "H"],
        positions_B=[[4.0, 0.0, 0.0], [4.0, 0.76, -0.47], [4.0, -0.76, -0.47]],
        ref_energy=-0.65,
    )

    return str(root)


class TestLoadNENCISubdirs:
    """Tests for NENCI loading with subdirectory-based subsets."""

    def test_loads_all_subdirs(self, nenci_subdir):
        ds = load_nenci(nenci_subdir)
        assert len(ds.reactions) == 3
        assert {r.subset for r in ds.reactions} == {"S66", "IonPi"}

    def test_subset_labels(self, nenci_subdir):
        ds = load_nenci(nenci_subdir)
        s66_rxns = [r for r in ds.reactions if r.subset == "S66"]
        ionpi_rxns = [r for r in ds.reactions if r.subset == "IonPi"]
        assert len(s66_rxns) == 2
        assert len(ionpi_rxns) == 1

    def test_filter_single_subset(self, nenci_subdir):
        ds = load_nenci(nenci_subdir, subsets=["S66"])
        assert len(ds.reactions) == 2
        assert all(r.subset == "S66" for r in ds.reactions)

    def test_filter_excludes_other(self, nenci_subdir):
        ds = load_nenci(nenci_subdir, subsets=["IonPi"])
        assert len(ds.reactions) == 1
        assert ds.reactions[0].subset == "IonPi"

    def test_dataset_name_single_subset(self, nenci_subdir):
        ds = load_nenci(nenci_subdir, subsets=["S66"])
        assert ds.name == "S66"

    def test_dataset_name_multiple_subsets(self, nenci_subdir):
        ds = load_nenci(nenci_subdir)
        assert ds.name == "NENCI_IonPi+S66"

    def test_resolve_dataset_with_subset_filter(self, nenci_subdir):
        ds = resolve_dataset(f"nenci:{nenci_subdir}:S66")
        assert len(ds.reactions) == 2
        assert all(r.subset == "S66" for r in ds.reactions)

    def test_resolve_dataset_multiple_subsets(self, nenci_subdir):
        ds = resolve_dataset(f"nenci:{nenci_subdir}:S66,IonPi")
        assert len(ds.reactions) == 3

    def test_empty_subset_raises(self, nenci_subdir):
        with pytest.raises(ValueError, match="No valid NENCI"):
            load_nenci(nenci_subdir, subsets=["Nonexistent"])


# ---------------------------------------------------------------------------
# MPCONF196 loading
# ---------------------------------------------------------------------------

def _write_mpconf_xyz(path, elements, positions):
    """Write one MPCONF196-format XYZ file (title line, no atom count)."""
    stem = os.path.splitext(os.path.basename(path))[0]
    lines = [f"{stem}.xyz"]
    for el, pos in zip(elements, positions):
        lines.append(f"{el}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")
    lines.append("")
    with open(path, "w") as f:
        f.write("\n".join(lines))


@pytest.fixture
def mpconf_dir(tmp_path):
    """Create a temp MPCONF196-style directory with 2 molecule groups."""
    d = tmp_path / "mpconf"
    d.mkdir()

    # Molecule A: 3 conformers (water at different orientations)
    for name, positions in [
        ("MolA_I", [[0.0, 0.0, 0.12], [0.0, 0.76, -0.47], [0.0, -0.76, -0.47]]),
        ("MolA_a", [[0.0, 0.0, 0.12], [0.0, 0.80, -0.47], [0.0, -0.80, -0.47]]),
        ("MolA_b", [[0.0, 0.0, 0.12], [0.0, 0.70, -0.47], [0.0, -0.70, -0.47]]),
    ]:
        _write_mpconf_xyz(str(d / f"{name}.xyz"), ["O", "H", "H"], positions)

    # Molecule B: 2 conformers (methane-like)
    for name, positions in [
        ("MolB01", [
            [0.0, 0.0, 0.0], [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63],
        ]),
        ("MolB02", [
            [0.0, 0.0, 0.0], [0.65, 0.65, 0.65],
            [-0.65, -0.65, 0.65], [-0.65, 0.65, -0.65], [0.65, -0.65, -0.65],
        ]),
    ]:
        _write_mpconf_xyz(str(d / f"{name}.xyz"), ["C", "H", "H", "H", "H"], positions)

    # Reference energies CSV
    import csv
    csv_path = str(d / "reference_energies.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["conformer", "molecule", "energy_kcal"])
        writer.writerow(["MolA_I", "MolA", "0.000"])
        writer.writerow(["MolA_a", "MolA", "1.500"])
        writer.writerow(["MolA_b", "MolA", "-0.800"])
        writer.writerow(["MolB01", "MolB", "0.300"])
        writer.writerow(["MolB02", "MolB", "-0.300"])

    return str(d)


class TestLoadMPCONF196:
    """Tests for loading MPCONF196 conformational energy benchmark."""

    def test_reaction_count(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir)
        # MolA: 3 conformers -> 3 reactions; MolB: 2 conformers -> 2 reactions
        assert len(ds.reactions) == 5

    def test_species_count(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir)
        assert len(ds.species) == 5

    def test_subset_labels(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir)
        subsets = {r.subset for r in ds.reactions}
        assert subsets == {"MolA", "MolB"}

    def test_mean_relative_stoichiometry(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir)
        # Each reaction: target conf gets (N-1)/N, others get -1/N
        # Coefficients should sum to 0.0
        for rxn in ds.reactions:
            assert sum(rxn.stoichiometry.values()) == pytest.approx(0.0)

    def test_reference_energies(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir)
        # Reference energies come directly from the CSV (mean-relative)
        mola_rxns = sorted(
            [r for r in ds.reactions if r.subset == "MolA"],
            key=lambda r: r.reference_energy,
        )
        # MolA_I=0.0, MolA_a=1.5, MolA_b=-0.8
        assert mola_rxns[0].reference_energy == pytest.approx(-0.800)
        assert mola_rxns[1].reference_energy == pytest.approx(0.000)
        assert mola_rxns[2].reference_energy == pytest.approx(1.500)

    def test_molecule_filter(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir, molecules=["MolA"])
        assert len(ds.reactions) == 3
        assert all(r.subset == "MolA" for r in ds.reactions)
        assert len(ds.species) == 3

    def test_dataset_name_full(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir)
        # molecules=None loads all -> name is "MPCONF196"
        assert ds.name == "MPCONF196"

    def test_dataset_name_single_molecule(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir, molecules=["MolB"])
        assert ds.name == "MolB"

    def test_d3_computable(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir)
        for sp in ds.species.values():
            im = precompute_d3(sp)
            energy = d3_energy_bj(im, 1.0, 1.9889, 0.3981, 4.4211)
            assert np.isfinite(energy)

    def test_missing_dir_raises(self):
        with pytest.raises(FileNotFoundError):
            load_mpconf196("/nonexistent/path")

    def test_missing_csv_raises(self, tmp_path):
        d = tmp_path / "no_csv"
        d.mkdir()
        _write_mpconf_xyz(str(d / "conf01.xyz"), ["H", "H"], [[0, 0, 0], [0, 0, 1]])
        with pytest.raises(FileNotFoundError, match="reference_energies.csv"):
            load_mpconf196(str(d))

    def test_empty_dir_raises(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        # Write CSV but no matching XYZ files
        import csv
        with open(str(d / "reference_energies.csv"), "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["conformer", "molecule", "energy_kcal"])
        with pytest.raises(ValueError, match="No valid MPCONF196"):
            load_mpconf196(str(d))

    def test_resolve_dataset_mpconf196(self, mpconf_dir):
        ds = resolve_dataset(f"mpconf196:{mpconf_dir}")
        assert len(ds.reactions) == 5

    def test_resolve_dataset_with_molecule_filter(self, mpconf_dir):
        ds = resolve_dataset(f"mpconf196:{mpconf_dir}:MolA")
        assert len(ds.reactions) == 3
        assert all(r.subset == "MolA" for r in ds.reactions)

    def test_element_symbols_valid(self, mpconf_dir):
        ds = load_mpconf196(mpconf_dir)
        for sp in ds.species.values():
            for el in sp.elements:
                assert el in _ELEMENT_TO_Z


# ---------------------------------------------------------------------------
# GMTKN55 loading (requires gmtkn package)
# ---------------------------------------------------------------------------

gmtkn = pytest.importorskip("gmtkn")

from dftd3.benchmark import load_gmtkn55  # noqa: E402


class TestLoadGMTKN55:
    """Tests for loading GMTKN55 subsets via the gmtkn package."""

    def test_load_single_subset(self):
        ds = load_gmtkn55(subsets=["S22"])
        assert len(ds.reactions) == 22
        assert all(r.subset == "S22" for r in ds.reactions)

    def test_load_multiple_subsets(self):
        ds = load_gmtkn55(subsets=["S22", "S66"])
        subsets = {r.subset for r in ds.reactions}
        assert subsets == {"S22", "S66"}
        assert len(ds.reactions) == 22 + 66

    def test_qualified_species_names(self):
        ds = load_gmtkn55(subsets=["S22"])
        for name in ds.species:
            assert name.startswith("S22-")

    def test_species_have_valid_atomnos(self):
        ds = load_gmtkn55(subsets=["S22"])
        for sp in ds.species.values():
            assert all(z > 0 for z in sp.atomnos)

    def test_element_symbols_normalized(self):
        """Ensure element symbols are title-cased and in _ELEMENT_TO_Z."""
        ds = load_gmtkn55(subsets=["G21EA"])
        for sp in ds.species.values():
            for el in sp.elements:
                assert el in _ELEMENT_TO_Z, f"{el} not in _ELEMENT_TO_Z"

    def test_stoichiometry_references_valid_species(self):
        ds = load_gmtkn55(subsets=["S22"])
        for rxn in ds.reactions:
            for sp_name in rxn.stoichiometry:
                assert sp_name in ds.species

    def test_weights_default_to_one(self):
        ds = load_gmtkn55(subsets=["S22"])
        assert all(r.weight == 1.0 for r in ds.reactions)

    def test_unknown_subset_raises(self):
        with pytest.raises(ValueError, match="Unknown"):
            load_gmtkn55(subsets=["NONEXISTENT"])

    def test_case_insensitive_subset_names(self):
        ds = load_gmtkn55(subsets=["s22"])
        assert len(ds.reactions) == 22

    def test_load_all_defaults(self):
        ds = load_gmtkn55()
        assert len(ds.reactions) > 1000
        subsets = {r.subset for r in ds.reactions}
        assert "S22" in subsets
        assert "S66" in subsets
        assert "BH76" in subsets


class TestResolveDataset:
    """Tests for resolve_dataset() dispatch."""

    def test_yaml_path(self, mini_yaml):
        ds = resolve_dataset(mini_yaml)
        assert len(ds.reactions) == 3

    def test_gmtkn55_single(self):
        ds = resolve_dataset("gmtkn55:S22")
        assert len(ds.reactions) == 22

    def test_gmtkn55_multiple(self):
        ds = resolve_dataset("gmtkn55:S22,S66")
        assert {r.subset for r in ds.reactions} == {"S22", "S66"}

    def test_gmtkn55_all(self):
        ds = resolve_dataset("gmtkn55")
        assert len(ds.reactions) > 1000

    def test_precompute_d3_on_gmtkn_species(self):
        """Ensure precompute_d3 works on gmtkn-loaded species."""
        ds = load_gmtkn55(subsets=["S22"])
        for sp in list(ds.species.values())[:5]:
            im = precompute_d3(sp)
            energy = d3_energy_bj(im, 1.0, 1.9889, 0.3981, 4.4211)
            assert np.isfinite(energy)
