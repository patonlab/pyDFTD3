"""Tests for D3 parameter optimization (benchmark + optimize modules)."""

import os

import pytest

from dftd3.benchmark import (
    Dataset,
    Reaction,
    Species,
    generate_orca_inputs,
    load_dataset,
    read_energies_csv,
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
