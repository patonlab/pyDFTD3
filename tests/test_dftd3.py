"""Tests for pyDFTD3.

Reference values for CH3F2TS are taken directly from Grimme's original
Fortran DFTD3 V2.1 output (examples/CH3F2TS_reald3).
README examples provide additional reference values for formic acid dimer.
"""

import os

import pytest
from cclib.io import ccread

from dftd3.dftd3 import (
    AUTOANG,
    AUTOKCAL,
    MAX_ELEM,
    MAXC,
    PERIODIC_TABLE,
    CalcD3,
    _c6ab,
    _element_index,
    _getc6,
    _ncoord,
    _r,
    _r2r4,
    read_file,
)
from dftd3.pars import bj_parms, d4_parms, elements, resolve_functional, zero_parms

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "..", "examples")


def _example(filename):
    return os.path.join(EXAMPLES_DIR, filename)


def _parse_and_prepare(filename):
    """Parse a file and return (atomtype, xco, yco, zco, mxc, cn)."""
    data = read_file(_example(filename))
    atom_nums = data.atomnos.tolist()
    atomtype = [PERIODIC_TABLE[atno] for atno in atom_nums]
    cartesians = data.atomcoords[-1].tolist()
    natom = len(atomtype)

    xco = [at[0] for at in cartesians]
    yco = [at[1] for at in cartesians]
    zco = [at[2] for at in cartesians]

    mxc = [0] * MAX_ELEM
    for j in range(MAX_ELEM):
        for k in range(natom):
            if atomtype[k] == elements[j]:
                for ll in range(MAXC):
                    entry = _c6ab[j][j][ll][ll]
                    if isinstance(entry, list) and entry[0] > 0:
                        mxc[j] += 1
                break

    cn = _ncoord(natom, atomtype, xco, yco, zco)
    return data, atomtype, xco, yco, zco, mxc, cn


# ---------------------------------------------------------------------------
# Parameter / alias tests
# ---------------------------------------------------------------------------

class TestFunctionalAliases:
    def test_pbe1pbe_resolves_to_pbe0(self):
        assert resolve_functional("PBE1PBE") == "PBE0"

    def test_pbepbe_resolves_to_pbe(self):
        assert resolve_functional("PBEPBE") == "PBE"

    def test_bhandhlyp_resolves_to_bhlyp(self):
        assert resolve_functional("BHANDHLYP") == "BHLYP"

    def test_cam_b3lyp_resolves(self):
        assert resolve_functional("CAM-B3LYP") == "CAMB3LYP"

    def test_case_insensitive(self):
        assert resolve_functional("pbe1pbe") == "PBE0"
        assert resolve_functional("b3lyp") == "B3LYP"

    def test_unknown_functional_returns_uppercased(self):
        assert resolve_functional("myfunc") == "MYFUNC"


class TestParameterDicts:
    def test_zero_parms_is_dict(self):
        assert isinstance(zero_parms, dict)

    def test_bj_parms_is_dict(self):
        assert isinstance(bj_parms, dict)

    def test_b3lyp_zero_params(self):
        s6, rs6, s8 = zero_parms["B3LYP"]
        assert s6 == 1.0
        assert rs6 == 1.261
        assert s8 == 1.703

    def test_b3lyp_bj_params(self):
        s6, a1, s8, a2 = bj_parms["B3LYP"]
        assert s6 == 1.0
        assert a1 == 0.3981
        assert s8 == 1.9889
        assert a2 == 4.4211

    def test_alias_lookup_works_with_dict(self):
        canonical = resolve_functional("PBE1PBE")
        assert canonical in bj_parms

    def test_new_functionals_present(self):
        assert "HF" in zero_parms
        assert "SCAN" in zero_parms
        assert "R2SCAN" in bj_parms
        assert "RSCAN" in bj_parms

    def test_d4_parms_is_dict(self):
        assert isinstance(d4_parms, dict)

    def test_d4_b3lyp_params(self):
        """D4 B3LYP params from dftd4 parameters.toml."""
        s6, s8, a1, a2 = d4_parms["B3LYP"]
        assert s6 == 1.0
        assert s8 == pytest.approx(2.02929367, abs=1e-6)
        assert a1 == pytest.approx(0.40868035, abs=1e-6)
        assert a2 == pytest.approx(4.53807137, abs=1e-6)

    def test_d4_pbe_params(self):
        """D4 PBE params from dftd4 parameters.toml."""
        s6, s8, a1, a2 = d4_parms["PBE"]
        assert s6 == 1.0
        assert s8 == pytest.approx(0.95948085, abs=1e-6)
        assert a1 == pytest.approx(0.38574991, abs=1e-6)
        assert a2 == pytest.approx(4.80688534, abs=1e-6)

    def test_d4_scan_family_present(self):
        for func in ("SCAN", "RSCAN", "R2SCAN", "R2SCAN0", "R2SCAN50"):
            assert func in d4_parms, f"{func} missing from d4_parms"

    def test_d4_tuple_format(self):
        """All D4 entries should be (s6, s8, a1, a2) tuples."""
        for name, params in d4_parms.items():
            assert len(params) == 4, f"{name} has {len(params)} params, expected 4"


# ---------------------------------------------------------------------------
# Benchmark against Grimme's Fortran DFTD3 output (CH3F2TS_reald3)
#
# Reference (B3LYP, zero-damping):
#   CN:      C=4.182, H=0.998, F=0.606
#   R0(AA):  C=1.455, H=1.091, F=1.150
#   C6(AA):  C=18.2, H=3.1, F=7.8
#   C8(AA):  C=527.4, H=37.4, F=134.2
#   molC6:   223.33 au
#   E6:      -0.0777 kcal/mol
#   E8:      -1.3790 kcal/mol
#   E6(ABC): +0.000352 kcal/mol
#   Edisp:   -1.4563 kcal/mol  /  -0.00232076 au
# ---------------------------------------------------------------------------

class TestGrimmeReference:
    """Test against Grimme's original Fortran DFTD3 V2.1 output for CH3F2TS."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.data, self.atomtype, self.xco, self.yco, self.zco, self.mxc, self.cn = (
            _parse_and_prepare("CH3F2TS.log")
        )

    # -- Coordination numbers --

    def test_cn_carbon(self):
        """Atom 1 (C): CN = 4.182"""
        assert self.cn[0] == pytest.approx(4.182, abs=0.001)

    def test_cn_hydrogen(self):
        """Atoms 2-4 (H): CN = 0.998"""
        for i in (1, 2, 3):
            assert self.cn[i] == pytest.approx(0.998, abs=0.001)

    def test_cn_fluorine(self):
        """Atoms 5-6 (F): CN = 0.606"""
        for i in (4, 5):
            assert self.cn[i] == pytest.approx(0.606, abs=0.001)

    # -- R0(AA) diatomic cutoff radii --

    def test_r0_carbon(self):
        z = _element_index("C")
        r0 = 0.5 * AUTOANG * _r[z][z]
        assert r0 == pytest.approx(1.455, abs=0.001)

    def test_r0_hydrogen(self):
        z = _element_index("H")
        r0 = 0.5 * AUTOANG * _r[z][z]
        assert r0 == pytest.approx(1.091, abs=0.001)

    def test_r0_fluorine(self):
        z = _element_index("F")
        r0 = 0.5 * AUTOANG * _r[z][z]
        assert r0 == pytest.approx(1.150, abs=0.001)

    # -- C6(AA) coefficients --

    def test_c6_carbon(self):
        """C6(AA) for C at CN=4.182 should be 18.2 au."""
        C6 = _getc6(self.mxc, self.atomtype, self.cn, 0, 0)
        assert C6 == pytest.approx(18.2, abs=0.1)

    def test_c6_hydrogen(self):
        """C6(AA) for H at CN=0.998 should be 3.1 au."""
        C6 = _getc6(self.mxc, self.atomtype, self.cn, 1, 1)
        assert C6 == pytest.approx(3.1, abs=0.1)

    def test_c6_fluorine(self):
        """C6(AA) for F at CN=0.606 should be 7.8 au."""
        C6 = _getc6(self.mxc, self.atomtype, self.cn, 4, 4)
        assert C6 == pytest.approx(7.8, abs=0.1)

    # -- C8(AA) coefficients --

    def test_c8_carbon(self):
        """C8(AA) for C should be 527.4 au."""
        C6 = _getc6(self.mxc, self.atomtype, self.cn, 0, 0)
        z = _element_index("C")
        C8 = 3.0 * C6 * _r2r4[z] ** 2
        assert C8 == pytest.approx(527.4, abs=0.1)

    def test_c8_hydrogen(self):
        """C8(AA) for H should be 37.4 au."""
        C6 = _getc6(self.mxc, self.atomtype, self.cn, 1, 1)
        z = _element_index("H")
        C8 = 3.0 * C6 * _r2r4[z] ** 2
        assert C8 == pytest.approx(37.4, abs=0.1)

    def test_c8_fluorine(self):
        """C8(AA) for F should be 134.2 au."""
        C6 = _getc6(self.mxc, self.atomtype, self.cn, 4, 4)
        z = _element_index("F")
        C8 = 3.0 * C6 * _r2r4[z] ** 2
        assert C8 == pytest.approx(134.2, abs=0.1)

    # -- Molecular C6 --

    def test_molecular_c6(self):
        """Sum of all C6(ij) should be 223.33 au."""
        natom = len(self.atomtype)
        mol_c6 = 0.0
        for j in range(natom):
            for k in range(natom):
                mol_c6 += _getc6(self.mxc, self.atomtype, self.cn, j, k)
        assert mol_c6 == pytest.approx(223.33, abs=0.01)

    # -- D3 energy components (zero-damping, B3LYP) --

    def test_e6_kcal(self):
        """E6 should be -0.0777 kcal/mol."""
        result = CalcD3(self.data, "B3LYP", damp="zero")
        assert result.attractive_r6_vdw == pytest.approx(-0.0777, abs=0.0001)

    def test_e8_kcal(self):
        """E8 should be -1.3790 kcal/mol."""
        result = CalcD3(self.data, "B3LYP", damp="zero")
        assert result.attractive_r8_vdw == pytest.approx(-1.3790, abs=0.0001)

    def test_edisp_total_kcal(self):
        """Edisp (E6+E8) should be -1.4563 kcal/mol (within rounding)."""
        result = CalcD3(self.data, "B3LYP", damp="zero")
        total = result.attractive_r6_vdw + result.attractive_r8_vdw
        assert total == pytest.approx(-1.4563, abs=0.001)

    def test_edisp_total_au(self):
        """Edisp should be -0.00232076 au (within rounding)."""
        result = CalcD3(self.data, "B3LYP", damp="zero")
        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        assert total_au == pytest.approx(-0.00232076, abs=1e-6)

    # -- 3-body ATM term --

    def test_e6abc(self):
        """E6(ABC) should be +0.000352 kcal/mol."""
        result = CalcD3(self.data, "B3LYP", damp="zero", abc=True)
        assert result.repulsive_abc == pytest.approx(0.000352, abs=0.000001)


# ---------------------------------------------------------------------------
# D3(BJ) - formic acid dimer (README reference values)
# ---------------------------------------------------------------------------

class TestD3BJ:
    """Test D3(BJ) corrections against README examples."""

    def test_formic_acid_dimer_log_bj(self):
        """README example 2: formic_acid_dimer.log with BJ damping."""
        data = ccread(_example("formic_acid_dimer.log"))
        result = CalcD3(data, "B3LYP", damp="bj")

        c6_au = result.attractive_r6_vdw / AUTOKCAL
        c8_au = result.attractive_r8_vdw / AUTOKCAL
        total_au = c6_au + c8_au

        assert c6_au == pytest.approx(-0.00455241, abs=1e-7)
        assert c8_au == pytest.approx(-0.00457708, abs=1e-7)
        assert total_au == pytest.approx(-0.00912948, abs=1e-6)

    def test_formic_acid_dimer_xyz_bj(self):
        """README example 5: formic_acid_dimer.xyz with BJ damping + B3LYP."""
        data = ccread(_example("formic_acid_dimer.xyz"))
        result = CalcD3(data, "B3LYP", damp="bj")

        c6_au = result.attractive_r6_vdw / AUTOKCAL
        c8_au = result.attractive_r8_vdw / AUTOKCAL
        total_au = c6_au + c8_au

        assert c6_au == pytest.approx(-0.00455241, abs=1e-7)
        assert c8_au == pytest.approx(-0.00457708, abs=1e-7)
        assert total_au == pytest.approx(-0.00912949, abs=1e-6)

    def test_formic_acid_dimer_sdf_bj(self):
        """SDF file should give same result as XYZ."""
        data = ccread(_example("formic_acid_dimer.sdf"))
        result = CalcD3(data, "B3LYP", damp="bj")

        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        assert total_au == pytest.approx(-0.00912949, abs=1e-5)

    def test_formic_acid_dimer_pdb_bj(self):
        """PDB file should give same result as XYZ."""
        data = ccread(_example("formic_acid_dimer.pdb"))
        result = CalcD3(data, "B3LYP", damp="bj")

        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        assert total_au == pytest.approx(-0.00912949, abs=1e-5)


# ---------------------------------------------------------------------------
# D3 zero-damping - formic acid dimer
# ---------------------------------------------------------------------------

class TestD3Zero:
    def test_formic_acid_dimer_log_zero(self):
        data = ccread(_example("formic_acid_dimer.log"))
        result = CalcD3(data, "B3LYP", damp="zero")

        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        assert total_au < 0


# ---------------------------------------------------------------------------
# Functional alias resolution end-to-end
# ---------------------------------------------------------------------------

class TestAliasIntegration:
    def test_pbe1pbe_gives_same_as_pbe0(self):
        data = ccread(_example("formic_acid_dimer.xyz"))
        result_pbe0 = CalcD3(data, "PBE0", damp="bj")
        result_alias = CalcD3(data, "PBE1PBE", damp="bj")

        assert result_pbe0.attractive_r6_vdw == pytest.approx(
            result_alias.attractive_r6_vdw, abs=1e-10)
        assert result_pbe0.attractive_r8_vdw == pytest.approx(
            result_alias.attractive_r8_vdw, abs=1e-10)

    def test_pbepbe_gives_same_as_pbe(self):
        data = ccread(_example("formic_acid_dimer.xyz"))
        result_pbe = CalcD3(data, "PBE", damp="bj")
        result_alias = CalcD3(data, "PBEPBE", damp="bj")

        assert result_pbe.attractive_r6_vdw == pytest.approx(
            result_alias.attractive_r6_vdw, abs=1e-10)


# ---------------------------------------------------------------------------
# 3-body term
# ---------------------------------------------------------------------------

class TestThreeBody:
    def test_three_body_term_computed(self):
        data = ccread(_example("formic_acid_dimer.log"))
        result = CalcD3(data, "B3LYP", damp="bj", abc=True)
        assert result.repulsive_abc != 0.0


# ---------------------------------------------------------------------------
# Ibuprofen D3 — cross-validated against ORCA, Q-Chem, and Gaussian
#
# All ibuprofen examples use the same C13H18O2 geometry (33 atoms).
# ORCA reference values include E6/E8 breakdown and molecular C6.
# Q-Chem reference values give total D3 energy (no E6/E8 breakdown).
# Gaussian reference values are obtained by SCF(with D3) - SCF(without D3).
# ---------------------------------------------------------------------------


class TestOrcaFallbackParser:
    """Verify the ORCA fallback parser produces valid data."""

    def test_parses_b3lyp_orca(self):
        data = read_file(_example("ibuprofen_b3lyp_orca.out"))
        assert data is not None
        assert len(data.atomnos) == 33
        assert data.atomcoords.shape == (1, 33, 3)

    def test_parses_pbe_orca(self):
        data = read_file(_example("ibuprofen_pbe_orca.out"))
        assert data is not None
        assert len(data.atomnos) == 33

    def test_detects_functional_b3lyp(self):
        data = read_file(_example("ibuprofen_b3lyp_orca.out"))
        assert data.metadata.get("functional") == "B3LYP"

    def test_detects_functional_pbe(self):
        data = read_file(_example("ibuprofen_pbe_orca.out"))
        assert data.metadata.get("functional") == "PBE"


class TestIbuprofenOrcaB3lypBJ:
    """D3(BJ)/B3LYP for ibuprofen parsed from ORCA output.

    ORCA reference (ibuprofen_b3lyp_orca.out):
      Edisp = -0.063582584493 au, E6 = -19.144204 kcal, E8 = -20.754474 kcal
      molecular C6(AA) = 9785.118961 au
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        self.data = read_file(_example("ibuprofen_b3lyp_orca.out"))
        self.result = CalcD3(self.data, "B3LYP", damp="bj")

    def test_total_dispersion_au(self):
        total_au = (self.result.attractive_r6_vdw + self.result.attractive_r8_vdw) / AUTOKCAL
        assert total_au == pytest.approx(-0.063582584, abs=1e-5)

    def test_e6_kcal(self):
        assert self.result.attractive_r6_vdw == pytest.approx(-19.144204, abs=0.01)

    def test_e8_kcal(self):
        assert self.result.attractive_r8_vdw == pytest.approx(-20.754474, abs=0.01)

    def test_molecular_c6(self):
        """molecular C6(AA) = 9785.12 au (ORCA reference)."""
        _, atomtype, xco, yco, zco, mxc, cn = _parse_and_prepare(
            "ibuprofen_b3lyp_orca.out"
        )
        natom = len(atomtype)
        mol_c6 = 0.0
        for j in range(natom):
            for k in range(natom):
                mol_c6 += _getc6(mxc, atomtype, cn, j, k)
        assert mol_c6 == pytest.approx(9785.12, abs=0.01)


class TestIbuprofenOrcaPbeZero:
    """D3(zero)/PBE for ibuprofen parsed from ORCA output.

    ORCA reference (ibuprofen_pbe_orca.out, run with D3ZERO ABC):
      E6 = -5.330527 kcal, E8 = -6.828818 kcal
      E6(ABC) = 0.111980 kcal
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        self.data = read_file(_example("ibuprofen_pbe_orca.out"))

    def test_e6_kcal(self):
        result = CalcD3(self.data, "PBE", damp="zero")
        assert result.attractive_r6_vdw == pytest.approx(-5.330527, abs=0.01)

    def test_e8_kcal(self):
        result = CalcD3(self.data, "PBE", damp="zero")
        assert result.attractive_r8_vdw == pytest.approx(-6.828818, abs=0.01)

    def test_abc_term(self):
        result = CalcD3(self.data, "PBE", damp="zero", abc=True)
        assert result.repulsive_abc == pytest.approx(0.111980, abs=0.001)


class TestIbuprofenQchemB3lypBJ:
    """D3(BJ)/B3LYP for ibuprofen parsed from Q-Chem output.

    Q-Chem reference (ibuprofen_b3lyp_qchem.out):
      D3(BJ) energy = -0.0635825846 hartrees
    """

    def test_total_dispersion_au(self):
        data = ccread(_example("ibuprofen_b3lyp_qchem.out"))
        result = CalcD3(data, "B3LYP", damp="bj")
        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        assert total_au == pytest.approx(-0.0635825846, abs=1e-5)


class TestIbuprofenQchemM062x:
    """D3(zero)/M062X for ibuprofen parsed from Q-Chem output.

    Q-Chem (ibuprofen_pbe_qchem.out — actually M062X despite filename):
      D3(zero) energy = -0.0016601747 hartrees
    """

    def test_total_dispersion_au(self):
        data = ccread(_example("ibuprofen_pbe_qchem.out"))
        result = CalcD3(data, "M062X", damp="zero")
        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        assert total_au == pytest.approx(-0.0016601747, abs=1e-5)


class TestIbuprofenGaussian:
    """D3 for ibuprofen validated against Gaussian SCF energy differences.

    The Gaussian log files contain 3 linked jobs each:
      Job 1: TPSS/def2TZVP  (with/without emp=GD3BJ)
      Job 2: M06/def2TZVP   (with/without emp=GD3)
      Job 3: PBE/def2TZVP   (with/without emp=GD3BJ)

    Reference D3 energies (SCF with D3 minus SCF without D3):
      TPSS-D3BJ:  -657.090997955 - (-657.042971105) = -0.048026850 au
      M06-D3zero: -656.478723887 - (-656.473417627) = -0.005306260 au
      PBE-D3BJ:   -656.139645483 - (-656.101381251) = -0.038264232 au
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        self.data = ccread(_example("ibuprofen_d3.log"))

    def test_tpss_d3bj(self):
        result = CalcD3(self.data, "TPSS", damp="bj")
        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        gaussian_ref = -657.090997955 - (-657.042971105)
        assert total_au == pytest.approx(gaussian_ref, abs=1e-5)

    def test_m06_d3zero(self):
        result = CalcD3(self.data, "M06", damp="zero")
        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        gaussian_ref = -656.478723887 - (-656.473417627)
        assert total_au == pytest.approx(gaussian_ref, abs=1e-5)

    def test_pbe_d3bj(self):
        result = CalcD3(self.data, "PBE", damp="bj")
        total_au = (result.attractive_r6_vdw + result.attractive_r8_vdw) / AUTOKCAL
        gaussian_ref = -656.139645483 - (-656.101381251)
        assert total_au == pytest.approx(gaussian_ref, abs=1e-5)


class TestIbuprofenCrossProgram:
    """Cross-program consistency: same geometry and functional gives same D3."""

    def test_orca_qchem_b3lyp_bj_consistency(self):
        """B3LYP-D3BJ should match across ORCA and Q-Chem parsed geometries."""
        data_orca = read_file(_example("ibuprofen_b3lyp_orca.out"))
        data_qchem = ccread(_example("ibuprofen_b3lyp_qchem.out"))
        result_orca = CalcD3(data_orca, "B3LYP", damp="bj")
        result_qchem = CalcD3(data_qchem, "B3LYP", damp="bj")
        orca_total = (result_orca.attractive_r6_vdw + result_orca.attractive_r8_vdw) / AUTOKCAL
        qchem_total = (result_qchem.attractive_r6_vdw + result_qchem.attractive_r8_vdw) / AUTOKCAL
        assert orca_total == pytest.approx(qchem_total, abs=1e-6)

    def test_orca_gaussian_b3lyp_bj_consistency(self):
        """B3LYP-D3BJ should match across ORCA and Gaussian parsed geometries."""
        data_orca = read_file(_example("ibuprofen_b3lyp_orca.out"))
        data_gauss = ccread(_example("ibuprofen_d3.log"))
        result_orca = CalcD3(data_orca, "B3LYP", damp="bj")
        result_gauss = CalcD3(data_gauss, "B3LYP", damp="bj")
        orca_total = (result_orca.attractive_r6_vdw + result_orca.attractive_r8_vdw) / AUTOKCAL
        gauss_total = (result_gauss.attractive_r6_vdw + result_gauss.attractive_r8_vdw) / AUTOKCAL
        assert orca_total == pytest.approx(gauss_total, abs=1e-6)
