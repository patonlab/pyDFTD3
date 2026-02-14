![pyDFTD3 logo](pyDFTD3_banner.png)

[![DOI](https://zenodo.org/badge/54939983.svg)](https://zenodo.org/badge/latestdoi/54939983)

pyDFTD3 computes Grimme's D3 dispersion energy corrections for molecular geometries. It implements both zero-damping and Becke-Johnson (BJ) damping schemes with optional 3-body Axilrod-Teller-Muto terms.

Supported input formats: Gaussian/ORCA output files (`.log`, `.out`), XYZ, PDB, and SDF. For Gaussian and ORCA output files the density functional is detected automatically; for other formats use `--func` to specify the functional.

## Installation

```bash
pip install -e .              # editable install (creates `pydftd3` CLI)
pip install -e ".[test]"      # with test dependencies (pytest)
```

## Usage

```bash
# As a CLI tool (after install)
pydftd3 <file> --func <functional> --damp <zero|bj>

# As a module
python -m dftd3 <file> --func <functional> --damp <zero|bj>
```

### Options

| Flag | Description |
|------|-------------|
| `--damp {zero,bj}` | Damping function (default: zero) |
| `--func NAME` | Density functional for default parameters |
| `--abc` | Include repulsive 3-body (ATM) term |
| `--pw` | Print pairwise dispersion breakdown |
| `--im "1-5:6-10"` | Compute only intermolecular dispersion |
| `--cutoff N` | Distance cutoff in Angstrom (default: no cutoff) |
| `--kcal` | Print energies in kcal/mol |
| `--cite` | Print citation information |
| `-v` | Verbose output |

For zero-damping, manual parameters are `--s6`, `--rs6`, `--s8`. For BJ-damping: `--s6`, `--s8`, `--a1`, `--a2`.

## Examples

Structure files are available in the `examples/` directory.

1. **D3 zero-damping** from a Gaussian output file. The functional (B3LYP) is detected automatically.

```
$ python -m dftd3 examples/formic_acid_dimer.log

   Species                                           D3(R6)     D3(R8)        ABC  Etot (Hartree)
   ----------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.log                 -0.000913  -0.004347                  -0.005259
   ----------------------------------------------------------------------------------------------
```

2. **D3(BJ) damping** from a Gaussian output file.

```
$ python -m dftd3 examples/formic_acid_dimer.log --damp bj

   Species                                           D3(R6)     D3(R8)        ABC  Etot (Hartree)
   ----------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.log                 -0.004552  -0.004577                  -0.009129
   ----------------------------------------------------------------------------------------------
```

3. **D3(BJ) with 3-body term** enabled via `--abc`.

```
$ python -m dftd3 examples/formic_acid_dimer.log --damp bj --abc

   Species                                           D3(R6)     D3(R8)        ABC  Etot (Hartree)
   ----------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.log                 -0.004552  -0.004577  -0.000000       -0.009129
   ----------------------------------------------------------------------------------------------
```

4. **Output in kcal/mol** using `--kcal`.

```
$ python -m dftd3 examples/formic_acid_dimer.log --damp bj --kcal

   Species                                           D3(R6)     D3(R8)        ABC Etot (kcal/mol)
   ----------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.log                     -2.86      -2.87                      -5.73
   ----------------------------------------------------------------------------------------------
```

5. **XYZ input** with explicit functional. For file formats without embedded DFT metadata (XYZ, PDB, SDF), the functional must be specified with `--func`.

```
$ python -m dftd3 examples/formic_acid_dimer.xyz --func b3lyp --damp bj

   Species                                           D3(R6)     D3(R8)        ABC  Etot (Hartree)
   ----------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.xyz                 -0.004552  -0.004577                  -0.009129
   ----------------------------------------------------------------------------------------------
```

## Cutoff Radius Benchmark

By default, all pairwise interactions are included (no cutoff). A distance cutoff can be applied with `--cutoff` for large systems to reduce computation time. The tables below show D3(BJ)/B3LYP dispersion energy convergence with respect to cutoff radius.

**Maitotoxin** (285 atoms):

| Cutoff | Edisp (kcal/mol) | Error (kcal/mol) | Time (s) |
|-------:|-----------------:|-----------------:|---------:|
| 6 | -522.718 | +5.337 | 0.9 |
| 9 | -527.363 | +0.690 | 1.4 |
| 12 | -527.863 | +0.191 | 2.0 |
| 15 | -527.978 | +0.076 | 2.7 |
| 20 | -528.036 | +0.018 | 4.7 |
| 30 | -528.053 | +0.001 | 8.5 |
| None | -528.054 | reference | 11.5 |

**Human Insulin A chain, PDB: 3I40** (446 atoms):

| Cutoff | Edisp (kcal/mol) | Error (kcal/mol) | Time (s) |
|-------:|-----------------:|-----------------:|---------:|
| 6 | -975.557 | +46.468 | 2.9 |
| 9 | -1014.332 | +7.693 | 7.3 |
| 12 | -1020.264 | +1.757 | 14.4 |
| 15 | -1021.600 | +0.421 | 23.7 |
| 20 | -1021.992 | +0.029 | 38.1 |
| 30 | -1022.021 | +0.000 | 45.4 |
| None | -1022.021 | reference | 44.2 |

A cutoff of 15-20 Angstrom recovers >99.99% of the full dispersion energy. The benchmarking script is available at `examples/benchmark_cutoff.py`.

## Testing

```bash
pytest                  # run all tests
pytest tests/ -v        # verbose
```

Tests verify D3 energies against Grimme's original Fortran DFTD3 V2.1 reference output and against the examples above.

## Citations

If you use pyDFTD3, please cite:

1. Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H. *J. Chem. Phys.* **2010**, *132*, 154104.
2. Grimme, S.; Ehrlich, S.; Goerigk, L. *J. Comput. Chem.* **2011**, *32*, 1456-1465.

---
License: [MIT](https://opensource.org/licenses/MIT)
