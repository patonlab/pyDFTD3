![pyDFTD3 logo](pyDFTD3_banner.png)

[![DOI](https://zenodo.org/badge/54939983.svg)](https://zenodo.org/badge/latestdoi/54939983)

pyDFTD3 computes Grimme's D3 dispersion energy corrections for molecular geometries. It implements both zero-damping and Becke-Johnson (BJ) damping schemes with optional 3-body Axilrod-Teller-Muto terms.

Supported input formats: Gaussian/ORCA/Q-Chem output files (`.log`, `.out`), XYZ, PDB, and SDF. For computational chemistry output files the density functional is detected automatically; for other formats use `--func` to specify the functional.

## Installation

```bash
pip install -e .              # editable install (creates `pydftd3` CLI)
pip install -e ".[test]"      # with test dependencies (pytest)
pip install -e ".[optimize]"  # with optimization dependencies (scipy, pyyaml)
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
| `--im auto` | Auto-detect fragments from covalent connectivity |
| `--cutoff N` | Distance cutoff in Angstrom (default: no cutoff) |
| `--kcal` | Print energies in kcal/mol |
| `--cite` | Print citation information |
| `-v` | Verbose output |

For zero-damping, manual parameters are `--s6`, `--rs6`, `--s8`. For BJ-damping: `--s6`, `--s8`, `--a1`, `--a2`.

## Examples

Structure files are available in the `examples/` directory.

1. **D3 zero-damping** from a Gaussian output file. The functional (B3LYP) is detected automatically.

```text
$ python -m dftd3 examples/formic_acid_dimer.log

   Species                                        D3(R6)     D3(R8)        ABC  Etot (Hartree)
   -------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.log              -0.000913  -0.004347                  -0.005259
   -------------------------------------------------------------------------------------------
```

2. **D3(BJ) damping** from a Gaussian output file.

```text
$ python -m dftd3 examples/formic_acid_dimer.log --damp bj

   Species                                        D3(R6)     D3(R8)        ABC  Etot (Hartree)
   -------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.log              -0.004552  -0.004577                  -0.009129
   -------------------------------------------------------------------------------------------
```

3. **D3(BJ) with 3-body term** enabled via `--abc`.

```text
$ python -m dftd3 examples/formic_acid_dimer.log --damp bj --abc

   Species                                        D3(R6)     D3(R8)        ABC  Etot (Hartree)
   -------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.log              -0.004552  -0.004577  -0.000000       -0.009129
   -------------------------------------------------------------------------------------------
```

4. **Output in kcal/mol** using `--kcal`.

```text
$ python -m dftd3 examples/formic_acid_dimer.log --damp bj --kcal

   Species                                        D3(R6)     D3(R8)        ABC Etot (kcal/mol)
   -------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.log                  -2.86      -2.87                      -5.73
   -------------------------------------------------------------------------------------------
```

5. **XYZ input** with explicit functional. For file formats without embedded DFT metadata (XYZ, PDB, SDF), the functional must be specified with `--func`.

```text
$ python -m dftd3 examples/formic_acid_dimer.xyz --func b3lyp --damp bj

   Species                                        D3(R6)     D3(R8)        ABC  Etot (Hartree)
   -------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.xyz              -0.004552  -0.004577                  -0.009129
   -------------------------------------------------------------------------------------------
```

6. **Automatic intermolecular mode** using `--im auto`. Fragments are detected automatically from covalent radii-based connectivity. Use `-v` to see fragment details.

```text
$ python -m dftd3 examples/formic_acid_dimer.xyz --func b3lyp --damp bj --im auto

   Caution: fragments detected automatically from covalent connectivity.
   Only intermolecular dispersion interactions are included.

   Species                                        D3(R6)     D3(R8)        ABC  Etot (Hartree)
   -------------------------------------------------------------------------------------------
   examples/formic_acid_dimer.xyz              -0.001521  -0.001449                  -0.002970
   -------------------------------------------------------------------------------------------
```

## Cutoff Radius Benchmark

By default, all pairwise interactions are included (no cutoff). A distance cutoff can be applied with `--cutoff` for large systems to reduce computation time. The tables below show B3LYP-D3(BJ) dispersion energy convergence with respect to cutoff radius.

**Maitotoxin** (285 atoms):

| Cutoff | Edisp (kcal/mol) | Error (kcal/mol) | Time (s) |
|-------:|-----------------:|-----------------:|---------:|
| 10 | -527.626 | +0.419 | 0.03 |
| 15 | -527.969 | +0.076 | 0.02 |
| 25 | -528.041 | +0.004 | 0.02 |
| None | -528.045 | reference | 0.02 |

**Maitotoxin with 3-body ATM term** (`--abc`). The ABC term is O(N³) and dominates the computation time, though its energy contribution is negligible. A cutoff is recommended when using `--abc` on large systems.

| Cutoff | Edisp (kcal/mol) | Time without ABC (s) | Time with ABC (s) |
|-------:|-----------------:|---------------------:|-------------------:|
| 10 | -527.626 | 0.02 | 1.7 |
| 15 | -527.969 | 0.02 | 3.0 |
| 25 | -528.041 | 0.02 | 9.1 |
| None | -528.045 | 0.02 | 17.9 |

**Human Insulin A chain, PDB: 3I40** (446 atoms):

| Cutoff | Edisp (kcal/mol) | Error (kcal/mol) | Time (s) |
|-------:|-----------------:|-----------------:|---------:|
| 10 | -1017.432 | +4.595 | 0.04 |
| 15 | -1021.606 | +0.421 | 0.04 |
| 25 | -1022.026 | +0.001 | 0.04 |
| None | -1022.027 | reference | 0.04 |

**Titin Z1Z2–Telethonin, PDB: 1YA5** (3930 atoms):

| Cutoff | Edisp (kcal/mol) | Error (kcal/mol) | Time (s) |
|-------:|-----------------:|-----------------:|---------:|
| 10 | -8918.611 | +65.417 | 4.4 |
| 15 | -8972.563 | +11.465 | 4.7 |
| 25 | -8983.077 | +0.951 | 4.9 |
| None | -8984.027 | reference | 4.3 |

A cutoff of 25 Angstrom recovers >99.99% of the full dispersion energy. The benchmarking script is available at `examples/benchmark_cutoff.py`.

## Parameter Optimization

pyDFTD3 can fit D3 damping parameters for functional/basis set combinations that lack published values, or to refit parameters against a custom benchmark. The workflow has two steps:

### 1. Generate DFT inputs (`prep`)

Create dispersion-free DFT input files for all species in a benchmark dataset:

```bash
pydftd3 optimize prep <dataset> -o <output_dir> --func <functional>
```

| Flag | Description |
|------|-------------|
| `<dataset>` | Dataset specifier (see below) |
| `-o, --output-dir` | Directory for generated input files |
| `--func` | DFT functional name |
| `--basis` | Basis set (default: def2-TZVP) |
| `--nprocs` | ORCA parallel processes (default: 8) |
| `--maxcore` | Memory per core in MB (default: 4000) |
| `--extra` | Additional ORCA keywords |

Run the generated ORCA calculations externally, then proceed to fitting.

### 2. Fit parameters (`fit`)

Optimize D3 parameters against benchmark reference energies using `scipy.optimize.differential_evolution`:

```bash
pydftd3 optimize fit <dataset> --energies <dir_or_csv> --damp bj
```

| Flag | Description |
|------|-------------|
| `<dataset>` | Dataset specifier (see below) |
| `--energies` | CSV file or directory of ORCA output files |
| `--damp {zero,bj}` | Damping scheme (default: bj) |
| `--test-frac` | Fraction for random test set (default: 0.2, 0=no split) |
| `--test-subsets` | Hold out entire subsets for validation (e.g. `S66,WATER27`) |
| `--folds` | k-fold cross-validation (0=single split, default: 0) |
| `--val-dataset` | External validation dataset (repeatable) |
| `--val-energies` | DFT energies for validation dataset (repeatable, paired with `--val-dataset`) |
| `--skip-missing` | Skip reactions with missing DFT energies |
| `--save-csv` | Save per-reaction errors to CSV |
| `-v` | Print per-subset breakdown |

### Dataset specifiers

| Format | Example | Description |
|--------|---------|-------------|
| YAML path | `AllElements_100.yaml` | DietGMTKN55 YAML file |
| `gmtkn55` | `gmtkn55` | All 55 GMTKN55 subsets (requires [gmtkn](https://github.com/obackhouse/gmtkn)) |
| `gmtkn55:subsets` | `gmtkn55:S22,S66,BH76` | Specific GMTKN55 subsets |
| `nenci:path` | `nenci:nenci2021/` | NENCI-2021 dimer XYZ files |

### Examples

```bash
# Prepare ORCA inputs for S22 and S66 subsets
pydftd3 optimize prep gmtkn55:S22,S66 -o s22_s66/ --func PBE

# Fit BJ parameters with 5-fold cross-validation
pydftd3 optimize fit gmtkn55:S22 --energies s22_dft/ --damp bj --folds 5

# Fit with an external validation set
pydftd3 optimize fit data.yaml --energies dft.csv \
    --val-dataset nenci:nenci2021/ --val-energies nenci_dft.csv

# Hold out entire subsets for testing
pydftd3 optimize fit gmtkn55:S22,S66,WATER27 --energies dft.csv \
    --test-subsets S66,WATER27 --damp bj -v
```

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
