"""
benchmark.py - Load benchmark datasets and generate ORCA inputs.

Provides tools for loading benchmark reaction datasets from GMTKN55 YAML
files, GMTKN55 (via the gmtkn package), NENCI-2021 XYZ files, or
MPCONF196 conformational energy data. Also
generates ORCA input files for dispersion-free DFT calculations and parses
ORCA outputs to extract single-point energies.
"""

import csv
import logging
import os
import re
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)

# Element symbol -> atomic number mapping
_ELEMENT_TO_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54,
    "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61,
    "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68,
    "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75,
    "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
    "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89,
    "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94,
}


@dataclass
class Species:
    """A single molecular species in the benchmark dataset."""

    name: str
    charge: int
    uhf: int
    elements: list
    positions: list

    @property
    def atomnos(self):
        """Atomic numbers as numpy array (compatible with CalcD3)."""
        return np.array([_ELEMENT_TO_Z[el] for el in self.elements], dtype=int)

    @property
    def atomcoords(self):
        """Coordinates as numpy array shaped (1, natom, 3) for CalcD3."""
        return np.array([self.positions])

    @property
    def natom(self):
        return len(self.elements)

    @property
    def multiplicity(self):
        """Spin multiplicity from UHF count."""
        return self.uhf + 1


@dataclass
class Reaction:
    """A benchmark reaction with reference energy and stoichiometry."""

    subset: str
    index: int
    reference_energy: float
    weight: float
    stoichiometry: dict


@dataclass
class Dataset:
    """A full benchmark dataset."""

    name: str
    reactions: list
    species: dict


def load_dataset(yaml_path):
    """Load a GMTKN55 YAML file into a Dataset object.

    Parameters
    ----------
    yaml_path : str
        Path to a GMTKN55 AllElements YAML file
        (e.g. AllElements_100.yaml).

    Returns
    -------
    Dataset
        Parsed dataset with reactions and unique species.
    """
    import yaml

    with open(yaml_path) as f:
        raw = yaml.safe_load(f)

    reactions = []
    species = {}

    for subset_name, subset_data in raw.items():
        if not isinstance(subset_data, dict):
            continue

        for rxn_index, rxn_data in subset_data.items():
            if not isinstance(rxn_data, dict):
                continue

            energy = rxn_data.get("Energy")
            weight = rxn_data.get("Weight")
            if energy is None or weight is None:
                logger.warning("Skipping %s-%s: missing Energy or Weight", subset_name, rxn_index)
                continue

            species_data = rxn_data.get("Species", {})
            stoichiometry = {}

            for sp_name, sp_info in species_data.items():
                count = sp_info.get("Count", 1)
                # Use subset-qualified name to avoid collisions
                qualified_name = f"{subset_name}-{sp_name}"
                stoichiometry[qualified_name] = count

                if qualified_name not in species:
                    elements = sp_info.get("Elements", [])
                    positions = sp_info.get("Positions", [])
                    charge = sp_info.get("Charge", 0)
                    uhf = sp_info.get("UHF", 0)

                    species[qualified_name] = Species(
                        name=qualified_name,
                        charge=charge,
                        uhf=uhf,
                        elements=elements,
                        positions=positions,
                    )

            reactions.append(Reaction(
                subset=subset_name,
                index=int(rxn_index),
                reference_energy=float(energy),
                weight=float(weight),
                stoichiometry=stoichiometry,
            ))

    dataset_name = os.path.splitext(os.path.basename(yaml_path))[0]
    logger.info("\n   Loaded %d reactions, %d unique species from %s",
                len(reactions), len(species), yaml_path)

    return Dataset(name=dataset_name, reactions=reactions, species=species)


# The 55 canonical GMTKN55 subset names (Goerigk et al., PCCP 2017)
GMTKN55_SUBSETS = frozenset([
    "W4_11", "G21EA", "G21IP", "DIPCS10", "PA26", "SIE4x4", "ALKBDE10",
    "YBDE18", "AL2X6", "HEAVYSB11", "NBPRC", "ALK8", "RC21", "G2RC",
    "FH51", "TAUT15", "DC13", "MB16_43", "DARC", "RSE43", "BSR36",
    "CDIE20", "ISO34", "ISOL24", "C60ISO", "PArel", "BH76", "BHPERI",
    "BHDIV10", "INV24", "BHROT27", "PX13", "WCPT18", "RG18", "ADIM6",
    "S22", "S66", "HEAVY28", "WATER27", "CARBHB12", "PNICO23", "HAL59",
    "AHB21", "CHB6", "IL16", "IDISP", "ICONF", "ACONF", "Amino20x4",
    "PCONF21", "MCONF", "SCONF", "UPU23", "BUT14DIOL",
])


def load_gmtkn55(subsets=None):
    """Load GMTKN55 benchmark data from the gmtkn package.

    Parameters
    ----------
    subsets : list of str, optional
        Subset names to load (e.g., ["S22", "S66", "BH76"]).
        If None, loads all 55 canonical GMTKN55 subsets.

    Returns
    -------
    Dataset

    Raises
    ------
    ImportError
        If the gmtkn package is not installed.
    ValueError
        If a requested subset name is not found.
    """
    try:
        import gmtkn
    except ImportError:
        raise ImportError(
            "The 'gmtkn' package is required for GMTKN55 support. "
            "Install it from: https://github.com/obackhouse/gmtkn"
        )

    available = set(gmtkn.sets.keys())
    # Case-insensitive lookup table
    available_lower = {k.lower(): k for k in available}

    if subsets is None:
        subset_names = sorted(GMTKN55_SUBSETS & available)
    else:
        subset_names = []
        for s in subsets:
            canonical = available_lower.get(s.lower())
            if canonical is None:
                raise ValueError(
                    f"Unknown GMTKN subset: '{s}'. "
                    f"Available: {sorted(available)}"
                )
            subset_names.append(canonical)

    reactions = []
    species = {}

    for subset_name in subset_names:
        subset_module = gmtkn.sets[subset_name]
        gmtkn_systems = subset_module.systems

        if not hasattr(subset_module, "reactions"):
            logger.warning("Subset %s has no reactions, skipping", subset_name)
            continue

        gmtkn_reactions = subset_module.reactions

        # Convert systems to Species objects
        for sys_name, sys_data in gmtkn_systems.items():
            qualified_name = f"{subset_name}-{sys_name}"
            if qualified_name not in species:
                species[qualified_name] = Species(
                    name=qualified_name,
                    charge=sys_data["charge"],
                    uhf=sys_data["spin"],
                    elements=[el.capitalize() for el in sys_data["atoms"]],
                    positions=sys_data["coords"],
                )

        # Convert reactions
        for rxn_idx, rxn_data in enumerate(gmtkn_reactions):
            sys_names = rxn_data["systems"]
            stoich_strs = rxn_data["stoichiometry"]
            reference = rxn_data["reference"]

            stoichiometry = {}
            skip = False
            for sys_name, coeff_str in zip(sys_names, stoich_strs):
                qualified_name = f"{subset_name}-{sys_name}"
                if qualified_name not in species:
                    logger.warning(
                        "Reaction %d in %s references unknown species '%s', skipping",
                        rxn_idx + 1, subset_name, sys_name,
                    )
                    skip = True
                    break
                stoichiometry[qualified_name] = int(float(coeff_str))

            if not skip:
                reactions.append(Reaction(
                    subset=subset_name,
                    index=rxn_idx + 1,
                    reference_energy=float(reference),
                    weight=1.0,
                    stoichiometry=stoichiometry,
                ))

    if subsets is None:
        name = "GMTKN55"
    elif len(subset_names) == 1:
        name = subset_names[0]
    else:
        name = "GMTKN55_" + "+".join(subset_names)

    logger.info(
        "\n   Loaded %d reactions, %d unique species from GMTKN (%d subsets)",
        len(reactions), len(species), len(subset_names),
    )

    return Dataset(name=name, reactions=reactions, species=species)


def load_nenci(xyz_dir, subsets=None):
    """Load NENCI-2021 dataset from a directory of dimer XYZ files.

    Each XYZ file contains one dimer configuration. The comment line (line 2)
    encodes monomer metadata and reference energies:

        dimer_charge dimer_mult monA_charge monA_mult monB_charge monB_mult
        natoms_A natoms_B  CCSD(T)/CBS  [other reference levels...]

    Monomers appear in order: the first natoms_A atoms are monomer A, the
    remaining natoms_B atoms are monomer B.

    Files can be organized flat or in subdirectories. When subdirectories
    are present (e.g., S66/, S101/, IonPi/), each subdirectory name is
    used as the subset label. Files in the root directory use "NENCI".

    Parameters
    ----------
    xyz_dir : str
        Directory containing .xyz files, optionally in subdirectories.
    subsets : list of str, optional
        Only load files from these subdirectories.

    Returns
    -------
    Dataset
    """
    if not os.path.isdir(xyz_dir):
        raise FileNotFoundError(f"NENCI directory not found: {xyz_dir}")

    # Collect (subset_label, filepath) pairs
    file_list = []
    for entry in sorted(os.listdir(xyz_dir)):
        entry_path = os.path.join(xyz_dir, entry)
        if os.path.isdir(entry_path):
            if subsets is not None and entry not in subsets:
                continue
            for fname in sorted(os.listdir(entry_path)):
                if fname.endswith(".xyz"):
                    file_list.append((entry, os.path.join(entry_path, fname)))
        elif entry.endswith(".xyz") and subsets is None:
            file_list.append(("NENCI", entry_path))

    species = {}
    reactions = []

    for subset_label, filepath in file_list:
        fname = os.path.basename(filepath)
        stem = os.path.splitext(fname)[0]

        with open(filepath) as f:
            natoms = int(f.readline().strip())
            comment_parts = f.readline().strip().split()

            if len(comment_parts) < 9:
                logger.warning("Skipping %s: comment line has fewer than 9 fields", fname)
                continue

            dimer_charge = int(comment_parts[0])
            dimer_mult = int(comment_parts[1])
            monA_charge = int(comment_parts[2])
            monA_mult = int(comment_parts[3])
            monB_charge = int(comment_parts[4])
            monB_mult = int(comment_parts[5])
            natoms_A = int(comment_parts[6])
            natoms_B = int(comment_parts[7])
            ref_energy = float(comment_parts[8])  # CCSD(T)/CBS in kcal/mol

            elements = []
            positions = []
            for _ in range(natoms):
                parts = f.readline().split()
                elements.append(parts[0].capitalize())
                positions.append([float(parts[1]), float(parts[2]), float(parts[3])])

        if natoms_A + natoms_B != natoms:
            logger.warning(
                "Skipping %s: natoms_A (%d) + natoms_B (%d) != natoms (%d)",
                fname, natoms_A, natoms_B, natoms,
            )
            continue

        # Dimer
        dimer_name = f"NENCI-{stem}"
        species[dimer_name] = Species(
            name=dimer_name,
            charge=dimer_charge,
            uhf=dimer_mult - 1,
            elements=elements,
            positions=positions,
        )

        # Monomer A (first natoms_A atoms)
        monA_name = f"NENCI-{stem}_monA"
        species[monA_name] = Species(
            name=monA_name,
            charge=monA_charge,
            uhf=monA_mult - 1,
            elements=elements[:natoms_A],
            positions=positions[:natoms_A],
        )

        # Monomer B (remaining atoms)
        monB_name = f"NENCI-{stem}_monB"
        species[monB_name] = Species(
            name=monB_name,
            charge=monB_charge,
            uhf=monB_mult - 1,
            elements=elements[natoms_A:],
            positions=positions[natoms_A:],
        )

        # Reaction: E_int = E_dimer - E_monA - E_monB
        reactions.append(Reaction(
            subset=subset_label,
            index=len(reactions) + 1,
            reference_energy=ref_energy,
            weight=1.0,
            stoichiometry={dimer_name: 1, monA_name: -1, monB_name: -1},
        ))

    if not reactions:
        raise ValueError(f"No valid NENCI XYZ files found in {xyz_dir}")

    # Dataset name reflects what was loaded
    loaded_subsets = sorted({r.subset for r in reactions})
    if loaded_subsets == ["NENCI"]:
        name = "NENCI-2021"
    elif len(loaded_subsets) == 1:
        name = loaded_subsets[0]
    else:
        name = "NENCI_" + "+".join(loaded_subsets)

    logger.info(
        "\n   Loaded %d reactions, %d unique species from NENCI-2021",
        len(reactions), len(species),
    )

    return Dataset(name=name, reactions=reactions, species=species)


def load_mpconf196(xyz_dir, molecules=None):
    """Load MPCONF196 conformational energy benchmark from XYZ files.

    MPCONF196 (Brauer et al., JCTC 2018) contains 196 conformers across
    13 peptide and macrocycle molecules. Reference energies are CCSD(T)/CBS
    or DLPNO-CCSD(T)/CBS relative conformational energies (kcal/mol),
    expressed relative to the per-molecule mean.

    Each XYZ file uses a non-standard format: line 1 is the filename (no
    atom count), followed by element/coordinate lines. A ``reference_energies.csv``
    file in the same directory provides benchmark energies with columns:
    ``conformer,molecule,energy_kcal``.

    Reactions measure each conformer's energy relative to the per-molecule
    mean: E(conformer_i) - mean(E_all), giving 196 reactions (one per conformer).

    Parameters
    ----------
    xyz_dir : str
        Directory containing .xyz files and ``reference_energies.csv``.
    molecules : list of str, optional
        Only load conformers for these molecule groups (e.g. ["FGG", "GFA"]).

    Returns
    -------
    Dataset
    """
    if not os.path.isdir(xyz_dir):
        raise FileNotFoundError(f"MPCONF196 directory not found: {xyz_dir}")

    # Read reference energies — check xyz_dir first, then parent directory
    csv_path = os.path.join(xyz_dir, "reference_energies.csv")
    if not os.path.isfile(csv_path):
        csv_path = os.path.join(os.path.dirname(os.path.normpath(xyz_dir)), "reference_energies.csv")
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(
            f"Reference energy file not found in {xyz_dir} or parent directory\n"
            "Expected a CSV with columns: conformer,molecule,energy_kcal"
        )

    ref_energies = {}  # conformer_name -> energy (kcal/mol)
    mol_groups = {}    # conformer_name -> molecule group
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["conformer"]
            mol = row["molecule"]
            if molecules is not None and mol not in molecules:
                continue
            ref_energies[name] = float(row["energy_kcal"])
            mol_groups[name] = mol

    # Read XYZ files (non-standard format: line 1 is title, no atom count)
    species = {}
    for fname in sorted(os.listdir(xyz_dir)):
        if not fname.endswith(".xyz"):
            continue
        stem = os.path.splitext(fname)[0]
        if stem not in ref_energies:
            continue

        filepath = os.path.join(xyz_dir, fname)
        elements = []
        positions = []
        with open(filepath) as f:
            f.readline()  # skip title line
            for line in f:
                parts = line.split()
                if len(parts) >= 4:
                    elements.append(parts[0].capitalize())
                    positions.append([float(parts[1]), float(parts[2]), float(parts[3])])

        if not elements:
            logger.warning("Skipping %s: no atoms found", fname)
            continue

        species[stem] = Species(
            name=stem, charge=0, uhf=0,
            elements=elements, positions=positions,
        )

    # Group conformers by molecule and create mean-relative reactions
    by_molecule = {}
    for name, mol in mol_groups.items():
        if name in species:
            by_molecule.setdefault(mol, []).append(name)

    reactions = []
    for mol in sorted(by_molecule):
        conformers = sorted(by_molecule[mol])
        n = len(conformers)
        for conf in conformers:
            # Reaction: E(conf) - mean(E_all) = conformer energy relative to group mean
            # Stoichiometry: conf gets +1, each conformer (including conf) gets -1/N
            stoich = {c: -1.0 / n for c in conformers}
            stoich[conf] = stoich[conf] + 1.0  # net: (N-1)/N for target conformer
            reactions.append(Reaction(
                subset=mol,
                index=len(reactions) + 1,
                reference_energy=ref_energies[conf],
                weight=1.0,
                stoichiometry=stoich,
            ))

    if not reactions:
        raise ValueError(f"No valid MPCONF196 conformers found in {xyz_dir}")

    # Dataset name
    loaded_molecules = sorted(by_molecule.keys())
    if molecules is None or len(loaded_molecules) == 13:
        name = "MPCONF196"
    elif len(loaded_molecules) == 1:
        name = loaded_molecules[0]
    else:
        name = "MPCONF_" + "+".join(loaded_molecules)

    logger.info(
        "\n   Loaded %d reactions, %d species from MPCONF196 (%d molecules)",
        len(reactions), len(species), len(loaded_molecules),
    )

    return Dataset(name=name, reactions=reactions, species=species)


def resolve_dataset(dataset_arg):
    """Resolve a dataset argument to a Dataset object.

    Supports:
    - A path to a YAML file (existing behavior)
    - "gmtkn55" to load all 55 GMTKN55 subsets
    - "gmtkn55:S22,S66,BH76" to load specific subsets
    - "nenci:/path/to/xyz/dir" to load all NENCI-2021 dimer XYZ files
    - "nenci:/path:S66,IonPi" to load specific NENCI subdirectories
    - "mpconf196:/path/to/xyz/dir" to load MPCONF196 conformer XYZ files
    - "mpconf196:/path:FGG,GFA" to load specific molecule groups

    Parameters
    ----------
    dataset_arg : str

    Returns
    -------
    Dataset
    """
    if dataset_arg.lower().startswith("gmtkn55"):
        parts = dataset_arg.split(":", 1)
        if len(parts) == 2 and parts[1].strip():
            subsets = [s.strip() for s in parts[1].split(",")]
        else:
            subsets = None
        return load_gmtkn55(subsets=subsets)
    elif dataset_arg.lower().startswith("mpconf196:") or dataset_arg.lower().startswith("mpconf:"):
        rest = dataset_arg.split(":", 1)[1]
        # Check for molecule filter: mpconf196:/path:FGG,GFA
        parts = rest.rsplit(":", 1)
        if len(parts) == 2 and "/" not in parts[1] and "\\" not in parts[1]:
            path, mol_str = parts
            molecules = [s.strip() for s in mol_str.split(",")]
        else:
            path = rest
            molecules = None
        return load_mpconf196(path, molecules=molecules)
    elif dataset_arg.lower().startswith("nenci:"):
        rest = dataset_arg[6:]  # strip "nenci:"
        # Check for subset filter: nenci:/path:S66,IonPi
        parts = rest.rsplit(":", 1)
        if len(parts) == 2 and "/" not in parts[1] and "\\" not in parts[1]:
            path, subset_str = parts
            subsets = [s.strip() for s in subset_str.split(",")]
        else:
            path = rest
            subsets = None
        return load_nenci(path, subsets=subsets)
    else:
        return load_dataset(dataset_arg)


def generate_orca_inputs(dataset, output_dir, functional, basis="def2-TZVP",
                         nprocs=8, maxcore=4000, extra_keywords=""):
    """Generate ORCA input files for all unique species in the dataset.

    Creates one .inp file per species for a single-point DFT calculation
    WITHOUT any D3 dispersion correction.

    Parameters
    ----------
    dataset : Dataset
    output_dir : str
        Directory to write .inp files.
    functional : str
        DFT functional name for ORCA (e.g. "PBE", "B3LYP").
    basis : str
        Basis set (default: "def2-TZVP").
    nprocs : int
        Number of parallel processes.
    maxcore : int
        Memory per core in MB.
    extra_keywords : str
        Additional ORCA keywords (e.g. "TightSCF RIJCOSX def2/J").

    Returns
    -------
    list of str
        Paths to generated .inp files.
    """
    os.makedirs(output_dir, exist_ok=True)
    paths = []
    name_map = {}  # safe_name -> original species name

    keywords = f"{functional} {basis}"
    if extra_keywords:
        keywords += f" {extra_keywords}"

    for sp_name, sp in sorted(dataset.species.items()):
        # Sanitize filename
        safe_name = re.sub(r'[^\w\-.]', '_', sp_name)
        inp_path = os.path.join(output_dir, f"{safe_name}.inp")
        name_map[safe_name] = sp_name

        lines = [
            f"! {keywords}",
            f"%pal nprocs {nprocs} end",
            f"%maxcore {maxcore}",
            f"* xyz {sp.charge} {sp.multiplicity}",
        ]
        for el, pos in zip(sp.elements, sp.positions):
            lines.append(f"  {el:2s}  {pos[0]:16.10f}  {pos[1]:16.10f}  {pos[2]:16.10f}")
        lines.append("*")
        lines.append("")

        with open(inp_path, "w") as f:
            f.write("\n".join(lines))

        paths.append(inp_path)

    # Write name mapping so parse_orca_outputs can map filenames back
    _write_name_map(name_map, os.path.join(output_dir, "species_map.csv"))

    logger.info("Generated %d ORCA input files in %s", len(paths), output_dir)
    return paths


def _write_name_map(name_map, csv_path):
    """Write filename-to-species name mapping CSV."""
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["safe_name", "species_name"])
        for safe, original in sorted(name_map.items()):
            writer.writerow([safe, original])


def _read_name_map(csv_path):
    """Read filename-to-species name mapping CSV."""
    name_map = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name_map[row["safe_name"]] = row["species_name"]
    return name_map


def parse_orca_outputs(output_dir):
    """Parse all ORCA .out files in a directory to extract DFT energies.

    Looks for the last occurrence of "FINAL SINGLE POINT ENERGY" in each
    .out file. Species name is derived from the filename stem.

    Parameters
    ----------
    output_dir : str
        Directory containing ORCA .out files.

    Returns
    -------
    dict of str to float
        {species_name: energy_in_hartree}
    """
    energies = {}
    pattern = re.compile(r"FINAL SINGLE POINT ENERGY\s+([-\d.]+)")
    func_pattern = re.compile(r"The (\S+) functional is recognized")
    basis_pattern = re.compile(r"Your calculation utilizes the basis:\s+(\S+)")
    input_echo_pattern = re.compile(r"\|\s+\d+>\s+!\s+(.*)")

    # Load name mapping if available (written by generate_orca_inputs)
    map_path = os.path.join(output_dir, "species_map.csv")
    name_map = _read_name_map(map_path) if os.path.exists(map_path) else {}

    detected_functional = None
    detected_basis = None
    input_echo_keywords = None

    for fname in sorted(os.listdir(output_dir)):
        if not fname.endswith(".out"):
            continue

        filepath = os.path.join(output_dir, fname)
        safe_name = os.path.splitext(fname)[0]
        sp_name = name_map.get(safe_name, safe_name)
        energy = None

        with open(filepath) as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    energy = float(match.group(1))
                if detected_functional is None:
                    fm = func_pattern.search(line)
                    if fm:
                        detected_functional = fm.group(1)
                if detected_basis is None:
                    bm = basis_pattern.search(line)
                    if bm:
                        detected_basis = bm.group(1)
                if input_echo_keywords is None:
                    em = input_echo_pattern.match(line)
                    if em:
                        input_echo_keywords = em.group(1).split()

        if energy is not None:
            energies[sp_name] = energy
        else:
            logger.warning("!  No FINAL SINGLE POINT ENERGY found in %s", fname)

    # Fallback: extract functional from input echo keywords
    if detected_functional is None and input_echo_keywords:
        from .pars import bj_parms, resolve_functional, zero_parms
        known = set(bj_parms) | set(zero_parms)
        for kw in input_echo_keywords:
            # Strip trailing /G (Gaussian compatibility notation)
            clean = kw.rstrip("/G").rstrip("/g")
            canonical = resolve_functional(clean)
            if canonical in known or clean.upper() in known:
                detected_functional = canonical if canonical in known else clean.upper()
                break

    logger.info("Parsed %d ORCA output files from %s", len(energies), output_dir)
    return energies, detected_functional, detected_basis


def write_energies_csv(energies, csv_path, functional=None, basis=None):
    """Write species DFT energies to a CSV file.

    Parameters
    ----------
    energies : dict of str to float
        {species_name: energy_in_hartree}
    csv_path : str
        Output CSV path.
    functional : str or None
        DFT functional name to store as metadata comment.
    basis : str or None
        Basis set name to store as metadata comment.
    """
    with open(csv_path, "w", newline="") as f:
        if functional:
            f.write(f"# functional: {functional}\n")
        if basis:
            f.write(f"# basis: {basis}\n")
        writer = csv.writer(f)
        writer.writerow(["species_name", "energy_hartree"])
        for name in sorted(energies):
            writer.writerow([name, f"{energies[name]:.12f}"])


def read_energies_csv(csv_path):
    """Read species DFT energies from a CSV file.

    Parameters
    ----------
    csv_path : str
        CSV file with columns: species_name, energy_hartree

    Returns
    -------
    (energies, functional, basis) : tuple
        energies is {species_name: energy_in_hartree},
        functional and basis are str or None (parsed from comment lines).
    """
    energies = {}
    functional = None
    basis = None
    with open(csv_path) as f:
        lines = f.readlines()

    data_lines = []
    for line in lines:
        if line.startswith("# functional:"):
            functional = line.split(":", 1)[1].strip()
        elif line.startswith("# basis:"):
            basis = line.split(":", 1)[1].strip()
        elif not line.startswith("#"):
            data_lines.append(line)

    reader = csv.DictReader(data_lines)
    for row in reader:
        energies[row["species_name"]] = float(row["energy_hartree"])
    return energies, functional, basis
