"""
benchmark.py - Load DietGMTKN55 benchmark datasets and generate ORCA inputs.

Provides tools for loading benchmark reaction datasets from DietGMTKN55 YAML
files, generating ORCA input files for dispersion-free DFT calculations,
and parsing ORCA outputs to extract single-point energies.
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
    """Load a DietGMTKN55 YAML file into a Dataset object.

    Parameters
    ----------
    yaml_path : str
        Path to a DietGMTKN55 AllElements YAML file
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
        from .pars import bj_parms, zero_parms, resolve_functional
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
