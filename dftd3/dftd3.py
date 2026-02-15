"""
dftd3.py - Python implementation of Grimme's D3 dispersion correction.

This is a translation of Grimme's D3 Fortran code into Python. It computes
D3-density independent dispersion terms for molecular geometries. Zero and
Becke-Johnson damping schemes are both implemented.

Written by Rob Paton and Kelvin Jackson.

References:
  [1] S. Grimme, J. Antony, S. Ehrlich, H. Krieg,
      J. Chem. Phys. 132, 154104 (2010).  DOI: 10.1063/1.3382344
  [2] S. Grimme, S. Ehrlich, L. Goerigk,
      J. Comput. Chem. 32, 1456-1465 (2011).  DOI: 10.1002/jcc.21759

License: MIT
"""

import logging
import math
import sys
from argparse import ArgumentParser
from datetime import datetime

import numpy as np
from cclib.io import ccread

try:
    from .pars import (
        bj_parms,
        elements,
        pars,
        r0ab,
        r2r4,
        rcov,
        resolve_functional,
        zero_parms,
    )
except ImportError:
    from pars import (
        bj_parms,
        elements,
        pars,
        r0ab,
        r2r4,
        rcov,
        resolve_functional,
        zero_parms,
    )

logger = logging.getLogger(__name__)

CITATION_ZERO = (
    "Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H. "
    "A Consistent and Accurate Ab Initio Parametrization of Density Functional "
    "Dispersion Correction (DFT-D) for the 94 Elements H-Pu. "
    "J. Chem. Phys. 2010, 132, 154104."
)
CITATION_ZERO_SHORT = "Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H. J. Chem. Phys. 2010, 132, 154104."
CITATION_BJ = (
    "Grimme, S.; Ehrlich, S.; Goerigk, L. "
    "Effect of the Damping Function in Dispersion Corrected Density Functional Theory. "
    "J. Comput. Chem. 2011, 32, 1456\u20131465."
)
CITATION_BJ_SHORT = "Grimme, S.; Ehrlich, S.; Goerigk, L. J. Comput. Chem. 2011, 32, 1456\u20131465."

SUPPORTED_EXTENSIONS = {"out", "log", "sdf", "xyz", "pdb"}

# Periodic table with empty string at index 0 so that atomic numbers map directly
PERIODIC_TABLE = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Uub", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
AUTOANG = 0.52917726
AUTOKCAL = 627.509541

# Exponents for distance-dependent damping factors (R6, R8, R10)
ALPHA6 = 14
ALPHA8 = ALPHA6 + 2

# Fractional connectivity constants
K1 = 16.0
K2 = 4.0 / 3.0
K3 = -4.0

# D3 is parameterized up to element 94
MAX_ELEM = 94
MAXC = 5

# ---------------------------------------------------------------------------
# Module-level pre-computed arrays (loaded once at import time)
# ---------------------------------------------------------------------------

def _build_r_matrix():
    """Convert diatomic cutoff radii from atomic units to Angstrom."""
    r = np.zeros((MAX_ELEM, MAX_ELEM))
    k = 0
    for i in range(MAX_ELEM):
        for j in range(i + 1):
            r[i, j] = r0ab[k] / AUTOANG
            r[j, i] = r0ab[k] / AUTOANG
            k += 1
    return r


def _build_r2r4():
    """Transform PBE0/def2-QZVP multipole coefficients."""
    r = np.array(r2r4[:MAX_ELEM])
    indices = np.arange(1, MAX_ELEM + 1, dtype=float)
    return np.sqrt(0.5 * r * np.sqrt(indices))


def _build_c6ab():
    """Load C6 reference data from the parameter array (legacy list format)."""
    nlines = 32385
    c6ab = [[None] * MAX_ELEM for _ in range(MAX_ELEM)]

    for iat in range(MAX_ELEM):
        for jat in range(MAX_ELEM):
            c6ab[iat][jat] = [[0] * MAXC for _ in range(MAXC)]

    for nn in range(nlines):
        kk = nn * 5
        iadr = 0
        jadr = 0
        iat = int(pars[kk + 1]) - 1
        jat = int(pars[kk + 2]) - 1

        while iat > 99:
            iadr += 1
            iat -= 100
        while jat > 99:
            jadr += 1
            jat -= 100

        c6ab[iat][jat][iadr][jadr] = [pars[kk], pars[kk + 3], pars[kk + 4]]
        c6ab[jat][iat][jadr][iadr] = [pars[kk], pars[kk + 4], pars[kk + 3]]

    return c6ab


def _build_c6ab_numpy():
    """Load C6 reference data into numpy arrays for vectorized computation.

    Returns (c6_ref, cn1_ref, cn2_ref, mxc) where each ref array has
    shape (94, 94, 5, 5) and mxc has shape (94,).
    """
    nlines = 32385
    c6_ref = np.zeros((MAX_ELEM, MAX_ELEM, MAXC, MAXC))
    cn1_ref = np.zeros((MAX_ELEM, MAX_ELEM, MAXC, MAXC))
    cn2_ref = np.zeros((MAX_ELEM, MAX_ELEM, MAXC, MAXC))

    for nn in range(nlines):
        kk = nn * 5
        iadr = 0
        jadr = 0
        iat = int(pars[kk + 1]) - 1
        jat = int(pars[kk + 2]) - 1

        while iat > 99:
            iadr += 1
            iat -= 100
        while jat > 99:
            jadr += 1
            jat -= 100

        c6_ref[iat, jat, iadr, jadr] = pars[kk]
        cn1_ref[iat, jat, iadr, jadr] = pars[kk + 3]
        cn2_ref[iat, jat, iadr, jadr] = pars[kk + 4]

        c6_ref[jat, iat, jadr, iadr] = pars[kk]
        cn1_ref[jat, iat, jadr, iadr] = pars[kk + 4]
        cn2_ref[jat, iat, jadr, iadr] = pars[kk + 3]

    # Max coordination reference count per element
    mxc = np.zeros(MAX_ELEM, dtype=int)
    for j in range(MAX_ELEM):
        for ll in range(MAXC):
            if c6_ref[j, j, ll, ll] > 0:
                mxc[j] += 1

    return c6_ref, cn1_ref, cn2_ref, mxc


_r = _build_r_matrix()
_r2r4 = _build_r2r4()
_c6ab = _build_c6ab()
_c6_ref, _cn1_ref, _cn2_ref, _mxc_arr = _build_c6ab_numpy()
_rcov = np.array(rcov[:MAX_ELEM])


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def parse_int_set(input_str=""):
    """Expand a string like '1-4, 7-9' into {1, 2, 3, 4, 7, 8, 9}."""
    selection = set()
    invalid = set()
    for token in (t.strip() for t in input_str.split(",")):
        if not token:
            continue
        if token.startswith("<"):
            token = f"1-{token[1:]}"
        try:
            selection.add(int(token))
        except ValueError:
            try:
                bounds = sorted(int(x.strip()) for x in token.split("-"))
                selection.update(range(bounds[0], bounds[-1] + 1))
            except ValueError:
                invalid.add(token)
    if invalid:
        logger.warning("Invalid set tokens: %s", invalid)
    return selection


def _format_atom_range(atoms):
    """Format a sorted list of integers into a compact range string, e.g. [1,2,3,5] -> '1-3, 5'."""
    if not atoms:
        return ""
    sorted_atoms = sorted(atoms)
    ranges = []
    start = end = sorted_atoms[0]
    for n in sorted_atoms[1:]:
        if n == end + 1:
            end = n
        else:
            ranges.append(f"{start}-{end}" if end > start else str(start))
            start = end = n
    ranges.append(f"{start}-{end}" if end > start else str(start))
    return ", ".join(ranges)


def _find_fragments(atomtype, xco, yco, zco, scale=1.3):
    """Identify separate molecules from Cartesian coordinates using covalent radii.

    Atoms are considered bonded if their distance is less than
    ``scale * (rcov[A] + rcov[B])``.  Connected components of the
    resulting graph are returned as a list of lists of 1-indexed atom
    numbers.
    """
    natom = len(atomtype)

    # Build adjacency list via covalent radii
    adj = [[] for _ in range(natom)]
    for i in range(natom):
        zi = _element_index(atomtype[i])
        for j in range(i + 1, natom):
            zj = _element_index(atomtype[j])
            dist = math.sqrt(
                (xco[i] - xco[j]) ** 2 + (yco[i] - yco[j]) ** 2 + (zco[i] - zco[j]) ** 2
            )
            if dist < scale * (rcov[zi] + rcov[zj]):
                adj[i].append(j)
                adj[j].append(i)

    # BFS to find connected components
    visited = [False] * natom
    fragments = []
    for start in range(natom):
        if visited[start]:
            continue
        frag = []
        queue = [start]
        visited[start] = True
        while queue:
            node = queue.pop(0)
            frag.append(node + 1)  # 1-indexed
            for neighbor in adj[node]:
                if not visited[neighbor]:
                    visited[neighbor] = True
                    queue.append(neighbor)
        fragments.append(frag)
    return fragments


def _element_index(symbol):
    """Return the 0-based index of an element symbol in the elements array."""
    for i, el in enumerate(elements):
        if el == symbol:
            return i
    raise ValueError(f"Unknown element: {symbol}")


def _getc6(mxc, atomtype, cn, a, b):
    """Obtain the C6 coefficient for interaction between atoms a and b."""
    iat = _element_index(atomtype[a])
    jat = _element_index(atomtype[b])

    c6mem = -1.0e99
    rsum = 0.0
    csum = 0.0
    c6 = 0.0

    for i in range(mxc[iat]):
        for j in range(mxc[jat]):
            entry = _c6ab[iat][jat][i][j]
            if isinstance(entry, list):
                c6 = entry[0]
                if c6 > 0:
                    c6mem = c6
                    cn1 = entry[1]
                    cn2 = entry[2]
                    r_val = (cn1 - cn[a]) ** 2 + (cn2 - cn[b]) ** 2
                    tmp1 = math.exp(K3 * r_val)
                    rsum += tmp1
                    csum += tmp1 * c6

    return csum / rsum if rsum > 0 else c6mem


def _getc6_all_pairs(elem_idx, cn):
    """Compute C6 coefficients for all atom pairs (vectorized).

    Groups atoms by element type and uses numpy broadcasting over the
    (max 5x5) reference coordination-number grid for each element pair.

    Parameters
    ----------
    elem_idx : np.ndarray, shape (natom,), dtype int
        0-based element indices.
    cn : np.ndarray, shape (natom,)
        Coordination numbers.

    Returns
    -------
    c6_matrix : np.ndarray, shape (natom, natom)
    """
    natom = len(elem_idx)
    c6_matrix = np.empty((natom, natom))

    unique_elems = np.unique(elem_idx)
    MAX_BLOCK = 20_000_000  # max entries in the 4D broadcast tensor

    for ei in unique_elems:
        atoms_i = np.where(elem_idx == ei)[0]
        cn_i = cn[atoms_i]
        ni = len(atoms_i)
        mi = int(_mxc_arr[ei])

        for ej in unique_elems:
            if ej < ei:
                continue

            atoms_j = np.where(elem_idx == ej)[0]
            cn_j = cn[atoms_j]
            nj = len(atoms_j)
            mj = int(_mxc_arr[ej])

            if mi == 0 or mj == 0:
                c6_matrix[np.ix_(atoms_i, atoms_j)] = 0.0
                if ei != ej:
                    c6_matrix[np.ix_(atoms_j, atoms_i)] = 0.0
                continue

            c6_block = _c6_ref[ei, ej, :mi, :mj]    # (mi, mj)
            cn1_block = _cn1_ref[ei, ej, :mi, :mj]   # (mi, mj)
            cn2_block = _cn2_ref[ei, ej, :mi, :mj]   # (mi, mj)
            valid = c6_block > 0
            c6mem = float(np.max(c6_block[valid])) if np.any(valid) else 0.0

            # Chunk over atoms_i to cap memory
            chunk_ni = max(1, MAX_BLOCK // max(1, nj * mi * mj))
            results = []

            for a_start in range(0, ni, chunk_ni):
                a_end = min(a_start + chunk_ni, ni)
                cn_i_chunk = cn_i[a_start:a_end]

                # Broadcasting: (chunk, 1, mi, mj) and (1, nj, mi, mj)
                dcn1 = (cn1_block[np.newaxis, np.newaxis, :, :]
                        - cn_i_chunk[:, np.newaxis, np.newaxis, np.newaxis])
                dcn2 = (cn2_block[np.newaxis, np.newaxis, :, :]
                        - cn_j[np.newaxis, :, np.newaxis, np.newaxis])

                r_val = dcn1 ** 2 + dcn2 ** 2  # (chunk, nj, mi, mj)
                weights = np.exp(K3 * r_val)
                weights *= valid[np.newaxis, np.newaxis, :, :]

                csum = np.sum(weights * c6_block[np.newaxis, np.newaxis, :, :],
                              axis=(2, 3))
                rsum = np.sum(weights, axis=(2, 3))

                results.append(np.where(rsum > 0, csum / rsum, c6mem))

            c6_result = np.vstack(results)  # (ni, nj)
            c6_matrix[np.ix_(atoms_i, atoms_j)] = c6_result
            if ei != ej:
                c6_matrix[np.ix_(atoms_j, atoms_i)] = c6_result.T

    return c6_matrix


def _ncoord_vectorized(coords, elem_idx):
    """Calculate fractional coordination numbers (vectorized).

    Parameters
    ----------
    coords : np.ndarray, shape (natom, 3)
    elem_idx : np.ndarray, shape (natom,), dtype int
        0-based element indices (atomic number - 1).

    Returns
    -------
    cn : np.ndarray, shape (natom,)
    dist_matrix : np.ndarray, shape (natom, natom)
        Pairwise distance matrix in Angstrom (reusable by caller).
    """
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]  # (N, N, 3)
    dist = np.sqrt(np.sum(diff ** 2, axis=2))  # (N, N)

    rcov_arr = _rcov[elem_idx]  # (natom,)
    rco = (rcov_arr[:, np.newaxis] + rcov_arr[np.newaxis, :]) * K2  # (N, N)

    safe_dist = np.where(dist < 1e-8, np.inf, dist)
    rr = rco / safe_dist

    cn_contrib = 1.0 / (1.0 + np.exp(-K1 * (rr - 1.0)))
    np.fill_diagonal(cn_contrib, 0.0)

    return np.sum(cn_contrib, axis=1), dist


def _ncoord(natom, atomtype, xco, yco, zco):
    """Calculate fractional coordination numbers for all atoms."""
    coords = np.column_stack([xco, yco, zco])
    elem_idx = np.array([_element_index(s) for s in atomtype])
    cn, _ = _ncoord_vectorized(coords, elem_idx)
    return cn.tolist()


def _lin(i1, i2):
    """Triangular index for a pair (i1, i2)."""
    imax = max(i1, i2)
    imin = min(i1, i2)
    return imin + imax * (imax - 1) // 2


def _lookup_functional(name, parm_dict):
    """Look up a functional in a parameter dictionary, resolving aliases.

    Returns (canonical_name, params_tuple) or (None, None).
    """
    canonical = resolve_functional(name)
    if canonical in parm_dict:
        return canonical, parm_dict[canonical]
    # Also try the raw upper-cased name (in case the alias mapping is incomplete)
    upper = name.upper()
    if upper in parm_dict:
        return upper, parm_dict[upper]
    return None, None


# ---------------------------------------------------------------------------
# ORCA output parser (fallback when cclib fails)
# ---------------------------------------------------------------------------

class _OrcaData:
    """Lightweight container mimicking cclib parsed data for ORCA outputs."""

    def __init__(self, atomnos, atomcoords, functional=None):
        self.atomnos = np.array(atomnos, dtype=int)
        self.atomcoords = np.array([atomcoords])
        self.metadata = {"functional": functional} if functional else {}


def _parse_orca(filepath):
    """Parse an ORCA output file for coordinates and functional name.

    Extracts the last 'CARTESIAN COORDINATES (ANGSTROEM)' block and the
    functional name from the 'The <name> functional is recognized' line.
    Returns an _OrcaData object compatible with CalcD3, or None on failure.
    """
    with open(filepath) as f:
        lines = f.readlines()

    # --- Parse coordinates (use last block in case of geometry optimization) ---
    coord_start = None
    for i, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            coord_start = i + 2  # skip header + dashes

    if coord_start is None:
        return None

    atomnos = []
    coords = []
    for line in lines[coord_start:]:
        parts = line.split()
        if len(parts) != 4:
            break
        sym = parts[0]
        if sym not in PERIODIC_TABLE:
            break
        atomnos.append(PERIODIC_TABLE.index(sym))
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    if not atomnos:
        return None

    # --- Parse functional name ---
    functional = None
    for line in lines:
        if "functional is recognized" in line:
            # "The B3LYP functional is recognized"
            parts = line.split()
            idx = parts.index("functional")
            if idx > 0:
                functional = parts[idx - 1]
            break

    return _OrcaData(atomnos, coords, functional)


def read_file(filepath):
    """Parse a molecular structure file, with ORCA fallback.

    Tries cclib first.  If cclib fails or returns None (common with newer
    ORCA 6 outputs), falls back to a lightweight manual ORCA parser.
    """
    try:
        logging.getLogger("cclib").setLevel(logging.ERROR)
        try:
            from openbabel import openbabel
            openbabel.obErrorLog.SetOutputLevel(0)
        except ImportError:
            pass
        data = ccread(filepath)
        if data is not None:
            return data
    except Exception:
        pass

    # Fallback: try manual ORCA parser
    data = _parse_orca(filepath)
    if data is not None:
        logger.info("Parsed %s using ORCA fallback parser", filepath)
    return data


# ---------------------------------------------------------------------------
# Main D3 calculator
# ---------------------------------------------------------------------------

class CalcD3:
    """Compute the DFT-D3 dispersion correction for a molecular geometry.

    Parameters
    ----------
    file_data : cclib parsed data
        Parsed structure from cclib (must have atomnos and atomcoords).
    functional : str or None
        Name of the density functional.
    damp : str
        Damping scheme: 'zero' or 'bj'.
    s6, rs6, s8, a1, a2 : float
        Manual damping parameters (override defaults if all non-zero).
    abc : bool
        Include the 3-body Axilrod-Teller-Muto term.
    intermolecular : str or False
        Fragment specification for intermolecular-only calculation.
    pairwise : bool
        Print pairwise energy breakdown.

    Attributes
    ----------
    attractive_r6_vdw : float
        R^-6 contribution in kcal/mol.
    attractive_r8_vdw : float
        R^-8 contribution in kcal/mol.
    repulsive_abc : float
        3-body ATM contribution in kcal/mol.
    """

    def __init__(self, file_data, functional, damp="zero",
                 s6=None, rs6=None, s8=None, a1=None, a2=None,
                 abc=False, intermolecular=False, pairwise=False,
                 cutoff=None):

        atom_nums = file_data.atomnos
        coords = file_data.atomcoords[-1]  # (natom, 3) numpy array
        natom = len(atom_nums)
        elem_idx = (atom_nums - 1).astype(int)  # 0-based element indices
        atomtype = [PERIODIC_TABLE[int(atno)] for atno in atom_nums]

        # For _find_fragments compatibility
        xco = coords[:, 0].tolist()
        yco = coords[:, 1].tolist()
        zco = coords[:, 2].tolist()

        self.attractive_r6_vdw = 0.0
        self.attractive_r8_vdw = 0.0
        self.repulsive_abc = 0.0
        self.pairwise_terms = []

        # rs8 is always 1.0 for the standard D3 scheme
        rs8 = 1.0

        # -------------------------------------------------------------------
        # Resolve damping parameters
        # -------------------------------------------------------------------
        if damp == "zero":
            if s6 is None or rs6 is None or s8 is None:
                if functional is not None:
                    _, prm = _lookup_functional(functional, zero_parms)
                    if prm is not None:
                        s6, rs6, s8 = prm
                else:
                    logger.warning("No functional information could be read!")

        elif damp == "bj":
            if s6 is None or s8 is None or a1 is None or a2 is None:
                if functional is not None:
                    _, prm = _lookup_functional(functional, bj_parms)
                    if prm is not None:
                        s6, a1, s8, a2 = prm
                else:
                    logger.warning("No functional information could be read!")

        # -------------------------------------------------------------------
        # Intermolecular fragments
        # -------------------------------------------------------------------
        mols = None
        self.fragments = None
        if intermolecular:
            if intermolecular == "auto":
                fragments = _find_fragments(atomtype, xco, yco, zco)
                if len(fragments) < 2:
                    logger.warning("Only 1 molecule detected — cannot compute intermolecular dispersion")
                    intermolecular = False
                else:
                    self.fragments = fragments
                    logger.info("  Auto-detected %d fragments:", len(fragments))
                    mols = [0] * natom
                    for frag_idx, frag in enumerate(fragments):
                        for at in frag:
                            mols[at - 1] = frag_idx
            else:
                mols = [-1] * natom
                all_atoms = set(range(1, natom + 1))
                assigned = set()
                for frag_idx, frag_str in enumerate(intermolecular.split(":")):
                    frag_atoms = parse_int_set(frag_str)
                    overlap = assigned & frag_atoms
                    if overlap:
                        logger.warning("Atom(s) %s assigned to multiple fragments!", sorted(overlap))
                    for at in frag_atoms:
                        if at < 1 or at > natom:
                            logger.warning("Atom %d is out of range (1-%d), ignoring", at, natom)
                            continue
                        mols[int(at) - 1] = frag_idx
                    assigned |= frag_atoms
                missing = all_atoms - assigned
                if missing:
                    logger.warning("Atom(s) %s not assigned to any fragment — "
                                   "their interactions will be included", sorted(missing))
                    for at in missing:
                        mols[at - 1] = -1  # unique "fragment" so no interactions are skipped

        # -------------------------------------------------------------------
        # Vectorized computation
        # -------------------------------------------------------------------

        # 1. Coordination numbers + distance matrix (reused below)
        cn, dist_ang = _ncoord_vectorized(coords, elem_idx)

        # 2. C6 coefficients for all pairs
        c6_matrix = _getc6_all_pairs(elem_idx, cn)

        # 3. C8 = 3 * C6 * r2r4[A] * r2r4[B]
        r2r4_atoms = _r2r4[elem_idx]  # (natom,)
        c8_matrix = 3.0 * c6_matrix * np.outer(r2r4_atoms, r2r4_atoms)

        # 4. Extract upper-triangle pairs (j < k)
        j_idx, k_idx = np.triu_indices(natom, k=1)
        dist_ang_pairs = dist_ang[j_idx, k_idx]
        c6_pairs = c6_matrix[j_idx, k_idx]
        c8_pairs = c8_matrix[j_idx, k_idx]

        # 5. Apply cutoff mask
        if cutoff is not None:
            mask = dist_ang_pairs <= cutoff
            j_idx = j_idx[mask]
            k_idx = k_idx[mask]
            dist_ang_pairs = dist_ang_pairs[mask]
            c6_pairs = c6_pairs[mask]
            c8_pairs = c8_pairs[mask]

        # 6. Distance in atomic units
        dist_au = dist_ang_pairs / AUTOANG

        # 7. Scalefactor for intermolecular mode
        if mols is not None:
            mols_arr = np.array(mols)
            scalefactor = np.where(mols_arr[j_idx] == mols_arr[k_idx], 0.0, 1.0)
        else:
            scalefactor = np.ones(len(j_idx))

        # 8. Damping and energy terms
        atomA = elem_idx[j_idx]
        atomB = elem_idx[k_idx]

        if damp == "zero":
            rr = _r[atomA, atomB] / dist_au
            damp6 = 1.0 / (1.0 + 6.0 * (rs6 * rr) ** ALPHA6)
            damp8 = 1.0 / (1.0 + 6.0 * (rs8 * rr) ** ALPHA8)

            r6_terms = -s6 * c6_pairs * damp6 / dist_au ** 6 * AUTOKCAL * scalefactor
            r8_terms = -s8 * c8_pairs * damp8 / dist_au ** 8 * AUTOKCAL * scalefactor

        elif damp == "bj":
            rr = np.sqrt(c8_pairs / c6_pairs)
            tmp1 = a1 * rr + a2
            damp6 = tmp1 ** 6
            damp8 = tmp1 ** 8

            r6_terms = -s6 * c6_pairs / (dist_au ** 6 + damp6) * AUTOKCAL * scalefactor
            r8_terms = -s8 * c8_pairs / (dist_au ** 8 + damp8) * AUTOKCAL * scalefactor

        # 9. Accumulate totals
        self.attractive_r6_vdw = float(np.sum(r6_terms))
        self.attractive_r8_vdw = float(np.sum(r8_terms))

        # 10. Pairwise output
        if pairwise:
            nonzero = scalefactor != 0
            for idx in np.where(nonzero)[0]:
                self.pairwise_terms.append((
                    int(j_idx[idx]) + 1, int(k_idx[idx]) + 1,
                    float(r6_terms[idx]), float(r8_terms[idx])
                ))

        # -------------------------------------------------------------------
        # 3-body ATM term
        # -------------------------------------------------------------------
        if abc:
            npairs = natom * (natom - 1) // 2 + 1
            icomp = np.zeros(npairs, dtype=bool)
            cc6ab_arr = np.zeros(npairs)
            r2ab = np.zeros(npairs)
            dmp_arr = np.zeros(npairs)

            # Build pair data indexed by _lin(j, k)
            imax = np.maximum(j_idx, k_idx)
            imin = np.minimum(j_idx, k_idx)
            lin_idx = imin + imax * (imax - 1) // 2

            if damp == "zero":
                rr_3body = _r[atomA, atomB] / dist_au
                dmp_vals = (1.0 / rr_3body) ** (1.0 / 3.0)
            elif damp == "bj":
                rr_3body = np.sqrt(c8_pairs / c6_pairs)
                dmp_vals = (1.0 / rr_3body) ** (1.0 / 3.0)

            icomp[lin_idx] = True
            cc6ab_arr[lin_idx] = np.sqrt(c6_pairs)
            r2ab[lin_idx] = dist_au ** 2
            dmp_arr[lin_idx] = dmp_vals

            e63 = 0.0
            for iat in range(natom):
                for jat in range(natom):
                    ij = int(_lin(jat, iat))
                    if not icomp[ij]:
                        continue
                    for kat in range(jat, natom):
                        ik = int(_lin(kat, iat))
                        jk = int(_lin(kat, jat))

                        if kat > jat > iat and icomp[ik] and icomp[jk]:
                            if mols is not None and mols[iat] == mols[jat] == mols[kat]:
                                continue
                            rav = (4.0 / 3.0) / (dmp_arr[ik] * dmp_arr[jk] * dmp_arr[ij])
                            tmp = 1.0 / (1.0 + 6.0 * rav ** ALPHA6)

                            c9 = cc6ab_arr[ij] * cc6ab_arr[ik] * cc6ab_arr[jk]
                            d2 = [r2ab[ij], r2ab[jk], r2ab[ik]]
                            t1 = (d2[0] + d2[1] - d2[2]) / math.sqrt(d2[0] * d2[1])
                            t2 = (d2[0] + d2[2] - d2[1]) / math.sqrt(d2[0] * d2[2])
                            t3 = (d2[2] + d2[1] - d2[0]) / math.sqrt(d2[1] * d2[2])
                            ang = 0.375 * t1 * t2 * t3 + 1.0
                            e63 += tmp * c9 * ang / (d2[0] * d2[1] * d2[2]) ** 1.5

            self.repulsive_abc = s6 * e63 * AUTOKCAL


# Keep old name as an alias for backwards compatibility
calcD3 = CalcD3


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    # Dispatch to optimize subcommand if requested
    if len(sys.argv) > 1 and sys.argv[1] == "optimize":
        from dftd3.optimize import optimize_main
        return optimize_main(sys.argv[2:])

    parser = ArgumentParser(
        description="Compute Grimme's DFT-D3 dispersion correction.",
        epilog=f"{CITATION_ZERO}\n{CITATION_BJ}",
    )
    parser.add_argument("files", nargs="*", help="Input structure file(s)")
    parser.add_argument("-v", dest="verbose", action="store_true", default=False,
                        help="Turn on verbose printing")
    parser.add_argument("--damp", dest="damp", default="zero", type=str.lower,
                        choices=("zero", "bj"),
                        help="Type of D3-damping function (zero, bj)")
    parser.add_argument("--func", dest="functional", default=None,
                        help="Use default D3 parameters for this density functional")
    parser.add_argument("--s6", dest="s6", default=None, type=float,
                        help="s6 parameter (used in zero and bj damping)")
    parser.add_argument("--rs6", dest="rs6", default=None, type=float,
                        help="rs6 parameter used in zero damping")
    parser.add_argument("--s8", dest="s8", default=None, type=float,
                        help="s8 parameter used in zero damping")
    parser.add_argument("--a1", dest="a1", default=None, type=float,
                        help="a1 parameter used in bj damping")
    parser.add_argument("--a2", dest="a2", default=None, type=float,
                        help="a2 parameter used in bj damping")
    parser.add_argument("--kcal", dest="kcal", action="store_true", default=False,
                        help="Print energies in kcal/mol")
    parser.add_argument("--abc", dest="threebody", action="store_true", default=False,
                        help="Turn on repulsive 3-body (ABC) term")
    parser.add_argument("--pw", dest="pairwise", nargs="*", type=int, default=None,
                        help="Print pairwise dispersion terms (optionally specify atom indices, e.g. --pw 1 2)")
    parser.add_argument("--cutoff", dest="cutoff", type=float, default=None,
                        help="Distance cutoff in Angstrom (default: no cutoff)")
    parser.add_argument("--im", dest="intermolecular", type=str, default=False,
                        help="Compute only intermolecular dispersion terms")
    parser.add_argument("--cite", dest="cite", action="store_true", default=False,
                        help="Print citation information and exit")

    options, extra_args = parser.parse_known_args()

    # Configure logging
    logging.basicConfig(
        format="%(message)s",
        level=logging.INFO if options.verbose else logging.WARNING,
    )

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(rf"""    ____       __
   /\  _`\   /'__`\
   \ \ \/\ \/\_\L\ \
    \ \ \ \ \/_/_\_<_
     \ \ \_\ \/\ \L\ \
      \ \____/\ \____/
       \/___/  \/___/   ¯\_(ツ)_/¯  {now}
""")

    if options.cite:
        print("Please cite the following when using DFT-D3 corrections:\n")
        print(f"  [1] {CITATION_ZERO}")
        print(f"  [2] {CITATION_BJ}")
        print()
        return 0

    # Collect input files from positional args and any extra args
    files = []
    for arg in list(options.files) + extra_args:
        if "." in arg and arg.rsplit(".", 1)[-1] in SUPPORTED_EXTENSIONS:
            files.append(arg)
    # Fallback: scan sys.argv for files (backwards compatibility)
    if not files:
        for arg in sys.argv[1:]:
            if "." in arg and arg.rsplit(".", 1)[-1] in SUPPORTED_EXTENSIONS:
                files.append(arg)

    if not files:
        print("\nNo valid files found!\n")
        return 1

    # Parse all input files
    parm_dict = zero_parms if options.damp == "zero" else bj_parms
    parsed_files = []
    for filepath in files:
        data = read_file(filepath)
        if data is None:
            logger.warning("Could not parse file: %s", filepath)
            continue
        parsed_files.append((filepath, data))

    if not parsed_files:
        print("\nNo valid files could be parsed!\n")
        return 1

    # Resolve the functional name through aliases
    dft_functional = None
    if options.functional is not None:
        canonical = resolve_functional(options.functional)
        if canonical in parm_dict:
            dft_functional = canonical
        else:
            # Try the raw upper-cased name
            upper = options.functional.upper()
            if upper in parm_dict:
                dft_functional = upper
            else:
                print(f"\nUnable to match requested functional '{options.functional}' "
                      f"to stored {options.damp}-damping parameters!")
                print(f"Available functionals: {', '.join(sorted(parm_dict.keys()))}\n")
                return 1
    else:
        # Auto-detect functional from parsed files and check consistency
        detected = {}
        for filepath, data in parsed_files:
            try:
                parsed_func = data.metadata.get("functional")
                if parsed_func:
                    canonical = resolve_functional(parsed_func)
                    if canonical in parm_dict:
                        detected[filepath] = canonical
                    elif parsed_func.upper() in parm_dict:
                        detected[filepath] = parsed_func.upper()
            except (AttributeError, KeyError):
                pass

        unique_functionals = set(detected.values())
        if len(unique_functionals) > 1:
            print("\nInconsistent functionals detected across input files:")
            for fp, func in detected.items():
                print(f"  {fp}: {func}")
            print("\nPlease use --func to specify a single functional.\n")
            return 1
        elif len(unique_functionals) == 1:
            dft_functional = unique_functionals.pop()

    # Check that we have a functional or manual parameters
    manual_params = (options.s6 is not None and options.s8 is not None and
                     (options.damp == "zero" and options.rs6 is not None or
                      options.damp == "bj" and options.a1 is not None and options.a2 is not None))
    if dft_functional is None and not manual_params:
        print("\nError: No functional detected. When using XYZ/PDB/SDF files, you must specify")
        print("a functional with --func <name>, or provide damping parameters manually")
        print(f"(e.g. --s6 --s8 {'--rs6' if options.damp == 'zero' else '--a1 --a2'}).")
        print(f"\nAvailable {options.damp}-damping functionals: {', '.join(sorted(parm_dict.keys()))}\n")
        return 1

    # Table formatting constants
    name_w = 42  # width of the species column
    c1_w = 10   # D3(R6)
    c2_w = 10   # D3(R8)
    c3_w = 10   # ABC
    c4_w = 15   # Total
    fmt_dp = 2 if options.kcal else 6  # decimal places
    total_label = "Etot (kcal/mol)" if options.kcal else "Etot (Hartree)"
    do_pairwise = options.pairwise is not None
    if do_pairwise:
        species_label = f"{'Species':<35}{'i':>3} {'j':>2}"
    else:
        species_label = "Species"
    header = f"   {species_label:<{name_w}} {'D3(R6)':>{c1_w}} {'D3(R8)':>{c2_w}} {'ABC':>{c3_w}} {total_label:>{c4_w}}"
    rule = "   " + "-" * (name_w + c1_w + c2_w + c3_w + c4_w + 4)

    if options.damp == "bj":
        citation = CITATION_BJ_SHORT
        print(f"   D3(BJ): {citation}")
    else:
        citation = CITATION_ZERO_SHORT
        print(f"   D3(0): {citation}")

    if options.verbose:
        manual_params = (options.s6 is not None and options.s8 is not None and
                         (options.damp == "zero" and options.rs6 is not None or
                          options.damp == "bj" and options.a1 is not None and options.a2 is not None))
        if options.damp == "zero":
            print("\n   D3-dispersion correction with zero-damping")
            if manual_params:
                print("   Manual parameters have been defined")
                print(f"   Zero-damping parameters: s6 = {options.s6}  rs6 = {options.rs6}  s8 = {options.s8}")
            elif dft_functional is not None:
                _, prm = _lookup_functional(dft_functional, zero_parms)
                if prm is not None:
                    s6, rs6, s8 = prm
                    print(f"   Detected {dft_functional} functional - using default zero-damping parameters")
                    print(f"   Zero-damping parameters: s6 = {s6}  rs6 = {rs6}  s8 = {s8}")
        elif options.damp == "bj":
            print("\n   D3-dispersion correction with Becke-Johnson damping")
            if manual_params:
                print("   Manual parameters have been defined")
                print(f"   BJ-damping parameters: s6 = {options.s6}  s8 = {options.s8}  "
                      f"a1 = {options.a1}  a2 = {options.a2}")
            elif dft_functional is not None:
                _, prm = _lookup_functional(dft_functional, bj_parms)
                if prm is not None:
                    s6, a1, s8, a2 = prm
                    print(f"   Detected {dft_functional} functional - using default BJ-damping parameters")
                    print(f"   BJ-damping parameters: s6 = {s6}  s8 = {s8}  a1 = {a1}  a2 = {a2}")
        if options.threebody:
            print("   Including the Axilrod-Teller-Muto repulsive 3-body dispersion term")

    if options.intermolecular == "auto":
        print("\n   Caution: fragments detected automatically from covalent connectivity.")
        print("   Only intermolecular dispersion interactions are included.")

    print()
    print(header)
    print(rule)

    exit_code = 0
    for filepath, data in parsed_files:
        try:
            pw_atoms = set(options.pairwise) if do_pairwise and options.pairwise else None

            result = CalcD3(
                data, dft_functional, options.damp,
                options.s6, options.rs6, options.s8,
                options.a1, options.a2,
                options.threebody, options.intermolecular, do_pairwise,
                options.cutoff,
            )

            unit_factor = 1.0 if options.kcal else 1.0 / AUTOKCAL

            # Print pairwise breakdown (same column layout as summary)
            if do_pairwise and result.pairwise_terms:
                abc_blank = " " * c3_w
                for at1, at2, pw_r6, pw_r8 in result.pairwise_terms:
                    if pw_atoms and not (at1 in pw_atoms and at2 in pw_atoms):
                        continue
                    pw_r6_u = pw_r6 * unit_factor
                    pw_r8_u = pw_r8 * unit_factor
                    pw_total = pw_r6_u + pw_r8_u
                    pw_label = f"{filepath:<38}{at1:>3d} {at2:>2d}"
                    print(f"   {pw_label:<{name_w}} {pw_r6_u:>{c1_w}.{fmt_dp}f} "
                          f"{pw_r8_u:>{c2_w}.{fmt_dp}f} {abc_blank:>{c3_w}} "
                          f"{pw_total:>{c4_w}.{fmt_dp}f}")

            # Print fragment details in verbose mode
            if options.verbose and result.fragments is not None:
                for i, frag in enumerate(result.fragments):
                    print(f"   Fragment {i + 1}: atoms {_format_atom_range(frag)} ({len(frag)} atoms)")

            c6_term = result.attractive_r6_vdw * unit_factor
            c8_term = result.attractive_r8_vdw * unit_factor
            threebody_term = result.repulsive_abc * unit_factor

            total_vdw = c6_term + c8_term
            if options.threebody:
                total_vdw += threebody_term
                abc_str = f"{threebody_term:{c3_w}.{fmt_dp}f}"
            else:
                abc_str = " " * c3_w

            print(f"   {filepath:<{name_w}} {c6_term:>{c1_w}.{fmt_dp}f} "
                  f"{c8_term:>{c2_w}.{fmt_dp}f} {abc_str:>{c3_w}} "
                  f"{total_vdw:>{c4_w}.{fmt_dp}f}")

        except Exception as e:
            logger.error("Error processing %s: %s", filepath, e)
            exit_code = 1

    print(rule)
    print()
    return exit_code

if __name__ == "__main__":
    main()
