"""
Atom classification for RABDAM 3 structure preparation.
"""

from collections.abc import Iterable

from input.reader import AtomRecord
from structure.filters import is_hydrogen
from structure.models import PreparedAtom


AMINO_ACIDS = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    }
)


NUCLEOTIDES = frozenset(
    {
        "A",
        "C",
        "G",
        "U",
        "T",
        "DA",
        "DC",
        "DG",
        "DT",
        "DU",
    }
)


SOLVENTS = frozenset(
    {
        "HOH",
        "WAT",
        "DOD",
    }
)


def classify_atoms(atoms: Iterable[AtomRecord]) -> tuple[PreparedAtom, ...]:
    """
    Classify atom records for RABDAM preparation.
    """

    return tuple(classify_atom(atom) for atom in atoms)


def classify_atom(atom: AtomRecord) -> PreparedAtom:
    """
    Classify one atom record for RABDAM preparation.
    """

    component_name = atom.residue_name.strip().upper()
    record_type = atom.record_type.strip().upper()

    return PreparedAtom(
        record=atom,
        is_hydrogen=is_hydrogen(atom),
        is_protein=component_name in AMINO_ACIDS,
        is_nucleic_acid=component_name in NUCLEOTIDES,
        is_solvent=component_name in SOLVENTS,
        is_hetatm=record_type == "HETATM",
    )
