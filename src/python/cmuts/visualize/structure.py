import os
import subprocess
from typing import Union

import numpy as np
from Bio import PDB, Align


def _seq_from_cif(filename: str, chain_id: Union[str, None] = None) -> str:
    if chain_id is None:
        chain_id = "A"

    parser: Union[PDB.MMCIFParser, PDB.PDBParser]
    if filename.endswith(".cif"):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        parser = PDB.PDBParser(QUIET=True)

    rna_dict = {
        "A": "A",
        "ADE": "A",
        "U": "U",
        "URA": "U",
        "URI": "U",
        "G": "G",
        "GUA": "G",
        "GTP": "G",
        "C": "C",
        "CYT": "C",
        "PSU": "U",
        "5MU": "U",
        "H2U": "U",
        "M2G": "G",
        "M7G": "G",
        "1MA": "A",
        "5MC": "C",
    }

    try:
        structure = parser.get_structure("rna", filename)

        # Get the specified chain
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                sequence = ""

                for residue in chain:
                    resname = residue.get_resname().strip()

                    # Check if it's an RNA nucleotide
                    if resname in rna_dict:
                        sequence += rna_dict[resname]
                    elif PDB.is_aa(residue):
                        # Skip amino acids
                        continue
                    else:
                        # Unknown nucleotide, use 'N'
                        sequence += "N"

                return sequence.replace("U", "T")

        print(f"Chain {chain_id} not found in structure")
        print(f"The valid chains are {[chain.id for chain in model]}")
        return ""

    except Exception as e:
        print(f"Error parsing file: {e}")
        return ""


def _seq_align(seq1: str, seq2: str) -> tuple[str, str]:
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]

    return best_alignment[0], best_alignment[1]


def _data_aln(data: np.ndarray, aln1: str, aln2: str) -> np.ndarray:
    """
    Map data values to align with seq2, using zeros for gaps.
    """

    chars = len([c for c in aln1 if c != "-"])
    if data.shape[0] != chars:
        raise ValueError(
            f"Data length ({data.shape[0]}) must match number of non-gap characters in aln1 ({chars})"
        )

    mapped_data: list[float] = []
    data_idx = 0

    for _, (char1, char2) in enumerate(zip(aln1, aln2)):
        if char2 == "-":
            # Gap in seq2, skip this position entirely
            if char1 != "-":
                data_idx += 1  # Still advance data index for seq1
            continue

        if char1 == "-":
            # Gap in seq1, insert zero for seq2
            mapped_data.append(0.0)
        else:
            # Both sequences have characters, use data value
            mapped_data.append(float(data[data_idx]))
            data_idx += 1

    result: np.ndarray = np.array(mapped_data)
    return result


def _to_defattr(
    values: np.ndarray,
    out: str,
    start: int = 1,
    attr_name: str = "value",
    pad5: int = 0,
    pad3: int = 0,
    chain: Union[str, None] = None,
):
    """
    Convert numpy array to ChimeraX defattr format.
    """

    values = np.concatenate([np.zeros(pad5), values, np.zeros(pad3)])

    with open(out, "w") as f:
        f.write(f"attribute: {attr_name}\n")
        f.write("recipient: residues\n\n")

        for i, value in enumerate(values):
            resnum = start + i
            if chain is not None:
                f.write(f"\t/{chain}:{resnum}\t{value}\n")
            else:
                f.write(f"\t:{resnum}\t{value}\n")


def _to_defattr_atom(
    values: np.ndarray,
    out: str,
    atoms: list[str],
    sizes: np.ndarray,
    start: int = 1,
    attr_name: str = "value",
    pad5: int = 0,
    pad3: int = 0,
    chain: Union[str, None] = None,
):
    """
    Convert a numpy array to ChimeraX defattr format, with one value per atom.
    """

    values = np.concatenate([np.zeros(pad5), values, np.zeros(pad3)])
    if len(atoms) != values.shape[0]:
        raise ValueError("The number of atoms must match the number of reactivity values.")
    if len(atoms) != sizes.sum():
        raise ValueError("The number of atoms must match the sum of all residue sizes.")

    ix = 0
    with open(out, "w") as f:
        f.write(f"attribute: {attr_name}\n")
        f.write("recipient: atoms\n\n")

        for jx in range(sizes.shape[0]):
            for _ in range(sizes[jx]):
                resnum = start + jx
                value = values[ix]
                atom = atoms[ix].replace("p", "'")
                if chain is not None:
                    f.write(f"\t/{chain}:{resnum}@{atom}\t{value}\n")
                else:
                    f.write(f"\t:{resnum}@{atom}\t{value}\n")
                ix += 1


def chimerax_command(
    cif: str,
    defattr: str,
    color: str = "indianred",
    max_value: float = 1.0,
    chain: Union[str, None] = None,
) -> str:
    """Build the ChimeraX command string that visualizes ``defattr`` on ``cif``.

    Use the basename of ``cif`` and ``defattr`` (rather than absolute paths)
    when generating a command for users to run locally — they will typically
    keep both files in the same directory.
    """
    return (
        f"open {cif}; "
        + "close #1.2-999; "
        + ("" if chain is None else f"del ~/{chain}; ")
        + "hide pseudobonds; "
        "color grey; "
        + "graphics quality 5; "
        + "renumber start 1 relative false; "
        + f"open {defattr}; "
        + f"color byattribute value palette white:{color} range 0,{max_value}; "
        + "hide cartoons; "
        + "nucleotides atoms; "
        + "style sphere; "
        + "lighting soft; "
        + "lighting ambientIntensity 1.3"
    )


def _color_by_defattr(
    bin: str,
    cif: str,
    defattr: str,
    color: str,
    max: float = 1,
    chain: Union[str, None] = None,
) -> None:
    """Color a .cif file with the given .defattr file by launching ChimeraX."""
    chmx_cmd = chimerax_command(cif, defattr, color, max, chain)
    subprocess.run([bin, "--cmd", chmx_cmd])


def _color_by_reactivity(
    bin: str,
    cif: str,
    reactivity: np.ndarray,
    chain: Union[str, None],
    color: str,
) -> None:
    """
    Color a .cif file with the given reactivity.
    """

    defattr = ".cmuts-visualize.defattr"
    _to_defattr(reactivity, defattr, chain=chain)
    _color_by_defattr(bin, cif, defattr, color, reactivity.max(), chain)


def _color_atoms_by_value(
    bin: str,
    cif: str,
    reactivity: np.ndarray,
    atoms: list[str],
    sizes: np.ndarray,
    chain: Union[str, None],
    color: str,
) -> None:
    """
    Color a .cif file with the given reactivity.
    """

    defattr = ".cmuts-visualize.defattr"
    _to_defattr_atom(reactivity, defattr, atoms, sizes, chain=chain)
    _color_by_defattr(bin, cif, defattr, color)


def make_defattr(
    reactivity: np.ndarray,
    seq: str,
    cif: str,
    defattr_path: str,
    chain: Union[str, None] = None,
) -> float:
    """Align ``reactivity`` against the ``cif`` sequence and write a defattr file.

    Args:
        reactivity: Per-position reactivity values for ``seq``.
        seq: Reference sequence the reactivities correspond to (T/U interchangeable).
        cif: Path to the structure file used to extract the chain sequence for alignment.
        defattr_path: Output path for the .defattr file.
        chain: Chain id to align against (default: A).

    Returns:
        The maximum (post-alignment, NaN-zeroed) reactivity value, suitable for
        passing as ``max_value`` to :func:`chimerax_command`.
    """
    if not os.path.exists(cif):
        raise FileNotFoundError(f"Structure file not found: {cif}")
    if len(reactivity) != len(seq):
        raise ValueError(
            f"Reactivity length ({len(reactivity)}) must match sequence length ({len(seq)})"
        )

    cif_seq = _seq_from_cif(cif, chain)
    aln1, aln2 = _seq_align(seq, cif_seq)
    aln_data = _data_aln(reactivity, aln1, aln2)
    aln_data = np.nan_to_num(aln_data, nan=0.0)
    _to_defattr(aln_data, defattr_path, chain=chain)
    return float(aln_data.max()) if aln_data.size else 1.0


def visualize_structure(
    reactivity: np.ndarray,
    seq: str,
    cif: str,
    color: str = "indianred",
    chain: Union[str, None] = None,
    bin: str = "ChimeraX",
    defattr_out: Union[str, None] = None,
    run_chimerax: bool = True,
) -> str:
    """Visualize reactivity data on a 3D molecular structure.

    Writes a defattr file and (by default) launches ChimeraX to display the
    coloured structure. Set ``run_chimerax=False`` to write only the defattr
    and return the ChimeraX command — useful for server-side use where the
    user runs ChimeraX on their own machine.

    Args:
        defattr_out: Path to write the defattr file. Defaults to a hidden file
            in the working directory.
        run_chimerax: When True (default), launch ChimeraX to render the result.
            When False, the defattr is written and the ChimeraX command is
            returned without launching anything.

    Returns:
        The ChimeraX command string used to (or that would) render the structure.
    """
    defattr = defattr_out or ".cmuts-visualize.defattr"
    max_value = make_defattr(reactivity, seq, cif, defattr, chain)
    cmd = chimerax_command(cif, defattr, color, max_value, chain)
    if run_chimerax:
        subprocess.run([bin, "--cmd", cmd])
    return cmd


def visualize_structure_atoms(
    reactivity: np.ndarray,
    atoms: list[str],
    sizes: np.ndarray,
    seq: str,
    cif: str,
    color: str = "indianred",
    chain: Union[str, None] = None,
    bin: str = "ChimeraX",
) -> None:
    """Visualize reactivity data on specific atoms of a 3D structure."""
    if not os.path.exists(cif):
        raise FileNotFoundError(f"Structure file not found: {cif}")

    cif_seq = _seq_from_cif(cif)
    aln1, aln2 = _seq_align(seq, cif_seq)
    sizes = _data_aln(sizes, aln1, aln2)

    _color_atoms_by_value(bin, cif, reactivity, atoms, sizes, chain, color)


def main():
    """CLI entry point for cmuts-visualize."""
    import argparse

    import h5py

    from cmuts.internal import ProbingData

    parser = argparse.ArgumentParser(
        prog="cmuts-visualize", description="Visualize reactivity on a 3D structure"
    )
    parser.add_argument("file", help="HDF5 file with reactivity data")
    parser.add_argument("cif", help="CIF/PDB structure file")
    parser.add_argument("--group", default="", help="Group name in HDF5 file")
    parser.add_argument("--index", type=int, default=0, help="Sequence index to visualize")
    parser.add_argument("--chain", help="Chain ID (default: A)")
    parser.add_argument("--color", default="indianred", help="Color for high reactivity")
    parser.add_argument("--chimerax", default="ChimeraX", help="Path to ChimeraX executable")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        raise FileNotFoundError(f"File not found: {args.file}")

    with h5py.File(args.file, "r") as f:
        data = ProbingData.load(args.group, f)

    # Get reactivity for specified index
    if data.reactivity.ndim > 1:
        reactivity = data.reactivity[args.index]
    else:
        reactivity = data.reactivity

    # Get sequence if available
    if data.sequences is not None and len(data.sequences) > args.index:
        # Decode sequence tokens (A=0, C=1, G=2, U=3)
        token_map = {0: "A", 1: "C", 2: "G", 3: "U"}
        seq = "".join(token_map.get(int(t), "N") for t in data.sequences[args.index])
    else:
        # Use placeholder sequence matching reactivity length
        seq = "N" * len(reactivity)

    visualize_structure(reactivity, seq, args.cif, args.color, args.chain, args.chimerax)


if __name__ == "__main__":
    main()
