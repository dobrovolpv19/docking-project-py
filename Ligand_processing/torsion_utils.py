# torsion_utils.py
from pathlib import Path
import logging
import csv


class PDBQTTorsionProcessor:
    """
    A class to process PDBQT files for torsions, active torsion counts,
    and warnings about torsion limits.
    """

    @staticmethod
    def inactivate_CX_torsions(pdbqt_file: Path) -> None:
        """
        Checks torsional bonds in the PDBQT file and inactivates non-C-C bonds
        only if a warning about exceeded torsions exists. Changes are made
        directly to the original file.

        Args:
            pdbqt_file (Path): Path to the PDBQT file.
        """
        warning_found = False
        with pdbqt_file.open("r") as infile:
            lines = infile.readlines()

        for line in lines:
            if "REMARK WARNING: 32 MAX_TORS EXCEEDED!!!" in line:
                warning_found = True
                break

        if not warning_found:
            logging.info(f"No torsion warning found in {pdbqt_file}. No action taken.")
            return

        with pdbqt_file.open("w") as outfile:
            for line in lines:
                if line.startswith("REMARK") and "between atoms:" in line and " A " in line:
                    parts = line.split()
                    atom1, atom2 = parts[-3], parts[-1]
                    element1 = atom1.split("_")[0][0]
                    element2 = atom2.split("_")[0][0]
                    if not (element1 == "C" and element2 == "C"):
                        line = line.replace(" A ", " I ")
                outfile.write(line)

    @staticmethod
    def update_active_torsions(pdbqt_file: Path) -> None:
        """
        Updates the count of active torsions in the PDBQT file and removes
        warnings if active torsions are fewer than 32.

        Args:
            pdbqt_file (Path): Path to the PDBQT file.
        """
        warning_found = False
        with pdbqt_file.open("r") as infile:
            lines = infile.readlines()

        for line in lines:
            if "REMARK WARNING: 32 MAX_TORS EXCEEDED!!!" in line:
                warning_found = True
                break

        if not warning_found:
            logging.info(f"No torsion warning found in {pdbqt_file}. No action taken.")
            return

        active_torsions_count = 0
        updated_lines = []
        warning_removed = False

        for line in lines:
            if line.startswith("REMARK") and "between atoms:" in line and " A " in line:
                active_torsions_count += 1

        for line in lines:
            if line.startswith("REMARK WARNING: 32 MAX_TORS EXCEEDED!!!"):
                if active_torsions_count < 32:
                    warning_removed = True
                    continue

            if line.startswith("REMARK") and "active torsions:" in line:
                line = f"REMARK  {active_torsions_count} active torsions:\n"

            if line.startswith("TORSDOF"):
                line = f"TORSDOF {active_torsions_count}"

            updated_lines.append(line)

        with pdbqt_file.open("w") as outfile:
            outfile.writelines(updated_lines)

        logging.info(f"File {pdbqt_file} updated with {active_torsions_count} active torsions.")
        if warning_removed:
            logging.info("Warning removed as active torsions are fewer than 32.")

    @staticmethod
    def check_torsion_warning(pdbqt_file: Path) -> bool:
        """
        Checks the PDBQT file for the presence of a torsion warning.

        Args:
            pdbqt_file (Path): Path to the PDBQT file.

        Returns:
            bool: True if the torsion warning is found, otherwise False.
        """
        with pdbqt_file.open("r") as file:
            for line in file:
                if "REMARK WARNING: 32 MAX_TORS EXCEEDED!!!" in line:
                    logging.warning(f"Torsion warning found in file: {pdbqt_file}")
                    return True
        return False


def process_torsions(pdbqt_file: Path, mol_id: int) -> bool:
    """
    Processes torsions for a given PDBQT file.

    Args:
        pdbqt_file (Path): Path to the PDBQT file.
        mol_id (int): Unique identifier for the molecule.
    """
    try:
        PDBQTTorsionProcessor.inactivate_CX_torsions(pdbqt_file)
        PDBQTTorsionProcessor.update_active_torsions(pdbqt_file)
        PDBQTTorsionProcessor.check_torsion_warning(pdbqt_file)
        return True
    except Exception as e:
        logging.error(f"Error processing torsions for {pdbqt_file}: {e}")
        return False
