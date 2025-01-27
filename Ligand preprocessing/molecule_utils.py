import logging
import csv
import os
from typing import Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from modeller import Environ
from modeller.scripts import complete_pdb



def validate_molecule(mol: Optional[Chem.Mol]) -> bool:
    """
    Проверяет валидность молекулы.
    
    Args:
        mol (Optional[Chem.Mol]): Объект молекулы RDKit.
    
    Returns:
        bool: True, если молекула валидна, иначе False.
    """
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        return True
    except Exception as e:
        logging.warning(f"Молекула некорректна: {e}")
        return False

def replace_atoms_with_carbon(mol: Chem.Mol, mol_id: int, pdb_file: str) -> {bool, Chem.Mol}:
    """
    Заменяет атомы кремния (Si) и бора (B) в молекуле на углерод (C) и сохраняет изменения в PDB-файл.
    
    Args:
        mol (Chem.Mol): Объект молекулы RDKit.
        mol_id (int): Идентификатор молекулы.
        pdb_file (str): Путь для сохранения PDB-файла.
    
    Returns:
        True/False
        Chem.Mol: Молекула с заменёнными атомами.
    """
    try:
        replaced_atoms = 0
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ["Si", "B"]:
                logging.info(f"Молекула {mol_id}: замена {atom.GetSymbol()} на C (атом ID {atom.GetIdx()})")
                atom.SetAtomicNum(6)
                replaced_atoms += 1
        Chem.MolToPDBFile(mol, pdb_file)
        return True, mol
    except Exception as e:
        logging.error(f"Ошибка при замене атомов в молекуле: {e}")
        return False, mol

def optim_mol_geometry(mol: Chem.Mol, mol_id: int, pdb_file: str) -> {bool, Chem.Mol}:
    """
    Оптимизирует геометрию молекулы с использованием силовых полей.
    
    Args:
        mol (Chem.Mol): Объект молекулы RDKit.
        mol_id (int): Идентификатор молекулы.
        pdb_file (str): Путь к файлу PDB.
    """
    mol = AllChem.RemoveHs(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    force_fields = {
        "MMFF94": AllChem.MMFFOptimizeMolecule,
        "UFF": AllChem.UFFOptimizeMolecule
    }
    for force_field_name, optimization_function in force_fields.items():
        try:
            AllChem.EmbedMolecule(mol, maxAttempts=100000, randomSeed=42)
            if optimization_function(mol, maxIters=100000) == 0:
                Chem.MolToPDBFile(mol, pdb_file)
                return True, mol
        except Exception as e:
            logging.warning(f"Оптимизация {force_field_name} не удалась: {e}")
            return False, mol
    logging.error(f"Молекула {mol_id}: невозможно оптимизировать геометрию")
    return False, mol

def fill_structure_gaps(mol: Chem.Mol, pdb_file: str) -> {bool, Chem.Mol}:
    """
    Заполняет пропуски в структуре молекулы.
    """
    try:
        Chem.SanitizeMol(mol)
        mol = AllChem.RemoveHs(mol)
        mol = Chem.AddHs(mol, addCoords=True)
        if is_protein(pdb_file):
            return fill_gaps_in_protein(mol, pdb_file)
        return fill_gaps_smallmol_rdkit(mol, pdb_file)
    except Exception as e:
        logging.error(f"Ошибка при заполнении пропусков: {e}")
        return False, mol

def fill_gaps_smallmol_rdkit(mol: Chem.Mol, pdb_file: str) -> {bool, Chem.Mol}:
    """
    Восстанавливает молекулу RDKit:
    - Заменяет неопределенные атомы (`[*]`, `X`, `?`) на `C`.
    - Соединяет разорванные фрагменты (если возможно).
    - Оптимизирует структуру с MMFF.
    - Сохраняет результат в PDB.

    :param mol: RDKit Mol-объект
    :param output_pdb: путь к файлу PDB для сохранения
    :return: True, если успешно, иначе False
    """
    try:
        # Заменяем неопределенные атомы на углерод (C)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:  # Неопределенный атом ([*])
                atom.SetAtomicNum(6)  # C (углерод)

        # Проверяем разорванные фрагменты
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) > 1:
            editable_mol = Chem.EditableMol(frags[0])
            for i in range(1, len(frags)):
                editable_mol.AddBond(0, i, Chem.BondType.SINGLE)

            mol = editable_mol.GetMol()

        # Сохранение в PDB
        Chem.MolToPDBFile(mol, pdb_file)

        return True, mol

    except Exception as e:
        print(f"Ошибка: {e}")
        return False, mol

def fill_gaps_in_protein(mol: Chem.Mol, pdb_file: str) -> {bool, Chem.Mol}:
    """
    Заполняет пропуски и достраивает атомы белка через Modeller.
    """
    try:
        env = Environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        Chem.MolToPDBFile(mol, pdb_file)
        complete_pdb(env, pdb_file)
        return True, mol
    except Exception as e:
        logging.error(f"Ошибка заполнения пропусков в белке: {e}")
        return False, mol

def is_protein(pdb_file: str) -> bool:
    """
    Определяет, содержит ли PDB-файл аминокислотные остатки.
    
    Args:
        pdb_file (str): Путь к файлу PDB.
    
    Returns:
        bool: True, если файл содержит аминокислоты, иначе False.
    """
    standard_aa = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                   "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                   "TYR", "VAL"}
    try:
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() in standard_aa:
                    return True
    except Exception as e:
        logging.error(f"Ошибка при обработке файла {pdb_file}: {e}")
    return False

def log_to_csv(data: list, csv_file: str) -> None:
    """
    Записывает данные в CSV-файл.
    
    Args:
        data (list): Данные для записи в CSV.
        csv_file (str): Путь к файлу CSV.
    """
    file_exists = os.path.exists(csv_file)
    with open(csv_file, mode='a', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        if not file_exists:
            writer.writerow([
                "Molecule ID",  
                "Is Valid", 
                "Filling gaps",
                "Atom replacing",
                "Geometry optimization",
                "Torsion handling",
                "Status",
                ])
        writer.writerow(data)