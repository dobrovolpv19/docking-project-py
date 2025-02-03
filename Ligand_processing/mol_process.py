import logging
import subprocess
import json
import os
from typing import Tuple, Optional
from rdkit import Chem
from multiprocessing import Pool, TimeoutError
from molecule_utils import (
    validate_molecule, optim_mol_geometry, replace_atoms_with_carbon, 
    fill_structure_gaps, log_to_csv
)
from torsion_utils import PDBQTTorsionProcessor
from pathlib import Path
import sys
from pathlib import Path

# Путь к текущей директории, где находится mol_process.py
BASE_DIR = Path(__file__).parent

# Путь к конфигу относительно скрипта
CONFIG_PATH = BASE_DIR / "config.json"

def get_path(path_str: str) -> Path:
    """Преобразует строку пути в абсолютный путь."""
    path = Path(path_str)
    return path if path.is_absolute() else BASE_DIR.parent / path

# Читаем конфиг
import json
with open(CONFIG_PATH, "r", encoding="utf-8") as f:
    config = json.load(f)
def process_molecule(mol_data: Tuple[Optional[Chem.Mol], int], config: dict) -> None:
    """
    Обрабатывает одну молекулу, вызывая функции из molecule_utils и torsion_utils.
    """
    mol, mol_id = mol_data
    log_data = [mol_id]

    if not validate_molecule(mol):
        logging.warning(f"Молекула {mol_id} некорректна")
        log_data.append("Failed")
        log_to_csv(log_data, get_path(config["log_csv_file"]))
        return
    else:
        log_data.append("Success")
        
    try:
        original_dir = os.getcwd()
        os.chdir(get_path(config["output_folder"]))

        mol_name = f"ligand_{mol_id}"
        pdb_file = f"{get_path(config['output_folder'])}/{mol_name}.pdb"
        pdbqt_file = f"{get_path(config['output_folder'])}/{mol_name}.pdbqt"

        results = []
        Chem.MolToPDBFile(mol, pdb_file)

        success, mol = fill_structure_gaps(mol, pdb_file)
        results.append(success)

        success, mol = optim_mol_geometry(mol, mol_id, pdb_file)
        results.append(success)

        success, mol = replace_atoms_with_carbon(mol, mol_id, pdb_file)
        results.append(success)

        log_data.extend(["Success" if res else "Failed" for res in results])
        Chem.MolToPDBFile(mol, pdb_file)

        cmd = [
            config["pythonsh_path"],
            config["prepare_ligand_script"],
            "-l", pdb_file,
            "-o", pdbqt_file,
            "-A", "hydrogens",
            "-U", "nphs"
        ]
        subprocess.run(cmd, check=True)

        pdbqt_path = Path(pdbqt_file)
        if PDBQTTorsionProcessor.check_torsion_warning(pdbqt_path):
            logging.info(f"Молекула {mol_id} превышает лимит торсионов, повторная генерация (-Z)")
            cmd_z = [
                config["pythonsh_path"],
                config["prepare_ligand_script"],
                "-l", pdb_file,
                "-o", pdbqt_file,
                "-A", "hydrogens",
                "-U", "nphs",
                "-Z"
            ]
            subprocess.run(cmd_z, check=True)

        log_data.append("Success")
        log_to_csv(log_data, get_path(config["log_csv_file"]))

    except subprocess.CalledProcessError as e:
        logging.error(f"Ошибка MGLTools для молекулы {mol_id}: {e}")
        log_to_csv([mol_id, "Failed", f"MGLTools error: {e}"], get_path(config["log_csv_file"]))
    except Exception as e:
        logging.error(f"Общая ошибка для молекулы {mol_id}: {e}")
        log_to_csv([mol_id, "Failed", str(e)], get_path(config["log_csv_file"]))
    finally:
        os.chdir(original_dir)

def process_molecule_with_timeout(data, config, timeout=600):
    """Обёртка для выполнения process_molecule с таймаутом."""
    try:
        process_molecule(data, config)
    except Exception as e:
        logging.error(f"Ошибка обработки молекулы {data[1]}: {e}")

def main() -> None:
    """Главная функция программы."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    
    config_path = sys.argv[1]

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            config = json.load(f)

        logging.info(f"Чтение входного SDF файла: {get_path(config['input_sdf'])}")
        supplier = Chem.SDMolSupplier(get_path(config["input_sdf"]))
        mol_data = [
            (mol, i + config.get("start_from", 1))
            for i, mol in enumerate(supplier) 
            if (mol is not None and i + 1 >= config.get("start_from", 1))
        ]

        num_processes = config.get("num_process", 4)
        with Pool(num_processes) as pool:
            results = [pool.apply_async(process_molecule_with_timeout, (data, config, 600)) for data in mol_data]
            
            for res in results:
                try:
                    res.get(timeout=600)
                except TimeoutError:
                    logging.warning(f"Превышено время обработки молекулы")
                except Exception as e:
                    logging.error(f"Ошибка обработки: {e}")

        logging.info("Обработка завершена")
    except Exception as e:
        logging.critical(f"Критическая ошибка: {e}")

if __name__ == "__main__":
    main()
