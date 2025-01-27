import logging
import subprocess
import json
import os
from typing import Tuple, Optional
from rdkit import Chem
from multiprocessing import Pool, TimeoutError
from molecule_utils import (validate_molecule, optim_mol_geometry, replace_atoms_with_carbon, fill_structure_gaps, log_to_csv)
from torsion_utils import PDBQTTorsionProcessor
from pathlib import Path

def process_molecule(mol_data: Tuple[Optional[Chem.Mol], int], config: dict) -> None:
    """
    Обрабатывает одну молекулу, вызывая функции из molecule_utils и torsion_utils.

    Args:
        mol_data (Tuple[Optional[Chem.Mol], int]): Кортеж с молекулой RDKit и её ID.
        config (dict): Словарь с настройками из config.json.
    """
    mol, mol_id = mol_data
    log_data = [mol_id]

    if not validate_molecule(mol):
        logging.warning(f"Молекула {mol_id} некорректна")
        log_data.append("Failed")
        log_to_csv(log_data, config["log_csv_file"])
        return
    else:
        log_data.append("Success")
        

    try:
        original_dir = os.getcwd()
        os.chdir(config["output_folder"])

        mol_name = f"ligand_{mol_id}"
        pdb_file = f"{config['output_folder']}/{mol_name}.pdb"
        pdbqt_file = f"{config['output_folder']}/{mol_name}.pdbqt"

        results = []
        Chem.MolToPDBFile(mol, pdb_file)

        # Заполнение пропусков в структуре
        success, mol = fill_structure_gaps(mol, pdb_file)
        results.append(success)

        # Оптимизация геометрии молекулы
        success, mol = optim_mol_geometry(mol, mol_id, pdb_file)
        results.append(success)

        # Замена атомов кремния и бора на углерод
        success, mol = replace_atoms_with_carbon(mol, mol_id, pdb_file)
        results.append(success)

        # Запись промежуточных результатов
        log_data.extend(["Success" if res else "Failed" for res in results])
        Chem.MolToPDBFile(mol, pdb_file)

        # Запуск MGLTools для подготовки первичного PDBQT
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

        # Проверка превышения лимита торсионов
        if PDBQTTorsionProcessor.check_torsion_warning(pdbqt_path):
            logging.info(f"Молекула {mol_id} превышает лимит торсионов, повторная генерация с отколючением торсионов (-Z)")

            # Повторное создание PDBQT с отключением торсионов (-Z)
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
        log_to_csv(log_data, config["log_csv_file"])

    except subprocess.CalledProcessError as e:
        logging.error(f"Ошибка при выполнении MGLTools для молекулы {mol_id}: {e}")
        log_to_csv([mol_id, "Failed", f"MGLTools error: {e}"], config["log_csv_file"])
    except Exception as e:
        logging.error(f"Общая ошибка для молекулы {mol_id}: {e}")
        log_to_csv([mol_id, "Failed", str(e)], config["log_csv_file"])
    finally:
        os.chdir(original_dir)

def process_molecule_with_timeout(data, config, timeout=600):
    """
    Обёртка для выполнения process_molecule с таймаутом.
    """
    try:
        process_molecule(data, config)
    except Exception as e:
        logging.error(f"Ошибка обработки молекулы {data[1]}: {e}")

def main() -> None:
    """
    Главная функция программы. Читает входной файл, обрабатывает молекулы и записывает результаты.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    try:
        with open("config.json", "r", encoding="utf-8") as f:
            config = json.load(f)

        logging.info(f"Чтение входного SDF файла: {config['input_sdf']}")
        supplier = Chem.SDMolSupplier(config["input_sdf"])
        mol_data = [
            (mol, i + config.get("start_from", 1))
            for i, mol in enumerate(supplier) 
                if (mol is not None and i + 1 >= config.get("start_from", 1))
        ]

        def process_with_limit(data):
                    num_processes = config.get("num_process")
                    with Pool(num_processes) as pool:
                        async_result = pool.apply_async(process_molecule, (data, config))
                        try:
                            async_result.get(timeout=600)  # 600 секунд = 10 минут
                        except TimeoutError:
                            logging.warning(f"Время обработки молекулы {data[1]} истекло.")
                        except Exception as e:
                            logging.error(f"Ошибка обработки молекулы {data[1]}: {e}")

        for data in mol_data:
            process_with_limit(data)

        logging.info("Обработка завершена")
    except Exception as e:
        logging.critical(f"Критическая ошибка: {e}")

if __name__ == "__main__":
    main()

