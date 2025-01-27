import os
import json
import pandas as pd
import subprocess

def load_config(config_path):
    """Загружает конфигурацию из JSON-файла."""
    with open(config_path, 'r', encoding='utf-8') as file:
        return json.load(file)

def run_binana(binana_script, receptor_pdbqt, ligand_pdbqt, output_json, output_csv):
    """Запускает BINANA и анализирует взаимодействия."""
    cmd = [
        "python", binana_script,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--output_json", output_json,
        "--output_csv", output_csv
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"BINANA output for {ligand_pdbqt}:\n{result.stdout}")
        return output_json
    except subprocess.CalledProcessError as e:
        print(f"Ошибка при запуске BINANA для {ligand_pdbqt}:\n{e.stderr}")
        return None

def extract_interactions(json_file):
    """Извлекает все взаимодействия из JSON-файла BINANA."""
    if not os.path.exists(json_file):
        return None  # JSON не был создан
    
    with open(json_file, 'r', encoding='utf-8') as file:
        data = json.load(file)

    return {
        "File Name": os.path.basename(json_file).replace('.json', '.pdbqt'),
        "Hydrogen Bonds": len(data.get("hydrogenBonds", [])),
        "Salt Bridges": len(data.get("saltBridges", [])),
        "pi-pi Stacking": len(data.get("piPiStacking", [])),
        "T-shaped Interactions": len(data.get("tStacking", [])),
        "Halogen Bonds": len(data.get("halogenBonds", [])),
        "Cation-Pi Interactions": len(data.get("cationPiInteractions", [])),
        "Hydrophobic Interactions": len(data.get("hydrophobicContacts", [])),
        "Electrostatic Interactions": len(data.get("electrostaticInteractions", [])),
        "Metal Coordination": len(data.get("metalCoordinations", [])),
        "Close Contacts": len(data.get("closeContacts", []))
    }

def process_binana(binana_script, receptor_pdbqt, ligands_dir, output_directory, output_csv):
    """Обрабатывает все PDBQT-лиганды с BINANA и записывает результаты в CSV."""
    data = []
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for file_name in os.listdir(ligands_dir):
        if file_name.endswith('.pdbqt'):
            ligand_path = os.path.join(ligands_dir, file_name)
            json_path = os.path.join(output_directory, file_name.replace('.pdbqt', '.json'))
            csv_path = os.path.join(output_directory, file_name.replace('.pdbqt', '.csv'))
            
            run_binana(binana_script, receptor_pdbqt, ligand_path, json_path, csv_path)

            if os.path.exists(json_path):
                interaction_data = extract_interactions(json_path)
                if interaction_data:
                    data.append(interaction_data)
                else:
                    data.append({"File Name": file_name, "Error": "BINANA output empty"})
            else:
                data.append({"File Name": file_name, "Error": "BINANA failed"})

    df = pd.DataFrame(data)
    df.to_csv(os.path.join(output_directory, output_csv), index=False, encoding='utf-8')
    print(f"Результаты сохранены в {os.path.join(output_directory, output_csv)}")

def main():
    """Основная функция."""
    config_path = "C:/Users/Asus/Desktop/Molecules/config3.json"
    config = load_config(config_path)
    binana_script = config.get("binana_path")
    receptor_pdbqt = config.get("receptor_pdbqt")
    ligands_dir = config.get("ligands_directory")
    output_directory = config.get("output_directory")
    output_csv = config.get("output_csv")

    if not os.path.exists(binana_script):
        print(f"Ошибка: BINANA не найден по пути {binana_script}")
        return
    if not os.path.exists(receptor_pdbqt):
        print(f"Ошибка: файл рецептора {receptor_pdbqt} не найден")
        return
    if not os.path.exists(ligands_dir):
        print(f"Ошибка: директория с лигандами {ligands_dir} не найдена")
        return

    process_binana(binana_script, receptor_pdbqt, ligands_dir, output_directory, output_csv)
    
if __name__ == "__main__":
    main()
