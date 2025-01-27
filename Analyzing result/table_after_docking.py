import os
import json
import pandas as pd

def load_config(config_path):
    """Загружает конфигурацию из JSON-файла."""
    with open(config_path, 'r', encoding='utf-8') as file:
        return json.load(file)

def extract_binding_energy_from_dlg(file_path):
    """Извлекает энергию связывания из файла DLG."""
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            if "DOCKED: USER    Estimated Free Energy of Binding" in line:
                try:
                    return float(line.split('=')[1].split()[0])
                except (IndexError, ValueError):
                    return None
    return None

def extract_best_conformation(file_path):
    """Извлекает лучшую конформацию из DLG."""
    best_conf = []
    capturing = False
    found_model = False
    found_atoms = False

    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            if "DOCKED: MODEL" in line and not found_model:
                capturing = True
                found_model = True
                best_conf = []
            elif "DOCKED: ENDMDL" in line and capturing:
                capturing = False
                break
            elif capturing:
                if line.startswith("DOCKED: ATOM") or line.startswith("DOCKED: HETATM"):
                    found_atoms = True
                    best_conf.append(line.replace("DOCKED: ", ""))
    
    if not found_model:
        print(f"Ошибка: MODEL не найден в {file_path}")
    if not found_atoms:
        print(f"Ошибка: ATOM не найден в {file_path}")

    return best_conf

def convert_dlg_to_pdbqt(file_path, output_directory):
    """Создаёт PDBQT-файл из DLG."""
    pdbqt_content = extract_best_conformation(file_path)
    output_file_path = os.path.join(output_directory, os.path.basename(file_path).replace('.dlg', '.pdbqt'))

    if pdbqt_content:
        with open(output_file_path, 'w', encoding='utf-8') as output_file:
            output_file.writelines(pdbqt_content)
        print(f"Успешно создан PDBQT: {output_file_path}")
        return output_file_path
    else:
        print(f"Ошибка: не удалось создать PDBQT для {file_path}")
        return None

def process_single_dlg(file_path, output_directory, results_csv):
    """Обрабатывает один DLG-файл и записывает результат в CSV."""
    file_name = os.path.basename(file_path)
    pdbqt_path = convert_dlg_to_pdbqt(file_path, output_directory)
    binding_energy = extract_binding_energy_from_dlg(file_path)

    errors = []
    if binding_energy is None:
        errors.append("Ошибка: не удалось извлечь энергию связывания")
    if pdbqt_path is None:
        errors.append("Ошибка: не удалось создать PDBQT")

    df = pd.DataFrame([{  
        'File Name': file_name,
        'Binding Energy': binding_energy,
        'Error': "; ".join(errors) if errors else None
    }])

    df.to_csv(results_csv, mode='a', header=not os.path.exists(results_csv), index=False, encoding='utf-8')
    
    return file_name, binding_energy, "; ".join(errors) if errors else None

def main():
    """Основная функция, загружает конфиг и обрабатывает файлы."""
    config_path = "C:/Users/Asus/Desktop/Molecules/config2.json"
    config = load_config(config_path)
    dlg_directory = config.get("input_directory")
    pdbqt_directory = config.get("output_directory")
    results_csv = os.path.join(pdbqt_directory, config.get("output_csv", "binding_energies.csv"))

    if not os.path.exists(dlg_directory):
        print(f"Ошибка: директория {dlg_directory} не найдена")
        return
    if not os.path.exists(pdbqt_directory):
        os.makedirs(pdbqt_directory)

    for file_name in os.listdir(dlg_directory):
        if file_name.endswith('.dlg'):
            file_path = os.path.join(dlg_directory, file_name)
            file_name, energy, error = process_single_dlg(file_path, pdbqt_directory, results_csv)
            print(f"Обработан {file_name}: Энергия = {energy}, Ошибки = {error}")

    print(f"Все результаты сохранены в {results_csv}")

if __name__ == "__main__":
    main()
