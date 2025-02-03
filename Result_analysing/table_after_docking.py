import os
import re
import pandas as pd

def extract_binding_energy_from_dlg(file_path):
    """
    Извлекает энергию связывания из файла DLG.
    """
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
            for line in file:
                if "DOCKED: USER    Estimated Free Energy of Binding" in line:
                    try:
                        energy = line.split('=')[1].split()[0]
                        return float(energy)
                    except (IndexError, ValueError):
                        return None
    except FileNotFoundError:
        print(f"Ошибка: файл {file_path} не найден.")
    except Exception as e:
        print(f"Ошибка при обработке {file_path}: {e}")
    return None

def convert_dlg_to_pdbqt(file_path, output_directory):
    """
    Преобразует файл DLG в формат PDBQT, извлекая лучшую конформацию лиганда.
    """
    pdbqt_content = []
    capture = False
    base_name = os.path.basename(file_path)
    root, _ = os.path.splitext(base_name)
    output_file_path = os.path.join(output_directory, f"{root}.pdbqt")
    
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
            for line in file:
                stripped_line = line.strip()
                parts = stripped_line.split()
                # Проверяем начало блока модели
                if len(parts) >= 3 and parts[0] == 'DOCKED:' and parts[1] == 'MODEL' and parts[2] == '1':
                    capture = True
                elif stripped_line == 'DOCKED: ENDMDL':
                    capture = False
                    break
                elif capture:
                    if line.startswith("DOCKED:"):
                        pdbqt_content.append(line.replace("DOCKED: ", ""))
        
        if pdbqt_content:
            with open(output_file_path, 'w', encoding='utf-8') as output_file:
                output_file.writelines(pdbqt_content)
        else:
            print(f"Ошибка: не удалось извлечь PDBQT из {file_path}")
    except FileNotFoundError:
        print(f"Ошибка: файл {file_path} не найден.")
    except Exception as e:
        print(f"Ошибка при обработке {file_path}: {e}")
    
    return output_file_path if pdbqt_content else None

def process_dlg_files(directory, output_directory):
    """
    Обрабатывает все файлы DLG в указанной директории, извлекает данные и сохраняет в PDBQT.
    """
    data = []
    
    if not os.path.exists(directory):
        print(f"Ошибка: директория {directory} не существует.")
        return pd.DataFrame()
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for file_name in os.listdir(directory):
        if file_name.lower().endswith('.dlg'):
            file_path = os.path.join(directory, file_name)
            binding_energy = extract_binding_energy_from_dlg(file_path)
            pdbqt_path = convert_dlg_to_pdbqt(file_path, output_directory)
            data.append({
                'File Name': file_name,
                'Binding Energy': binding_energy,
                'PDBQT Path': pdbqt_path if pdbqt_path else 'Ошибка'
            })
    
    return pd.DataFrame(data)

if __name__ == "__main__":
    dlg_directory = "C:/Users/yhalias/Downloads/docking-project-py/Docking/docking_results"
    pdbqt_directory = "C:/Users/yhalias/Downloads/docking-project-py/Docking/pdqbt_files"

    results_df = process_dlg_files(dlg_directory, pdbqt_directory)

    if not results_df.empty:
        output_csv = "binding_energies_and_pdbqt.csv"
        results_df.to_csv(output_csv, index=False, encoding='utf-8')
        print(f"Результаты сохранены в файл: {output_csv}")
    else:
        print("Ошибка: таблица пуста. Проверьте файлы DLG.")