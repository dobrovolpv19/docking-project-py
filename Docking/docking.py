import os
import subprocess
from pathlib import Path
import logging
from datetime import datetime
import sys
import platform
import threading
import queue
import time
import json
from multiprocessing import Pool, cpu_count
from functools import partial
import shutil

from pathlib import Path

# Путь к текущей директории, где находится mol_process.py
BASE_DIR = Path(__file__).parent

# Путь к конфигу относительно скрипта
CONFIG_PATH = BASE_DIR / "config.json"

def get_path(path_str: str) -> Path:
    """Преобразует строку пути в абсолютный путь."""
    path = Path(path_str)
    return path if path.is_absolute() else BASE_DIR.parent / path


def prepare_grid_parameter_file(receptor_path, grid_center, grid_size, spacing=0.375):
    """
    Создание GPF файла для AutoGrid
    Args:
        receptor_path: путь к файлу рецептора
        grid_center: центр грида (x, y, z)
        grid_size: размер грида в точках (nx, ny, nz)
        spacing: расстояние между точками грида
    """
    receptor_name = os.path.splitext(os.path.basename(receptor_path))[0]
    gpf_content = f"""npts {grid_size[0]} {grid_size[1]} {grid_size[2]}
gridfld {receptor_name}.maps.fld
spacing {spacing}
receptor_types A C HD N OA SA NA S
ligand_types A C HD N OA SA NA S Cl F Br P I
receptor {receptor_path}
gridcenter {grid_center[0]} {grid_center[1]} {grid_center[2]}
smooth 0.5
map {receptor_name}.A.map
map {receptor_name}.C.map
map {receptor_name}.HD.map
map {receptor_name}.N.map
map {receptor_name}.NA.map
map {receptor_name}.OA.map
map {receptor_name}.SA.map
map {receptor_name}.S.map
map {receptor_name}.Cl.map
map {receptor_name}.F.map
map {receptor_name}.Br.map
map {receptor_name}.P.map
map {receptor_name}.I.map
elecmap {receptor_name}.e.map
dsolvmap {receptor_name}.d.map
"""
    with open(f"{receptor_name}.gpf", 'w', encoding='utf-8') as gpf_file:
        gpf_file.write(gpf_content)
    return f"{receptor_name}.gpf"

class StreamReader(threading.Thread):
    """Класс для асинхронного чтения потоков вывода процесса"""
    def __init__(self, stream, queue, stream_name):
        threading.Thread.__init__(self, daemon=True)
        self.stream = stream
        self.queue = queue
        self.stream_name = stream_name

    def run(self):
        while True:
            line = self.stream.readline()
            if not line:
                break
            self.queue.put((self.stream_name, line.strip()))
        self.stream.close()


def run_autogrid(gpf_file, timeout=3600):
    """
    Запуск AutoGrid с асинхронным чтением вывода
    Args:
        gpf_file: путь к GPF файлу
        timeout: таймаут в секундах
    """
    logging.info(f"Запуск AutoGrid для файла: {gpf_file}")
    
    if not os.path.exists(gpf_file):
        logging.error(f"GPF файл не найден: {gpf_file}")
        return False

    try:
        # Создаём очередь для сбора вывода
        output_queue = queue.Queue()
        
        # Запускаем процесс
        process = subprocess.Popen(
            ['autogrid4', '-p', gpf_file, '-l', f'{gpf_file}.glg'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8',
            bufsize=1,
            universal_newlines=True
        )
        
        logging.info("Процесс AutoGrid запущен")

        # Создаем и запускаем потоки для чтения stdout и stderr
        stdout_reader = StreamReader(process.stdout, output_queue, "stdout")
        stderr_reader = StreamReader(process.stderr, output_queue, "stderr")
        
        stdout_reader.start()
        stderr_reader.start()

        # Устанавливаем время начала
        start_time = time.time()
        
        # Основной цикл мониторинга
        while True:
            # Проверяем таймаут
            if time.time() - start_time > timeout:
                process.kill()
                logging.error(f"Превышено время ожидания ({timeout} секунд)")
                return False

            # Проверяем, завершился ли процесс
            return_code = process.poll()
            
            # Читаем доступный вывод
            try:
                while True:  # Читаем все доступные сообщения
                    stream_name, line = output_queue.get_nowait()
                    if stream_name == "stdout":
                        logging.info(f"AutoGrid: {line}")
                    else:
                        logging.error(f"AutoGrid ошибка: {line}")
            except queue.Empty:
                pass  # Нет больше сообщений в очереди

            # Если процесс завершился
            if return_code is not None:
                # Даем время на получение оставшегося вывода
                time.sleep(0.1)
                
                # Читаем оставшиеся сообщения
                try:
                    while True:
                        stream_name, line = output_queue.get_nowait()
                        if stream_name == "stdout":
                            logging.info(f"AutoGrid: {line}")
                        else:
                            logging.error(f"AutoGrid ошибка: {line}")
                except queue.Empty:
                    pass

                if return_code == 0:
                    logging.info("AutoGrid успешно завершил работу")
                    return True
                else:
                    logging.error(f"AutoGrid завершился с ошибкой (код {return_code})")
                    return False
            
            # Небольшая пауза перед следующей итерацией
            time.sleep(0.1)

    except FileNotFoundError:
        logging.error("Исполняемый файл autogrid4 не найден. Проверьте установку и PATH")
        return False
    except PermissionError:
        logging.error("Отказано в доступе при запуске autogrid4")
        return False
    except Exception as e:
        logging.error(f"Неожиданная ошибка при выполнении AutoGrid: {str(e)}")
        if 'process' in locals() and process.poll() is None:
            process.kill()
        return False
    finally:
        # Убеждаемся, что процесс завершен
        if 'process' in locals() and process.poll() is None:
            process.kill()

def check_grid_files(receptor_name):
    """
    Проверка создания необходимых файлов после выполнения AutoGrid
    Args:
        receptor_name: имя рецептора (без расширения)
    """
    required_files = [
        f"{receptor_name}.maps.fld",
        f"{receptor_name}.A.map",
        f"{receptor_name}.C.map",
        f"{receptor_name}.HD.map",
        f"{receptor_name}.N.map",
        f"{receptor_name}.NA.map",
        f"{receptor_name}.N.map",
        f"{receptor_name}.OA.map",
        f"{receptor_name}.SA.map",
        f"{receptor_name}.e.map",
        f"{receptor_name}.d.map"
    ]
    
    missing_files = []
    for file in required_files:
        if not os.path.exists(file):
            missing_files.append(file)
    
    if missing_files:
        logging.error(f"Отсутствуют следующие файлы после выполнения AutoGrid: {', '.join(missing_files)}")
        return False
    
    logging.info("Все необходимые файлы grid map успешно созданы")
    return True

def check_autodock():
    """Проверка наличия установленного AutoDock и AutoGrid"""
    autodock_found = False
    autogrid_found = False
    
    # Проверяем разные возможные имена исполняемых файлов
    possible_names = {
        'autodock4': ['autodock4', 'autodock4.exe'],
        'autogrid4': ['autogrid4', 'autogrid4.exe']
    }
    
    for program, names in possible_names.items():
        for name in names:
            # Проверяем наличие в PATH
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, name)
                if os.path.isfile(exe_file):
                    if program == 'autodock4':
                        autodock_found = True
                    else:
                        autogrid_found = True
                    break
    
    if not autodock_found:
        logging.error("AutoDock4 не найден в системе. Убедитесь, что он установлен и добавлен в PATH")
    if not autogrid_found:
        logging.error("AutoGrid4 не найден в системе. Убедитесь, что он установлен и добавлен в PATH")
    
    return autodock_found and autogrid_found

def prepare_docking_parameter_file(receptor_name, ligand_path, ga_runs=10):
    """
    Create DPF file with optimized parameters to prevent energy calculation issues
    """
    ligand_name = os.path.splitext(os.path.basename(ligand_path))[0]
    dpf_content = f"""autodock_parameter_version 4.2
outlev 1
set_psw1 none
intelec                                      # calculate internal electrostatics
seed pid time                                # seeds for random generator
ligand_types A C HD N OA SA NA S Cl F Br P I             # atoms types in ligand
fld {receptor_name}.maps.fld                 # grid_data_file
map {receptor_name}.A.map
map {receptor_name}.C.map
map {receptor_name}.HD.map
map {receptor_name}.N.map
map {receptor_name}.NA.map
map {receptor_name}.OA.map
map {receptor_name}.SA.map
map {receptor_name}.S.map
map {receptor_name}.Cl.map
map {receptor_name}.F.map
map {receptor_name}.Br.map
map {receptor_name}.P.map
map {receptor_name}.I.map
elecmap {receptor_name}.e.map                # electrostatics map
desolvmap {receptor_name}.d.map              # desolvation map
move {ligand_path}                           # small molecule
about 0.0 0.0 0.0                           # small molecule center
tran0 random                                 # initial coordinates/random
quaternion0 random                           # initial orientation/random
dihe0 random                                 # initial dihedrals/random
torsdof 4                                    # torsional degrees of freedom
rmstol 2.0                                   # cluster_tolerance/A
extnrg 50.0                               # external grid energy
e0max 0.0 10000                             # max initial energy; max number of retries
ga_pop_size 50                             # number of individuals in population
ga_num_evals 500000                        # maximum number of energy evaluations
ga_num_generations 5000                     # maximum number of generations
ga_elitism 1                                # number of top individuals to survive to next generation
ga_mutation_rate 0.02                       # rate of gene mutation
ga_crossover_rate 0.8                       # rate of crossover
ga_window_size 10                           # 
ga_cauchy_alpha 0.0                         # Alpha parameter of Cauchy distribution
ga_cauchy_beta 1.0                          # Beta parameter Cauchy distribution
set_ga                                      # set the above parameters for GA or LGA
sw_max_its 100                              # iterations of Solis & Wets local search
sw_max_succ 4                               # consecutive successes before changing rho
sw_max_fail 4                               # consecutive failures before changing rho
sw_rho 1.0                                  # size of local search space to sample
sw_lb_rho 0.1                             # lower bound on rho
ls_search_freq 0.06                         # probability of performing local search on individual
unbound_model bound                         # state of unbound ligand
ga_run {ga_runs}                           # do this many hybrid GA-LS runs
analysis                                    # perform a ranked cluster analysis
"""
    dpf_file = f"{ligand_name}_{receptor_name}.dpf"
    with open(dpf_file, 'w', encoding='utf-8') as f:
        f.write(dpf_content)
    return dpf_file

def run_autodock(dpf_file, timeout=7200):
    """Run AutoDock with proper monitoring"""
    try:
        output_queue = queue.Queue()
        process = subprocess.Popen(
            ['autodock4', '-p', dpf_file, '-l', f'{dpf_file}.dlg'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8',
            bufsize=1,
            universal_newlines=True
        )
        
        # Setup async output readers
        stdout_reader = StreamReader(process.stdout, output_queue, "stdout")
        stderr_reader = StreamReader(process.stderr, output_queue, "stderr")
        stdout_reader.start()
        stderr_reader.start()

        start_time = time.time()
        error_count = 0
        while True:
            if time.time() - start_time > timeout:
                process.kill()
                logging.error(f"Timeout exceeded ({timeout} seconds)")
                return False

            return_code = process.poll()
            
            try:
                while True:
                    stream_name, line = output_queue.get_nowait()
                    if "All energies are equal in population" in line:
                        error_count += 1
                        if error_count > 5:
                            process.kill()
                            logging.error("Too many energy calculation errors")
                            return False
                    logging.info(f"AutoDock: {line}" if stream_name == "stdout" else f"AutoDock Error: {line}")
            except queue.Empty:
                pass

            if return_code is not None:
                # Check final DLG file
                if os.path.exists(f"{dpf_file}.dlg"):
                    with open(f"{dpf_file}.dlg", 'r', encoding='utf-8') as f:
                        content = f.read()
                        if "Estimated Free Energy of Binding    =   +0.00 kcal/mol" in content:
                            logging.error("Invalid binding energy calculation")
                            return False
                
                return return_code == 0
            
            time.sleep(0.1)

    except Exception as e:
        logging.error(f"Error running AutoDock: {str(e)}")
        if 'process' in locals() and process.poll() is None:
            process.kill()
        return False

def setup_logging():
    """Настройка логирования с поддержкой UTF-8"""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    if platform.system() == 'Windows':
        sys.stdout.reconfigure(encoding='utf-8')
        sys.stderr.reconfigure(encoding='utf-8')
    
    logging.basicConfig(
        filename=f'autodock_batch_{timestamp}.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        encoding='utf-8'
    )

def process_single_ligand(ligand_path, args, receptor_name):
    """Обработка одного лиганда в отдельном процессе"""
    try:
        # Проверяем существование файла результата
        ligand_name = os.path.splitext(os.path.basename(ligand_path))[0]
        result_file_name = f"{ligand_name}_{receptor_name}.dpf.dlg"
        result_file_path = os.path.join(args["output_dir"], result_file_name)
        
        if os.path.exists(result_file_path):
            # Проверяем, что файл не пустой и содержит результаты
            if os.path.getsize(result_file_path) > 0:
                with open(result_file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                    if "FINAL DOCKED STATE" in content:
                        logging.info(f"Skipping {ligand_path} - results already exist")
                        return True
                    
        # Если файл не существует или пустой/поврежден, продолжаем обработку
        process_dir = f"temp_process_{os.getpid()}"
        os.makedirs(process_dir, exist_ok=True)
        
        # Копируем необходимые файлы во временную директорию
        required_files = [
            f"{receptor_name}.maps.fld",
            f"{receptor_name}.A.map",
            f"{receptor_name}.C.map",
            f"{receptor_name}.HD.map",
            f"{receptor_name}.N.map",
            f"{receptor_name}.NA.map",
            f"{receptor_name}.OA.map",
            f"{receptor_name}.SA.map",
            f"{receptor_name}.S.map",
            f"{receptor_name}.Cl.map",
            f"{receptor_name}.F.map",
            f"{receptor_name}.Br.map",
            f"{receptor_name}.P.map",
            f"{receptor_name}.I.map",
            f"{receptor_name}.e.map",
            f"{receptor_name}.d.map"
        ]
        
        for file in required_files:
            if os.path.exists(file):
                shutil.copy2(file, process_dir)
        
        # Копируем файл лиганда
        temp_ligand_path = os.path.join(process_dir, os.path.basename(ligand_path))
        shutil.copy2(ligand_path, temp_ligand_path)
        
        # Переходим во временную директорию
        original_dir = os.getcwd()
        os.chdir(process_dir)
        
        # Создаем DPF файл
        dpf_file = prepare_docking_parameter_file(receptor_name, os.path.basename(temp_ligand_path), args["ga_runs"])
        
        # Запускаем докинг
        success = run_autodock(dpf_file)
        
        if success:
            # Копируем результаты в основную директорию с результатами
            result_file = f"{dpf_file}.dlg"
            if os.path.exists(result_file):
                output_path = os.path.join(original_dir, args["output_dir"], result_file)
                shutil.copy2(result_file, output_path)
                logging.info(f"Results saved for {ligand_path}: {output_path}")
        
        # Возвращаемся в исходную директорию
        os.chdir(original_dir)
        
        # Очищаем временную директорию
        shutil.rmtree(process_dir, ignore_errors=True)
        
        return success
        
    except Exception as e:
        logging.error(f"Error processing ligand {ligand_path}: {str(e)}")
        return False
    
def main():

    config_path = sys.argv[1]

    with open(config_path, "r", encoding="utf-8") as f:
        config = json.load(f)

    args = {
        "receptor": get_path(config["receptor"]),
        "ligands_dir": get_path(config["ligands_dir"]),
        "grid_center": config["grid_center"],
        "grid_size": config["grid_size"],
        "output_dir": get_path(config["output_dir"]),
        "ga_runs": config["ga_runs"]
    }
    # Настройка и проверки
    os.makedirs(args["output_dir"], exist_ok=True)
    setup_logging()
    if not check_autodock():
        return
    
    # Проверка файлов
    if not os.path.isfile(args["receptor"]) or not os.path.isdir(args["ligands_dir"]):
        logging.error("Receptor file or ligands directory not found")
        return
    
    try:
        # Подготовка грида (выполняется один раз)
        grid_center = tuple(map(float, args["grid_center"].split(',')))
        grid_size = tuple(map(int, args["grid_size"].split(',')))
        receptor_name = os.path.splitext(os.path.basename(args["receptor"]))[0]
        
        gpf_file = prepare_grid_parameter_file(args["receptor"], grid_center, grid_size)
        if not run_autogrid(gpf_file) or not check_grid_files(receptor_name):
            logging.error("Grid preparation failed")
            return
    
        # Получаем список всех лигандов
        ligand_files = list(Path(args["ligands_dir"]).glob("*.pdbqt"))
        if not ligand_files:
            logging.error(f"No ligands found in {args['ligands_dir']}")
            return
        
        # Проверяем сколько лигандов уже обработано
        existing_results = 0
        total_ligands = len(ligand_files)
        for ligand_path in ligand_files:
            ligand_name = os.path.splitext(os.path.basename(ligand_path))[0]
            result_file = os.path.join(args["output_dir"], f"{ligand_name}_{receptor_name}.dpf.dlg")
            if os.path.exists(result_file) and os.path.getsize(result_file) > 0:
                existing_results += 1
        
        logging.info(f"Found {existing_results} existing results out of {total_ligands} total ligands")
        if existing_results == total_ligands:
            logging.info(f"All ligands are docked")
            return

        # Определяем количество процессов
        num_processes = 1 # Оставляем один ядро свободным
        logging.info(f"Starting parallel processing with {num_processes} processes")
        
        # Создаем частичную функцию с фиксированными аргументами
        process_ligand = partial(process_single_ligand, args=args, receptor_name=receptor_name)
        
        # Запускаем параллельную обработку
        with Pool(processes=num_processes) as pool:
            results = pool.map(process_ligand, ligand_files)
        
        # Подсчет результатов
        successful = sum(1 for r in results if r)
        failed = len(results) - successful
        logging.info(f"Processing completed. Successful: {successful}, Failed: {failed}")
        
    except Exception as e:
        logging.error(f"Error in main process: {str(e)}")
        raise


if __name__ == "__main__":
    main()