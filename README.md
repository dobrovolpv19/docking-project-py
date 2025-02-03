# docking-project-py

## Автоматизированный молекулярный докинг

### Описание
Автоматизированный молекулярный докинг — это проект для проведения молекулярного докинга с использованием Modeller, RDKit и MGLTools. Проект предназначен для работы на локальной машине и выполняет автоматизированный анализ молекулярных взаимодействий.
Рекомендуется прочесть docx файл, для корректного запуска работы.

### Зависимости
Перед запуском убедитесь, что установлены следующие зависимости:

- **Python**: 3.10.13
- **Modeller**: 10.6
- **MGLTools**: 1.5.7
- **RDKit**: 2022.03.5
- **AutoDock**: 4.2.6

#### Команды для установки:
```sh
pip install rdkit==2022.03.5
```

Modeller, MGLTools и AutoDock 4.2.6 не устанавливаются через `pip`, их необходимо загрузить и установить вручную с официальных сайтов:

- [Modeller](https://salilab.org/modeller/)
- [MGLTools](http://mgltools.scripps.edu/downloads)
- [AutoDock 4.2.6](http://autodock.scripps.edu/downloads)

Скачайте и установите **Python 3.10.13** с официального сайта:
[Python 3.10.13](https://www.python.org/downloads/release/python-31013/)

После установки убедитесь, что используется нужная версия:
```sh
python --version
```

### Установка и запуск

1. **Клонирование репозитория**
```sh
git clone https://github.com/dobrovolpv19/docking-project-py.git
cd docking-project-py
```
2. **Запуск подготовки лигандов**
```sh
python Ligand_processing/mol_process.py Ligand_processing/config.json
```
3. **Запуск молекулярного докинга**
```sh
python Docking/docking.py Docking/config2.json
```
4. **Запуск сбора результатов**
```sh
python Result_analysing/table_after_docking.py
```
To do:
5. **Экспериментальный процесс, ещё предстоит оптимизировать**
```sh
python Result_analysing/binding_analyzing.py Result_analysing/config3.json
```