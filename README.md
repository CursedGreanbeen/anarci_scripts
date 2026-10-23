# filtering

Набор скриптов для фильтрации результатов после BLAST для датасета антител. Он делит попадания BLAST на релизы до 2024 года и с 2024 года по сводке SAbDab. Затем оставляет «свежие» запросы: те, у которых нет совпадений с антителами до 2024 года или максимальная идентичность H3 меньше 90%. После этого объединяет оба набора в итоговый список запросов и PDB. Ещё несколько автономных инструментов (поиск растворителей и лигандов, дедупликация H+L) дополняют набор.

Все скрипты используют только стандартную библиотеку (без numpy и pandas) и запускаются через `/opt/conda/bin/python3.10`. Большинство принимает каталог результатов явным аргументом или использует абсолютные пути, поэтому текущий каталог не важен. Исключения отмечены отдельно. Родительский скрипт `scripts/run_filtering.py` относится к другому конвейеру (`antibody_dataset.filtering`) и эти файлы не использует.

## Конвейер (для каждого подкаталога в `/mnt/1857e392-.../results/`, например `fasta-fv-only-levenshtein`)

```
combined_filtered_hits.csv  (вход, из пакетов BLAST-search)
  → filter_pre2024_antibody_cli.py   → pre_2024_antibody_hits.csv
  → extract_2024_only_queries.py     → queries_2024_only.txt / _details.csv
  → extract_2024_only_pdbs.py        → 2024_only_pdbs.txt / 2024_only_hits.csv
  → extract_unique_queries.py        → pre_2024_unique.csv   (запросы, у которых все попадания имеют идентичность H3 < 90%)
  → merge_unique_with_2024.py        → unique_fresh_queries.csv
  → extract_unique_pdbs.py           → unique_pdbs.txt
```

## Запуск

### filter_pre2024_antibody_cli.py

Оставляет попадания, у которых PDB есть в SAbDab, дата релиза раньше 2024-01-01, а описание содержит ключевые слова антител.

```
cd /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6
/opt/conda/bin/python3.10 scripts/filtering/filter_pre2024_antibody_cli.py \
    results/fasta-fv-only-levenshtein sabdab_summary_all.tsv
```

Позиционные аргументы: `<results_dir>` (в нём должен лежать `combined_filtered_hits.csv`) и необязательный `[sabdab_file]` (по умолчанию `./sabdab_summary_all.tsv`, поэтому важен текущий каталог). Скрипт записывает в тот же каталог `pre_2024_antibody_hits.csv` (с добавленными столбцами `release_date` и `pdb_id`) и `pre_2024_antibody_summary.txt`.

### extract_2024_only_queries.py / extract_2024_only_pdbs.py

В каталоге результатов должны лежать `combined_filtered_hits.csv` и `pre_2024_antibody_hits.csv` (результат первого этапа). Скрипты находят запросы, у которых есть только попадания 2024+, и PDB ID, которых нет в наборе до 2024 года.

```
/opt/conda/bin/python3.10 .../filtering/extract_2024_only_queries.py <results_dir>
/opt/conda/bin/python3.10 .../filtering/extract_2024_only_pdbs.py    <results_dir>
```

Один позиционный аргумент, флагов нет. Результаты: `queries_2024_only{.txt,_details.csv}`, `2024_only_pdbs.txt` и `2024_only_hits.csv`.

### extract_unique_queries.py, merge_unique_with_2024.py, extract_unique_pdbs.py

Аргументов командной строки нет. Пути жёстко заданы: `RESULTS_DIR=/mnt/1857e392-.../results` и три подкаталога `fasta-fv-only-0-100`, `fasta-fv-only-cdr-80-cutoff`, `fasta-fv-only-levenshtein`. Чтобы сменить цель, отредактируйте константы `SUBDIRS` и `RESULTS_DIR`. Запускайте в порядке конвейера. Каждый скрипт выводит, что он записал (`pre_2024_unique.csv` → `unique_fresh_queries.csv` → `unique_pdbs.txt`).

```
/opt/conda/bin/python3.10 .../filtering/extract_unique_queries.py
/opt/conda/bin/python3.10 .../filtering/merge_unique_with_2024.py
/opt/conda/bin/python3.10 .../filtering/extract_unique_pdbs.py
```

### filter_hits_by_date.py

Разовый вариант первого этапа: фильтрует только по дате, без проверки ключевых слов антител. Все пути жёстко заданы (SAbDab лежит в `/mnt/1857e392-.../sabdab_summary_all.tsv`, каталог с попаданиями — подкаталог levenshtein). Аргументов нет. Записывает `pre_2024_hits.csv` и `date_filter_summary.txt`.

### count_unique_pdbs.py

Небольшая утилита: считает уникальные значения в столбце TSV-файла.

```
/opt/conda/bin/python3.10 .../filtering/count_unique_pdbs.py <file.tsv> [--column N]
```

`--column` нумеруется с нуля (по умолчанию 0). Скрипт выводит заголовок и число уникальных значений.

### detect_solvent_ligands.py

Сканирует файлы `*.cif` и относит записи `_chem_comp` к воде, ионам или малым молекулам. Результат записывается в `solvent_ligands_report.tsv`. Перед запуском отредактируйте `cif_folder` и `results_folder`. Аргументов нет.

### deduplicate_antibodies.py

**Деструктивный скрипт с устаревшими путями.** Группирует комплексы антител по точной паре последовательностей (H, L) из `sequence_occurrence_table.tsv`. Для повторяющихся пар оставляет первый по алфавиту комплекс, а у остальных удаляет файлы `.cif` и `.fasta`. Сначала отредактируйте `table_file` и `base_dir`. Аргументов нет.

## Примечания

- `filter_pre2024_antibody.py` (без суффикса `_cli`) — непараметризованный предшественник первого этапа. Он использует относительные пути `sabdab_summary_all.tsv` и `results/fasta_fv_only_h3_summary/...`, поэтому запускается только из родительского каталога `scripts/`. Лучше использовать версию `_cli`.
- Парсер дат SAbDab ожидает формат MM/DD/YY. Попадания, PDB ID которых нет в SAbDab, пропускаются.