# results/

Výstupy experimentů, jedna složka na experiment (`NNN_short_name/`, stejné jméno jako v
`experiments/`). Obsah a pravidla: `experiments/README.md`.

Commitují se jen malé soubory (CSV, MD, JSON, PNG; < 5 MB). Velké mezivýsledky (`*.parquet`,
`*.npy`, `*.pkl`) jsou v `.gitignore`.

Paralelní session zapisují do `NNN_short_name/partK/` a sloučený výsledek vzniká
`voytools.merge_runs([...], "results/NNN_short_name")`.
