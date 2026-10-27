# experiments/ – protokol

Každý experiment má vlastní složku `NNN_short_name/` s těmito soubory:

| Soubor | Obsah |
|---|---|
| `spec.md` | zadání: otázka, hypotézy a co by ukázal smysluplný vs. generovaný text, data a filtry, metriky, nulové modely, počet seedů/bootstrapů, výstupy, kritéria přijetí |
| `run.py` | kód; používá `voytools` a `voytools.experiment.run_seeds`; parametry přes argparse (`--seeds`, `--out`, `--chunk`), výchozí hodnoty = hodnoty ze spec |
| `notes.md` | (volitelné) deník session: rozhodnutí, odchylky od spec, slepé uličky |

Výsledky jdou do `results/NNN_short_name/` (stejné číslo a jméno):

| Soubor | Obsah |
|---|---|
| `runs.csv` | jeden řádek na seed/běh (zapisuje `run_seeds`, průběžně, navazovatelně) |
| `summary.csv` | agregace (`summarize_runs`) nebo vlastní tabulka efektů s intervaly |
| `report.md` | čísla (jen z `results/`), interpretace ve 2–3 větách, co výsledek **neukazuje**, další krok |
| `*.png` | grafy (matplotlib), malé |
| `meta.json` | parametry, časy, počet seedů (zapisuje `run_seeds`) |

Pravidla:

1. **Nulové modely ve stejném běhu.** Každá metrika se počítá pro reálný text i pro každý nulový
   model ve stejné funkci a stejném filtru; výsledky v jedné tabulce (sloupec `source`).
2. **ZL i GC.** Pokud spec neříká jinak, všechno dvakrát (`data/ZL3b-n.txt`, `data/GC2a-n.txt`).
3. **Tři politiky mezer** tam, kde jde o tokeny (`split`, `join`, `strict`); u glyfů preset
   `extended` a případně BPE jednotky.
4. **Seedy a bloky.** `run_seeds` zapisuje po každém seedu; bloky volte tak, aby jeden příkaz
   skončil do 10 minut. Paralelní session = disjunktní `seed_range`, podsložky `partK/`, nakonec
   `merge_runs`.
5. **Commit a push** po každém dokončeném bloku a nejméně každých ~15 minut během dlouhého běhu.
   Větev `exp/NNN-short-name`.
6. **Report píše ten, kdo počítal.** Bez čísel z `results/` žádný závěr. Interpretace opatrná:
   co by vysvětlil generátor, co pozice, co je artefakt abecedy.
7. **Nikdy** nenavrhovat čtení, zvukové hodnoty ani „řešení“; nepracovat s okrajovou literaturou
   (viz `CLAUDE.md`).

Šablona: `_template/` (`spec.md` + `run.py`). Zkopírujte ji a doplňte.

## Seznam experimentů

| # | Název | Stav | Otázka |
|---|---|---|---|
| 001 | autocopy_calibration | hotovo (kritérium 1 nesplněno, viz report) | kalibrace autokopírovacího generátoru na ZL; které vlastnosti textu zůstanou nevysvětlené |
| 002 | paradigm_coherence | zadáno | chovají se kmeny a afixy jako paradigma, nebo jako šablona? |
