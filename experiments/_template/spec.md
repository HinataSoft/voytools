# NNN – krátký název

## Otázka
Jedna věta. Co se dozvíme, když výsledek dopadne tak, nebo onak.

## Hypotézy a predikce
- **Smysluplný text:** co očekáváme (směr, velikost efektu).
- **Generovaný text / blábol:** co by ukázal autokopírovací generátor, míchání, slotový generátor.
- **Artefakt:** co by znamenalo, že efekt zmizí v GC (v101) nebo při jiné politice mezer.

## Data a filtry
- Soubory: `data/ZL3b-n.txt` (primární), `data/GC2a-n.txt` (kontrola).
- Tokeny: `filter_tokens(tables.tokens, generic_type="P", clean=True)`; politiky mezer `split|join|strict`.
- Podskupiny: `language` A/B, `illustration`, `hand`; co se vyhodnocuje zvlášť.

## Metriky
Přesná definice každé metriky (funkce z `voytools`, parametry), jednotky, očekávaný rozsah.

## Nulové modely
Které (`shuffle_within_lines`, `shuffle_within_pages`, `autocopy_like(...)` s parametry z 001, jiné),
kolik opakování, jaký seed.

## Nejistota
Bootstrap po stránkách/kvaternech (kolik), permutační testy, počet seedů; co je „jeden běh“ v `runs.csv`.

## Výstupy
`results/NNN_name/`: `runs.csv`, `summary.csv`, `report.md`, grafy (`*.png`) – které.

## Kritéria přijetí
Kdy je experiment hotový; co musí report obsahovat; jaké rozsahy hodnot by byly podezřelé (bug).

## Rozpočet
Odhad času na jeden seed, počet seedů, doporučené rozdělení mezi session (`seed_range`).
