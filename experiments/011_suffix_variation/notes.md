# 011 – deník session

## Rozhodnutí a interpretace spec

* **Kód.** Estimátory v `voytools/variation.py` (`conditional_entropy_decomposition`, `crossval_explained`,
  `neighbour_conditional_mi`, `mode_consistency`, `add_covariates`, `affix_by_covariate`), testy
  `tests/test_variation.py` (syntetická data se známou odpovědí: sufix = f(kmen, pozice) → podíl 1; sufix
  nezávislý → 0; soused kopíruje sufix → MI > 0; stránkový sufix → konzistence > očekávání).
* **Zdroje.** Reál a kontroly přes `experiments/009_register_profile/run.py` (`Ctx`, `control_table`),
  tj. přes `voytools.controls` se stejnými loadery a layoutem ZL split jako v 009. Klon
  `lrozanova/voynich-units` v `../lrozanova/voynich-units` (výchozí cesta 009). Nully: generátory 009
  (`ts2020` = config 21, `tmpl` = šablona 005 bez `suffix_by`, `ac86` = 003 config 86), seedy 0..29.
  ts2020 jen pro Eva (ZL split, ZL join), šablona a ac86 pro ZL split a GC split (jako 009).
* **Kontroly.** Celsus, Plinius (`pliny_prose`, voynich-units), Caesar, Culpeper (`culpeper_prose`),
  Linnaeus – oba layouty (`linnaeus_prose`, `linnaeus_items`). U kontrol jsou `language` a `hand` převzaté
  z layoutu ZL (nemají s textem nic společného), sekce = kniha/oddíl kontroly.
* **Dekompozice (b)** = naučené afixy + poziční glyfy podle 008 (`POS_AFFIXES`, Eva i v101); jen pro
  transliterace a nully, ne pro kontroly (latinka).
* **Kovariáty**: poziční třída 007 (22 obsazených tříd v ZL), `rel_pos` 5 binů, místo v řádku (4),
  `line_in_par` 0/1/2+, jazyk, písař, sekce, stránka; kombinace pozice+jazyk, pozice+písař,
  pozice+jazyk+sekce, vše kromě stránky, vše.
* **Korekce** podle spec: permutace X uvnitř kmene (100× bod, 20× subsamply a nully – rozpočet; 007 mělo
  50/10). Vysvětlený podíl = (H_perm − H(S|Z,X)) / H(S|Z).
  - Vlastnost, kterou spec nezmiňuje: pro X s mnoha hodnotami (stránka 207, „vše" ~4 000 buněk kmen × X)
    je tato korekce **konzervativní** – pod nulou permutace zaplní řídké buňky „náhodnou" redukcí, kterou
    skutečný vztah částečně už vyčerpal; korigovaný podíl je dolní mez a není monotónní (vše < vše kromě
    stránky). Ověřeno na syntetice (test: plný vztah přes 1 000 buněk → korigováno 0,78, surově 1,0).
  - Proto navíc **křížově validovaný** podíl (`crossval_explained`, 5 foldů po tokenech, vyhlazení
    s backoffem na p(s|kmen), α z mřížky 1–100 vybrané na held-out foldech – mírně optimistické).
    Slouží jako citlivostní kontrola; hlavní číslo je podle spec permutační.
  - Při psaní testu se ukázalo, že syntetika i foldy braly stejný proud náhodných čísel (seed 0 → fold =
    kmen); test používá jiný seed.
* **Zbytková struktura** (jen sufix, po podmínění kmen + poziční třída + jazyk + písař):
  - soused: I(S_t; S_{t−1} | ·) a I(S_t; P_{t+1} | ·), páry **uvnitř řádku** (010: vazba končí na konci
    řádku). Korekce: **podmíněná permutace** souseda uvnitř podmiňujících buněk cíle (kmen × pozice × jazyk
    × písař) – zjemnění „permutace uvnitř poziční třídy" ze spec. První verze permutovala uvnitř poziční
    třídy souseda; dávala záporné hodnoty (−0,006, permutace rozbije i vazbu souseda na kmen cíle a nadsadí
    bias), proto změna.
  - konzistence: pro kmeny ≥ 20 tokenů a buňky (kmen, stránka/odstavec) s ≥ 2 tokeny podíl tokenů s modálním
    sufixem buňky vs očekávání po permutaci stránek mezi tokeny téhož kmene **se stejnou poziční třídou,
    jazykem a písařem** (100× / 20×).
  - sekce: I(S; sekce | kmen, pozice, jazyk, písař) toutéž funkcí jako rozklad (permutace sekce uvnitř
    buněk podmínění); navíc stránka.
* **Intervaly** jako 007/009: bod + 2 × kvantily odchylek 80% subsamplů stránek (spec „SD × 2"); kontroly
  30 subsamplů (rozpočet), nully kvantily přes 30 seedů.
* **Deskriptivně** (`top_stems.csv`): p(afix | kmen, místo v řádku) a p(afix | kmen, řádek odstavce) pro
  10 nejčastějších neprázdných kmenů (ZL/GC split, obě dekompozice). Celá poziční třída (22) by byla
  nečitelná.

## Průběh

* Bod ZL split ≈ 40 s, subsample ≈ 10 s; vše na pozadí (`xargs -P 4`), průběžný auto-commit každých 15 min.
* Kód 012 (`voytools/registers.py`) jsem psal lokálně během výpočtu 011; commit a spuštění až po uzavření
  011 (pořadí podle zadání).

## Výsledek a uzavření

* Všech 16 zdrojů doběhlo bez chyb (4 reálné sady × 101, 6 kontrol × 31, 6 nullů × 30). Report:
  `results/011_suffix_variation/report.md`.
* Nečekané: kopírovací nully (ts2020, ac86) mají silnější vazbu sufixu na stránku a větší konzistenci
  uvnitř stránky než Voynich; Voynich má naopak výrazně silnější vazbu sufix → prefix následujícího tokenu.
* Stav: **hotovo** (všechna kritéria přijetí splněna; žádná neřešitelná překážka).
