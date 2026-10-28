# 008 – deník session

## Stav vstupů
* `main` při startu = eb25723 (004 i 005 sloučené). Větev `exp/008-ts2020-recheck`.
* Prostředí: Python 3.13, pandas připnut na 2.3.3 (jako u vlastníka), numpy 2.5; bez scipy.

## Změny kódu (bez změny výstupů pro dosavadní zdroje)
* `experiments/007_position_vs_content/run.py`: nové zdroje `selfcitation` a `template` (`--source`),
  volba `--out-root` (výchozí `results/007_…` beze změny). Vzorky ts2020 se cachují v
  `scratch/007_cache/<sada>/ts_<seed>.txt` jako u autocopy. Do `meta.json` nových běhů jdou parametry nullu.
* `experiments/006_topic_at_distance/run.py`: nové zdroje `selfcitation_21` a `template_005` mimo `SOURCES`
  (takže `--sources all` a `summarize` 006 se nemění), volba `--out-root`.
* `experiments/008_ts2020_recheck/run.py`: kroky 2, 4, 5, kontrola reprodukce a souhrny; importuje skripty
  005/006/007 (jejich funkce beze změny). V balíku `voytools/` se nic nemění.

## Nully
* **ts2020** = `selfcitation_like(real, seed, p_copy=0,1, replace_weights=(20,50,30), p_same_position=14)`
  (řádek `ts2020_calB` v `results/004_selfcitation_port/best_configs.csv`, config 21). Seedy 2100 + rep =
  seedy fáze B 004 → vzorky jsou totožné s 004 (ověřeno: seed 2100 ZL split 7 876 typů, hapax 0,5924 jako
  v `results/004_selfcitation_port/phaseB/ZL_split/runs.csv`). V 006 a 007 30 seedů (2100–2129), v syntéze
  (řádky 1–6) 50 seedů fáze B 004. Jen ZL (split, join): generátor vyrábí Eva; u GC null chybí.
* **tmpl** = plný model 005 (`ParadigmTemplate`; ZL: 8 000 kmenů, zipf 1,4; GC: 4 000, zipf 1,3; fixed,
  bez vazby prefixu, α = 5, λ = 0,3; konfigurace z `meta.json` ablace `best` 005). Seedy 100000 + rep,
  `pool_seed = rep` = vzorky ablace `best` 005 (ověřeno: seed 100000 ZL split 7 561 typů, GC 9 352, jako
  v 005). U `join` se použije konfigurace ze `split` (ablace 005 běžela jen na split), naučená na join textu.
  30 seedů.
* Pozor (zjištěno v kódu portu, `voytools/selfcitation.py`, `choose_source`): ts2020 **nekopíruje jen uvnitř
  stránky**. První řádek stránky bere zdroj z posledního řádku předchozí stránky, počáteční řádky odstavců
  z kteréhokoli dřívějšího počátečního řádku odstavce v celém textu a „suggestion“ z globálních četností.
  Předpoklad spec „plochá mezistránková křivka“ tedy není z konstrukce zaručený.

## Výklad spec
* **„9 řádků oddílu 1 syntézy“**: tabulka oddílu 1 má 6 řádků. Doplnil jsem tři nálezy 006/007, na které se
  ptá otázka spec a kritérium 2: (7) sekční informace kmenů nezávislá na písaři I(kmen; sekce | písař),
  (8) strmý pokles sdílení se vzdáleností (koef. `log_dist`, cos, K = 30, model spec+size),
  (9) „žádná MI za d = 2“ (plocha d = 2…10, kmen, korekce `pc`).
* **Reprodukuje** = |z| < 3, z = (reál − průměr nullu) / SD nullu přes seedy (jako SYNTHESIS). Verdikt řádku
  „ano“ jen když platí ve všech sadách, kde null existuje (ts2020: ZL split i join; tmpl: ZL i GC).
* **Krok 2 (b)**: spec píše prefixy „cp, cf, ct, ck“. V Eva (preset `basic`) jsou bench gallows jeden glyf
  `cph cfh cth ckh`; prefix `cp` by v parseru nebyl glyfovou hranicí. Použil jsem `cph cfh cth ckh`.
  Pro GC (v101) jsou ekvivalenty odvozené zarovnáním lokusů ZL/GC se stejným počtem tokenů
  (`parser_affix_map.csv`): p = g/j, f = f/u, t = k, k = h, cph = G/J, cfh = F/U, cth = K, ckh = H;
  m = p/P, g = *, am = ap, om = op. Parser má jediný prefixový slot (greedy longest match), takže gallows
  za `qo-`/`o-` zůstávají v kmeni i ve variantě (b).
* **Krok 2 (c)**: jednotky = BPE 64 sloučení (`learn_bpe` z 007, preset `extended` u Eva, znaky u v101);
  kmen = spojení prostředních jednotek; token s ≤ 2 jednotkami má prázdný kmen (vlastní kategorie „“;
  podíl v `parser_variants.csv`, `share_empty_stem`).
* **Krok 2 metriky**: stejný odhad jako 007 (plug-in MI, korekce permutací Y uvnitř vrstev, top 200 kmenů +
  `<other>`), bod (50 permutací) + 100 subsamplů (10 permutací; seedy 0–99 = prvních 100 subsamplů 007).
  Navíc `I_place` (jen místo v řádku: first/mid/last/single, 4 hodnoty) jako kontrola kapacity proměnné:
  22 pozičních tříd proti 7 sekcím. Párové rozdíly (I_pos − I_sec …) počítané v témže subsamplu.
* **Krok 4**: JSD sufixů (podpora = všechny sufixy, jako `ab_contrasts` 005) A vs B celkem a A3 vs B3.
  Podmíněná hodnota = vážený průměr JSD uvnitř 22 pozičních tříd 007, váha = počet tokenů obou skupin
  v třídě; třída se započte jen při ≥ 20 tokenech v každé skupině (pokrytí ve sloupci `coverage`).
  Vnitřní tokeny = ani první, ani poslední v řádku (nepodmíněné i podmíněné uvnitř vnitřních tříd).
  Null: permutace štítku jazyka mezi stránkami **uvnitř písaře** (200×, spec) a navíc mezi všemi stránkami
  (200×, jako 005). Uvnitř písaře se u A vs B celkem mění jen písař 3 (ostatní písaři píší jeden jazyk),
  takže tento null je u „A vs B celkem“ z konstrukce blízký reálu. JSD menších skupin je vychýlená nahoru,
  proto rozhoduje přebytek nad permutací; interval reálu jako v 005 (SD subsamplu · 2, 200 subsamplů
  80 % stránek, seedy 1–200) a interval přebytku ± 1,96·√(SD_reál² + SD_perm²).
* **Krok 5, řádky 1–3**: šablona 005 tyto metriky nepočítala → `run.py profile` (stejná volání jako
  `stats_for` 004: `edge_glyph_mi` 20 míchání, `positional_profile`, `bpe_dependence_curve` po řádcích),
  30 seedů ZL i GC split; reál GC spočítán stejným kódem.

## Průběh
* Fáze 1+3: 12 úloh (`scratch/008/job.sh`, 4 procesy): 007 `decomp`+`mi` pro ts2020 (ZL split/join) a šablonu
  (4 sady), 006 `pages` pro ts2020 (ZL) a šablonu (4 kombinace). ts2020 ~20 s na seed, šablona ZL ~1–2 min
  (zásoba 8 000 umělých kmenů se tvoří pro každý `pool_seed` znovu).
* Cache vzorků šablony v 007 přidána během běhu (`tpl_<seed>.txt`); první verze čtení cache skládala tabulku
  přes `as_table` bez presetu, což u v101 mění `len_glyphs` (sbalení jako Eva). Opraveno předáním presetu
  reálného textu dřív, než se cache poprvé četla; ověřeno, že vzorek z cache = čerstvý vzorek ve všech
  sloupcích (ZL i GC, šablona i ts2020).
* Vedlejší pozorování v 007 (neopraveno, mimo rozsah): `Context.autocopy` skládá vzorek `as_table(lines, real)`
  s výchozím presetem Eva i pro GC, takže cíl „délka“ u GC autocopy počítá v101 znaky po sbalení `ch/sh/…`.
  Týká se jen řádků `length` autocopy u GC v `results/007_…`.
* Krok 4: první verze podmíněné JSD (průměr JSD uvnitř 22 tříd) je u písaře 3 (A jen 733 tokenů) silně
  vychýlená malými buňkami (permutace uvnitř písaře dává 0,10 proti reálu 0,12, pokrytí tříd 46 %). Přidána
  **přímá standardizace** (`std`, `mid_std`): rozdělení sufixů A i B převážená na společný mix pozičních tříd
  (třídy s ≥ 5 tokeny v obou skupinách), pak jedna JSD. Krok 4 přepočítán celý (výstup první verze smazán).
  Hlavní číslo kroku 4 je `std`, průměr uvnitř tříd (`cond`) je uveden pro úplnost.
* Slepá ulička: `pkill -f "run.py ab"` a `pgrep -f "[0]08…"` zabily vlastní shell (vzor se shodoval s příkazem,
  stejně jako v 006). Řešeno spouštěcím skriptem `scratch/008/ab.sh`. Parser a `check` tím zasaženy nebyly.
* Kontrola reprodukce (`check_reproduction.json`): 007 rozklad i MI(d) přesně (max. rozdíl 1·10⁻¹⁶). 006 `real`:
  matice podobností jsou totožné s `pairs.csv` (rozdíl 0), ale koeficienty se liší až o 0,0015 (K = 0),
  0,0007 (K = 30) a 0,012 (K = 100, GC; jen R²). Příčiny: (1) 5 stránek bez jazyka (f65v, f67r2, f68v2, f68v3,
  f70r2) dává pandas 2.3.3 jako `None`, ne `NaN` jako v session 006 → `same_language` se liší u 10 z 20 100
  dvojic; (2) u K = 100 pořadí kmenů se shodnou četností (v GC se tím mění množina vyřazených kmenů).
  Nully 008 se srovnávají s reálem ze souborů 006, jak je spec žádá; rozdíl je o řád menší než rozdíly reál–null.
* Souhrny: `run.py summarize` → `decomp_ts2020.csv`, `mi_distance_ts2020.csv`, `pages_ts2020.csv`,
  `parser_variants.csv`, `ab_by_position.csv`, `synthesis_update.csv`. Řádky 7–9 syntézy mají po několika
  metrikách (verdikt po řádku × metrice × nullu); u řádku 9 je primární plocha s korekcí `line` a MI při d = 3,
  plocha s korekcí `pc` měří podle 007 hlavně složení řádku.
* `results/SYNTHESIS.md` neupravuji (udržuje Cowork); podklad je `synthesis_update.csv` a report.
