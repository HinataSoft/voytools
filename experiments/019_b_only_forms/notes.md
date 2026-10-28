# 019 – deník session

## Rozhodnutí a interpretace spec

* **Kód.** Nový modul `voytools/bforms.py` (používá `registers.py`: `_deletes`, `edit_script`, `levenshtein`, `jsd`):
  `neighbours_within` (všichni kandidáti do vzdálenosti d, index symetrických smazání jako v 016), `is_d_op`,
  `classify_forms`, `composition`, `match_by_count`, `form_profile` (+ `bigram_model`, `bigram_loglik`, `gries_dp`,
  `entropy_bits`), `suffix_pairs`, `pair_choice`, `both_forms`. Testy `tests/test_bforms.py` (sousedé proti hrubé
  síle, třídy a skladba na ručním příkladu, vyrovnání četností, profil, pravidlo). Žádná změna existujících signatur.
* **Reprezentace.** Nosná je **celý token** (016: kmeny dekompozice jsou pro Levenshtein příliš krátké); jednotky =
  glyfy presetu `extended` (ZL), znaky v101 (GC), písmena (latina). **Kmen** tvaru = dekompozice (a) s afixy naučenými
  na celém reálném textu (jako 012/016; u nullů seznamy reálného textu téže sady). **Nový kmen** = kmen bez jediného
  tokenu na druhé straně (v témže losu; v seznamech v plném textu).
* **Třídy tvarů X** (B hlavní, A symetricky): `shared` (≥ 1 token v Y), `d_variant` (X-only, tj. ≥ 2 v X a 0 v Y,
  s protějškem v Y (≥ 2 tokeny) dosažitelným **jednou** Currierovou operací: vložení/smazání složené jen z `e`, `ee`,
  `d`, `q` nebo záměna `ch`↔`sh`; sousední vložení sloučená do jedné operace, proto hledání do 2 editací, takže
  `ins:ed` i `ins:eed` jsou jedna operace), `other_le1` (X-only, jiný protějšek do 1 editace), `no_counterpart`
  (X-only, nic do 1 editace), `hapax` (1 token v X, 0 v Y). Mezi více kandidáty vyhrává ten s více tokeny v Y;
  d-varianta má přednost před jiným protějškem. **v101:** d-jednotky `c` (= e), `C` (= ee), `8` a `7` (obě podoby
  Eva `d`), `4` (= q); záměna `1`↔`2` (= ch↔sh). Latina: žádné d-operace, všechny protějšky ≤ 1 jdou do `other_le1`.
* **Skladba** (metrika 1): podíly tokenů tříd na všech tokenech X (`share_*`) i jen na X-only tokenech (`only_*`,
  bez hapaxů; hapaxy zvlášť `share_hapax`); `only_new_stem` = podíl X-only tokenů s novým kmenem.
* **Vyrovnání četností** (metrika 2): ke každému tvaru třídy se losuje sdílený tvar se stejným počtem tokenů v losu
  (bez vracení; když dojdou, nejbližší počet). Hapaxy se vyrovnávají sdílenými tvary s 1 tokenem v X (a ≥ 1 v Y).
  Profil třídy i vyrovnané základny týmž kódem; rozdíl v témže losu (`diff`), 200 losů → interval.
* **Profil.** `paradigm_forms` = počet různých tvarů vzorku se stejným kmenem (průměr přes tvary třídy);
  `H_suffix_given_stem` přes tokeny vzorku, jejichž kmen patří tvarům třídy; disperze: `pages_per_token`,
  Griesův `gries_dp` (stránky losu jako části), `top_page_share`; `section_entropy`; popisky z **plného textu**
  (`label_share` = L / (L + P) tokenů tvaru v jazyce stránky, `label_any` = podíl tvarů s aspoň jedním popiskem);
  kraje řádků (`first_share`, `last_share`, `par_first_share`, sdruženě přes tokeny); `cooc_ratio` = dvojice tokenů
  třídy v témž řádku / průměr 20 promíchání příslušnosti mezi tokeny vzorku; `bigram_ll` = průměrná log2
  věrohodnost přechodu (s okraji, add-0,5) pod bigramovým modelem **sdílených tvarů mimo vyrovnaný vzorek**
  (váženo tokeny); `seam_jsd` = JSD rozdělení prvního glyfu následujícího tokenu po tokenech třídy proti
  rozdělení po všech sdílených tokenech (vyrovnaná základna: totéž pro vyrovnaný vzorek – ukazuje úroveň šumu
  při stejné velikosti).
* **Srovnání v losu:** `AB y_only` = B-only vs A (hlavní), `AB x_only` = A-only vs B, `BB` = B1-only vs B2 (šum),
  `AA` = A1-only vs A2. Velikost N = polovina menšího jazyka (celé stránky), 200 losů; jen výřez celkem.
* **Pravidlo `-edy`** (metrika 3) na plném textu P: dvojice `-ey`/`-edy` (token končí `edy` → volba 1, jinak
  končí `ey` → 0; základ = zbytek, neprázdný), srovnávací dvojice `-y`/`-dy` (základy končící `e` vyloučeny, aby se
  nepočítala táž dvojice) a `-ol`/`-or`. **v101** (ověřeno zarovnáním ZL/GC po locích: `chey`→`1c9`, `chedy`→`1c89`
  / `1c79`, `ol`→`oe`, `or`→`oy`): `-c9`/`-c89`, `-c9`/`-c79` (GC rozlišuje dvě podoby `d`), `-9`/`-89` (základy na
  `c` vyloučeny), `-oe`/`-oy`. Kontexty: jednotka před sufixem, prefix dekompozice (a), pozice v řádku
  (první/uvnitř/poslední), poslední jednotka předchozího tokenu, jednotka+prefix. Entropie volby bez a s kontextem,
  200 permutací volby mezi tokeny, `share_determined` = podíl tokenů v kontextech (≥ 5 tokenů) s P ≤ 0,1 nebo ≥ 0,9,
  `share_middle` = P v [0,3; 0,7]; kmeny (= základy ≥ 2 tokeny) s oběma tvary a z nich s oběma na téže stránce.
  Latina: `-us`/`-um`, `-is`/`-es` jen jako hrubé měřítko.
* **Nully:** ts2020 config 21 (ZL split/join) a šablona 005 (ZL/GC split), 30 seedů, A/B podle layoutu, cache
  `scratch/019_cache/` (gitignored, 016 cache nebyla v repu). **Latina:** Plinius 16/18/21 vs 23–30 (téma) a Caesar
  vs Celsus (registr) přes `M9.control_table` jako v 016 (klon `lrozanova/voynich-units` v `../lrozanova/`).
* **Seznamy** (metrika 5): `classes.csv` = všechny X-only tvary a hapaxy v plném textu (B vs plné A, A vs plné B)
  s protějškem, operací, kmenem, stránkami, sekcemi, popisky; `top_new_stems.csv` = 30 nejčastějších nových kmenů
  B a A. Doplněk `new_stems_ed.csv`: podíl tokenů s novým kmenem, jejichž kmen obsahuje `ed` (v101 `c8`/`c7`) –
  test, zda „nový kmen“ není jen `d` uvnitř kmene, které parser neodtrhl.

## Průběh

* Jeden los ZL split ≈ 2,8 s (4 srovnání × skladba + 4 třídy × 2 profily); 200 losů ≈ 10 min na sadu samostatně,
  ~25 min při čtyřech paralelních bězích. Latina 200 losů, generátory 30 seedů; vše `nohup` na pozadí, logy
  `scratch/019_logs/`.
* Prostředí: Python 3.13 / pandas 3.0.5 (cloud); kód drží 3.10 / pandas 2.3.

## Slepá ulička

* První běh vyrovnával hapaxy nejbližším počtem (sdílených tvarů s 1 tokenem je v losu jen ~180 na ~1 070 hapaxů),
  takže „vyrovnaná základna“ hapaxů měla 2–5 tokenů a bigramový model se učil z prázdné množiny. Opraveno:
  `match_by_count(tol=0.25)` – přesná shoda do 3 tokenů, do ±25 % výš, nespárované tvary se vynechají (profil
  hapaxů tedy stojí na ~190–210 hapaxech na los). Všechny losy spočítány znovu.

## Výsledek a uzavření

* Report: `results/019_b_only_forms/report.md`. Pětina B = 59–63 % d-variant, 25–27 % `l`-variant (jiný protějšek ≤ 1),
  12–15 % bez protějšku (kombinace obou); nové kmeny 6,5–7,7 % tokenů pětiny, pod polovinami B (7,6–10,9 %). Profil
  nových tvarů = sdílené (disperze, sekce, popisky, kraje, souvýskyt, bigramy v intervalech přes nulu), liší se jen šev
  (+0,04–0,075) a konec řádku u d-variant (−0,04); latina (téma i registr) dává shluky, cizí bigramy a souvýskyt.
  `-ey`/`-edy` v B je volba tvaru (P 0,3–0,8, redukce 0,02–0,05 bitu, 79–80 % kmenů s oběma tvary), v A neexistuje.
* Soubory: `runs/<zdroj>__<sada>/runs.csv` (0,6–1,9 MB), `composition.csv`, `class_profiles.csv`, `latin_reference.csv`,
  `edy_rule.csv`, `edy_contexts.csv`, `classes.csv`, `top_new_stems.csv`, `new_stems_ed.csv`, `summary.csv`, `meta.json`,
  tři grafy; `extras/<zdroj>__<sada>/` po sadách.
* Stav: **hotovo**.
