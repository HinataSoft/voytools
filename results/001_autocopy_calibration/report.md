# 001 – Kalibrace autokopírovacího generátoru: co zůstane nevysvětleno

Kód: `experiments/001_autocopy_calibration/run.py`; deník a výklad spec: `experiments/001_autocopy_calibration/notes.md`.
Všechna čísla níže pocházejí ze souborů v této složce (`grid_summary.csv`, `best_configs.csv`,
`phaseB_configs.csv`, `profile.csv`, `diag_summary.csv`, `real/*/point.csv`, `summary_info.json`).

## Co bylo spočítáno

| Krok | Rozsah | Soubor |
|---|---|---|
| Reálný text (seed 0) + bootstrap po stránkách (200×) | ZL3b, GC2a × `split`, `join`; P text, clean | `real/<sada>/point.csv`, `real/<sada>/bootstrap/` |
| Fáze A: 576 konfigurací × 3 seedy, cílové metriky + edge/token MI | ZL split | `runs.csv` (sloučeno z `phaseA/part0..3`), `grid_summary.csv` |
| Fáze B: 10 nejlepších z A × 50 seedů, všechny metriky | ZL split, ZL join, GC split, GC join | `phaseB/<sada>/runs.csv`, `phaseB_configs.csv`, `best_configs.csv` |
| Čisté míchání v řádcích, 50 seedů, všechny metriky | všechny 4 sady | `shuffle/<sada>/runs.csv` |
| **Mimo spec:** diagnostika `p_new` / seed slovník / bez editací (16 bodů × 3 seedy) | ZL split | `diag/`, `diag_summary.csv` |
| **Mimo spec:** nejlépe kalibrovaný diagnostický bod × 50 seedů, všechny metriky | všechny 4 sady | `diag_full/<sada>/runs.csv` |

Generátor: `autocopy_like(real, allowed_from=<všechna reálná slova>, seed, **config)`. Bigramový model
legálních editací je odhadnut ze všech reálných slov a seed tokeny jsou prvních 100 tokenů. Na GC
(v101) substituuje mezi libovolnými glyfy, omezeno jen bigramovým modelem (prázdná `edit_table`).
H2 je na ZL počítáno s presetem `extended`. Na GC je preset `None`, protože v101 znaky jsou glyfy a
Eva kompozity by je chybně slučovaly. Referenční hodnoty Rozanova & Temerev tab. 7 v repu nejsou a
spec nedovoluje je stahovat, proto v `profile.csv` chybí.

## Čísla

**Reálný text** (`real/*/point.csv`; ± = 1,96·SD page bootstrapu):

| metrika | ZL split | ZL join | GC split | GC join |
|---|---|---|---|---|
| share_near | 0,1726 ± 0,008 | 0,1524 | 0,2023 | 0,1789 |
| near_ratio | 1,183 | 1,197 | 1,139 | 1,160 |
| n_types | 6 796 ± 374 | 7 729 | 8 351 | 9 011 |
| hapax_share | 0,676 ± 0,070 | 0,711 | 0,702 | 0,722 |
| mean_len | 4,988 | 5,363 | 3,874 | 4,126 |
| H2 | 2,689 ± 0,044 | 2,688 | 2,897 | 2,897 |
| edge_mi_excess | 0,185 | 0,152 | 0,226 | 0,182 |
| token_mi_excess_share | 0,0089 | 0,0071 | 0,0129 | 0,0101 |
| bpe_gap64 | 1,058 | 1,065 | 1,126 | 1,138 |
| pos_jsd_first_mid | 0,186 | 0,183 | 0,195 | 0,191 |
| cmi_excess | 0,088 | 0,078 | 0,089 | 0,077 |
| frac_stems_multi_suffix | 0,912 | 0,907 | 0,839 | 0,871 |
| mean_suffixes_per_stem | 5,25 | 5,09 | 4,72 | 4,84 |

Hodnoty ZL split odpovídají baseline ze spec (0,173 / 0,676 / ~6 800 / 4,99 / 2,69 / 0,185 / 0,9 % / 1,05).

**Kalibrace (kritérium 1 – nesplněno).** Nejlepší konfigurace mřížky (fáze B, 50 seedů) je
`config_id 573` (window 32, p_copy_exact 0,5, p_new 0,15, p_two_edits 0,5, recency 0,
edit_type_probs 0,4/0,3/0,3) s **loss 0,978** a 3 netrefenými složkami: share_near (0,233),
n_types (0,189) a H2 (0,117); hapax 0,044, mean_len 0,0003. Žádná z 576 konfigurací nemá loss
< 0,15; nejnižší loss ve fázi A je 0,964 a každá konfigurace má aspoň 2 netrefené složky. Na ostatních
sadách vyšla nejlépe `config_id 536` (loss: ZL join 0,890, GC split 0,649, GC join 0,604) a
H2 je netrefené všude. Všech 10 nejlepších konfigurací leží na okraji mřížky (window 32, p_new 0,15).

Proč to nejde (`grid_summary.csv`, `diag_summary.csv`):

1. **H2 je v mřížce nedosažitelné.** Přes 576 konfigurací je průměrné H2 3,52–3,69 bitu, reálně
   2,69 (minimální relativní odchylka 31 %, složka ≥ 0,095). Diagnostika ukazuje příčinu: bez
   editací (`p_copy_exact = 1`) je H2 2,66–2,76, tedy na úrovni reálného textu. Každá editace ale
   vytváří slovo, které je sice složené z doložených bigramů, ale celkově je méně předvídatelné než
   reálná slova. Editace jsou přitom jediný zdroj nových typů a hapaxů. Bez editací generátor
   vytvoří málo typů (87 typů při seedu 100 tokenů; 1 765–4 542 při seedu z celého slovníku), takže
   hapax_share a n_types nejdou trefit.
2. **share_near a n_types/hapax jdou proti sobě.** Nízké share_near vyžaduje hodně editovaných
   kopií, které rozšiřují slovník. Malý slovník vyžaduje přesné kopie, které zvyšují share_near.
   Žádná konfigurace mřížky netrefí share_near a n_types zároveň (0 z 576). Jediná páka, která snižuje
   obojí, je `p_new` (čerpání ze 100 seed tokenů), a ta je v mřížce omezena na ≤ 0,15.
3. Mimo mřížku (p_new 0,5, p_copy_exact 0,5, jinak config 573) má generátor na 50 seedech loss
   **0,202** (ZL split; ZL join 0,213, GC split 0,175, GC join 0,176). share_near, n_types, hapax a
   mean_len jsou trefené, jedinou netrefenou složkou zůstává H2 (0,063 / 0,064 / 0,053 / 0,046; na GC
   join tedy těsně pod prahem). ZL split: H2 3,365 vs 2,689 reálně. I tento bod nesplňuje loss < 0,15.

**Necílové metriky** (`profile.csv`; z = (reál − průměr generátoru) / SD generátoru přes 50 seedů).
`best` = nejlepší konfigurace mřížky pro danou sadu, `diag` = diagnostický bod z odst. 3:

| metrika | ZL split best / diag | ZL join best / diag | GC split best / diag | GC join best / diag | reál vs. generátor |
|---|---|---|---|---|---|
| edge_mi_excess | 228 / 237 | 148 / 185 | 165 / 177 | 164 / 168 | generátor ≈ 0 (např. 0,0001 vs 0,185) |
| token_mi_excess_share | 7,6 / 7,7 | 4,4 / 6,5 | 10,2 / 13,2 | 5,9 / 10,6 | generátor ≈ 0 (|x| ≤ 0,0006 vs 0,007–0,013) |
| pos_jsd_first_mid | 403 / 432 | 461 / 480 | 317 / 365 | 309 / 362 | generátor 0,001–0,003 vs 0,18–0,20 |
| cmi_excess | −20,5 / −24,6 | −14,9 / −24,3 | −15,7 / −27,1 | −15,4 / −24,5 | generátor **vyšší** (0,19–0,24 vs 0,08–0,09) |
| frac_stems_multi_suffix | 23,9 / 12,6 | 17,0 / 15,8 | 8,5 / 5,5 | 13,9 / 10,1 | generátor 0,51–0,74 vs 0,84–0,91 |
| mean_suffixes_per_stem | 77,0 / 52,3 | 44,5 / 50,2 | 30,9 / 25,1 | 34,5 / 29,1 | generátor 2,3–3,3 vs 4,7–5,3 |
| bpe_gap64 | 37,2 / **2,7** | 28,0 / **0,6** | 51,5 / 12,9 | 55,8 / 10,0 | ZL diag: 0,975 vs 1,058; GC diag 0,93 vs 1,13 |
| near_ratio (doplňková) | −16,4 / −13,2 | −17,2 / −12,5 | −16,3 / −14,3 | −16,9 / −14,4 | generátor 1,31–1,69 vs 1,14–1,20 |

Míchání v řádcích (`shuffle_lines`) dává edge/token MI ≈ 0 a pos JSD 0,003–0,005, stejně jako
generátor. Typové metriky (n_types, hapax, mean_len, paradigmata) jsou při míchání beze změny, jak
mají být; z-skóre tam proto není definováno (NaN).

**Kritérium 2.** Necílové metriky, u kterých je reálný text > 3 SD od nejlepšího generátoru
ve **všech čtyřech** sadách a pro oba generátory (best i diag): `edge_mi_excess`, `token_mi_excess_share`,
`pos_jsd_first_mid`, `cmi_excess` (reál níž), `frac_stems_multi_suffix`, `mean_suffixes_per_stem`.
Generátor reprodukuje: `bpe_gap64`, ale jen na ZL a jen u lépe kalibrovaného diagnostického bodu
(z 2,7 a 0,6). Na GC zůstává reál o 10–13 SD výš.

**Kritérium 4 (podezřelé hodnoty).** 18 běhů fáze A má share_near > 0,9 (konfigurace 72, 74,
108–111: window 4, p_new 0,02, p_two_edits 0,1). Jde o konfigurace s malým oknem a málo novými
slovy, kde vysoký podíl blízkých sousedů je očekávaný, nikoli chyba. H2 mimo 1,5–4 ani záporný hapax
se nevyskytly. Fáze B, míchání i diagnostika mají 0 podezřelých běhů (`summary_info.json`).

**Kritérium 3.** `python -m pytest -q`: 166 passed, 1 skipped. Všechny příkazy `run.py` jsou
navazovatelné (`run_seeds`) a dopočítají jen chybějící seedy.

## Interpretace

Jednoduchý generátor „kopíruj a obměň“ nejde v zadané mřížce nakalibrovat na všech pět cílových
vlastností najednou. Jeho editace vytvářejí slova s vyšší glyfovou entropií (H2 o 0,6–1,0 bitu), než
má reálný text, a zároveň nejde současně udržet nízký podíl blízkých sousedů a reálnou velikost
slovníku. Vazba okrajových glyfů, pořadí tokenů a poziční profil prvního glyfu řádku jsou v
generátoru nulové (stejně jako po míchání). Reálný text se navíc od generátoru liší v paradigmatické
struktuře, a to ve všech čtyřech sadách a nezávisle na kalibraci. Kmeny reálného textu se pojí s
více sufixy (5,3 vs 2,3–2,8 na kmen, ZL split) a podmíněná závislost prefixu a sufixu daného kmene je
u reálného textu *nižší* než u generátoru (0,088 vs 0,22–0,23). Hloubka BPE žlabu se naopak ukázala
citlivá na kalibraci: když se generátor přiblíží cílovým vlastnostem, rozdíl na ZL zmizí. BPE žlab
proto podle pravidla „platí jen v ZL i GC“ nebereme jako indicii.

## Co výsledek neukazuje

* Neukazuje, že text nese význam. Ukazuje jen, že **tento** jednoduchý generátor neprodukuje vazbu na
  hranici slov, poziční efekty ani paradigmatickou strukturu. Publikovaný algoritmus Timm & Schinner
  (2020) má podrobnější pravidla (poziční, řádkové) a mohl by část těchto vlastností vytvořit.
  Poziční profil a edge MI nemá generátor jak vytvořit, protože nemá žádný poziční ani mezislovní
  mechanismus. Velká z-skóre tam tedy testují spíš konstrukci generátoru než text.
* z-skóre se počítají vůči SD generátoru. Nejistota reálného textu (page bootstrap) je větší, např.
  H2 ± 0,044, hapax ± 0,07. Rozdíly v H2, edge MI, pos JSD a paradigmatech jsou ale mnohonásobně
  větší než obě nejistoty. Bootstrap pro paradigmata, BPE a MI metriky se nepočítal (spec ho
  vyžaduje jen u hapax, n_types, H2 a share_near).
* Generátor netrefuje H2 a v mřížce ani n_types. Část rozdílů v necílových metrikách (hlavně
  paradigmata a BPE) může proto plynout z nekalibrace, ne z „obsahu“. Diagnostický bod tuto obavu
  pro paradigmata zmenšuje (rozdíly zůstávají), pro BPE ji potvrzuje.
* Afixy jsou naučené greedy rozkladem na reálném textu (`decompose_table`). Paradigmatické metriky
  závisí na této segmentaci, kterou spec bere jako danou.
* Diagnostika a `diag_full` jsou mimo zadanou mřížku. `best_configs.csv` obsahuje podle spec
  10 nejlepších konfigurací mřížky.

## Doporučený další krok

1. Rozhodnout, který nulový model mají číst další experimenty. `best_configs.csv` (config 573,
   loss 0,98) odpovídá spec, ale diagnostický bod (p_new 0,5, loss 0,20) je lépe kalibrovaný.
   Navrhuji rozšířit mřížku o `p_new ∈ {0,3; 0,5; 0,8}` a seed slovník.
2. Snížit H2 generovaných slov: editace vážené trigramovým (nebo slotovým, Zattera 2022) modelem místo
   nezávislých bigramů. Ověřit, zda pak zůstane rozdíl v paradigmatech (`frac_stems_multi_suffix`,
   `mean_suffixes_per_stem`, `cmi_excess`). To je přímý vstup pro experiment 002.
3. Přidat do generátoru poziční mechanismus (řádkové prefixy podle Feaster/Currier) a zjistit, kolik
   z `pos_jsd_first_mid` a `edge_mi_excess` tím vysvětlí. Pak teprve tyto metriky brát jako indicie.
