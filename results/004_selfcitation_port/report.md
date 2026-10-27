# 004 – Port generátoru Timm & Schinner (self-citation) jako kanonický kopírovací null

Kód: `voytools/selfcitation.py` (port), `experiments/004_selfcitation_port/run.py` (experiment), deník a výklad
spec: `experiments/004_selfcitation_port/notes.md`. Všechna čísla níže jsou ze souborů v této složce
(`validation.csv`, `validation/identity.csv`, `grid_A.csv`, `grid_B.csv`, `best_configs.csv`, `profile.csv`,
`vocab_matching.csv`, `summary_info.json`, `real/*/point.csv`). Reference 001/002/003 jsou převzaté z jejich
složek v `results/`, trigram 86 z 003 je tu dopočítán s plnou sadou metrik (`trigram86/`).

## Co bylo spočítáno

| Krok | Rozsah | Soubor |
|---|---|---|
| Port | `SelfCitation` + tabulky z `lrozanova/voynich-units` (commit 956a7c4, MIT), logika beze změny; API `selfcitation_generate` / `selfcitation_template` / `selfcitation_like` | `voytools/selfcitation.py`, `tests/test_selfcitation.py` |
| Reál | ZL3b, P, clean, `split` (34 268 tokenů) a `join`; všechny metriky, seed 0 | `real/ZL_*/point.csv` (shoda s 003 na 2e-16, s 002 přesná) |
| 1. Validace | shoda token po tokenu s jejich vzorky; metriky jejich 3 vzorků vs portu (5 seedů, R&T parametry i Java výchozí) vs tabulky 7 | `validation/identity.csv`, `validation.csv` |
| 2. Kalibrace A a B | 60 konfigurací (spec mřížka 48 + mimo spec R&T váhy 0/0/100), 220 vzorků, ZL split | `grid/`, `grid_A.csv`, `grid_B.csv` |
| 2. Fáze B | default (= calA), calB, R&T body × 50 seedů × split/join, všechny metriky vč. 002 | `phaseB/`, `runs.csv`, `summary.csv`, `best_configs.csv` |
| 3. Velikost slovníku | calB: reál vs generátor na šabloně škálované na n_stems_k reálu, 20 seedů | `vocab/`, `vocab_matching.csv` |
| 4. Profil | 23 metrik × 2 sady × 10 zdrojů, z-skóre, tři seznamy | `profile.csv`, `profile.png`, `summary_info.json` |

GC2a se nepoužívá: tabulky generátoru jsou Eva-specifické. Robustnost místo „ZL i GC“ = `split` i `join`.

## Čísla

**1. Validace (kritérium 1 – splněno).** Klíčové zjištění: jejich publikované vzorky
`selfcite_seed{1,2,3}.txt` nejsou vyrobené s Java výchozími parametry (náhrady 1×/2×/3× s vahami 30/50/20), ale
s výsledkem jejich kalibrace **p_copy 0, replace_weights 0/0/100** (vždy trojí náhrada). S těmito parametry a
s jejich šablonou (jejich ZL3b, jejich čištění, prodloužení na 37 500 tokenů) náš port reprodukuje **všechny tři
vzorky token po tokenu** (4 653/4 653 řádků, `validation/identity.csv`). Port tedy není jen „v toleranci“, je
identický s jejich kódem.

Metriky na prvních 34 268 tokenech (`validation.csv`; jejich vzorky: rozsah 3 seedů; port: průměr 5 seedů na
naší ZL split šabloně):

| metrika | jejich vzorky | port, R&T param. | port, Java default | tab. 7 (self-citation) | reál ZL split | tolerance |
|---|---|---|---|---|---|---|
| H2 (bit) | 2,565–2,701 | 2,672 | 2,647 | 2,6–2,7 | 2,689 | ±0,1 |
| hapax | 0,585–0,602 | 0,585 | 0,590 | 0,59–0,60 | 0,676 | ±0,05 |
| D64 (`bpe_gap64`) | 0,660–0,854 | 0,787 | 0,711 | 0,65–0,86 | 1,058 | ±0,1 |
| edge MI excess (bit) | 0,0025–0,0066 | 0,0069 | 0,0057 | 0,006–0,009 | 0,185 | ±0,01 |
| pořadí tokenů (excess/H) | 0,57–0,61 % | 0,52 % | 0,56 % | 0,5–0,6 % | 0,89 % | ±0,3 p. b. |
| near_ratio (naše def.) | 0,98–1,02 | 1,01 | 1,01 | „edit-1“ 1,06–1,23 | 1,18 | – |

Všech pět metrik s tolerancí leží pro oba porty v rozmezí jejich vzorků i tabulky 7 (`within_*_tol` = 1).
Jedinou odchylkou je `near_ratio`: naše definice (okno 5, Levenshtein ≤ 1, míchání v řádcích) dává na jejich
vlastních vzorcích 1,00, tabulka 7 uvádí 1,06–1,23. Je to rozdíl definic metrik, ne chyba portu.

**2. Kalibrace (kritérium 2 – splněno, `provisional = 0`).**

* Kalibrace A (ztráta R&T: délka/0,25, typy/500, hapax/0,03; `grid_A.csv`): nejlepší bod spec mřížky je
  **config 1 = Java výchozí parametry** (p_copy 0, 30/50/20, p_same_position 28), loss_A 10,6. Ztráta je velká
  kvůli hapaxu (složka 8,2, tj. 0,590 proti 0,676) – žádný bod mřížky hapax nedosáhne (nejvýš ~0,59). R&T
  bod 0/0/100 má loss_A 14,4 (šablona 34 268 tokenů; R&T kalibrovali na 32 747 tokenech). calA = default,
  ve fázi B jde o tytéž vzorky.
* Kalibrace B (naše ztráta z 001; `grid_B.csv`): nejlepší **config 21** – p_copy 0,1, replace_weights 20/50/30,
  p_same_position 14 – loss_B 0,039 na 3 seedech.
* Fáze B, 50 seedů, ZL split (`best_configs.csv`):

| zdroj | config | loss_B | share_near | hapax | H2 | n_types | mean_len |
|---|---|---|---|---|---|---|---|
| reál | – | – | 0,173 | 0,676 | 2,689 | 6 796 | 4,99 |
| **ts2020_calB** (doporučený null) | 21 | **0,052** | 0,179 | 0,579 | 2,618 | 6 530 | 5,36 |
| ts2020_calA = ts2020_default | 1 | 0,053 | 0,177 | 0,580 | 2,635 | 6 238 | 5,29 |
| ts2020_RT (mimo spec mřížku) | 13 | 0,044 | 0,173 | 0,588 | 2,672 | 6 183 | 5,24 |
| trigram 86 (003) | – | 0,269 | 0,183 | 0,557 | 3,329 | 8 785 | 4,86 |

calB splňuje kritérium (loss < 0,15, všechny složky < 0,05; největší je hapax 0,020), proto je v prvním řádku
s `provisional = 0`. Splňují ho ale **všechny** tři body generátoru; rozdíly mezi nimi jsou menší než rozptyl
mezi seedy. Nejnižší loss_B má R&T bod (0,044), ten ale není ve spec mřížce. Proti 003 (0,269) je to zlepšení
o řád a poprvé je H2 trefené (2,62 vs 2,69).
**Hapax zůstává mimo**: 0,58 proti 0,68 (z = 8,3 split, 12,5 join). Pod prahem složky je jen proto, že práh 0,05
odpovídá relativní odchylce 22 %.

**Necílové metriky** (`profile.csv`; z = (reál − průměr generátoru) / SD generátoru přes 50 seedů):

| metrika | reál split / join | calB split / join | z calB split / join | trigram 86 split | bigram 573 split |
|---|---|---|---|---|---|
| edge_mi_excess | 0,185 / 0,152 | 0,007 / 0,009 | 40 / 24 | 0,000 | 0,000 |
| token_mi_excess_share | 0,0089 / 0,0071 | 0,0044 / 0,0039 | 2,5 / 1,4 | 0,0001 | −0,0003 |
| bpe_gap64 | 1,058 / 1,065 | 0,692 / 0,714 | 5,7 / 5,1 | 0,759 | 0,601 |
| pos_jsd_first_mid | 0,186 / 0,183 | 0,045 / 0,044 | 22 / 18 | 0,001 | 0,002 |
| cmi_excess | 0,088 / 0,078 | 0,115 / 0,112 | −1,4 / −1,2 | 0,209 | 0,229 |
| frac_stems_multi_suffix | 0,912 / 0,908 | 0,815 / 0,795 | 3,0 / 3,5 | 0,620 | 0,507 |
| mean_suffixes_per_stem | 5,25 / 5,09 | 3,79 / 3,72 | 10,2 / 10,7 | 2,71 | 2,28 |
| frac_multi_suffix_5to9 | 0,832 / 0,832 | 0,711 / 0,685 | 2,7 / 3,4 | 0,479 | 0,349 |
| H_suffix_given_stem | 1,817 / 1,888 | 1,649 / 1,596 | 1,9 / 2,3 | 1,490 | 1,363 |
| mi_prefix_suffix | 0,205 / 0,183 | 0,187 / 0,197 | 0,7 / −0,4 | 0,106 | 0,088 |
| productivity_slope | 0,309 / 0,328 | 0,305 / 0,310 | 0,2 / 0,9 | 0,411 | 0,480 |
| stem_dispersion | 0,459 / 0,474 | 0,331 / 0,339 | 2,7 / 2,0 | 0,309 | 0,005 |
| jaccard_stem_sections | 0,325 / 0,296 | 0,241 / 0,246 | 4,5 / 1,9 | 0,184 | 0,115 |
| jsd_stem_sections | 0,227 / 0,263 | 0,243 / 0,244 | −0,4 / 0,4 | 0,360 | 0,493 |
| jsd_suffix_AB | 0,088 / 0,093 | 0,004 / 0,004 | 33 / 35 | 0,001 | 0,003 |
| n_stems_k (kontrola) | 308 / 335 | 503 / 472 | −3,1 / −1,7 | 736 | 947 |

**Kritérium 3 – tři seznamy** (calB, |z| < 3; `summary_info.json` → `lists`; „v obou“ = split i join):

* (a) **reprodukované v obou:** `token_mi_excess_share`, `cmi_excess`, `H_suffix_given_stem`, `mi_prefix_suffix`,
  `productivity_slope`, `stem_dispersion`, `jsd_stem_sections`.
* (b) **nereprodukované v obou:** `edge_mi_excess`, `bpe_gap64`, `pos_jsd_first_mid`, `mean_suffixes_per_stem`,
  `jsd_suffix_AB`.
* (c) **nekonzistentní split vs. join:** `frac_stems_multi_suffix` (z 3,0 / 3,5), `frac_multi_suffix_5to9`
  (2,7 / 3,4), `jaccard_stem_sections` (4,5 / 1,9) – všechny na hraně prahu.

Pro default a R&T bod vycházejí seznamy skoro stejně (`summary_info.json`). Pro trigram 86 z 003 je seznam (a)
prázdný a všech 15 metrik je v (b).

Pozor na seznam (a): z se počítá vůči SD generátoru a ta je u tohoto „rich-get-richer“ procesu velká.
Například `stem_dispersion` má SD 0,048 (trigram 86: 0,015). Průměr generátoru je přitom pořád jasně pod reálem
(0,33 vs 0,46), podobně `token_mi_excess_share` (polovina reálu). „Reprodukováno“ tu znamená „reál leží
v rozptylu jednotlivých běhů generátoru“, ne „generátor má v průměru stejnou hodnotu“. Bez výhrad sedí
`cmi_excess`, `mi_prefix_suffix`, `productivity_slope` a `jsd_stem_sections`.

**Velikost slovníku (krok 3, `vocab_matching.csv`).** Generátor má **víc** kmenů s ≥ 5 tokeny než reál
(503 vs 308), ne méně. Varianta (a) ze spec (zmenšit reál podvzorkováním stránek na n_stems_k generátoru) proto
nejde: podvzorkování n_stems_k jen dál snižuje. Byla vynechána (sloupce `real_sub*` = NaN). Varianta (b) je
obrácená: šablona zkrácená na prvních 46 % tokenů (`gen_scaled_scale` 0,462) dává n_stems_k 295 ± 31 (reál 308).
Rozdíly reál − generátor se tím **nezmenší**: mean_suffixes_per_stem 1,46 → 1,53, frac_stems_multi_suffix
0,097 → 0,109, cmi_excess −0,028 → −0,030, H_suffix_given_stem 0,17 → 0,21. U sekčních metrik
(stem_dispersion 0,13 → 0,18, jsd_stem_sections −0,02 → −0,10) zkrácení mění i skladbu sekcí (prvních 46 %
stránek je hlavně herbář), takže tam srovnání neplatí. Paradigmatické rozdíly z 002 tedy u self-citation
generátoru **nejsou** artefaktem velikosti slovníku. Platí to i pro 001/003: tam má generátor 2–3× víc kmenů a
rozdíly jsou větší, ale stejného směru.

**Kritérium 4.** `python -m pytest -q`: 195 passed, 1 skipped na Pythonu 3.13 / pandas 3.0 i 3.10.20 /
pandas 2.3.3. Regresní test `autocopy_generate` je v `tests/test_selfcitation.py::test_autocopy_generate_is_untouched`
(a dál v `test_wordgrammar.py`). Bez scipy.

**Kritérium 5 (podezřelé).** Validace je v toleranci u všech pěti metrik. Žádný běh fáze B nemá H2 < 2,0 ani
share_near > 0,9 (`summary_info.json`).

## Interpretace

Port je s kódem Rozanové a Temereva identický a reprodukuje jejich čísla. Mechanismus Timma a Schinnera trefí
H2, share_near, počet typů i délku tokenů najednou, což naše generátory 001/003 neuměly. Ani on ale nedá dost
hapaxů (0,58 vs 0,68) a nevytváří vazbu přes hranici slov (edge MI 0,007 vs 0,185), A/B rozdíl v sufixech, BPE
žlab ani plnou šířku paradigmat (3,8 vs 5,3 sufixu na kmen).

Pro výzkumnou otázku je podstatné, že **predikce ze spec se u paradigmat nepotvrdila**. Self-citation generátor
reprodukuje závislost prefix–sufix daná kmenem (`cmi_excess`), MI prefix–sufix i sklon produktivity. `stem_dispersion` je 0,33 (reál 0,46). Bigramové kopírování
z 001 dává 0,005, trigram 86 z 003 podobných 0,31, ale ten paradigmata nereprodukuje. Tyto vlastnosti 002 tedy samy
o sobě neodlišují text od kopírovacího procesu, který má bohaté glyfové náhrady a pravidla pro začátky řádků.
Zbývající rozdíly (sufixů na kmen, edge MI, A/B, poziční JSD 0,045 vs 0,186) jsou nyní nejsilnější kandidáti.

## Co výsledek neukazuje

* Neukazuje, že text nese význam, ani že nebyl vygenerován. Ukazuje jen, které vlastnosti **tento** generátor
  s **těmito** parametry nereprodukuje. Edge MI, A/B rozdíl a poziční profil nemá generátor jak vytvořit
  (žádný mechanismus mezislovní vazby, žádné rozlišení A/B, pravidla řádků jen pro začátek). Velká z tam testují
  hlavně jeho konstrukci.
* „Reprodukováno“ (|z| < 3) je vůči velkému rozptylu mezi seedy. Průměry se u `stem_dispersion` a
  `token_mi_excess_share` od reálu stále liší.
* Kalibrace pokrývá jen p_copy, váhy náhrad a p_same_position. Ostatní parametry (p_add_remove,
  p_combine_split, tabulky podobnosti glyfů) zůstaly na Java hodnotách. Hapax by možná šel trefit jinými
  parametry nebo delší šablonou (R&T kalibrovali na 32 747 tokenech).
* Generátor je Eva-only. Robustnost vůči abecedě (GC) testovat nejde, místo ní je tu jen split vs join.
* Tabulky glyfových podobností a gramatika „curve/line“ jsou ručně odvozené z rukopisu (Timm & Schinner).
  Generátor tedy do slov vnáší znalost tvarů glyfů Voynichu. To je pro null výhoda, ne handicap.
* Krok 3 je jen varianta (b), (a) nešla (viz výše). Zkrácená šablona mění skladbu sekcí.

## Doporučený další krok

1. Nahradit provizorní null z 001/003 konfigurací **ts2020_calB** (`best_configs.csv`, config 21) a přepočítat
   paradigmatické srovnání 002 a experimenty 005–007 proti ní. Vždy ale uvádět, že hapax netrefuje.
2. Hledat paradigmatické metriky, které self-citation nevysvětlí: počet sufixů na kmen při stejné frekvenci kmene,
   A/B rozdíly v afixech uvnitř kmene, vazba konce slova na začátek dalšího (edge MI) podmíněná kmenem.
3. Rozšířit generátor o rozlišení A/B (samostatné zdrojové řádky / seed line podle jazyka) a zjistit, kolik
   z `jsd_suffix_AB` tím vznikne. Jde o test, zda je A/B rozdíl vysvětlitelný jen odlišným „startem“ kopírování.
