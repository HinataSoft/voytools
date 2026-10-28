# 018 – Lokalita slovníku: drift, nebo bloky? A fyzické pořadí stránek

Větev `exp/018-locality-drift-vs-blocks`; kód `experiments/018_locality_drift_vs_blocks/run.py` + `summarize.py`,
funkce `voytools/locality.py` (testy `tests/test_locality.py`), rozšíření `voytools/selfcitation.py` o paměť přes
stránky (`p_prev_page`, `page_memory`; výchozí hodnoty dávají bitově shodný výstup, test
`test_memory_defaults_are_bit_identical`). Volby při výkladu spec jsou v `experiments/018_locality_drift_vs_blocks/notes.md`.
Čísla níže jsou z `similarity_profile.csv`, `memory_calibration.csv`, `segmentation_gain.csv`, `crossval_blocks.csv`,
`changepoints.csv` / `changepoints_list.csv`, `physical_order.csv`, `neighbours.csv` / `neighbours_list.csv`,
`quire_matrix.csv` a `summary.csv` (3 580 běhů, `meta.json`).

## Co se měřilo

Stránka = folio s čistým `P` textem (ZL i GC 207 stránek), vzdálenost = pořadí v souboru. Hlavní míra podobnosti
(a) = sdílení **vzácných kmenů** (2–5 výskytů, dekompozice a): pro dvojici stránek pozorovaný počet dvojic výskytů
téhož kmene a jeho očekávání při permutaci tokenů na stránky s vahou ∝ velikosti (015); podobnost množiny dvojic
= Σ pozorované / Σ očekávané (= poměr k náhodě z 015). Míry (b) kosinus log(1+f) bez 30 funkčních kmenů (006) a
(c) Jaccard kmenů ≥ 2 výskytů jsou doplňkové (nesou závislost na velikosti stránky, kterou sdílejí všechny nully
v témž layoutu). Hladký model f(d): sdružený poměr po binech vzdálenosti 1…40 (+ 41+), izotonicky nerostoucí;
reziduum R = N − f(d)·W. Bloky: optimální segmentace na K = 1…30 souvislých bloků (DP) podle jednostranné deviance
bloku (Poisson pro (a), Gauss pro (b, c); doslovný součet reziduí je degenerovaný – viz notes). **Zisk bloků** =
sdružené reziduum uvnitř bloků, křížově ověřený: bloky na jedné polovině kmenů, zisk na druhé polovině proti
jejímu vlastnímu binovanému profilu (očekávání 0 pro bloky nezávislé na druhé polovině). Pro (a) je zisk
v jednotkách poměru k náhodě (0,10 = uvnitř bloků o 10 % víc sdílených dvojic, než dá vzdálenost). Reál: 50 rozdělení
(řádek = rozdělení; interval = kvantily přes rozdělení, `se` = směrodatná chyba průměru), nully a generátory:
5 rozdělení na vzorek, interval přes vzorky (200 permutací stránek, 30 generátorů). Výřezy `all`, `within` (jen
dvojice v téže sekci), `herbalA` (H & jazyk A, 95 stránek); K_ref = počet souvislých úseků sekce (24; herbář A:
kvaternů, 9).

## 1. Profil vzdálenosti (`similarity_profile.csv`, `distance_profile.png`)

Sdílení vzácných kmenů, poměr pozorované/očekávané podle vzdálenosti d (výřez `all`):

| zdroj | 1 | 2 | 3 | 5 | 8 | 12 | 16 | 24 | 32 | 41+ |
|---|---|---|---|---|---|---|---|---|---|---|
| **ZL split, reál** | 1,97 | 1,75 | 1,44 | 1,55 | 0,70 | 0,93 | 0,99 | 1,31 | 0,74 | 0,86 |
| **GC split, reál** | 1,73 | 1,57 | 1,30 | 1,23 | 1,60 | 1,55 | 1,29 | 1,24 | 0,88 | 0,85 |
| ZL join / GC join, reál | 2,25 / 1,66 | 1,86 / 1,19 | 1,41 / 1,30 | 1,55 / 1,46 | 0,96 / 1,75 | 0,89 / 1,40 | 1,16 / 1,55 | 1,07 / 1,38 | 0,81 / 0,53 | 0,86 / 0,86 |
| permutace stránek globálně (ZL) | 0,97 | 1,00 | 1,01 | 0,98 | 1,00 | 0,99 | 1,03 | 0,96 | 0,98 | 0,98 |
| permutace stránek v sekci (ZL) | 1,36 | 1,35 | 1,36 | 1,28 | 1,23 | 1,18 | 1,12 | 1,01 | 0,93 | 0,90 |
| ts2020 cfg 21 (ZL) | 1,56 | 0,81 | 0,76 | 0,72 | 0,73 | 0,74 | 0,79 | 0,71 | 0,73 | 0,68 |
| **ts2020-mem** (0,3 / 8; ZL) | 1,92 | 1,60 | 1,49 | 1,34 | 1,25 | 1,08 | 0,99 | 0,88 | 0,71 | 0,71 |
| šablona 005 s tématem (ZL) | 0,78 | 0,76 | 0,79 | 0,86 | 0,79 | 0,71 | 0,81 | 0,81 | 0,78 | 0,79 |
| Plinius (knihy) | 3,90 | 2,23 | 1,86 | 1,71 | 1,23 | 1,05 | 0,76 | 0,82 | 0,94 | 0,75 |
| Celsus (knihy) | 4,11 | 2,54 | 2,35 | 1,96 | 1,65 | 1,00 | 0,93 | 0,82 | 0,93 | 0,66 |
| Linnaeus po položkách | 4,67 | 1,70 | 1,04 | 1,09 | 0,89 | 0,64 | 0,46 | 0,91 | 0,81 | 0,55 |
| blábol (60 stránek) | 4,54 | 1,03 | 0,61 | 0,58 | 0,25 | 0,39 | 0,52 | 0,22 | 0,16 | 0,28 |

Uvnitř sekcí (`within`): ZL split reál 2,13 / 1,84 / 1,48 / 1,72 / 0,77 / 1,21 / 0,98 (d = 1, 2, 3, 5, 8, 12, 16),
GC split 1,79 / 1,57 / 1,26 / 1,21 / 1,47 / 2,12 / 1,66; ts2020-mem 1,90 / 1,54 / 1,50 / 1,34 / 1,30 / 1,07 / 1,14;
permutace v sekci plochá 1,41 / 1,43 / 1,43 / 1,38 / 1,37 / 1,42 / 1,41. Herbář A: ZL split 2,23 / 1,34 / 1,67 /
1,44 / 0,71 / 1,02 / 0,99, GC join 1,77 / 1,36 / 1,40 / 1,10 / 1,07 / 1,47 / 2,10; ts2020-mem 1,98 / 1,55 / 1,35 /
1,33 / 1,24 / 0,93 / 0,97; permutace uvnitř sekce 1,00 / 0,98 / 0,99 / 0,96 / 0,97 / 1,01 / 1,04.

Výklad: profil 015 se reprodukuje (d = 1: 1,97; pokles do ~8–15 stránek), platí i uvnitř sekcí a uvnitř herbáře A
(kde permutace uvnitř sekce dá přesně 1). ts2020 s pamětí přes stránky dává stejný hladký pokles ve všech třech
výřezech – z konstrukce, bez znalosti sekcí.

## 2. ts2020 s pamětí přes stránky (`memory_calibration.csv`)

Mřížka p_prev_page ∈ {0,1; 0,2; 0,3; 0,5} × page_memory ∈ {2, 4, 8, 16}, 5 seedů, ztráta = Σ (poměr binu − reál)²
přes biny 1 / 2–3 / 4–7 / 8–15 profilu 015. Vítěz **p_prev_page = 0,3, page_memory = 8** (ztráta 0,036; druhý 0,5/16:
0,21; cfg 21 bez paměti: 1,17). Vítěz, 30 seedů [2,5–97,5 %]:

| zdroj | 0 | 1 | 2–3 | 4–7 | 8–15 | 16–31 | 32+ | `x_far_nl` | přebytek ≥ 16 v sekci | v písaři | tvar 0 / 1–3 / 4–15 / 16+ | loss B (004) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| reál ZL split (015) | 2,90 | 1,97 | 1,60 | 1,38 | 1,24 | 1,02 | 0,87 | −0,053 | **+0,081** | +0,068 | −0,02 / +0,01 / +0,01 / 0,00 | 0,0 |
| ts2020 cfg 21 | 35,2 | 1,56 [0,82; 2,17] | 0,79 | 0,72 | 0,75 [0,59; 0,84] | 0,75 | 0,68 | −0,008 | +0,014 [−0,023; 0,056] | +0,004 | +0,08 / −0,02 / −0,04 / −0,04 | 0,049 |
| **ts2020-mem** | 18,6 | 1,93 [1,27; 2,35] | 1,55 | 1,36 | 1,14 [0,94; 1,31] | 0,92 | 0,72 | −0,074 | +0,008 [−0,029; 0,075] | +0,002 | +0,09 / +0,03 / +0,01 / −0,03 | **0,114** (share_near 0,141 vs 0,173; hapax 0,58 vs 0,68) |

Výklad: paměť přes stránky dožene **celý pokles 1–15 stránek** (biny v intervalu reálu) a lokalitu uvnitř sekcí
i herbáře A (bod 1), ale **ne přebytek dálkových dvojic v téže sekci** (+0,008 proti +0,081, interval reál
nedosahuje) ani tvar opakování (kopie na téže stránce +0,09 nad paradigma jako cfg 21); vrchol na stránce zůstává
6× vyšší než reál (18,6 vs 2,9). Ztráta 004 vzroste z 0,049 na 0,114 (> 0,1): kopírování z dřívějších stránek
sníží podíl blízkých sousedů v řádku (share_near 0,141 proti 0,173 reálu) – generátor s pamětí už není kalibrovaný
na cíle 004 (uvedeno podle spec).

## 3. Zisk bloků nad model vzdálenosti (`segmentation_gain.csv`, `crossval_blocks.csv`, `segmentation_gain.png`)

Křížově ověřený zisk bloků, vzácné kmeny; průměr [interval], u reálu navíc `se` (`crossval_blocks.csv`):

| výřez / zdroj | K = 2 | K = 8 | K = 24 (K_ref) | K = 30 |
|---|---|---|---|---|
| **all, ZL split, reál** | 0,063 [−0,00; 0,12] se 0,005 | 0,093 [−0,07; 0,32] se 0,015 | **0,118** [−0,09; 0,38] se 0,021 | 0,150 se 0,026 |
| all, ZL join / GC split / GC join, reál | 0,073 / 0,024 / 0,055 | 0,134 / 0,071 / 0,066 | 0,131 / **0,324** / 0,133 | 0,152 / 0,459 / 0,157 |
| all, permutace globálně (ZL) | −0,002 [−0,04; 0,04] | −0,006 [−0,08; 0,06] | 0,011 [−0,12; 0,15] | 0,020 [−0,12; 0,17] |
| all, permutace v sekci (ZL) | 0,036 [0,00; 0,08] | 0,167 [0,01; 0,31] | **0,215** [−0,02; 0,51] | 0,210 |
| all, ts2020 cfg 21 (ZL) | 0,002 | −0,004 | 0,013 [−0,07; 0,13] | 0,033 |
| all, **ts2020-mem** (ZL) | 0,003 | 0,013 | **0,049** [−0,12; 0,40] | 0,067 [−0,11; 0,44] |
| all, šablona 005 s tématem (ZL / GC) | 0,013 / 0,004 | 0,012 / 0,018 | 0,064 / 0,051 | 0,086 / 0,085 |
| all, Plinius (K_ref = 7) | 0,028 | 0,096 [0,03; 0,25] | 0,211 | 0,272 |
| all, Plinius, permutace v knize | 0,028 | 0,267 | 0,591 [0,41; 0,77] | 0,658 |
| all, Celsus / Linnaeus | 0,010 / 0,014 | 0,375 / −0,014 | 0,450 / −0,011 | 0,478 / 0,004 |
| **within, ZL split, reál** | 0,004 | 0,000 | 0,032 [−0,21; 0,30] se 0,022 | 0,071 se 0,029 |
| within, ZL join / GC split / GC join, reál | −0,006 / −0,019 / 0,001 | −0,022 / 0,024 / −0,029 | −0,037 / 0,226 / 0,005 | 0,032 / 0,311 / 0,006 |
| within, permutace globálně / v sekci (ZL) | −0,002 / 0,005 | −0,001 / 0,031 | 0,011 [−0,15; 0,16] / 0,065 | 0,021 / 0,089 |
| within, ts2020-mem / šablona (ZL) | −0,002 / −0,001 | −0,015 / 0,006 | 0,006 / 0,066 | 0,049 / 0,086 |
| **herbalA, ZL split, reál** (K = 2 / 9 / 18 / 27) | 0,088 [0,01; 0,17] | 0,191 se 0,029 | **0,432** [−0,24; 1,02] se 0,059 | 0,535 se 0,064 |
| herbalA, ZL join / GC split / GC join, reál | 0,073 / 0,050 / 0,053 | 0,315 / 0,066 / 0,322 | 0,646 / 0,079 / 0,682 | 0,755 / 0,072 / 0,813 |
| herbalA, permutace globálně / v sekci (ZL) | −0,001 / −0,000 | 0,004 / −0,003 | 0,024 [−0,16; 0,25] / −0,003 | 0,035 / 0,004 |
| herbalA, ts2020 / ts2020-mem / šablona (ZL) | 0,004 / 0,005 / −0,001 | −0,006 / 0,010 / 0,001 | −0,022 / 0,015 / 0,032 | −0,004 / 0,026 / 0,042 |

Kosinus (výřez all, K = 24): reál ZL 0,019 [−0,00; 0,04], permutace globálně 0,002 [−0,004; 0,008], **permutace
v sekci 0,019**, ts2020 0,022, ts2020-mem 0,028, šablona 0,020, Plinius 0,020, Linnaeus 0,032 – u kosinu a Jaccardu
dává stejný zisk každý zdroj v témž layoutu (velikost stránek a jazyk A/B podle pozice), nerozlišují nic; dál jen (a).

Výklad: (i) V celém rukopisu je zisk bloků nad vzdálenost kladný a nad permutaci stránek (ZL split 0,12 ± 0,02
proti 0,01 [−0,12; 0,15] při K = 24; ve všech čtyřech kombinacích; při K = 2 je první hranice f70r2 | f75r, tj.
konec herbáře/astro a začátek balneo – bloky = sekce a jazyk). ts2020-mem dá 0,05, šablona s tématem 0,06–0,09,
ts2020 0,01. Permutace stránek uvnitř sekcí má zisk **vyšší než reál** (0,22): když se zruší pozvolný pokles uvnitř
sekce, všechna sekční podobnost zůstane blokům; v reálu část sekční podobnosti pohltí f(d). Plinius s knihami je
stejný případ (0,21 proti 0,59 po permutaci v knize). (ii) **Uvnitř sekcí** zisk nad permutaci není (ZL split 0,03
± 0,02, ZL join −0,04, GC join 0,01 při K = 24; jen GC split 0,23 ± 0,03 – nepotvrzeno v join). (iii) **V herbáři A**
je zisk kladný a 10× nad všemi nully a generátory ve třech ze čtyř kombinací (ZL split 0,43 ± 0,06, ZL join 0,65,
GC join 0,68 při K = 18; permutace 0,02 [−0,16; 0,25], ts2020-mem 0,02, šablona 0,03); GC split dá jen 0,08 ± 0,03.
Hodnoty jsou v jednotkách poměru: uvnitř nalezených bloků herbáře A je o 40–70 % víc sdílených vzácných kmenů,
než dá vzdálenost.

## 4. Body změny (`changepoints.csv`, `changepoints_list.csv`)

Přesnost hranic (± 1 stránka) při K = K_ref, výřez all, vzácné kmeny; náhodná základna = náhodné hranice:

| zdroj | sekce: přesnost [perm. interval] / náhoda | úplnost | kvatern: přesnost / náhoda | jazyk (K = 2 K_ref): přesnost / náhoda |
|---|---|---|---|---|
| **ZL split, reál** | 0,30 / 0,21 | 0,43 / 0,30 | 0,35 / 0,21 | 0,41 / 0,36 |
| ZL join / GC split / GC join, reál | 0,35 / 0,35 / 0,26 | 0,48 / 0,39 / 0,35 | 0,30 / 0,39 / 0,22 | 0,34 / 0,48 / 0,41 |
| permutace globálně (ZL) | 0,20 [0,04; 0,35] | 0,29 | 0,21 [0,04; 0,35] | 0,35 |
| permutace v sekci (ZL) | 0,35 [0,17; 0,57] | 0,45 | 0,23 | 0,27 |
| ts2020 / ts2020-mem / šablona (ZL) | 0,20 / 0,24 / 0,24 | 0,28 / 0,36 / 0,36 | 0,22 / 0,23 / 0,25 | 0,36 / 0,37 / 0,32 |
| kosinus, ZL split reál | **0,48** [perm 0,04; 0,35] | 0,52 | 0,30 | 0,52 |
| Plinius (knihy, K = 7) | **0,33** [perm 0,00; 0,33] / 0,09 | 0,33 | 0,17 / 0,22 | – |
| Celsus (knihy, K = 5) | **0,50** [0,00; 0,25] / 0,05 | 0,50 | 0,00 | – |
| Linnaeus (řády, K = 4) | 0,00 / 0,04 | 0,00 | 0,33 | – |
| blábol (pisatelé, K = 30) | 1,00 / 0,86 | 0,97 | 0,66 / 0,64 | – |

Herbář A (K = 9, hranice kvaternů): ZL split a ZL join 0,50 (úplnost 0,38) proti náhodě 0,27 a permutaci
[0,00; 0,50]; GC 0,12; ts2020-mem 0,28, šablona 0,28. Hranice reálu v ZL split (vzácné kmeny, K = 24): 7 z 23 sedí na
hranici sekce (f68v2|f69r, f70r2|f75r, f75v|f76r, f86v3|f87r, f87v|f88r, f89v1|f90r1, f96v|f99r), 8 na hranici
kvaternu, zbytek uvnitř herbáře A (f9r|f9v, f13r|f13v, f38r|f38v, f49r|f49v …) a uvnitř Q20 (f104v|f105r,
f107v|f108r, f114r|f114v); kosinus: 11 z 23 na sekci (f56v|f57r, f67r2|f67v2, f70r2|f75r, f84v|f85r1, f102v1|f103r …).

Výklad: body změny vzácných kmenů sedí na hranice sekcí jen na horním okraji permutačního intervalu (0,30–0,35
proti [0,04; 0,35]); kosinus v ZL split je nad ním (0,48), ale permutace v sekci dává totéž (0,48) – jde o skok
sekce/jazyka, ne o bloky uvnitř. Pozitivní kontrola funguje: knihy Plinia a Celsa se najdou (0,33 a 0,50 proti
náhodě 0,05–0,09), Linnaeus po položkách ne. Uvnitř herbáře A sedí hranice bloků na kvaterny v ZL (0,50 proti 0,27),
v GC ne.

## 5. Fyzické pořadí (`physical_order.csv`, `quire_matrix.csv`, `quire_matrix.png`)

Přebytek rezidua dvojic na témž bifoliu / v témž kvaternu proti ostatním dvojicím na stejné vzdálenosti (d ≤ 15;
(a) v jednotkách poměru, kosinus v jednotkách kosinu):

| zdroj | bifolio, (a) | bifolio, kosinus | kvatern, (a) | kvatern, kosinus | pořadí kvaternů: nejpodobnější = soused | ρ spektrální seriace vs čtení |
|---|---|---|---|---|---|---|
| **ZL split, reál** (337 / 1 430 dvojic) | +0,36 | **+0,046** | +0,35 | +0,062 | 0,00 (kos.) / 0,13 (a) | **0,76** (kos.) / 0,44 (a) |
| ZL join / GC split / GC join, reál | +0,35 / +0,45 / +0,70 | +0,043 / +0,029 / +0,033 | +0,33 / +0,10 / +0,29 | +0,054 / +0,067 / +0,062 | 0,06 / 0,00 / 0,00 | 0,72 / 0,71 / 0,63 |
| permutace globálně (ZL) | +0,01 [−0,49; 0,46] | −0,000 [−0,008; 0,010] | −0,00 [−0,53; 0,40] | −0,000 [−0,014; 0,012] | 0,11 [0,00; 0,25] | 0,19 [0,01; 0,54] |
| permutace v sekci (ZL) | +0,09 [−0,35; 0,53] | +0,009 [0,001; 0,019] | +0,21 [−0,21; 0,52] | +0,048 [0,037; 0,058] | 0,03 | 0,28 [0,02; 0,61] |
| ts2020 / ts2020-mem / šablona (ZL) | −0,03 / +0,05 / +0,07 | −0,004 / −0,007 / +0,004 | −0,14 / −0,20 / −0,16 | +0,026 / +0,031 / +0,024 | 0,06 / 0,10 / 0,06 | 0,31 / 0,30 / 0,28 |
| Plinius / Celsus / Linnaeus | −0,00 / −0,19 / −0,03 | −0,001 / −0,010 / −0,015 | −0,10 / −0,17 / −0,44 | +0,015 / +0,014 / +0,009 | 0,50 / 0,13 / 0,31 | 0,03 / 0,58 / 0,43 |

Výklad: **dvojice stránek na témž bifoliu jsou si podobnější než jiné dvojice ve stejné vzdálenosti čtení**:
kosinus +0,029 až +0,046 ve všech čtyřech kombinacích proti permutaci globálně [−0,008; +0,010], permutaci v sekci
[+0,001; +0,019] a generátorům (≤ +0,004); u vzácných kmenů stejný směr (+0,35 až +0,70), ale interval permutace je
široký (významně jen GC join). Latinské kontroly v témž layoutu dávají 0 (bifolio v lineárním textu nic neznamená).
Kvatern: kosinus +0,05 až +0,07 jen těsně nad permutací v sekci (+0,048 [0,037; 0,058]) – skok na hranici kvaternu
nad rámec sekce je nanejvýš malý. Spektrální pořadí kvaternů (z matice reziduí) koreluje s pořadím čtení
(|ρ| 0,63–0,76 kosinus proti permutacím [0,01; 0,61]); u (a) 0,32–0,57 uvnitř intervalu.

## 6. Sousedé (`neighbours.csv`, `neighbours_list.csv`)

Tři nejpodobnější stránky téhož jazyka; podíl sousedů ve vzdálenosti ≤ 5 / > 16, počet stránek se všemi třemi
sousedy > 8 stránek daleko (ze 197 stránek s ≥ 3 partnery):

| zdroj | reziduum (a): ≤ 5 / > 16 / anomálie | reziduum kosinus | surový kosinus: ≤ 5 / > 16 / anomálie |
|---|---|---|---|
| **ZL split, reál** | 0,074 / 0,83 / 139 | 0,090 / 0,73 / 140 | **0,220 / 0,57 / 88** |
| ZL join / GC split / GC join, reál | 0,083 / 0,095 / 0,112 | 0,086 / 0,076 / 0,085 | 0,227 / 0,223 / 0,245 (anomálie 84 / 83 / 76) |
| permutace globálně (ZL) | 0,067 [0,05; 0,09] / 0,82 / 144 [131; 156] | 0,072 [0,05; 0,10] / 0,80 / 141 | 0,076 [0,05; 0,10] / 0,80 / 138 [123; 153] |
| ts2020-mem (ZL) | 0,086 / 0,77 / 135 | 0,092 / 0,75 / 132 | – |
| Plinius / Celsus / Linnaeus | 0,178 / 0,179 / 0,146 | 0,162 / 0,166 / 0,279 | 0,443 / 0,425 / 0,699 (anomálie 30 / 45 / 8) |

Výklad: po odečtení f(d) nejsou nejpodobnější stránky blíž než při permutaci (0,07–0,11 proti [0,05; 0,10]) –
lokalita je v hladké složce, ne v jednotlivých „párech“ stránek. Surovou podobností je 22 % sousedů do 5 stránek
(permutace 8 %), ale **88 stránek ze 197 má všechny tři sousedy > 8 stránek daleko** (permutace 138 [123; 153];
Plinius 30, Celsus 45, Linnaeus 8). Seznam (`neighbours_list.csv`, `cos_raw`, `anomalous = 1`): převážně herbář A
v Q1–Q7 (f1v–f55v; sousedé jiné stránky herbáře A o 15–100 stránek dál), herbářové B listy f26, f31, f33, f39,
f40, f41v, f43, f46r, f48, f50, f55 (sousedé v Q13–Q20: f111r, f108v, f106v, f115v), f66r, f75r, f84v, f85–f86 (astro/kosmo), f94–f95 (herbář B v Q17), recepty f103r–f106v (sousedé f111–
f115 10–18 stránek dál) a f116r. Herbářové B stránky uprostřed herbáře A mají sousedy v B sekcích na konci
rukopisu – to je známá hypotéza o přeskládání (herbář B v herbáři A), ale zde jde jen o jazyk A/B: sousedé se
hledají v témž jazyce, a B stránky v Q4–Q7 nemají poblíž žádnou jinou B stránku.

## Kontroly a kritéria přijetí

* Kritérium 1: vše pro ZL i GC, split i join (všechny body), tři výřezy; ts2020 a ts2020-mem jen ZL; šablona split.
* Kritérium 3 (podezřelé hodnoty): zisk bloků permutace stránek globálně je ≈ 0 ve všech výřezech (−0,006 až
  +0,02, intervaly obsahují 0; `segmentation_gain.csv`); ts2020-mem s p_prev_page = 0 je bitově shodný s ts2020
  (test); Plinius má body změny na hranicích knih (0,33 proti náhodě 0,09; Celsus 0,50 proti 0,05).
* Kosinus a Jaccard: zisk bloků stejný u reálu i u všech nullů v témž layoutu (velikost stránek) – použity jen pro
  fyzické pořadí, kde permutace dávají 0.
* Neurčitost zisku (a): intervaly přes rozdělení kmenů jsou široké (jedna polovina vzácných kmenů ≈ 240 kmenů);
  závěry stojí na průměru 50 rozdělení (se 0,02–0,06) proti rozdělení 200 permutací.

## Odpovědi (kritérium 2)

| otázka | ZL | GC |
|---|---|---|
| (i) reprodukuje ts2020-mem profil 1–15 a zároveň přebytek v téže sekci a lokalitu uvnitř sekcí? | profil 1–15 **ano** (1,93 / 1,55 / 1,36 / 1,14 proti 1,97 / 1,60 / 1,38 / 1,24), lokalita uvnitř sekcí a herbáře A **ano**; přebytek ≥ 16 v sekci **ne** (+0,008 proti +0,081); tvar opakování ne; loss 004 0,11 (> 0,1) | – (jen Eva) |
| (ii) křížově ověřený zisk bloků nad ts2020-mem? | celý rukopis **ano, malý** (0,12 ± 0,02 proti 0,05; Plinius s knihami 0,21); uvnitř sekcí **ne** (0,03 ± 0,02); herbář A **ano, velký** (0,43 ± 0,06 proti 0,02) | celý 0,32 / 0,13 (split / join); uvnitř sekcí 0,23 / 0,01; herbář A 0,08 / 0,68 |
| (iii) sedí body změny na hranice sekcí/kvaternů nad náhodu? | sekce: jen na okraji permutačního intervalu (0,30–0,35 vs [0,04; 0,35]; kosinus 0,48 = skok sekce/jazyka jako permutace v sekci); kvaterny uvnitř herbáře A **ano** (0,50 vs 0,27, interval [0; 0,50]) | sekce 0,26–0,35; kvaterny v herbáři A ne (0,12) |
| (iv) sleduje podobnost pořadí čtení, nebo bifolia/kvaterny; anomálie? | **bifolio nad pořadí čtení** (kosinus +0,046 vs perm. ≤ +0,019); kvatern ≈ sekce; pořadí kvaternů sleduje čtení (ρ 0,76); 88 stránek s cizími sousedy (herbář A v Q1–Q7, herbářové B listy, recepty) | bifolio +0,029 / +0,033; ρ 0,71 / 0,63; 83 / 76 anomálií |

## Shrnutí

**Čísla:** profil vzácných kmenů 1,97 / 1,75 / 1,44 (d = 1, 2, 3) → 0,86 (41+), uvnitř sekcí 2,13 / 1,84 / 1,48,
v herbáři A 2,23 / 1,34 / 1,67; ts2020-mem (0,3 / 8) 1,92 / 1,60 / 1,49 → 0,71, biny 015 1,93 / 1,55 / 1,36 / 1,14
(reál 1,97 / 1,60 / 1,38 / 1,24), přebytek v sekci +0,008 (reál +0,081), loss 004 0,11. Zisk bloků (K = 24): celý
rukopis ZL 0,12 ± 0,02 (GC 0,32 / 0,13), permutace 0,01 [−0,12; 0,15], permutace v sekci 0,22, ts2020-mem 0,05,
šablona 0,06, Plinius 0,21; uvnitř sekcí 0,03 ± 0,02 (GC 0,23 / 0,01); herbář A (K = 18) 0,43 ± 0,06 (ZL join
0,65, GC join 0,68, GC split 0,08), nully ≤ 0,04. Body změny na sekcích 0,30–0,35 vs permutace [0,04; 0,35], Plinius
0,33 vs 0,09; kvaterny v herbáři A 0,50 vs 0,27. Bifolio: kosinus +0,029 až +0,046 vs permutace [−0,008; +0,010] a
v sekci [+0,001; +0,019]; ρ pořadí kvaternů 0,63–0,76 vs [0,01; 0,61]; anomální stránky 76–88 vs 138 [123; 153].

**Interpretace:** Lokalita slovníku má dvě složky. Hladký pokles přes 1–15 stránek, včetně poklesu uvnitř sekcí
a uvnitř herbáře A, dá kopírování s pamětí přes stránky z konstrukce (p_prev_page 0,3, střední dosah 8 stránek) –
sám o sobě drift nerozliší od tématu. Nad tento model zůstávají **bloky**: v celém rukopisu malé a odpovídající
sekcím a jazyku A/B (zisk 0,12, body změny f70r2|f75r, f86v3|f87r, f96v|f99r …), uvnitř sekcí obecně žádné, ale
**uvnitř herbáře A velké** (o 40–70 % víc sdílených vzácných kmenů v blocích než dá vzdálenost, 10× nad všemi nully,
ve 3 ze 4 kombinací), s hranicemi sedícími na kvaterny v ZL. K tomu je podobnost vázána na **fyzický list**: stránky
téhož bifolia si jsou podobnější než stránky ve stejné vzdálenosti čtení (ve všech čtyřech kombinacích, u latinských
kontrol v témž layoutu nula), což paměťové kopírování nedává. Věta: **drift i bloky** – hladká složka slučitelná
s pamětí přes stránky, blokovou složku (sekce/jazyk globálně, kvaternové bloky v herbáři A, vazba na bifolio) paměť
nedá; podobnost sleduje spíš fyzickou jednotku (bifolio) než čistě pořadí čtení.

**Co výsledek neukazuje:** nerozliší, zda bloky herbáře A jsou témata (skupiny rostlin), dávky psaní nebo
samostatné sešity (kvaterny) – zisk a hranice na kvaternech jsou s obojím slučitelné; GC split dá v herbáři A
jen 0,08 (tři ze čtyř kombinací nestačí na „platí v obou transliteracích“ bez výhrady); přebytek dálkových dvojic
v téže sekci (015) zůstal nevysvětlený (ts2020-mem +0,008) a nebyl zde dál rozložen; kosinus a Jaccard nic
nerozlišují kvůli velikosti stránek; intervaly zisku (a) jsou široké a závěry stojí na průměrech přes rozdělení kmenů.
Anomálie sousedů jsou z velké části důsledek jazyka A/B (B listy v herbáři A nemají poblíž žádnou B stránku) a
neopravňují k závěru o přeskládání.

**Další krok:** bloky herbáře A rozložit podle obsahu: párovat stránky podle rostliny/ilustrace (totéž jako návrh
015 – herbář vs farma) a podle bifolia uvnitř kvaternu; ověřit, zda bloky v herbáři A sedí na bifolia (zde jen
kvaterny) a zda je dá ts2020-mem s pamětí definovanou přes fyzický list místo pořadí čtení; GC split prověřit na
úrovni dekompozice (proč jediný dává v herbáři A nulu). Pro 020 (odstavec jako jednotka) vzít stejný rámec
(f(d) + bloky) na úrovni řádků.
