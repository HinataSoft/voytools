# 014 – Inventář místních pravidel: kolik konvencí dělá pozici a šev slov

Kód: `voytools/placerules.py` (rozklad MI na položky, řetězový rozklad pozice, test švu, stabilita, post-procesor
R1–R3; testy `tests/test_placerules.py`), `experiments/014_place_rules_inventory/run.py`, `summarize.py`.
Rozhodnutí a interpretace spec: `experiments/014_place_rules_inventory/notes.md`. Všechna čísla níže jsou
z `results/014_place_rules_inventory/` (soubor v závorce u tabulky).

## 1. Rozklad poziční informace na položky (`coverage.csv`, `position_contributions.csv`, `position_chain.csv`)

Pointwise členy I(Y; poziční třída 007) sečtené po hodnotě Y (afix / glyf) a po třídě; počet položek na 50 / 80 / 95 %
celku. Dekompozice (a) naučené afixy, (b) + poziční glyfy (`p f t k…` / `m g am om`). Plná I (plug-in) a po korekci
na bias (permutace, 50×):

| Y (ZL split) | I plug-in / korig. (bit) | položek 50/80/95 % | tříd 50/80/95 % | párů 80 % | top položky (podíl I) |
|---|---|---|---|---|---|
| prefix (a) | 0,138 / 0,133 | **5 / 8 / 11** z 13 | 2 / 8 / 17 z 22 | 8 | s 0,13, ∅ 0,12, d 0,11, y 0,11, t 0,11, sh 0,10 |
| prefix (b) | 0,190 / 0,181 | 5 / 10 / 16 z 20 | 2 / 9 / 18 | 10 | **p 0,23**, s 0,09, d 0,08, y 0,08, t 0,08 |
| sufix (a) | 0,054 / 0,048 | **3 / 7 / 12** z 15 | 4 / 11 / 18 | 7 | **am 0,38**, ∅ 0,11, or 0,09, dy 0,07 |
| sufix (b) | 0,071 / 0,064 | 3 / 9 / 14 z 18 | 4 / 11 / 18 | 9 | am 0,28, m 0,17, ∅ 0,07, or 0,07, g 0,06 |
| první glyf | 0,191 / 0,180 | **4 / 9 / 15** z 30 | 2 / 8 / 18 | 10 | **p 0,23**, ch 0,11, y 0,09, d 0,08, t 0,08 |
| poslední glyf | 0,071 / 0,059 | **2 / 8 / 19** z 29 | 3 / 11 / 18 | 6 | **m 0,49**, r 0,08, g 0,07, o 0,05 |

GC split: prefix (a) 0,102 bit, 4 / 8 / 12 položek (9, k, 8, 1c…); sufix (a) 0,038, 4 / 9 / 13; první glyf 0,201,
5 / 12 / 32 z 63 (g 0,20); poslední glyf 0,084, 3 / 16 / 40 z 59 (p 0,39). ZL join a GC join se liší o ±1 položku.
Srovnání: ts2020 (3 seedy) prefix 0,062–0,069 bit, 6–7 položek na 80 %, ale **jen 3–4 třídy** (vs 8 u reálu);
katalog Linnaeus po položkách prefix 0,008 bit (8 položek), blábol ve vlastních řádcích 0,005 bit (10 položek,
16 tříd – difúzní).

Které třídy nesou informaci: prefix – `first|0|2` 0,28, `first|0|0` 0,24, `first|0|1` 0,12 (první token řádku
celkem 64 %), `last|4|2` 0,06, vnitřní třídy ≤ 0,04; sufix – `last|4|·` 0,38, první token 0,23; poslední glyf –
`last|4|·` 0,52. Nejsilnější páry: (∅ prefix, první token prvního řádku) 0,023 bit, (s-, první token) 0,021,
(y-, první token) 0,019, (t-, první token prvního řádku) 0,016; glyfy (p, první token prvního řádku) 0,047,
(m, poslední token) 0,027 + 0,009 + 0,009.

**Kraj řádku vs drift uvnitř** (řetězový rozklad, korigované bity; `position_chain.csv`):

| Y | I(pos) | místo v řádku (kraj) | první řádek odstavce \| místo | rel_pos \| místo, řádek (drift) |
|---|---|---|---|---|
| ZL split prefix (a) | 0,133 | **0,081** (61 %) | 0,041 (31 %) | 0,011 (8 %) |
| ZL split sufix (a) | 0,048 | 0,033 (70 %) | 0,011 (23 %) | 0,003 (6 %) |
| ZL split první glyf | 0,180 | 0,109 (61 %) | 0,058 (32 %) | 0,012 (7 %) |
| ZL split poslední glyf | 0,059 | 0,053 (88 %) | 0,005 | 0,002 |
| GC split prefix (a) | 0,097 | 0,061 (63 %) | 0,026 (27 %) | 0,010 (10 %) |
| GC split první glyf | 0,178 | 0,110 (62 %) | 0,056 (32 %) | 0,012 (7 %) |
| ts2020 prefix (a), 3 seedy | 0,062–0,069 | 0,013–0,015 | **0,048–0,054** | 0,000 |
| Linnaeus položky, prefix | 0,008 | 0,005 | 0,000 | 0,003 |

Shodně v ZL i GC, split i join: 60–70 % poziční informace afixů je **kraj řádku** (první/poslední token), 25–30 %
**první řádek odstavce** (gallows `p-`/`f-`, `t-`), jen 6–10 % je drift uvnitř řádku (Feaster). ts2020 má svou
polovinu informace téměř celou v prvním řádku odstavce (vlastní pravidlo generátoru), kraj řádku skoro nemá.

## 2. Rozklad švu (`seam_contributions.csv`, `coverage.csv`, `seam_test.csv`)

Okrajová MI (poslední glyf → první glyf dalšího tokenu uvnitř řádku), plug-in / korigovaná:

| zdroj | I (bit) | párů 50/80/95 % | po posl. glyfu 80 % | bez nejsilnějšího páru: podíl, párů 80 % | nejsilnější páry (podíl I) |
|---|---|---|---|---|---|
| ZL split | 0,207 / 0,186 | **2 / 4 / 6** z 340 | 6 | y→q **45 %**; 4 | y→q 0,45, r→a 0,22, s→a 0,10, iin→o 0,08, y→l 0,07, o→l 0,06 |
| ZL join | 0,174 / 0,152 | 2 / 3 / 4 | 5 | y→q 53 %; 3 | |
| GC split | 0,266 / 0,225 | 2 / 7 / 10 z 722 | 8 | 9→4 35 %; 8 | 9→4 0,35, y→a 0,17, s→a 0,09, m→o 0,06, 9→e 0,05 |
| GC join | 0,227 / 0,182 | 3 / 5 / 7 | 7 | 9→4 42 %; 5 | |
| ts2020 ZL (3 seedy) | 0,015 / 0,008 | 1–2 / 2 / 2 | 7–8 | r→a / l→o 44–57 % | difúzní, log-poměry ≤ 0,5 |
| Linnaeus položky | 0,337 / 0,297 | 2 / 9 / 16 | 10 | t→i 37 %; 12 | t→i, h→p, y→l, n→e (morfologie latiny) |
| blábol vlastní řádky | 0,106 / 0,013 | 17 / 33 / 42 | 22 | 5 %; 33 | difúzní (téměř celé je bias) |

Vyhýbání (záporné členy) v ZL: y→ch −12 % I, y→sh −8 %, y→o −8 %, r→q −6 %. Vazba sufix → prefix dalšího tokenu
(a): I 0,129 / 0,112 bit, **2 páry** na 80 % (dy→qo 0,50, ∅→∅ 0,32, y→qo 0,24; GC 89→4o 0,39, ∅→∅ 0,29, 9→4o 0,27).

**Je šev jen glyfový?** (korekce permutací tokenů uvnitř řádků, 50×)

| zdroj | I(posl. glyf; první glyf) | I(sufix; prefix) | I(sufix; prefix \| oba glyfy) | podíl zbytku | I(glyfy \| afixy) |
|---|---|---|---|---|---|
| ZL split (a) / (b) | 0,186 | 0,112 / 0,127 | **0,0010 / 0,0008** | 0,9 % / 0,7 % | 0,034 / 0,019 |
| ZL join (a) | 0,152 | 0,103 | 0,0008 | 0,8 % | 0,017 |
| GC split (a) / (b) | 0,225 | 0,108 / 0,121 | 0,0013 / 0,0004 | 1,2 % / 0,4 % | 0,048 / 0,033 |
| GC join (a) | 0,182 | 0,098 | 0,0008 | 0,9 % | 0,031 |
| Linnaeus položky | 0,297 | 0,033 | 0,0056 | 17 % | 0,189 |
| blábol | 0,013 | 0,005 | 0,0003 | 6 % | 0,006 |

Nad glyfovým párem nezbývá nic (≤ 0,001 bit, < 1,3 % afixové vazby) – v obou transliteracích i politikách mezer;
glyfy naopak nesou 0,02–0,05 bitu nad afixy. U latinského katalogu zbývá 17 % (shoda tvarů), tedy test rozlišuje.

## 3. Stabilita napříč A/B, písaři a sekcemi (`stability.csv`, `coverage_by_group.csv`, `stability.png`)

JSD mezi skupinami (na polovinách, stejné velikosti) / JSD mezi náhodnými polovinami téže skupiny (200×, řádky);
podíl pravidel (pointwise členy na 80 % Σ|příspěvek|) se stejným směrem v obou skupinách. ZL split, dekompozice (a):

| tabulka | A vs B: poměr; sdílené | písaři 1/2/3: medián poměru (rozsah); sdílené | sekce (≥ 2 000): medián (rozsah); sdílené |
|---|---|---|---|
| R1 prefix \| první token, první řádek odstavce | **0,67**; 2/2 | **0,96** (0,70–1,37); 1,00 | 2,8 (1,5–5,3); 1,00 |
| R1 prefix \| první token, jiný řádek | 10,1; 5/6 | 6,7 (5,5–11,0); 1,00 | 2,0 (1,0–11,1); 1,00 (min 0,83) |
| R1 prefix \| vnitřní token prvního řádku | 15,1; 2/3 | 7,3 (6,4–16,6); 1,00 | 3,0 (1,2–16,7); 1,00 (min 0,00 z 2) |
| první glyf \| první token | 6,5; 6/7 | 5,1 (4,7–6,7); 1,00 | 2,1 (1,4–8,0); 1,00 |
| R2 sufix \| poslední token | 8,1; 5/6 | 5,3 (2,8–9,5); 1,00 | 1,9 (1,0–7,4); 0,92 |
| poslední glyf \| poslední token | 5,5; 4/4 | 5,3 (2,4–6,1); 1,00 | 1,8 (1,3–7,4); 1,00 |
| R3 první glyf \| poslední glyf před | 10,8; 31/36 | 8,2 (2,9–11,7); 0,82 | 3,6 (1,8–10,3); 0,88 |
| R3 prefix \| poslední glyf před | 15,6; 28/34 | 12,4 (4,1–14,3); 0,82 | 4,4 (2,3–13,9); 0,93 |
| šev prefix \| sufix před | 9,4; 47/55 | 7,6 (3,0–8,9); 0,80 | 3,0 (1,6–8,2); 0,86 |

GC split dává totéž (A vs B 3,6–14; písaři 3,5–9,6; sekce 1,4–4; sdílené 0,75–1,00). Výjimka v obou: tabulka prvního
tokenu **prvního řádku odstavce** (gallows) je mezi A a B i mezi písaři nerozlišitelná od polovin (poměr 0,7–1,0).
Všude jinde se tabulky mezi A a B (a mezi písaři, kteří jsou s jazykem zaměněni) liší 5–15× víc než poloviny
(p < 0,005), ale **směr pravidel je sdílený**: 83–100 % položek u krajů řádku, 80–93 % u švu. Mezi sekcemi
(uvnitř jazyka: H vs P, S vs B…) je poměr 1–3, tj. blízko polovin. Po skupinách (`coverage_by_group.csv`): v A i B,
u každého písaře a každé sekce 7–9 prefixů a 2–7 tříd na 80 %; stejné top položky (∅/s, d, y, t, sh; am/or; p; m;
y→q, r→a, s→a) – liší se váhy (A: `d-`, `-or`, `y→k`, `l→d`; B: `s-`, `-am`, `r→a`, `s→a`).

## 4. Doplněný generátor ts2020 + R1 / R1+R2 / R1+R2+R3 (`augmented_profile.csv`, `augmented_lists.csv`, `augmented.png`)

Tabulky učené na ZL (afixy b), minimální přepis: R1 mění 3 548 tokenů (10 %), R2 1 157, R3 8 185 (24 %); 30 seedů,
celý profil 009 + 010 + 011 + kalibrace 004. z = (Voynich − null)/SD nullu; „dožene“ = |z| < 3 v split **i** join.

**Voynich-specifické vlastnosti (SYNTHESIS odd. 1), ZL split: reál / ts2020 (z) / +R1 (z) / +R1+R2 (z) / +R1+R2+R3 (z):**

| vlastnost | metrika | reál | ts2020 | +R1 | +R12 | +R123 | verdikt (split & join) |
|---|---|---|---|---|---|---|---|
| poziční vazba prefixu | prefix__I_pos (009) | 0,134 | 0,063 (12,8) | 0,160 (−1,0) | 0,164 (−1,1) | 0,115 (4,9) | R1, R12 **dožene**; R123 přiblíží |
| | prefix pos explained (011) | 0,046 | 0,026 (4,9) | 0,078 (−1,8) | 0,082 (−1,7) | 0,038 (2,2) | R123 dožene; R1/R12 jen split (join −3,1: přestřelí) |
| poziční vazba sufixu | suffix__I_pos | 0,048 | 0,029 (3,8) | 0,026 (6,7) | 0,068 (−2,5) | 0,068 (−2,5) | R12, R123 **dožene** |
| | suffix pos explained (011) | 0,026 | 0,015 (3,9) | 0,015 | 0,043 (−4,5) | 0,042 (−4,2) | stále ne (přestřelí) |
| začátek řádku | jsd_first_mid | 0,194 | 0,062 (21) | 0,207 (−0,5) | 0,212 (−0,7) | 0,177 (3,1) | R1, R12 **dožene**; R123 jen join |
| konec řádku | jsd_last_mid | 0,107 | 0,048 (7,2) | 0,050 | 0,151 (−1,8) | 0,150 (−1,7) | R12, R123 **dožene** |
| | end_mg_last (010) | 0,172 | 0,253 (−3,2) | 0,257 | 0,189 (−3,6) | 0,189 | jen join (split: SD nullu malá) |
| | end_enrich_max | 17,5× | 8,0× (5,9) | 7,8× | 6,5× (5,6) | 6,5× | **stále ne** (generátor má `-m` i uvnitř řádku: 0,035 vs 0,009) |
| šev slov | edge_mi_excess | 0,186 | 0,006 (70) | 0,006 | 0,005 | **0,126 (7,0)** | R123 přiblíží (2/3 hodnoty) |
| | bnd_edge_in (010) | 0,195 | 0,008 (69) | 0,007 | 0,008 | 0,131 (8,6) | přiblíží |
| | sufix → prefix dalšího (011) | 0,038 | 0,004 (9,7) | 0,004 | 0,005 | **0,024 (1,8)** | R123 **dožene** |
| uzavřený řádek | bnd_*_ratio_line | ≈ 0 | ≈ 0 (|z| < 1,1) | | | | beze změny (ts2020 má z konstrukce) |
| záhlaví řádku | I_first_rest, head_rest | 0,057; 0,012 | 0,029 (5,0); 0,007 (2,7) | 0,038 (1,9) | 0,039 (1,6) | 0,023 (6,3); 0,004 (7,4) | R1/R12 split dožene; **R123 pokazí** |

**Co se pokazí** (ts2020 mělo |z| < 3, varianta ≥ 3) a co se zhorší:

| metrika | reál | ts2020 | +R1 | +R12 | +R123 | poznámka |
|---|---|---|---|---|---|---|
| mi_prefix_suffix (009) | 0,205 | 0,182 (1,0) | 0,148 (3,3) | 0,144 (4,0) | **0,094 (11,9)** | přepis prefixu nezávisle na sufixu ruší vazbu prefix–sufix (005: ingredience šablony) |
| prefix__I_sec, I_sec\|pos | 0,085; 0,104 | 0,012 (20); 0,011 (25) | 0,010 | 0,010 (24) | 0,006 (42) | sekční vazba prefixu ještě klesne (R3 přepisuje 24 % prefixů podle globální tabulky) |
| stem__I_pos | 0,191 | 0,109 (7,2) | 0,098 (12) | 0,091 (12) | 0,083 (13) | zhorší: nové řezy (prefix + kmen) |
| tpl_gain_prev (010) | 0,088 | 0,035 (10) | 0,031 | 0,031 | 0,034 (19) | SD nullu klesne, hodnota ne |
| I_first_rest, head_rest | viz výše | | | | pokazí | R3 přepisuje první tokeny řádků? ne – přepisuje vnitřní, čímž ruší vazbu první ↔ zbytek z R1 |
| end_mg_other (010) | 0,009 | 0,035 (−2,8) | 0,035 (−3,6) | | | jen zúžení SD |
| a__cons_page__excess (011) | 0,023 | 0,054 (−2,4) | 0,056 (−3,1) | 0,053 | 0,053 | jen zúžení SD |

Kalibrační cíle 004 a slovník (`calibration.csv`): H2 glyfů (extended) 2,42 reál / 2,39 ts2020 / 2,43 / 2,46 / 2,53
(+R123, z −1,6); délka 4,06 / 4,26 / 4,31 / 4,29 / 4,39; typy 6 796 / 6 500 / 7 361 / 7 678 / **8 927** (+31 %;
z jen −1,6, protože ts2020 má přes seedy SD 1 126); hapax 0,676 / 0,578 (9,0) / 0,608 / 0,619 / **0,645 (2,0;
join 7,4)**. Paradigma: H(sufix|kmen) 1,82 / 1,61 / 1,63 / 1,76 / 1,78 (R12 i R123 blíž reálu), sufixů na kmen
beze změny (3,8 vs 5,25), cmi_excess beze změny. Hapax a počet typů tedy rostou přepisem (nové kombinace
prefix + kmen), ne kopírováním: „zlepšení“ hapaxu je artefakt, který se platí 2 000 typy navíc.

Souhrn verdiktů (65 metrik × varianta): R1 dožene 3 (+5 jen v jedné politice), pokazí 1; R12 dožene 6 (+9), pokazí 1,
zhorší 2; R123 dožene 4 (+10), přiblíží 4, pokazí 3, zhorší 4.

## 5. Odpovědi na kritéria

1. **Kolik položek nese 80 % poziční informace:** prefix 8 (ZL i GC, split i join; dekompozice b 10), sufix 7–9,
   první glyf 9 (GC 12), poslední glyf 8 (GC 16); páry (položka, třída) 6–10; tříd 8–11. Polovinu nese
   5 / 3 / 4 / 2 položek. 60–70 % informace je kraj řádku, 25–30 % první řádek odstavce, < 10 % drift uvnitř.
   Jde tedy o **malou sadu**: `s-/d-/y-/∅` a `t-/p-` na začátku řádku, `-am/-m/-g` (poslední glyf `m` sám 49 %)
   na konci, `p-` v prvním řádku odstavce.
2. **Šev:** 4 páry nesou 80 % okrajové MI (y→q 45 %, r→a 22 %, s→a 10 %, iin→o 8 %); bez y→q 4 páry z 55 %
   zbytku; GC 7 (9→4 35 %, y→a 17 %, s→a 9 %). Vazba sufix → prefix dalšího tokenu = 2 páry (dy→qo, y→qo).
   Nad glyfovým párem **nezbývá nic** (≤ 0,001 bit, < 1,3 %; u latinského katalogu 17 %).
3. **Sdílení:** směr pravidel je sdílený napříč A/B, písaři i sekcemi (83–100 % položek u krajů, 80–93 % u švu),
   a tabulka prvního tokenu prvního řádku odstavce je mezi A a B i písaři nerozlišitelná od polovin. Váhy ale
   sdílené nejsou: ostatní tabulky se mezi A a B liší 5–15× víc než poloviny téže skupiny (sekce uvnitř jazyka
   1–3×). Podezřelý případ „JSD mezi skupinami < JSD polovin“ nastal jen u gallows tabulky (poměr 0,67–0,97),
   kde jsou poloviny malé (≈ 700 tokenů) a tabulka má 2 pravidla – interpretujeme jako shodu, ne chybu.
4. **Doplněný ts2020:** R1 dožene poziční vazbu prefixu a začátek řádku (jsd_first_mid, I_first_rest); R2 navíc
   sufix a konec řádku (jsd_last_mid, end_mg_last v join), ne obohacení koncového glyfu (17× vs 6–8×, protože
   generátor píše `-m` i uvnitř řádku); R3 dožene vazbu sufix → prefix dalšího tokenu a dostane okrajovou MI
   na 2/3 (0,126 vs 0,186; zbytek je v tokenech bez prefixu, kde R3 první glyf nemění), ale **pokazí** vazbu
   prefix–sufix (0,094 vs 0,205), sekční vazbu prefixu, záhlaví řádku a přidá 2 000 typů (hapax 0,645 je pak
   artefakt přepisu). Nic z toho nemění paradigma (sufixů na kmen 3,8 vs 5,25), pořadí tokenů ani A/B.

## Interpretace

Poziční a švové vlastnosti Voynichu jsou **malá sada konvencí** – řádově deset glyfových pravidel (začátek řádku,
konec řádku, první řádek odstavce, čtyři švové páry), ne difúzní jev; hypotéza „konvence“ ze spec platí v obou
transliteracích a obou politikách mezer. Tato pravidla jsou stejná co do směru napříč jazyky A/B, písaři
i sekcemi (sdílený systém, ne idiosynkrazie pisatele), ale A a B je uplatňují s jinými vahami, což zapadá do 012
(dva slovníky). Generátor ts2020 s dodanými tabulkami dožene pozici a část švu, ale jen za cenu rozbití
vazby prefix–sufix a zvýšení počtu typů – kopírování plus lokální přepis nejsou mechanismus, který by dal
pozici, šev **a** vnitroslovní soudržnost zároveň.

## Co výsledek neukazuje

* Neříká, *proč* konvence existují (ozdobný zvyk, zápis, jazyk): katalog Linnaeus má podobně malou sadu švových
  párů (9 na 80 %) z morfologie latiny, blábol nemá žádnou; sama koncentrace nerozhoduje mezi jazykem a zvykem.
* Post-procesor je jedna konkrétní (minimální, nezávislá) realizace; přepis podmíněný i na sufixu nebo na
  prvním glyfu kmene by vazbu prefix–sufix neporušil a šev by dohnal celý – to jsme netestovali (spec: R1–R3).
* Stabilita měří tabulky, ne pravidla po odečtení základních frekvencí afixů v A/B; JSD mezi A a B tak částečně
  odráží rozdíl slovníků (012), nejen rozdíl pozičních pravidel. Podíl sdílených směrů tím netrpí.
* z-skóre podle SD 30 seedů je přísné u metrik s malým rozptylem (end_mg_other, cons_page: „pokazí“ jen zúžením SD).

## Další krok

Podmíněný přepis (prefix | sufix, poslední glyf před) jako R3′ a první glyf kmene pro tokeny bez prefixu: zjistit,
zda jde dohnat šev celý bez ztráty prefix–sufix vazby – a pokud ne, zda je ta vazba (H(sufix | kmen) + mi_prefix_suffix
+ pozice + šev dohromady) vlastnost, kterou žádný lokální mechanismus nedá. Dále 015 (vzácné kmeny) a 016 (A/B varianty).
