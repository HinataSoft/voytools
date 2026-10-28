# 019 – Tvary jen v B: jména, jiná morfologie, nebo slovník druhého systému?

Kód: `experiments/019_b_only_forms/run.py`; funkce `voytools/bforms.py` (třídy tvarů, vyrovnání četností, profily,
volba `-ey`/`-edy`; testy `tests/test_bforms.py`). Rozhodnutí a odchylky od spec: `experiments/019_b_only_forms/notes.md`.

## Co se počítalo

* ZL3b a GC2a, `P`, clean, split a join; **celé tokeny** jako n-tice glyfů (preset `extended` / znaky v101), kmen
  dekompozice (a). Stejná velikost jako 012/016 (N = polovina menšího jazyka: ZL 5 497 / 5 094, GC 5 692 / 5 389 tokenů
  na stranu, celé stránky), **200 losů**; v každém losu A × B oběma směry (**B-only vs A** hlavní, A-only vs B),
  dvě poloviny B (šum) a dvě poloviny A.
* Třídy tvarů X: `shared` (≥ 1 token v Y), **B-only** (≥ 2 v X, 0 v Y) = `d_variant` (protějšek ≤ 2 tokeny Y jednou
  Currierovou operací: vložení `d`/`e`/`ee`/`q`, `ch`↔`sh`; v101 `8`/`7`/`c`/`C`/`4`, `1`↔`2`), `other_le1` (jiný
  protějšek do 1 editace), `no_counterpart` (nic do 1 editace); `hapax` (1 v X, 0 v Y). **Nový kmen** = kmen (a) bez
  tokenu na druhé straně.
* Profil každé třídy proti **sdíleným tvarům vyrovnaným počtem tokenů** (stejný počet; hapaxy proti sdíleným tvarům
  s 1 tokenem; rozdíl v témže losu, interval 2,5–97,5 % přes losy): tvarů na kmen, H(sufix | kmen), stránek na výskyt,
  Griesův DP, podíl na nejčastější stránce, sekční entropie, popisky (plný text), kraje řádků, souvýskyt dvou tvarů
  třídy v řádku proti promíchání, bigramová log2-věrohodnost pod modelem sdílených tvarů, šev (JSD prvního glyfu
  následujícího tokenu proti sdíleným tokenům).
* Pravidlo `-ey`/`-edy` na plném textu: P(`-edy` | glyf před, prefix, pozice v řádku, poslední glyf předchozího
  tokenu); v101 `-c9`/`-c89` a `-c9`/`-c79` (GC rozlišuje dvě podoby `d`), měřítko `-y`/`-dy` (základy na `e`
  vyloučeny; v101 `-9`/`-89`) a `-ol`/`-or` (v101 `-oe`/`-oy`); 200 permutací volby.
* Kontroly: latina Plinius knihy 16/18/21 vs 23–30 (téma) a Caesar vs Celsus (registr), 200 losů; ts2020 config 21
  (ZL split/join) a šablona 005 (ZL/GC split), 30 seedů, A/B podle layoutu. Latina `-us`/`-um`, `-is`/`-es` jako hrubé
  měřítko kontextové vázanosti.

## (i) Skladba pětiny (`composition.csv`, `composition.png`, `new_stems_ed.csv`)

Podíl tokenů B připadajících na B-only tvary a jejich rozdělení do tříd (průměr [2,5 %; 97,5 %] přes losy;
„poloviny“ = tatáž třída mezi dvěma polovinami B):

| sada | B-only podíl tokenů B | poloviny | hapax | poloviny | d-varianty | jiný ≤ 1 | bez protějšku | nový kmen (z B-only) | poloviny | nový kmen (z hapaxů) | poloviny |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ZL split | **0,219** [0,180; 0,262] | 0,052 | 0,186 | 0,176 | **0,630** [0,537; 0,710] | 0,246 [0,181; 0,319] | 0,124 [0,082; 0,172] | **0,065** [0,041; 0,100] | 0,076 | 0,476 | 0,408 |
| GC split | **0,208** [0,176; 0,250] | 0,057 | 0,221 | 0,208 | **0,621** [0,542; 0,689] | 0,264 [0,199; 0,341] | 0,115 [0,079; 0,152] | **0,075** [0,052; 0,107] | 0,109 | 0,509 | 0,455 |
| ZL join | 0,232 [0,192; 0,276] | 0,057 | 0,230 | 0,220 | 0,589 [0,496; 0,674] | 0,258 [0,189; 0,332] | 0,154 [0,106; 0,207] | 0,071 [0,047; 0,104] | 0,084 | 0,528 | 0,466 |
| GC join | 0,223 [0,192; 0,261] | 0,060 | 0,258 | 0,244 | 0,590 [0,507; 0,664] | 0,271 [0,211; 0,345] | 0,140 [0,101; 0,179] | 0,077 [0,054; 0,100] | 0,107 | 0,528 | 0,467 |
| ts2020 ZL split | 0,131 [0,066; 0,164] | 0,129 | 0,167 | 0,167 | 0,229 | 0,506 | 0,265 | 0,295 | 0,299 | 0,480 | 0,483 |
| šablona ZL split | 0,075 [0,063; 0,092] | 0,075 | 0,209 | 0,208 | 0,325 | 0,472 | 0,203 | 0,086 | 0,089 | 0,285 | 0,294 |
| šablona GC split | 0,087 [0,069; 0,100] | 0,085 | 0,247 | 0,247 | 0,322 | 0,482 | 0,196 | 0,104 | 0,100 | 0,278 | 0,282 |
| Plinius (téma) | 0,189 [0,164; 0,221] | 0,102 | 0,315 | 0,294 | 0 | 0,195 | **0,805** | **0,708** | 0,638 | 0,830 | 0,795 |
| Caesar/Celsus (registr) | 0,301 [0,283; 0,319] | 0,068 | 0,157 | 0,125 | 0 | 0,202 | **0,798** | **0,769** | 0,585 | 0,830 | 0,715 |

A-only strana (A-only vs B; poloviny A): podíl 0,154 / 0,148 / 0,165 / 0,158 (ZL split, GC split, ZL join, GC join;
poloviny 0,056–0,061); d-varianty 0,26–0,30, jiný ≤ 1 **0,61–0,65**, bez protějšku 0,08–0,11; nový kmen 0,106–0,144
(poloviny 0,086–0,123).

Inventář operací v plném textu (`classes.csv`): u d-variant B nese `ins:ed@mid` 0,41 a `ins:d@mid` 0,30 tokenů (ZL;
GC `ins:c8` 0,25, `ins:8` 0,16, `ins:c7` 0,15, `ins:7` 0,13), u „jiného protějšku“ `sub:o>l@begin` 0,20 a `ins:l` 0,20
(ZL `lkaiin`/`okaiin`, `olkain`/`okain`; GC `sub:o>e@begin` 0,15, `ins:e` 0,14 – v101 `e` = Eva `l`). Tvary **bez
protějšku** jsou v plném textu z větší části kombinace obou (`lkedy`, `lkchedy`, `qolchedy`, `rchedy`; GC `e2c89`,
`ehc89`): dvě editace od A. Z tokenů B-only s **novým kmenem** má kmen `ed` v 0,40 (ZL) / `c8`,`c7` v 0,26 (GC) a celý
tvar v 0,58 / 0,39 (`ked`, `keed`, `ted`, `pched`, `ched`: kmeny A `ke`, `te`, `pche` s `d` uvnitř, které parser neodtrhl);
u A-only nových kmenů 0,00.

* **Odpověď (i):** pětina B (0,21–0,23 tokenů; poloviny 0,05–0,06) se skládá z **≈ 60 % d-variant sdílených tvarů**
  (0,59–0,63), **≈ 25 % jiných tvarů do jedné editace** (0,25–0,27, převážně `l-` za `o-`), **≈ 12–15 % tvarů bez
  protějšku do jedné editace** (většinou `l-` + `-ed-` zároveň) – a jen **6,5–7,7 % tokenů s novým kmenem**, což je
  **pod úrovní polovin B (7,6–10,9 %)**; po odečtení kmenů s `ed` uvnitř zbývají ≈ 4 % tokenů pětiny. Hapaxy B
  (0,19–0,26 tokenů) i podíl nových kmenů mezi nimi (0,48–0,53) jsou na úrovni polovin (0,18–0,24; 0,41–0,47).
  ZL i GC, split i join shodně. Latina s jiným tématem i registrem má opačnou skladbu: 0,80 bez protějšku a 0,71–0,77
  nových kmenů (poloviny 0,59–0,64); generátory dávají A × B = poloviny ve všem.

## (ii) Profil tříd (`class_profiles.csv`, `latin_reference.csv`, `class_profiles.png`)

Třída minus vyrovnané sdílené tvary [interval přes losy]; „poloviny B“ = táž třída definovaná mezi polovinami B.
Hlavní třída **bez protějšku** (ZL 45 / GC 48 tvarů na los):

| metrika | ZL B-only | ZL poloviny B | GC B-only | GC poloviny B | Plinius (téma) | Caesar/Celsus (registr) | ts2020 | šablona ZL |
|---|---|---|---|---|---|---|---|---|
| tvarů na kmen | −12,3 [−19,8; −5,1] | −15,3 [−35,5; −1,7] | −9,1 [−15,0; −2,5] | −15,3 [−35,0; −1,7] | −0,7 [−0,8; −0,4] | −0,6 [−0,7; −0,4] | −10,5 | −60,1 |
| H(sufix \| kmen) | 0,004 [−0,524; 0,455] | −0,462 [−1,828; 0,848] | 0,242 [−0,186; 0,575] | −0,409 [−1,861; 0,937] | −0,239 [−0,374; −0,045] | −0,289 [−0,399; −0,187] | −0,358 | −0,676 |
| stránek na výskyt | −0,012 [−0,066; 0,053] | −0,004 [−0,167; 0,146] | −0,017 [−0,096; 0,052] | −0,008 [−0,208; 0,167] | **−0,064 [−0,092; −0,036]** | **−0,070 [−0,094; −0,046]** | −0,273 | −0,172 |
| Gries DP | −0,005 [−0,026; 0,013] | −0,005 [−0,041; 0,035] | −0,004 [−0,026; 0,020] | 0,002 [−0,046; 0,051] | −0,000 [−0,010; 0,008] | 0,002 [−0,003; 0,005] | 0,047 | 0,023 |
| podíl na nejčastější stránce | 0,011 [−0,049; 0,071] | 0,005 [−0,146; 0,167] | 0,017 [−0,052; 0,092] | 0,008 [−0,167; 0,208] | **0,061 [0,030; 0,087]** | **0,056 [0,032; 0,076]** | 0,274 | 0,172 |
| sekční entropie | −0,084 [−0,295; 0,107] | −0,058 [−0,500; 0,376] | −0,075 [−0,273; 0,124] | 0,030 [−0,500; 0,501] | **−0,165 [−0,262; −0,067]** | **−0,303 [−0,408; −0,206]** | −0,399 | −0,227 |
| podíl tvarů v popiscích | −0,072 [−0,147; 0,000] | −0,061 [−0,301; 0,125] | −0,065 [−0,137; 0,001] | −0,040 [−0,333; 0,168] | – | – | – | – |
| první v řádku | 0,025 [−0,085; 0,144] | 0,028 [−0,223; 0,333] | 0,010 [−0,083; 0,113] | 0,090 [−0,200; 0,401] | −0,000 | −0,001 | −0,010 | 0,008 |
| poslední v řádku | −0,011 [−0,121; 0,083] | 0,073 [−0,223; 0,389] | −0,026 [−0,132; 0,072] | 0,070 [−0,333; 0,400] | 0,000 | −0,010 | −0,013 | 0,001 |
| první řádek odstavce | 0,081 [−0,041; 0,210] | 0,144 [−0,138; 0,500] | 0,055 [−0,088; 0,179] | 0,144 [−0,301; 0,530] | 0,025 [−0,013; 0,061] | 0,021 [0,002; 0,043] | −0,134 | −0,023 |
| souvýskyt v řádku / náhoda | 0,082 [−0,653; 0,949] | −1,910 [−18,3; 2,2] | 0,210 [−0,564; 0,954] | 0,176 [−3,307; 3,182] | **0,268 [0,069; 0,424]** | **0,174 [0,045; 0,341]** | 0,424 | 0,117 |
| bigramová log2-věrohodnost | 0,096 [−0,209; 0,422] | −0,072 [−1,313; 0,896] | −0,189 [−0,583; 0,240] | −0,117 [−1,037; 1,185] | **−0,124 [−0,202; −0,037]** | **−0,221 [−0,296; −0,140]** | −0,421 | −0,839 |
| šev (JSD) | 0,051 [−0,001; 0,113] | −0,004 [−0,165; 0,152] | **0,075 [0,021; 0,132]** | 0,021 [−0,179; 0,215] | 0,010 [0,003; 0,017] | 0,008 [0,003; 0,014] | 0,004 | 0,001 |

Absolutní hodnoty ZL (třída / vyrovnané sdílené): stránek na výskyt 0,866 / 0,879, nejčastější stránka 0,525 / 0,513,
sekční entropie 0,586 / 0,670, v popiscích 0,019 / 0,091 tvarů, bigram −2,79 / −2,89, šev 0,080 / 0,029.
Join (`ZL_join`, `GC_join`): stejné znaménko a velikost (stránek na výskyt −0,013 / −0,013, šev +0,037 [0,003; 0,083] /
+0,073 [0,035; 0,134], bigram +0,04 / −0,33 s intervalem přes nulu, popisky −0,068 / −0,062).

Ostatní třídy (ZL / GC split, B-only; poloviny B v závorce):

* **d-varianty** (111 / 135 tvarů): šev **+0,040 [0,020; 0,068] / +0,065 [0,044; 0,093]** (poloviny +0,004 / +0,008),
  poslední v řádku **−0,044 [−0,089; −0,006] / −0,045 [−0,084; −0,008]** (poloviny −0,040 / −0,036, interval přes nulu),
  H(sufix | kmen) +0,14 / **+0,29 [0,02; 0,58]**, první řádek odstavce +0,047 [−0,003; 0,106] / +0,026; disperze, sekce,
  popisky, bigram: 0 (např. stránek na výskyt +0,010 / +0,007, popisky 0,000 / −0,020, bigram +0,17 / +0,08 s intervaly
  přes nulu). Join stejně (šev +0,035 / +0,064, poslední −0,035 / −0,047).
* **jiný protějšek ≤ 1** (85 / 96 tvarů): bigram **−0,28 [−0,51; −0,06] / −0,38 [−0,65; −0,11]** – ale poloviny B
  totéž (−0,36 / −0,60), tvarů na kmen +0,5 / −3,4 (poloviny −5,3 / −8,0); disperze, sekce, popisky, kraje, šev
  (+0,007 / +0,011) na nule. Latina: bigram −0,04 / −0,15, stránek na výskyt −0,07 / −0,05, nejčastější stránka +0,06 /
  +0,05.
* **hapaxy** (193 / 206 tvarů, vyrovnány sdílenými s 1 tokenem): bigram −0,26 / −0,45 (poloviny −0,48 / −0,55), první
  řádek odstavce +0,10 [0,02; 0,20] / +0,09 (poloviny +0,09 / +0,09), popisky −0,038 / −0,035 (poloviny −0,038 /
  −0,036): hapaxy B se od hapaxů jedné poloviny B neliší v ničem.
* A-only vs B (ZL split): bez protějšku – stránek na výskyt +0,013, sekční entropie −0,147 [−0,379; 0,119], popisky
  −0,045, bigram −0,20 [−0,60; 0,22], šev +0,017 [−0,035; 0,069]: nic mimo interval; jiný ≤ 1 – sekční entropie −0,141
  [−0,227; −0,036], Gries DP +0,008 [0,001; 0,016] (malé, GC +0,007 [0,001; 0,014]).
* **Krit. 3 spec (podezřelé hodnoty):** třída z polovin B s profilem odlišným od sdílených nad interval – jen tam, kde
  totéž platí pro A × B (tvarů na kmen, bigram u jiných ≤ 1 a hapaxů, první řádek odstavce u hapaxů): není chyba
  vyrovnání, ale vlastnost řídkých tvarů obecně. Bigram třídy = sdílených na 3 místa: nenastalo (rozdíly −0,19 až +0,17
  s širokými intervaly). P(`-edy`) bez rozptylu přes kontexty: nenastalo (0,32–0,80, viz (iii)).
* **Odpověď (ii):** nové tvary B (bez protějšku i jiný ≤ 1) se od sdílených tvarů B při stejné četnosti **neliší v
  disperzi po stránkách, v sekčním rozložení, v popiscích, na krajích řádků, v souvýskytu ani v glyfové stavbě**
  (rozdíly −0,02 až +0,02 s intervaly přes nulu, ZL i GC, split i join) – tj. **ne vzor Plinius/Caesar**, kde slova jen
  v jiném tématu/registru jsou soustředěna na méně stránek (−0,06 až −0,07), v méně sekcích (−0,17 / −0,30), častěji
  spolu v řádku (+0,17 / +0,27) a mají cizí bigramy (−0,12 / −0,22), s intervaly mimo nulu. Jediné, čím se B-only
  tvary od sdílených liší, je **šev**: po d-variantách (+0,04 / +0,065) a po tvarech bez protějšku (+0,05 / +0,075)
  začíná následující token jinak než po sdílených tokenech stejné četnosti, a d-varianty se vyhýbají konci řádku
  (−0,04). Nové kmeny B jsou příliš řídké na samostatný profil; jejich podíl není nad polovinami.

## (iii) Pravidlo `-ey` → `-edy` (`edy_rule.csv`, `edy_contexts.csv`, `edy_rule.png`)

Plný text, kontext = glyf před sufixem; `determ.` = podíl tokenů v kontextech (≥ 5 tokenů) s P ≤ 0,1 nebo ≥ 0,9,
`pásmo` = P v [0,3; 0,7]; `oba tvary` = podíl kmenů (≥ 2 tokeny) s oběma tvary, `táž stránka` = z nich s oběma na jedné stránce.

| sada | dvojice | jazyk | n | P(2. tvar) | H | H \| glyf | redukce | perm. | z | determ. | pásmo | kmenů | oba tvary | táž stránka |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ZL split | `-ey`/`-edy` | A | 974 | 0,023 | 0,156 | 0,149 | 0,006 | 0,008 | −0,7 | 0,984 | 0,000 | 92 | 0,163 | 0,133 |
| ZL split | `-ey`/`-edy` | **B** | 6570 | **0,603** | 0,969 | 0,922 | **0,047** | 0,002 | 71 | **0,000** | **0,826** | 254 | **0,787** | **0,365** |
| ZL split | `-y`/`-dy` | A | 2333 | 0,249 | 0,809 | 0,301 | **0,508** | 0,008 | 227 | 0,652 | 0,108 | 245 | 0,290 | 0,310 |
| ZL split | `-y`/`-dy` | B | 3460 | 0,397 | 0,969 | 0,548 | **0,421** | 0,006 | 298 | 0,512 | 0,390 | 296 | 0,416 | 0,276 |
| ZL split | `-ol`/`-or` | A | 2482 | 0,406 | 0,974 | 0,958 | 0,016 | 0,008 | 4 | 0,000 | 0,974 | 138 | 0,804 | 0,351 |
| ZL split | `-ol`/`-or` | B | 1809 | 0,338 | 0,923 | 0,863 | 0,060 | 0,010 | 18 | 0,084 | 0,670 | 127 | 0,677 | 0,395 |
| GC split | `-c9`/`-c89` | A | 686 | 0,020 | 0,144 | 0,123 | 0,021 | 0,014 | 1 | 0,992 | 0,000 | 80 | 0,100 | 0,125 |
| GC split | `-c9`/`-c89` | **B** | 3796 | **0,575** | 0,984 | 0,961 | **0,023** | 0,008 | 12 | **0,000** | **0,919** | 183 | **0,803** | **0,381** |
| GC split | `-c9`/`-c79` | B | 2564 | 0,371 | 0,952 | 0,901 | 0,050 | 0,010 | 21 | 0,009 | 0,779 | 150 | 0,727 | 0,367 |
| GC split | `-9`/`-89` | A | 2596 | 0,194 | 0,710 | 0,256 | **0,454** | 0,014 | 199 | 0,720 | 0,086 | 265 | 0,249 | 0,197 |
| GC split | `-9`/`-89` | B | 6001 | 0,240 | 0,795 | 0,461 | **0,334** | 0,008 | 256 | 0,569 | 0,352 | 502 | 0,299 | 0,320 |
| GC split | `-oe`/`-oy` | B | 1754 | 0,309 | 0,892 | 0,826 | 0,066 | 0,018 | 15 | 0,101 | 0,463 | 128 | 0,633 | 0,333 |
| Plinius | `-us`/`-um` | A / B | 2017 / 3070 | 0,57 / 0,60 | 0,99 / 0,97 | 0,71 / 0,66 | 0,28 / 0,31 | 0,01 | 112 / 214 | 0,27 / 0,42 | 0,58 / 0,52 | 305 / 425 | 0,18 / 0,17 | 0,19 / 0,27 |
| Caesar/Celsus | `-us`/`-um` | A / B | 5333 / 4843 | 0,60 / 0,57 | 0,97 / 0,99 | 0,55 / 0,73 | 0,42 / 0,26 | 0,00 | 567 / 309 | 0,51 / 0,18 | 0,45 / 0,44 | 534 / 569 | 0,22 / 0,23 | 0,22 / 0,22 |

Join totéž (ZL `-ey`/`-edy` B: redukce 0,048, pásmo 0,824, oba tvary 0,775, táž stránka 0,375; GC `-c9`/`-c89` B: 0,024,
0,894, 0,784, 0,368). Ostatní kontexty v B (ZL / GC): prefix redukce 0,005 / 0,014, pozice v řádku 0,001 / 0,000,
poslední glyf předchozího tokenu 0,012 / 0,023, glyf + prefix 0,071 / 0,065 (77 / 70 kontextů, pásmo 0,72 / 0,88,
determ. 0,003 / 0,002). P(`-edy` | glyf před) v B ZL: po `e` 0,49 (n 2268), `ch` 0,67 (1770), `sh` 0,66 (965), `k` 0,75 (730),
`t` 0,80 (399), `ee` 0,32 (179), `ckh` 0,37 (111); v A 0,00–0,03 všude. Nully (30 seedů, A/B podle layoutu): šablona
`-ey`/`-edy` redukce 0,03–0,04, pásmo 0,91–0,93, oba tvary 0,66–0,68; ts2020 P(`-edy`) 0,91–0,94 v obou (generátor přebírá
B-ovitý slovník), `-ol`/`-or` redukce 0,01–0,02 – kontext nic nedává ani v generátorech.

* **Odpověď (iii): volba tvaru, ne deterministické pravidlo kontextu.** V B leží P(`-edy`) pro všechny hlavní kontexty
  v pásmu 0,3–0,8 (83 % tokenů ZL, 92 % GC v pásmu 0,3–0,7; 0 % v deterministických kontextech), kontext uvnitř slova
  sníží entropii volby jen o 0,047 bitu (ZL) / 0,023 (GC) z 0,97–0,98 (permutace 0,002–0,008, z 12–71: reálný, ale
  malý vliv), prefix, pozice v řádku a předchozí token nic; **79–80 % kmenů B má oba tvary a 37–38 % z nich na téže
  stránce**. To je stejný obraz jako u `-ol`/`-or` (redukce 0,06, oba tvary 0,68–0,80) a **méně** vázané než `-y`/`-dy`
  (redukce 0,42–0,51, 51–72 % tokenů v deterministických kontextech, oba tvary 0,29–0,42) i než latinské `-us`/`-um`
  (0,26–0,42). Deterministický je jen **jazyk**: v A je `-edy` 2 % (`-c89` 2 %, `-c79` 0,3 %) bez ohledu na kontext.

## (iv) Seznamy (`top_new_stems.csv`, `classes.csv`)

Nejčastější nové kmeny B (ZL, plný text; tokeny / tvarů / stránek / sekce): `ked` 59 / 24 / 34 (S 36, H 11, B 8, T 4; tvary
`qokedar`, `okedar`, `okedal`…), `keed` 36 / 22 / 20, `ted` 35 / 14 / 25 (8 v popiscích), `lche` 31 / 9 / 19 (B 18),
`pched` 30 / 10 / 20 (63 % na začátku řádku), `ched` 29 / 22 / 18 (62 % na začátku řádku), `keche` 21 / 10 / 13, `eed` 20,
`pchd` 18, `pshe` 16. GC: `hc8` 29 / 16 / 20, `5cc` 27 / 11 / 26, `c7` 26 / 15 / 22, `1c8` 23 / 19 / 11, `e1c` 22 / 9 / 14 –
tytéž kmeny ve v101 (`hc8` = `ked`, `1c8` = `ched`, `e1c` = `lche`). Jsou to kmeny A s `d` (v101 `8`/`7`) uvnitř nebo s `l-`
(`e`) před, široce rozprostřené (20–34 stránek, 3–5 sekcí), s mnoha tvary (9–24), téměř bez popisků (0–1; `ted` 8).
Nové kmeny A (ZL): `cthol` 9 / 5 / 9, `cto` 9 / 4 / 9, `chol` 7, `cthod` 7, `keom` 7, `kom` 7 – všechny ≤ 9 tokenů, herbář.

## Interpretace

Pětina B-tokenů, které A nemá, není nový slovník: ze tří pětin jsou to `-ed-`/`-d-` varianty sdílených tvarů, z jedné
čtvrtiny tvary s `l-` místo `o-` a zbytek kombinace obou; tvarů s novým kmenem je 6–8 % pětiny, ne víc než mezi
dvěma polovinami B, a i ty jsou z 40–60 % kmeny A s `ed` uvnitř. Nové tvary B se při stejné četnosti chovají jako sdílené
tvary B ve všem, co by odlišilo jména či položky (disperze, sekce, popisky, kraje řádků, souvýskyt, glyfová stavba) –
vzor latinského tématu i registru (soustředěné, cizí, spolu v řádku) chybí – a liší se jen na švu a v tom, že d-varianty
nestojí na konci řádku. Vložení `d` je v B volbou tvaru (P 0,3–0,8 ve všech kontextech, čtyři kmeny z pěti s oběma
tvary, třetina na téže stránce) a v A prostě neexistuje. Věta: **pětina B = 60 % d-variant + 25 % `l`-variant + 12–15 %
jejich kombinací, 7 % nových kmenů (= poloviny); profil nových tvarů = sdílené (druhý systém / táž konvence), ne jména
ani registr; `-ey → -edy` = volba tvaru vázaná na jazyk, ne na kontext.** To je slučitelné s A a B jako dvěma
konvencemi téhož systému (012, 014), které se liší v inventáři přípustných vložení (`ed`, `l-`), ne s B jako textem o
jiném obsahu nebo s B jako jinou morfologií sdílených kmenů; generátory bez jazyka nedávají nic.

## Co výsledek neukazuje

* Neříká, že A a B mají stejný obsah – jen že B-only tvary nenesou stopy jmen/položek ani tématu, které latina nese.
  Sekční test má v B malou sílu (čtyři sekce, B-only tvary jsou v 3–5 z nich jako sdílené).
* „Nový kmen“ závisí na parseru dekompozice (a): `ed` uvnitř kmene parser neodtrhne, proto jsou nové kmeny B většinou
  `-ed-` varianty; jiný parser by dal jinou skladbu, ne jiný závěr (profily jsou na celých tokenech).
* Šev je jediná odlišnost B-only tvarů a není zde vysvětlen: může jít o to, že `-edy`/`l-` tvary stojí v jiných
  pozicích vět/řádků (014: šev = glyfová dvojice), nebo o vazbu `-edy` → `qo-` apod.; test proti polovinám B je u
  tvarů bez protějšku slabý (13–15 tvarů na los, široké intervaly).
* Latinské třídy nemají d-operace (vše „bez protějšku“); srovnání profilů Voynich vs latina je po třídách, ne po
  stejném mechanismu. `-us`/`-um` je gramatická alternace, ne volné psaní – jako měřítko „běžné variace“ je to horní mez.
* Jazyk a písař zůstávají zmatené (012); zda „volbu“ `-edy` dělá písař nebo systém, odtud nelze říci.

## Další krok

* Šev u `-edy`/`l-` tvarů: rozložit JSD na glyfové dvojice (014) a ověřit, zda jde o vazbu `-dy` → `qo-`/`l-` (tj. tatáž
  místní pravidla, jen s jinými vahami v B) nebo o novou vazbu; ts2020+R (014) s B-vahami jako null.
* Nové kmeny B bez `ed`/`l` (≈ 4 % pětiny, `pai`, `5cc`): ručně projít stránky a sekce; pokud se shlukují, je to jediný
  kandidát na obsah v B – test na úrovni stránek (020).
