# SYNTHESIS – co která vlastnost textu vysvětluje a co ne

Průběžná tabulka; čísla z `results/NNN_*/`. **Stav po 001–016** (28. 10. 2026).

Nully a kontroly: **ts2020** (port Timm & Schinner, 004, config 21; Eva-only → jen ZL; zdroje kopií jen z předchozích
řádků téže stránky); **ts2020+R** (ts2020 s naučenými místními tabulkami R1 začátek řádku / R2 konec řádku / R3 šev,
014); **tmpl** (šablona s paradigmaty, tématem, lokalitou, 005); **ac-tri** (003); **pos-shuf**, **line-shuf**, **perm**
(permutace stránek/štítků, uvnitř sekce i globálně); **próza** = latinský narativ, lékařská a botanická latina,
anglický herbář, italština (009); **katalog** = Linnaeus po položkách (009); **blábol** = 38 pisatelů Gaskell & Bowern,
10 k tokenů, srovnáváno s Voynichem oříznutým (013); **Naibbe** = šifra latiny (Greshko), referenční bod, ne hypotéza
(013); **latina téma/registr** = Plinius knihy 16–21 vs 23–30, Caesar vs Celsus (016).

## 1. Vlastnosti, které má Voynich a nemá je žádný skutečný text, blábol, šifra ani generátor

| Vlastnost | Reál ZL split | nejbližší | Co o ní víme (014–016) |
|---|---|---|---|
| **Poziční vazby afixů** I(prefix; pozice) | 0,134 bit | ts2020 0,063; próza ≤ 0,003; katalog 0,042; blábol 0,006; Naibbe 0; **ts2020+R1 dožene** | **malá sada konvencí**: 8 prefixů (ZL i GC, split i join) nese 80 % informace, 5 polovinu (`s-/d-/y-/∅`, `t-/p-` na začátku řádku, `p-` v prvním řádku odstavce); 60–70 % informace je kraj řádku, 25–30 % první řádek odstavce, < 10 % drift uvnitř; směr pravidel sdílený A/B, písaři i sekcemi (83–100 % položek), ale **váhy** A vs B 5–15× dál než poloviny téže skupiny; vysvětlí jen 2–5 % variace tvaru při daném kmeni (011) |
| **Šev slov**: okrajová MI (poslední→první glyf) a I(sufix_t; prefix_{t+1} \| kmen, pozice, …) | 0,185 bit; 0,038 bit | ts2020 0,007 / 0,004; próza 0,02–0,07 / 0,006–0,010; blábol 0,013 / 0,011; Naibbe 0,005 / 0,001; **ts2020+R3 dožene 2/3 švu** (0,126) a sufix→prefix celý | **4 páry nesou 80 %** (y→q 45 %, r→a 22 %, s→a 10 %, iin→o 8 %; GC 7 párů), sufix→prefix = 2 páry (dy→qo, y→qo); **nad glyfovým párem nezbývá nic** (≤ 0,001 bit; latinský katalog 17 %) – šev je glyfová konvence, ne shoda tvarů; jen d = 1, končí na hranici řádku |
| **Konec řádku a uzavřený řádek**: `-m/-g` na konci 17,5×; vazba přes hranici řádku / uvnitř 0,02 | – | ts2020 7×; ts2020+R2 6–8× (píše `-m` i uvnitř); próza 1–2×, poměr 0,6–1,3; blábol 2×, 0,79; katalog 0,05–0,13 | poslední glyf `m` sám nese 49 % koncové informace; řádek jako jednotka (Currier); položkový rys bez položkového obsahu (010) |
| **A/B jako dva slovníky s ostrou hranicí** (012, 016) | tvary jen v B = **0,21–0,23 tokenů B** (poloviny téhož jazyka 0,05–0,06; herbář 0,24–0,27); jen v A 0,15 | ts2020 / tmpl: A×B = poloviny; blábol: pisatelé se liší 5× víc (každý jinak); **latina jiné téma 0,19, jiný registr 0,30** | B-only tvary **nejsou glyfové varianty A**: protějšek ≤ 1 editace 0,81–0,85 vs poloviny 0,90–0,94 (herbář −0,2), tj. dál než náhoda; jediná systematická operace **vložení `d`** (`-ey → -edy`, v101 `8`) nese 28–46 % protějšků, zbytek dlouhý ohon; přepis top-10 ani všemi ~300 operacemi sdílení na poloviny nevrátí (Jaccard 0,25 → 0,31 / 0,39 vs 0,42–0,44); A-only tvary vs B ≈ poloviny (B = A + něco navíc); latina s jiným tématem/registrem dává **stejný obraz** (B-only dál než poloviny, plochý inventář). Věta 016: obojí, převaha nového slovníku. Přežívá písaře (005), pozici (008); jazyk vs pisatelský okruh nerozlišeno |
| **Kombinace** pozice + šev + vazba prefix–sufix (0,205 bit) + šířka paradigmatu zároveň | – | ts2020+R1+R2+R3 dožene pozici a 2/3 švu, ale **pokazí** prefix–sufix (0,094), sekční vazbu prefixu, záhlaví řádku, přidá 2 000 typů; paradigma, pořadí ani A/B nemění | kopírování + lokální přepis nedá pozici, šev a vnitroslovní soudržnost najednou (014); přepis podmíněný i sufixem (R3′) netestován → 017 |

## 2. Vlastnosti, které Voynich sdílí s lidským blábolem (a ne s jazykem ani generátory)

| Vlastnost | Reál | blábol | próza | ts2020 / tmpl | Poznámka |
|---|---|---|---|---|---|
| Šířka a vyrovnanost paradigmat: sufixů na kmen; H(sufix \| kmen) | 5,25 (oříznuto 4,9); 1,82 | 3,96; **1,89** | ≤ 2,4; ≤ 0,88 | 3,8; 1,65 / 4,4; 2,8 | Naibbe také (4,5; 2,0) – vzniká i z více zápisů téhož písmene |
| >90 % variace tvaru při daném kmeni nevysvětlí pozice, registr, písař, sekce (011) | 92 % zbytek | zbytek ~65 % (pisatel vysvětlí 30–35 %) | latina: zbytek ~98 % (variace je flexe) | ts2020 98 % | zbytek má jen slabou strukturu: šev (výše) a malá konzistence uvnitř stránky |
| **Tvar opakování nezávislý na místě** (015): identický tvar téhož vzácného kmene na téže stránce / na dálku proti paradigmatu | −0,02 / −0,00 (vše v intervalu nuly) | **nablízko +0,07**, na dálku −0,12 | nablízko +0,04–0,05, na dálku −0,01 až −0,02 (pád) | ts2020 +0,08 / −0,04; tmpl 0 / 0 | Voynich volí tvar jako by nezávisle na místě – „nezávislejší“ než latina, blábol i kopírování; shoda s 002/011 (varianty rozptýlené); n na stránce malé (34–47) |
| Hapax | 0,676 | 0,85 | 0,47–0,67 | 0,58 / 0,61 | Naibbe opačně (0,41) |
| Slabé pořadí tokenů, MI kmenů na d = 1 | 0,9 %; 0,051 | ano (kmeny) | ≥ 2,2 %; ≥ 0,091 | 0,44 % / 0 | Naibbe také (0,012; 0,029) |

## 3. Vlastnosti, které Voynich sdílí s prózou a nemají je generátory, šablona ani blábol (015)

| Vlastnost | Reál (ZL / GC) | próza | generátory | blábol / katalog |
|---|---|---|---|---|
| **Pozvolný pokles opakování vzácných kmenů (2–5×) přes 1–15 stránek**: poměr pozorované/očekávané po binech 0 / 1 / 2–3 / 4–7 / 8–15 / 16–31 / 32+ | 2,9 / 2,0 / 1,6 / 1,4 / 1,2 / 1,0 / 0,87 (GC 3,1 / 1,7 / 1,4 / 1,45 / 1,3 / 1,1 / 0,85; kmeny 6–20 a dekompozice b totéž) | Plinius 7,8 / 3,9 / 2,1 / 1,6 / 1,2 / 0,86 / 0,76; Celsus podobně | ts2020 **35 / 1,6 / 0,8 / 0,7 / 0,75 / 0,75 / 0,7**; tmpl 26 / 0,8 / plochý; ac-tri 52 / 4,4 / 0,5 | blábol 13 / 4,5 / 0,8 / 0,5 …; katalog 37 / 4,7 / 1,4 / 1,0 / 0,7 – pokles jen do 1 stránky |
| **Dálkové dvojice (≥ 16 stránek) v téže sekci / jazyce / písaři nad náhodu** | +0,08 (0,236 vs 0,155) / +0,06 (0,207 vs 0,149) – všechny tři přebytky stejné | Celsus +0,07 (tematicky členěná próza) | tmpl s tématem +0,02; ts2020 +0,01 | – |
| Dálkový přebytek sám (podíl dvojic ≥ 16 − permutace) | −0,088 / −0,097 (záporný: opakují se **blíž** než náhoda) | −0,18 / −0,21 | −0,25 / −0,17 | −0,33 / −0,34 |
| Shluky vzácného kmene přes ≥ 3 kvaterny | 0,83 / 0,77 očekávání | ≈ tmpl / próza | – | – |

Permutace stránek uvnitř sekce dá 1,33 / 1,37 / 1,29 / 1,19 / 1,05 / 0,91 – reál má navíc **lokalitu uvnitř sekce**
(2,0 → 1,2). Profil 1–15 stránek je slučitelný s řazením podle tématu i s postupnou změnou slovníku psaní (drift);
kopírování s pamětí přes stránky by jej mohlo dát z konstrukce – ts2020 ji nemá (zdroje jen z téže stránky) → 018.
Popiskové kmeny (41–55) mají plošší profil, ale na závěr nestačí.

## 4. Vlastnosti vysvětlené nully (vyřazeny jako indicie)

`cmi_excess`, `mi_prefix_suffix`, produktivita, `jsd_stem_sections`, disperze kmenů (ts2020, tmpl) · sekční
informace kmenů v písaři (ts2020 80 %) · tvar MI(d) (všechny nully) · polovina poziční informace (ts2020), zbytek
R1/R2 (014) · vertikální podobnost (pos-shuf) · záhlaví a šablony řádků (pozice, próza) · sekční složka sdílení nad
vzdálenost, stopa popisků v textu (nic) · sekční profil kmene mezi A a B (nestabilní i u sdílených kmenů; test bez síly,
016) · Zipf, délky, typy, H2, blízcí sousedé (kalibrace).

## 5. Obraz po šestnácti experimentech

* **Místní pravidla jsou malý, sdílený systém.** Pozice v řádku, šev slov a konec řádku = řádově deset glyfových
  pravidel se stejným směrem u obou jazyků, všech písařů a sekcí; A a B je uplatňují s jinými vahami. Nad glyfovým
  párem na švu nezbývá nic – šev nese zvyk zápisu, ne shodu tvarů. Kopírovací generátor s těmito tabulkami dožene
  pozici a většinu švu, ale rozbije vazbu prefix–sufix: **kopírování + lokální přepis zatím nedá pozici, šev a
  soudržnost slova zároveň** (otevřeno: 017).
* **Vnitroslovní variace je jako u lidského blábolu** (široká paradigmata, hapax, tvar nevysvětlený ničím
  měřitelným) a navíc **nezávislá na místě**: ani kopie nablízko, ani jiný tvar na dálku – v tom je Voynich
  „nezávislejší“ než latina, blábol i kopírování.
* **Vzdálenostní profil je jako u prózy, ne jako u generátorů ani blábolu:** vzácné kmeny se vracejí s pozvolným
  poklesem přes sousední listy a kvatern, dálkové dvojice sedí do téže sekce/jazyka/písaře stejně jako u tematicky
  členěné latiny. Obsahové opakování na dálku (nad náhodu) ale **neexistuje**. Jediné, co odlišuje drift od tématu,
  je tvar lokality (hladký pokles vs bloky) – ještě neměřeno.
* **A a B jsou dva slovníky, ne dva zápisy:** pětina tokenů B jsou tvary, které A nemá a které nejsou jeho
  glyfovými variantami (až na jednu operaci, `-ey → -edy`). Totéž dává latina téhož jazyka s jiným tématem nebo
  registrem; generátory bez jazyka nic. Zda je to obsah (jiné věci), registr (jiný druh textu) nebo pisatelský
  okruh, data zatím nerozliší.
* Souhrnně: rukopis vypadá jako **systém s pevnými, sdílenými pravidly vázanými na místo, volnou obměnou tvarů uvnitř
  slov, dvěma slovníky a prózovitou lokalitou slovníku**. Nic z toho není důkaz obsahu; proti „bez obsahu“ ale stojí
  tři věci, které nedal žádný generátor ani blábol: lokalita přes 1–15 stránek, sekční vazba dálkových dvojic a nový
  slovník B – všechny tři jsou slučitelné i s driftem/registrem. Příští kolo míří přesně na ně.

## 6. Otevřené body (zadáno 017–020)

017 mechanismus: podmíněný přepis R3′ (prefix | sufix, předchozí glyf) – dá lokální mechanismus pozici, šev i vazbu
prefix–sufix zároveň? · 018 lokalita: drift (hladký pokles; ts2020 s pamětí přes stránky) vs bloky (segmentace,
hranice sekcí a odstavců na úrovni stránek), fyzické pořadí (kvatern, bifolio, anomálie sousedů) · 019 tvary jen v B:
`-edy` varianty vs zbytek – jména/položky (úzké paradigma, popisky, hapax) nebo jiná morfologie; pravidlo `-ey → -edy`
deterministické? · 020 odstavec jako jednotka obsahu: sdílení kmenů mezi řádky vs hranice odstavce a stránky ·
blábol jednoho pisatele na dlouhém textu (není korpus) · párování stránek podle ilustrace (chybí data) ·
průběžná zpráva (001–016).
