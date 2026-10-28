# 020 – Odstavec a stránka jako jednotky obsahu

Větev `exp/020-paragraph-as-unit`; kód `experiments/020_paragraph_as_unit/run.py` + `summarize.py`, modul
`voytools/paragraphs.py` (testy `tests/test_paragraphs.py`), kontroly s vlastními odstavci
`voytools/controls.py::paragraph_control` (testy v `tests/test_controls.py`). Volby při výkladu spec a popis
kontrol: `experiments/020_paragraph_as_unit/notes.md`. Čísla níže jsou z `boundary_step.csv`,
`line_pair_sharing.csv`, `within_par_decay.csv`, `first_line.csv`, `page_boundary.csv`, `form_by_boundary.csv`
a `summary.csv` (běhy po zdrojích v `runs/<combo>/<source>/runs.csv`, `meta.json`).

## Co se měřilo

Řádek = lokus `P` (clean); **sdílení** dvojice řádků téže stránky = počet různých kmenů společných oběma, dělený
očekáváním při délkách řádků z frekvencí kmenů **na dané stránce** (hlavní; „struktura nad slovníkem stránky“),
vedlejší z globálních frekvencí. Poměr s(d) = Σ pozorované / Σ očekávané pro dvojice ve vzdálenosti d řádků, zvlášť
v témž odstavci (`same`) a přes hranici (`diff`); **skok** Δ(d) = s_same(d) − s_diff(d) na stránkách s ≥ 2 odstavci
(ZL 159, GC 159). Kmeny = dekompozice (a) s ≥ 2 výskyty (varianta `st`; ZL split 752 kmenů, 32 675 tokenů), bez 50
nejčastějších (`st50`), celé tokeny (`tok`) a kmeny s odstraněnými počátečními gallows (`ng`). Navíc **Δ bez prvních
řádků** (dvojice, v nichž žádný řádek není první v odstavci; kontrola artefaktu prvního řádku).

Zdroje: reál (interval = bootstrap stránek 200×), null (a) přelosování pořadí délek odstavců na stránce (200×),
null (b) permutace řádků na stránce (200×), ts2020 config 21 (30, jen ZL), šablona 005 s tématem (30, ZL i GC);
kontroly v layoutu ZL split: Celsus a Plinius s **vlastními odstavci = kapitolami** (`celsus_chap`, `pliny_chap`; konec
kapitoly ukončí řádek), totéž s odstavci ze skupin vět ≥ 40 tokenů uvnitř kapitoly (`*_sent40`, slabá kontrola),
Linnaeus položka = odstavec (`linnaeus_par`), Plinius próza s hranicemi ZL (`pliny_prose`, hranice bez významu),
blábol s hranicemi ZL (`gibberish_pooled`) a s vlastními stránkami (`gibberish_native`, jen hranice stránky; přepisy
odstavce nenesou). p = podíl seedů nullu aspoň tak extrémních jako reál (nejmenší možné 0,005 u 200, 0,032 u 30).

Hodnoty níže: ZL split [95% interval]; GC split v závorce, join tamtéž, kde se liší.

## 1. Skok na hranici odstavce (`boundary_step.csv`, `line_pair_sharing.csv`, `boundary_step.png`)

Δ(d), kmeny `st`, očekávání ze stránky:

| zdroj | Δ(1) | Δ(2) | Δ(3) | Δ(4) | Δ(2) bez prvních řádků | Δ(3) bez prvních |
|---|---|---|---|---|---|---|
| **ZL split, reál** | **0,146 [0,099; 0,192]** | 0,052 [0,015; 0,089] | 0,012 [−0,027; 0,054] | 0,024 [−0,022; 0,063] | −0,013 [−0,058; 0,032] | 0,001 [−0,045; 0,043] |
| ZL join, reál | 0,175 [0,117; 0,222] | 0,047 [0,006; 0,095] | −0,006 | 0,016 | −0,030 [−0,084; 0,027] | −0,008 |
| **GC split, reál** | **0,115 [0,066; 0,171]** | 0,029 [−0,007; 0,072] | −0,013 | −0,010 | −0,012 [−0,052; 0,049] | −0,019 |
| GC join, reál | 0,141 [0,090; 0,200] | 0,045 [0,003; 0,097] | 0,008 | −0,002 | 0,018 [−0,028; 0,086] | 0,007 |
| (a) přelosování hranic, ZL | 0,035 [−0,003; 0,070] p 0,005 | 0,007 p 0,005 | −0,011 | 0,013 | 0,009 | 0,018 |
| (b) permutace řádků, ZL | 0,003 [−0,042; 0,045] p 0,005 | 0,005 p 0,015 | 0,009 | 0,013 | 0,011 | 0,015 |
| (b) permutace řádků, GC | −0,006 p 0,005 | −0,005 p 0,020 | −0,002 | 0,001 | −0,003 | 0,005 |
| **ts2020 cfg 21 (ZL)** | **0,232 [0,144; 0,305]** | 0,090 [0,058; 0,135] | 0,035 [0,009; 0,071] | 0,028 | **−0,026 [−0,076; 0,028]** | −0,020 |
| šablona 005 (ZL / GC) | −0,007 / −0,019 | −0,009 / −0,014 | −0,003 / −0,008 | 0,001 / −0,008 | −0,009 / −0,014 | −0,000 / −0,003 |
| **Celsus, kapitoly** | **0,385 [0,232; 0,589]** | **0,261 [0,113; 0,417]** | **0,179 [0,032; 0,347]** | 0,175 [0,078; 0,275] | 0,189 [−0,084; 0,446] | 0,180 [−0,046; 0,386] |
| **Plinius, kapitoly** | **0,212 [0,087; 0,347]** | **0,205 [0,117; 0,296]** | **0,146 [0,070; 0,222]** | 0,078 [0,005; 0,170] | **0,385 [0,232; 0,500]** | 0,136 [0,049; 0,218] |
| Celsus, věty ~40 | 0,067 [−0,029; 0,163] | 0,029 | −0,055 | −0,061 | 0,049 | −0,052 |
| Plinius, věty ~40 | 0,136 [0,062; 0,216] | 0,096 [0,031; 0,174] | 0,035 | −0,051 | 0,033 [−0,089; 0,146] | 0,009 |
| Linnaeus, položky | **−0,738 [−0,823; −0,659]** | −0,501 | −0,458 | −0,682 | −0,575 | −0,172 |
| Plinius próza, hranice ZL | −0,112 [−0,217; −0,001] | −0,071 | −0,029 | 0,004 | −0,090 | 0,007 |
| blábol, hranice ZL | −0,024 [−0,187; 0,155] | −0,174 | −0,100 | −0,145 | −0,236 | −0,219 |

Přelosování hranic (a) u kontrol s kapitolami dává Δ(1) 0,08–0,11 (p ≤ 0,02), permutace řádků 0,04–0,09; u Linnaea
−0,05 až −0,09.

Varianty (reál, Δ(1) / Δ(2)): `ng` (bez gallows) ZL **0,108 [0,060; 0,162]** / 0,025, GC 0,069 [0,016; 0,125] / 0,006
(ts2020 0,148 / 0,068); `tok` ZL 0,210 [0,126; 0,307] / 0,040, GC 0,133 / 0,027 (ts2020 0,268 / 0,106); `st50` ZL 0,185
[0,029; 0,337] / 0,057, GC 0,121 [−0,035; 0,257] / 0,131 (ts2020 0,592 / 0,230; intervaly široké, sdílených vzácných
kmenů je málo). Globální očekávání dává Δ(1) 0,20 (ZL) / 0,17 (GC), ale i permutace řádků 0,05 a přelosování 0,09 –
je zkreslené slovníkem stránky (viz `notes.md`), proto je hlavní očekávání stránkové.

Samotné s(d) (ZL split, reál, `same` / `diff`): d = 1 0,94 / 0,79, d = 2 0,92 / 0,87, d = 3 0,88 / 0,87, d = 4 0,88 /
0,87, d = 8 0,91 / 0,87; permutace řádků ≈ 0,84 všude. Pokles s(d) v témž odstavci mezi d = 1 a 3 je stejný jako u
dvojic přes hranici mezi d = 2 a 3; jediný velký rozdíl je `diff` při d = 1 (dvojice poslední řádek – první řádek
dalšího odstavce).

**Výklad.** Skok je v ZL i GC, split i join, nad oběma nully (p = 0,005), ale (i) je **jen při d = 1** (Δ(2) ≈ 0,03–0,05,
Δ(3), Δ(4) ≈ 0), (ii) **zmizí, když se vyřadí první řádky odstavců** (Δ bez prvních řádků −0,03 až +0,02 ve všech
čtyřech kombinacích), (iii) **ts2020 ho dává větší než reál** (0,23 vs 0,15; bez prvních řádků rovněž ≈ 0) – port
kopíruje první řádky odstavce ze zásobníku prvních řádků a přidává gallows (`choose_source`, režim
`paragraph_initial`), tj. skok z konstrukce přes první řádek, a (iv) odstranění gallows z kmenů ho sníží o třetinu
(0,146 → 0,108; GC 0,115 → 0,069). Próza s **vlastními kapitolami** se chová jinak: skok trvá přes d = 1…4 (Plinius
0,21 / 0,21 / 0,15 / 0,08, Celsus 0,39 / 0,26 / 0,18 / 0,18) a **zůstává bez prvních řádků** (Plinius 0,39 / 0,14,
Celsus 0,19 / 0,18) – to je stopa obsahové jednotky. Odstavce ze skupin vět (~ 40 tokenů, velikost odstavce Voynichu)
dávají jen slabý skok (Plinius 0,14 / 0,10 / 0,04, bez prvních řádků 0,03), tj. hranice nad kapitolu obsah nedělí.
Šablona 005 (téma na stránku) Δ ≈ 0, jak spec čeká. Skok Voynichu je tedy **hranový** (první řádek odstavce má jiné
kmeny, z větší části kvůli gallows), ne obsahový.

## 2. Uvnitř odstavce (`within_par_decay.csv`)

s(d) jen v témž odstavci, kmeny `st`, očekávání ze stránky; sklon na řádek (vážený počtem dvojic):

| zdroj | s(1) | s(2) | s(3) | s(4) | s(6) | sklon [interval] |
|---|---|---|---|---|---|---|
| **ZL split, reál** | 0,92 | 0,89 | 0,87 | 0,88 | 0,87 | **−0,012 [−0,020; −0,006]** |
| **GC split, reál** | 0,93 | 0,89 | 0,86 | 0,85 | 0,86 | **−0,016 [−0,024; −0,008]** |
| ZL join / GC join | 0,91 / 0,91 | 0,87 / 0,87 | 0,84 / 0,84 | 0,85 / 0,83 | 0,87 / 0,86 | −0,011 / −0,015 |
| permutace řádků (ZL / GC) | 0,84 / 0,85 | 0,84 / 0,85 | 0,84 / 0,85 | 0,85 / 0,85 | 0,86 / 0,85 | 0,003 / −0,000 (p 0,005) |
| ts2020 (ZL) | 0,98 | 0,96 | 0,94 | 0,93 | 0,93 | −0,012 [−0,018; −0,007] |
| šablona 005 (ZL / GC) | 0,94 / 0,93 | 0,94 / 0,92 | 0,94 / 0,92 | 0,95 / 0,92 | 0,95 / 0,93 | 0,002 / 0,000 |
| Celsus kapitoly / Plinius kapitoly | 0,74 / 0,64 | 0,68 / 0,69 | 0,70 / 0,64 | 0,67 / 0,65 | 0,65 / 0,64 | −0,017 [−0,024; −0,010] / −0,006 [−0,017; 0,004] |
| Celsus věty / Plinius věty | 0,73 / 0,61 | 0,69 / 0,68 | 0,63 / 0,63 | 0,62 / 0,57 | 0,45 / 0,48 | −0,048 / −0,022 |
| Plinius próza, hranice ZL | 0,61 | 0,62 | 0,63 | 0,60 | 0,59 | −0,003 [−0,015; 0,011] |

Bez gallows (`ng`) ZL −0,010 [−0,016; −0,003], GC −0,014; celé tokeny ZL −0,009 [−0,021; 0,002], GC −0,018.

**Výklad.** Uvnitř odstavce sdílení **klesá** se vzdáleností (≈ −0,012 až −0,016 na řádek, nad nulou permutace),
ne plocho jako pytel slov; stejný sklon má ts2020 (−0,012, kopírování z předchozích řádků) a Celsus (−0,017), šablona
(téma na stránku) 0. Pokles ale není vázaný na odstavec: tutéž hodnotu dává přelosování hranic (−0,011) a dvojice
přes hranici klesají s d stejně (bod 1) – jde o blízkost řádků na stránce, ne o odstavec.

## 3. První řádek odstavce (`first_line.csv`)

Sdílení prvního řádku s dalšími řádky vlastního odstavce (`first_own`) a s řádky následujícího odstavce
(`first_next`) v téže vzdálenosti d; totéž pro prostřední řádky (`mid_own`, `mid_next`); kmeny `st`:

| zdroj | d | first_own | first_next | mid_own | mid_next | „ohlašuje“ = (f_own − f_next) − (m_own − m_next) |
|---|---|---|---|---|---|---|
| **ZL split, reál** | 1 | **0,77 [0,73; 0,81]** (n 741) | – | **0,97 [0,94; 0,99]** (n 2 646) | – | – |
| | 2 | 0,78 [0,73; 0,82] | 0,79 [0,69; 0,92] (n 71) | 0,94 | 0,78 [0,74; 0,82] | −0,17 |
| | 4 | 0,70 [0,65; 0,76] | 0,75 [0,70; 0,79] (n 327) | 0,93 | 0,88 | **−0,09** |
| | 6 | 0,60 [0,53; 0,67] | 0,68 [0,61; 0,75] | 0,94 | 0,88 | −0,13 |
| **GC split, reál** | 1 | 0,79 [0,75; 0,84] | – | 0,97 [0,95; 0,99] | – | – |
| | 4 | 0,68 [0,62; 0,73] | 0,80 [0,75; 0,85] | 0,91 | 0,89 | **−0,15** |
| ZL join / GC join | 4 | 0,66 / 0,68 | 0,71 / 0,75 | 0,91 / 0,89 | 0,87 / 0,87 | −0,09 / −0,09 |
| ZL split, **bez gallows** (`ng`) | 1 | 0,87 [0,83; 0,90] | – | 0,95 [0,92; 0,97] | – | – |
| | 4 | 0,78 [0,72; 0,85] | 0,82 [0,78; 0,86] | 0,89 | 0,90 | −0,04 |
| GC split, bez gallows | 1 / 4 | 0,87 / 0,75 | – / 0,84 | 0,95 / 0,89 | – / 0,89 | – / −0,08 |
| ts2020 (ZL) | 1 / 4 | 0,87 / 0,78 | – / 0,81 | 1,02 / 0,98 | – / 0,95 | – / −0,05 |
| šablona 005 (ZL) | 1 / 4 | 0,94 / 0,93 | – / 0,95 | 0,94 / 0,95 | – / 0,95 | – / −0,02 |
| permutace řádků (ZL) | 1 / 4 | 0,84 / 0,83 | – / 0,84 | 0,84 / 0,85 | – / 0,86 | – / 0,00 |
| Plinius kapitoly | 4 | 0,62 [0,52; 0,72] | 0,47 [0,33; 0,66] (n 75) | 0,65 | 0,69 | +0,19 |
| Celsus věty / Plinius věty | 2 | 0,71 / 0,72 | 0,51 / 0,47 | 0,68 / 0,66 | 0,71 / 0,58 | +0,22 / +0,17 |
| | 4 | 0,66 / 0,60 | 0,60 / 0,50 | 0,58 / 0,55 | 0,71 / 0,68 | +0,20 / +0,23 |
| Plinius próza, hranice ZL | 4 | 0,62 | 0,71 | 0,59 | 0,61 | −0,07 |

**Výklad.** První řádek odstavce sdílí s vlastním odstavcem **méně** než prostřední řádky (0,77 vs 0,97 při d = 1; ZL i
GC), a s následujícím odstavcem **stejně nebo víc** než s vlastním (0,75 vs 0,70 při d = 4) – **neohlašuje** slovník
odstavce, je to kraj. Po odstranění gallows se rozdíl zmenší na polovinu (0,87 vs 0,95), ale nezmizí: první řádek má
jiné kmeny i mimo gallows (shoda s 014: `p-` v prvním řádku odstavce, `t-/p-` na začátku řádku, drift uvnitř řádku).
ts2020 dává totéž (0,87 vs 1,02) z konstrukce, šablona nic. Próza s vlastními odstavci má znaménko opačné: první řádek
sdílí s vlastním odstavcem víc než s dalším (+0,17 až +0,23 u vět, +0,19 u Pliniových kapitol při d = 4).

## 4. Hranice stránky (`page_boundary.csv`, `page_boundary.png`)

Poslední řádky stránky i vs první řádky stránky i + 1 (po sobě jdoucí stránky `P`; ZL 199 dvojic při d = 1), globální
očekávání, poměr k dvojicím uvnitř stránky v téže vzdálenosti; kmeny `st`:

| zdroj | přes/uvnitř d = 1 | d = 2 | d = 3 | i + 1 téhož jazyka a sekce, d = 1 (n 159) | jiného (n 40) | uvnitř stránky přes hranici odstavce / uvnitř odstavce, d = 1 |
|---|---|---|---|---|---|---|
| **ZL split, reál** | **0,74 [0,63; 0,84]** | 0,82 [0,75; 0,90] | 0,89 [0,84; 0,94] | 0,76 [0,63; 0,86] | 0,64 [0,43; 0,86] | 0,84 |
| ZL join | 0,70 [0,58; 0,82] | 0,80 | 0,87 | 0,71 | 0,65 | 0,81 |
| **GC split, reál** | **0,68 [0,58; 0,80]** | 0,84 [0,76; 0,91] | 0,90 [0,83; 0,96] | 0,71 [0,60; 0,82] | 0,54 [0,28; 0,83] | 0,87 |
| GC join | 0,61 [0,49; 0,74] | 0,79 | 0,88 | 0,62 | 0,53 | 0,84 |
| ZL split, bez gallows (`ng`) | 0,81 [0,71; 0,90] | 0,91 [0,83; 0,98] | 0,93 [0,87; 0,97] | 0,83 | 0,70 | 0,88 |
| GC split, bez gallows | 0,73 [0,64; 0,84] | 0,93 | 0,94 | 0,75 | 0,63 | 0,91 |
| ZL split, celé tokeny | 0,58 [0,38; 0,77] | 0,75 | 0,88 | 0,64 | 0,27 [0,00; 0,70] | 0,67 |
| permutace řádků (ZL / GC) | 0,93 [0,84; 1,03] p 0,005 / 0,93 p 0,005 | 0,93 p 0,005 / 0,93 | 0,93 p 0,055 / 0,93 p 0,090 | 0,96 / 0,96 | 0,80 / 0,80 | 0,96 / 0,97 |
| ts2020 (ZL) | 0,99 [0,88; 1,10] | 0,82 [0,72; 0,90] | 0,78 [0,71; 0,87] | 0,99 | 0,99 | 0,74 |
| šablona 005 (ZL / GC) | 0,95 / 0,93 | 0,94 / 0,92 | 0,94 / 0,92 | 0,96 / 0,94 | 0,87 / 0,86 | 0,99 / 1,01 |
| Celsus kapitoly / Plinius kapitoly | 0,99 [0,77; 1,22] / 1,02 [0,76; 1,25] | 1,06 / 0,98 | 1,01 / 1,10 | 1,03 / 0,91 | 0,83 / 1,52 | 0,57 / 0,72 |
| Celsus věty / Plinius věty | 0,96 / 1,12 | 0,87 / 0,95 | 0,98 / 1,22 | 0,98 / 1,15 | 0,87 / 0,97 | 0,90 / 0,77 |
| Linnaeus položky | 0,70 [0,45; 1,02] | 0,89 | 1,05 | 0,70 | 0,70 | 1,32 |
| Plinius próza, hranice ZL | 0,86 [0,62; 1,12] | 0,96 | 0,94 | 0,98 | 0,34 [0,00; 0,77] | 1,10 |
| blábol, vlastní stránky | 1,05 [0,65; 1,43] | 0,65 [0,45; 0,83] | 0,68 | 1,55 (n 50) | 0,28 (n 37) | – |

Stránky bez značky `<$>` na posledním odstavci: v ZL jen 2 (f68r1, f84v; `across_open` 0,00, n 2), v GC žádná –
podskupinu „pokračuje?“ nelze vyhodnotit. U kontrol s vlastními odstavci je takových stránek 184–203 a sdílení přes
hranici je tam **vyšší** (1,4–1,5 vs uvnitř stránky), jak má pokračující odstavec.

**Výklad.** Na hranici stránky sdílení **klesá** (0,74 ZL, 0,68 GC při d = 1; permutace řádků 0,93, p = 0,005) a
s odstupem se vrací (0,82, 0,89; při d = 3 už jen těsně pod nullem). Pokles je opět soustředěn do **prvního řádku
stránky** (= první řádek odstavce): bez gallows je 0,81 / 0,73 a při d = 2–3 už ≈ null (0,91–0,94). Hranice odstavce
uvnitř stránky dává při d = 1 pokles 0,84 (ZL) / 0,87 (GC), tj. menší než hranice stránky. Próza s vlastními odstavci
pokračuje přes stránku beze změny (0,96–1,12), stejně jako blábol s vlastními stránkami při d = 1 (1,05). ts2020 při
d = 1 nepadá (0,99; první řádek stránky kopíruje ze zásobníku prvních řádků, který není vázán na stránku) a padá při
d = 2–3 (0,82, 0,78; kopie jen z téže stránky); šablona 0,93–0,95. Voynich má tedy tvar **opačný než ts2020**: zlom v
prvním řádku, pak návrat; stránka i + 1 téhož jazyka a sekce 0,76 vs jiného 0,64 (n 40, intervaly se překrývají).

## 5. Tvar opakování na hranici (`form_by_boundary.csv`)

Podíl identických tokenů mezi dvojicemi výskytů společného kmene minus očekávání z paradigmatu, d ≤ 3:

| zdroj | týž odstavec | přes hranici | rozdíl | n dvojic výskytů |
|---|---|---|---|---|
| **ZL split, reál** | 0,074 [0,062; 0,087] | 0,047 [0,038; 0,058] | **0,027** | 24 730 / 9 522 |
| GC split, reál | 0,069 [0,058; 0,079] | 0,045 [0,036; 0,055] | 0,024 | 21 746 / 10 674 |
| ZL join / GC join | 0,073 / 0,069 | 0,048 / 0,051 | 0,025 / 0,019 | |
| přelosování hranic (ZL / GC) | 0,073 / 0,067 | 0,053 / 0,050 | 0,020 / 0,016 | |
| permutace řádků (ZL / GC) | 0,062 / 0,060 | 0,046 / 0,046 | 0,016 / 0,014 | |
| ts2020 (ZL) | 0,080 [0,049; 0,099] | 0,059 | 0,021 | |
| šablona 005 | −0,001 | −0,001 | 0,000 | |
| Celsus kapitoly / Plinius kapitoly | 0,064 / 0,056 | 0,027 / 0,011 | 0,037 / 0,044 | 3 289 / 168; 3 891 / 528 |
| Celsus věty / Plinius věty | 0,093 / 0,066 | 0,050 / 0,015 | 0,043 / 0,050 | |
| Plinius próza, hranice ZL | 0,064 | 0,059 | 0,005 | |

**Výklad.** Identický tvar je v témž odstavci o 0,02–0,03 častější než přes hranici, ale přelosování hranic dává
0,016–0,020 a permutace řádků 0,014–0,016 (tvarová shoda uvnitř stránky, 011/015) – rozdíl ≤ 0,01 nad nully, na
hranici intervalů; ts2020 totéž (0,021). Nic nad to, co dá stránka.

## Odpovědi (kritérium 2)

| otázka | ZL (split / join) | GC (split / join) |
|---|---|---|
| (i) Δ(d) > 0 nad přelosování a nad ts2020? | Δ(1) 0,146 / 0,175 nad oběma nully (p 0,005), ale **pod ts2020** (0,23 / 0,24) a **jen při d = 1**; bez prvních řádků Δ ≈ 0 (−0,01 / −0,03). Próza s kapitolami 0,21–0,39 přes d = 1…4 i bez prvních řádků (0,14–0,39); katalog položek **záporný** (−0,74) | Δ(1) 0,115 / 0,141, totéž: Δ(2) 0,03 / 0,05, Δ(3+) ≈ 0, bez prvních řádků −0,01 / +0,02; ts2020 pro GC není |
| (ii) uvnitř odstavce plocho, nebo pokles? | pokles −0,012 / −0,011 na řádek, jako ts2020 (−0,012) a Celsus (−0,017); stejný i přes hranici – blízkost, ne odstavec | −0,016 / −0,015 |
| (iii) ohlašuje první řádek? | **ne**: first_own 0,77 vs mid_own 0,97 (d = 1); first_own 0,70 < first_next 0,75 (d = 4); bez gallows 0,87 vs 0,95 (rozdíl na polovinu, nezmizí); próza s vlastními odstavci opačné znaménko (+0,2) | **ne**: 0,79 vs 0,97; 0,68 < 0,80; bez gallows 0,87 vs 0,95 |
| (iv) hranice stránky: pokračování, nebo zlom? | zlom při d = 1 (0,74 vs null 0,93), při d = 3 0,89 ≈ null; bez gallows 0,81 → 0,93; próza 0,96–1,12 (pokračuje), ts2020 0,99 → 0,78 (opačný tvar), šablona 0,95 | 0,68 → 0,90; bez gallows 0,73 → 0,94 |

**Věta:** odstavec a stránka jako jednotka obsahu – **ne; artefakt (konvence) prvního řádku.** Všechno, co se na
hranici odstavce i stránky děje, je soustředěno do prvního řádku (odstavce, stránky), z větší části do jeho gallows;
po vyřazení prvních řádků není mezi dvojicemi v témž odstavci a přes hranici rozdíl, zatímco latinská próza s vlastními
kapitolami takový rozdíl má (přes d = 1…4 i bez prvních řádků) a přes stránku pokračuje beze změny.

Podezřelé hodnoty (kritérium 3): Δ u permutace řádků ≈ 0 ve všech sadách (−0,006 až +0,003; v pořádku); **Δ u ts2020
≥ Δ reálu** (0,23 vs 0,15; ano – artefakt mechanismu prvního řádku, odečteno variantou bez prvních řádků a `ng`);
**Linnaeus bez skoku** – ano, dokonce **silně záporný** (−0,74): položky Linnaea jsou formulové (rodové jméno,
`habitat`, citace) a sousední položky téhož rodu sdílejí víc než řádky uvnitř jedné položky (diagnóza vs synonyma);
nully dávají −0,05 až −0,09, takže je to vlastnost kontroly, ne chyba výpočtu. Jako pozitivní kontrola skoku slouží
próza s vlastními kapitolami, ne katalog.

## Interpretace

Hranice odstavce a stránky nejsou v rukopise hranicemi obsahu: sdílení kmenů mezi řádky závisí na blízkosti a na tom,
zda je řádek první v odstavci, ne na tom, zda řádky leží v témž odstavci. První řádek je zvláštní (méně sdílí s čímkoli,
z poloviny kvůli gallows, z poloviny kvůli dalším tvarům), což je známá místní konvence (014), a kopírovací generátor s
pravidlem pro první řádky ji dává až nadměrně. Latinská próza vysázená do týchž řádků a stránek s vlastními kapitolami
ukazuje, jak by obsahová jednotka vypadala: skok trvající přes několik řádků, nezávislý na prvním řádku, a hladké
pokračování přes stránku – nic z toho Voynich nemá. Lokalita slovníku z 006/015 (pokles přes stránky) tedy není
členěna po odstavcích ani po stránkách; zbývá plynulý drift nebo členění ve větším měřítku (018).

## Co výsledek neukazuje

* Neukazuje, že odstavce *nemají* žádný význam – jen že se neprojevuje sdílením kmenů mezi řádky nad slovník stránky.
  Jednotka obsahu menší než řádek nebo větší než stránka tímto testem projít nemůže.
* Kmeny jsou z naučených afixů; varianta celých tokenů dává stejný obraz (Δ(1) 0,21 jen při d = 1, bez prvních řádků
  0,04), vzácné kmeny (`st50`) mají příliš široké intervaly (Δ(1) 0,19 [0,03; 0,34]).
* Kontroly s kapitolami mají jen 52 (Celsus) / 148 (Plinius) stránek s ≥ 2 odstavci a odstavce 2–3× delší než Voynich;
  varianta s větami má odstavce stejné velikosti, ale hranice uvnitř kapitol nejsou tematické (slabá kontrola, skok
  jen 0,07–0,14 při d = 1).
* Podskupina „poslední odstavec bez `<$>`“ má v ZL 2 stránky a v GC žádnou – otázka pokračování přes stránku se dá
  posoudit jen z poměru přes/uvnitř (zlom jen v prvním řádku), ne z transkripčních značek.
* ts2020 existuje jen pro ZL; pro GC je generátorový null jen šablona 005. Permutace řádků zachovává slovník stránky,
  proto u hranice stránky dává 0,93, ne 1 (lokalita stránek z 006); poměry jsou srovnatelné jen uvnitř sady.
* Očekávání je s vracením a z frekvencí celé stránky (včetně obou řádků); zkreslení je stejné pro reál i nully, ale
  absolutní s(d) < 1 nemá samostatný výklad.

## Další krok

* Oddělit „první řádek“ od „hranice“ přímo: profil kmenů prvních řádků odstavce (po odstranění gallows) proti ostatním
  řádkům – je to pevná sada tvarů (konvence jako `p-`), nebo jiný slovník (jména, nadpisy)? Souvisí s 019
  (popisky/jména).
* Lokalita přes stránky (006/015) bez odstavcové struktury: 018 (drift vs bloky) může hranice odstavců i stránek
  vypustit a pracovat jen s pořadím řádků.
* Kontrola s obsahovými odstavci stejné délky jako ve Voynichu (např. receptář nebo herbář po položkách s více řádky,
  Culpeper po odstavcích v `paragraph_control`), aby se skok obsahové jednotky změřil při 4 řádcích na odstavec.
