# 008 – Kanonický null ts2020 a šablona 005 v analýzách 006/007, „poziční kmen“ a A/B podle pozice

Spec: `experiments/008_ts2020_recheck/spec.md`. Kód: `experiments/007_position_vs_content/run.py` a
`experiments/006_topic_at_distance/run.py` (přidané zdroje, výstupy dosavadních zdrojů beze změny),
`experiments/008_ts2020_recheck/run.py` (kroky 2, 4, 5, kontrola, souhrny). Rozhodnutí a odchylky:
`experiments/008_ts2020_recheck/notes.md`. Všechna čísla níže jsou ze souborů v této složce.

**Nully.** **ts2020** = `selfcitation_like` s config 21 z 004 (p_copy 0,1, replace_weights 20/50/30,
p_same_position 14). V 006/007 běží 30 seedů 2100–2129, což jsou tytéž vzorky jako ve fázi B 004.
V řádcích 1–6 syntézy se používá 50 seedů z 004. Generátor vyrábí jen Eva, takže pro GC null chybí.
**tmpl** = plný model 005 (fixed paradigma, α = 5, λ = 0,3; ZL 8 000 kmenů / zipf 1,4, GC 4 000 / 1,3),
30 seedů 100000+ (vzorky ablace `best` 005), ZL i GC. Pro srovnání jsou v tabulkách i dřívější nully
(autocopy 86, permutace 006).
**Reprodukuje** = |z| < 3, kde z = (reál − průměr nullu) / SD nullu přes seedy.

**Důležité zjištění v kódu portu** (`voytools/selfcitation.py`, `choose_source`): ts2020 nekopíruje jen
uvnitř stránky, jak předpokládala spec. První řádek stránky bere zdroj z posledního řádku předchozí stránky.
Počáteční řádky odstavců berou zdroj z kteréhokoli dřívějšího počátečního řádku odstavce.

## Kritérium 3 – reprodukce dřívějších čísel

`check_reproduction.json`:

* **Rozklad 007 a MI(d) 007**: přepočet na reálném textu dává tytéž hodnoty (max. rozdíl 1·10⁻¹⁶).
* **A/B 005**: nepodmíněné JSD (`ab_by_position.csv`, `stat = all`) se shodují s `results/005_…/ab/*/real`
  přesně (rozdíl 0): ZL split 0,0883, písař 3: 0,0676.
* **Profil syntézy** (ZL split, `profile/`): 0,1846 / 0,1856 / 1,0579 = `results/004_…/real/ZL_split/point.csv`.
* **Parser, výchozí varianta**: oproti `I(kmen; ·)` v 007 max. rozdíl 0,0009 bitu. Rozdíl dělá jen jiné
  pořadí náhodných čísel v korekci mícháním.
* **006 `real`**: matice podobností jsou totožné s `results/006_…/pairs.csv`. Koeficienty se liší
  max. o 0,0007 (K = 30), 0,0015 (K = 0) a 0,012 (K = 100, GC, jen R²). Důvodem je jiná verze pandas
  (5 stránek bez jazyka → `same_language` u 10 z 20 100 dvojic) a pořadí kmenů se shodnou četností při
  K = 100 (`notes.md`). Srovnání s nully používá reál ze souborů 006.

## Krok 1 – 007 s ts2020 a šablonou (`decomp_ts2020.csv`, `mi_distance_ts2020.csv`)

Korigovaná MI v bitech. U nullu je uveden průměr přes 30 seedů, v závorce % reálu a z.
ZL split (ZL join se liší do ±5 procentních bodů, viz CSV):

| Y | člen | reál | autocopy 86 | **ts2020** | tmpl |
|---|---|---|---|---|---|
| prefix | I(Y; pozice) | 0,134 | 0,000 (0 %) | **0,063 (47 %, z 17)** | 0,000 (0 %) |
| prefix | I(Y; sekce \| pozice) | 0,104 | 0,004 (4 %) | 0,013 (12 %, z 12) | 0,001 (1 %) |
| sufix | I(Y; pozice) | 0,048 | 0,000 | **0,028 (59 %, z 4,1)** | 0,000 |
| sufix | I(Y; sekce \| pozice) | 0,080 | 0,003 (4 %) | 0,013 (16 %, z 14) | 0,002 (3 %) |
| kmen | I(Y; pozice) | 0,191 | 0,001 (0 %) | **0,109 (57 %, z 9,0)** | 0,002 (1 %) |
| kmen | I(Y; sekce) | 0,119 | 0,025 (21 %) | 0,053 (44 %, z 6,6) | 0,101 (85 %, z 0,9) |
| kmen | I(Y; sekce \| pozice) | 0,117 | 0,012 (11 %) | 0,035 (30 %, z 8,6) | 0,079 (68 %, z 1,8) |
| kmen | **I(Y; sekce \| písař)** | 0,055 | 0,025 (45 %) | **0,043 (80 %, z 1,6)** | **0,060 (109 %, z −0,4)** |
| délka | I(Y; pozice) | 0,024 | 0,000 | 0,009 (39 %, z 14) | 0,000 |

* **Poziční informace.** ts2020 jako první null dává poziční informaci. Je to asi polovina reálné
  (prefix 46–47 %, kmen 57–58 %, sufix 57–59 %; ZL split i join), přičemž generátor má pravidla jen pro
  začátky řádků a odstavců a pro volbu zdroje „stejná pozice v řádku nad“. Zbytek (z 4–17) nereprodukuje.
  Šablona nemá žádný poziční mechanismus, takže poziční informace je u ní 0.
* **Sekční informace.** I(kmen; sekce | písař) ts2020 **reprodukuje**: 80 % reálu, z 1,6 (split) a 1,3 (join).
  Nepodmíněnou sekční informaci kmene (44 %) a sekční informaci po podmínění pozicí (30 %) nereprodukuje.
  U afixů dává ts2020 12–18 % sekční informace a 37–62 % sekční informace uvnitř písaře (sufix z 2,9,
  na hraně). Šablona reprodukuje všechny sekční členy kmene (z ≤ 2,2), protože je má z konstrukce (α = 5).
  Afixy nereprodukuje, protože jejich frekvence nemá podle sekce.
* **MI(d), ts2020 má na d = 1 hranu**: korigovaná MI při d = 1 (korekce `pc`, ZL split) je hrana glyfů
  0,0053 (reál 0,143, tj. 4 %), BPE 0,037 (19 %), afixová třída 0,036 (29 %) a kmen 0,013 (26 %).
  Při d = 2 a dál je ts2020 nulové nebo mírně záporné (−0,003 až 0) se stejným tvarem jako reál:
  ostrý vrchol při d = 1, pak nic. Šablona má MI nulovou při všech d. Plocha d = 2…10 s korekcí `line`
  je u reálu i obou nullů nulová (kmen: reál −0,002, ts2020 −0,005, tmpl 0,000).

## Krok 2 – „poziční kmen“ a parser (`parser_variants.csv`, `parser_affix_map.csv`)

Bod a 95% interval ze 100 subsamplů, split. Varianta (b) přidává Eva gallows `p f t k cph cfh cth ckh`
a sufixy `m g am om`. U v101 jsou to ekvivalenty ze zarovnání lokusů. Varianta (c) bere jako kmen
prostřední BPE jednotky. `bpe_mid3` = (c) jen pro tokeny s ≥ 3 jednotkami: v (c) má 72–75 % tokenů
prázdný kmen.

| sada | varianta | I(kmen; pozice) | I(kmen; sekce) | I(kmen; sekce \| písař) | pozice − sekce [int.] | jen místo v řádku (4 hodnoty) − sekce |
|---|---|---|---|---|---|---|
| ZL | (a) výchozí | 0,191 | 0,119 | 0,054 | 0,072 [0,050; 0,097] | −0,013 [−0,032; 0,008] |
| ZL | (b) + poziční glyfy | 0,170 | 0,118 | 0,054 | 0,052 [0,031; 0,076] | −0,021 [−0,038; −0,000] |
| ZL | (c) BPE střed | 0,101 | 0,035 | 0,015 | 0,065 [0,054; 0,077] | 0,018 [0,010; 0,026] |
| ZL | (c) ≥ 3 jednotky | 0,254 | 0,111 | 0,053 | 0,143 [0,112; 0,172] | 0,043 [0,017; 0,072] |
| GC | (a) výchozí | 0,177 | 0,112 | 0,052 | 0,065 [0,042; 0,085] | −0,006 [−0,024; 0,014] |
| GC | (b) + poziční glyfy | 0,148 | 0,109 | 0,051 | 0,040 [0,018; 0,061] | −0,021 [−0,038; −0,003] |
| GC | (c) BPE střed | 0,101 | 0,040 | 0,019 | 0,061 [0,051; 0,072] | 0,012 [0,005; 0,020] |
| GC | (c) ≥ 3 jednotky | 0,248 | 0,124 | 0,059 | 0,124 [0,094; 0,155] | 0,019 [−0,005; 0,043] |

Poziční převaha kmene je **stabilní**. I(kmen; pozice) > I(kmen; sekce) i > I(kmen; sekce | písař)
vychází ve všech variantách, v ZL i GC, s intervalem rozdílu nad nulou. Totéž platí po podmínění
(I(pozice | sekce) − I(sekce | pozice) = 0,042–0,130, vše > 0). Odebrání pozičních glyfů do afixů ubere
z poziční informace kmene jen 11 % (ZL) a 16 % (GC). Sekční informace se přitom nezmění.
Je tu ale výhrada ke kapacitě proměnné. Samotné místo v řádku (first / mid / last / single) nese
u afixových parserů zhruba stejně informace jako 7 sekcí (rozdíl −0,021 až −0,006). Převaha
„pozice > sekce“ u afixových parserů tedy plyne hlavně z jemnějších pozičních tříd
(rel_pos × řádek v odstavci). U BPE kmenů vede pozice i při 4 hodnotách.

## Krok 3 – 006 s ts2020 a šablonou (`pages_ts2020.csv`)

Kosinus, K = 30, model spec+size, 30 seedů. U nullu průměr, z a p = podíl replik ≤ reál:

| | log_dist | same_section | kontrast po velikosti |
|---|---|---|---|
| reál ZL split | −0,0099 [−0,0139; −0,0057] | +0,0032 [−0,0098; +0,0194] | −0,060 [−0,092; −0,020] |
| ts2020 ZL split | −0,0028 (z −2,9; p 0,065) | +0,0039 (z −0,1) | −0,010 (z −4,5) |
| ts2020 ZL join (reál −0,0091) | −0,0021 (z −3,2; p 0,065) | +0,0065 (z −1,6) | −0,011 (z −3,6) |
| tmpl ZL split | −0,0009 (z −7,5) | +0,0060 (z −0,4) | −0,001 (z −4,7) |
| tmpl GC split (reál −0,0093 / +0,0060) | −0,0009 (z −7,8) | +0,0059 (z 0,0) | +0,003 (z −4,4) |
| autocopy 86 ZL split (006) | −0,0027 (z −9,5) | +0,0020 | −0,014 |

* **Mezistránková křivka ts2020 není plochá.** Stejná sekce klesá z 0,128 (vzdálenost 1) na 0,062 (32+),
  reál z 0,188 na 0,097 (ZL split, `table == curve`). Pokles nese řetězení přes začátky stránek
  (viz výše). Koeficient vzdálenosti je u ts2020 jen 23–28 % reálu, stejně jako u autocopy, ale ts2020
  má mezi seedy 3× větší rozptyl (SD 0,0022–0,0024). Proto z = −2,9 (split) a −3,2 (join) leží na prahu
  a reálu dosáhne 1 z 30 replik (p = (1 + 1) / 31 = 0,065). Kontrast „stejná sekce daleko − jiná sekce blízko“ ts2020
  nereprodukuje (z −4,5 / −3,6).
* **Šablona** má pokles se vzdáleností téměř nulový (λ = 0,3 je lokalita stránky, ne sousedství).
  Sekční složku má **stejně velkou jako reál** (ZL +0,0060 vs +0,0032, GC +0,0059 vs +0,0060). Obojí
  leží v šumu: interval reálu obsahuje 0 a šablona ji má z konstrukce (α = 5). Jazykovou složku šablona
  nemá (+0,0008 až +0,0034 vs reál +0,023 až +0,030).

## Krok 4 – posun sufixů A/B podmíněný pozicí (`ab_by_position.csv`)

`std` = přímá standardizace: rozdělení sufixů A i B převážené na společný mix 22 pozičních tříd 007,
pak jedna JSD. `cond` = vážený průměr JSD uvnitř tříd. U písaře 3 je `cond` vychýlený malými buňkami
(A3 má 733 tokenů, pokrytí tříd 39–47 %), proto je hlavním číslem `std`. Přebytek se počítá nad
permutací jazyka mezi stránkami. U A vs B celkem je to permutace přes všechny stránky, protože permutace
uvnitř písaře mění jen písaře 3 a je z konstrukce blízká reálu. U písaře 3 jde o permutaci uvnitř písaře
(spec).

| sada | A vs B: vše | A vs B: std (poměr přebytku) | A vs B: vnitřní tokeny | písař 3: vše (přebytek [int.]) | písař 3: std (přebytek, poměr, p) | písař 3: vnitřní std (poměr) |
|---|---|---|---|---|---|---|
| ZL split | 0,088 | 0,088 (1,00) | 0,102 | 0,068 (0,049 [0,023; 0,074]) | 0,065 (0,047; 0,96; p 0,005) | 0,071 (0,99) |
| ZL join | 0,093 | 0,092 (0,99) | 0,111 | 0,069 (0,048 [0,022; 0,074]) | 0,066 (0,046; 0,95; p 0,005) | 0,072 (0,97) |
| GC split | 0,065 | 0,064 (0,99) | 0,075 | 0,063 (0,045 [0,018; 0,072]) | 0,065 (0,047; 1,04; p 0,005) | 0,066 (0,99) |
| GC join | 0,075 | 0,074 (1,00) | 0,086 | 0,062 (0,042 [0,014; 0,071]) | 0,061 (0,041; 0,97; p 0,005) | 0,063 (0,95) |

Posun A/B **po podmínění pozicí zůstává celý** (poměr 0,95–1,04 ≥ ½), celkem i uvnitř písaře 3, v ZL i GC.
U vnitřních tokenů řádku je dokonce větší než celkem (0,102 vs 0,088), protože poslední slova řádku
se mezi A a B liší méně (`last` 0,043 vs `first` 0,097, ZL split). Interval přebytku `std` u písaře 3
z jackknife subsamplů je široký. Dolní mez je kladná jen v GC split (ZL split [−0,004; 0,097]), protože
subsamply často vyřadí část z mála stránek A3. Permutační test uvnitř písaře je ale ve všech čtyřech
sadách na minimu p = 0,005.

## Krok 5 – tabulka pro syntézu (`synthesis_update.csv`)

Řádek × metrika × sada × null, s průměrem, SD, z, `reproduces` a verdiktem. Verdikt je „ano“ jen tehdy,
když platí ve všech sadách, kde null existuje. Oddíl 1 syntézy má 6 řádků. Spec chce 9, takže jsem
doplnil tři nálezy 006/007, na které se ptá kritérium 2 (7–9; `notes.md`).

| # | vlastnost (ZL split) | reál | ts2020 (průměr, z) | reprodukuje ts2020 | tmpl ZL / GC (z) | reprodukuje tmpl |
|---|---|---|---|---|---|---|
| 1 | okrajová MI | 0,185 | 0,007 (40; join 23) | ne | 0,0002 (286 / 255) | ne |
| 2 | poziční JSD první glyf | 0,186 | 0,044 (22; join 18) | ne | 0,0015 (387 / 274) | ne |
| 3 | BPE žlab @64 | 1,058 | 0,692 (5,7; join 5,1) | ne | 0,756 (7,6 / 14,0) | ne |
| 4 | sufixů na kmen | 5,25 | 3,79 (10,2; join 10,7) | ne | 4,44 (7,0 / 3,9) | ne |
| 5 | hapax | 0,676 | 0,579 (8,3; join 12,5) | ne | 0,610 (10,6 / 12,4) | ne |
| 6 | jsd_suffix_AB | 0,088 | 0,004 (33; join 35) | ne | 0,0006 (269 / 165) | ne |
| 7 | I(kmen; sekce \| písař) | 0,055 | 0,043 (1,6; join 1,3) | **ano** | 0,060 (−0,4 / −0,6) | **ano** |
| 7 | I(kmen; sekce \| pozice) | 0,117 | 0,035 (8,6) | ne | 0,079 (1,8 / 2,1) | ano |
| 8 | log_dist sdílení | −0,0099 | −0,0028 (−2,9; join −3,2) | částečně (jen split) | −0,0009 (−7,5 / −7,8) | ne |
| 8 | kontrast po velikosti | −0,060 | −0,010 (−4,5; join −3,6) | ne | −0,001 (−4,7 / −4,4) | ne |
| 9 | plocha d = 2…10, korekce `line` | −0,002 | −0,005 (0,3; join −0,6) | **ano** | 0,000 (−0,3 / 0,0) | **ano** |
| 9 | MI d = 3, korekce `pc` | 0,0002 | −0,0014 (0,6) | ano | 0,000 (0,1 / 0,1) | ano |
| 9 | plocha d = 2…10, korekce `pc` | 0,031 | 0,001 (2,6; join 2,0) | ano | −0,000 (4,5 / 3,0) | ne |
| 9 | MI d = 1, kmen | 0,051 | 0,013 (10,1) | ne | 0,001 (27 / 22) | ne |

Řádky 1–6 jsou u ts2020 převzaté z 50 seedů fáze B 004 (čísla shodná se SYNTHESIS). U šablony jde
o řádky 4–6 a hapax z 50 seedů ablace `best` 005 a řádky 1–3 z 30 nových vzorků (`profile/`).

## Odpovědi (kritérium 2)

| otázka | ZL | GC |
|---|---|---|
| (a) přežije ts2020 sekční informace kmenů nezávislá na písaři? | **ne**: ts2020 dává I(kmen; sekce \| písař) 0,043 = 80 % reálu (z 1,6 / 1,3). Přežívá jen nepodmíněná a pozičně podmíněná sekční informace (ts2020 44 % a 30 %, z 6–9). Šablona reprodukuje všechny tři (z konstrukce). | ts2020 chybí. Šablona reprodukuje (z −0,6) → jako indicie smyslu **ne** |
| (b) přežije strmý pokles sdílení? | **ano, těsně**: ts2020 má průměr jen 23–28 % reálu, ale z −2,9 (split) / −3,2 (join), p 0,065. Kontrast po velikosti přežije jasně (z −4,5 / −3,6). Šablona pokles nedává (z −7,5) | ts2020 chybí. Šablona nedává (z −7,8 / −9,0) → **ano** |
| (c) „žádná MI za d = 2“ jako rozdíl od nullu? | **ne, u nullu je to stejné**: ts2020 i šablona mají při d ≥ 2 také nulu (plocha s korekcí `line` z 0,3 / −0,6; d = 3 z 0,6 / 1,2). Rozdíl je jen ve výšce vrcholu při d = 1 (kmen 0,051 vs 0,013, z 10). | totéž pro šablonu (z −0,05 / 0,6); ts2020 chybí |
| (d) přežije A/B uvnitř písaře podmínění pozicí? | **ano**: přebytek po standardizaci 0,047 / 0,046 = 95–96 % nepodmíněného, p 0,005 | **ano**: 0,047 / 0,041 (97–104 %), p 0,005 |
| (e) je poziční převaha kmene artefakt parseru? | **ne**: pozice − sekce > 0 ve všech 4 variantách (0,052–0,143), poziční glyfy ubírají 11 %. Výhrada: jen 4 místa v řádku nesou u afixových parserů tolik jako sekce | **ne**: 0,040–0,124, poziční glyfy ubírají 16 %; stejná výhrada |

## Interpretace

Kanonický kopírovací null mění obraz 006/007 ve dvou bodech. Za prvé vyrábí asi polovinu poziční informace
o tvaru tokenu a 80 % sekční informace kmenů uvnitř písaře. Ta tedy přestává být samostatnou indicií obsahu:
lokální kopírování v souvislých blocích téhož písaře stačí. Za druhé tvar MI(d) „vrchol při d = 1, pak nic“
mají všechny nully. To není rozdíl, je to vlastnost, kterou reál sdílí s generátorem (reál má jen vyšší vrchol).
Nepodmíněná sekční informace kmenů, strmý pokles sdílení se vzdáleností a posun sufixů A/B uvnitř písaře,
který nevysvětlí ani pozice v řádku, zůstávají nad ts2020 i nad šablonou. Poziční charakter kmene parser
nevyrábí.

## Co výsledek neukazuje

* ts2020 chybí pro GC (generátor zná jen Eva), takže verdikty „ts2020“ stojí jen na ZL. Kritérium
  „obě sady“ pro ts2020 znamená split + join.
* Šablona obsahuje téma a lokalitu z konstrukce (α, λ naladěné v 005 na `jsd_stem_sections` a
  `stem_dispersion`). Její reprodukce sekčních členů je důkazem existence parametrů, ne nezávislou predikcí.
* Verdikt (b) leží na prahu |z| = 3. Rozptyl ts2020 mezi seedy je velký, průměr je přitom jen čtvrtina reálu.
  „Reprodukuje“ u ts2020 tu znamená „reál leží na okraji rozptylu běhů“.
* U písaře 3 jde o 733 tokenů A na málo stránkách. Jackknife interval standardizovaného přebytku je široký
  (ve 3 ze 4 sad zasahuje pod 0), rozhoduje permutační test.
* Varianta (c) parseru má 72–75 % prázdných kmenů. Varianta `bpe_mid3` ji doplňuje na podmnožině tokenů,
  takže není přímo srovnatelná s (a)/(b).
* 30 seedů nullu na 006/007, jedna konfigurace každého generátoru. Jiné konfigurace ts2020 (např.
  p_same_position 28/42) by poziční a sekční podíl posunuly.

## Další krok

1. Do syntézy: řádek „sekční informace kmenů nezávislá na písaři“ přesunout do oddílu 2 (reprodukuje ts2020),
   „žádná MI za d ≥ 2“ vést jako vlastnost sdílenou s nully, „A/B uvnitř písaře“ ponechat v oddílu 1
   s poznámkou „není poziční“.
2. Rozložit poziční informaci reálu na podíl, který dává ts2020 (začátky řádků/odstavců, stejná pozice nad),
   a zbytek: například rozklad 007 jen na vnitřních tokenech řádků mimo řádek 0 odstavce.
3. Pokles sdílení se vzdáleností: ts2020 s delší pamětí přes stránky (zdroj z předchozích N stránek)
   jako cílený null, protože současný port má jen řetězení přes začátek stránky.
