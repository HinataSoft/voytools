# 015 – Opakování vzácných kmenů na dálku

Větev `exp/015-rare-stem-recurrence`; kód `experiments/015_rare_stem_recurrence/run.py` + `summarize.py`,
funkce `voytools/recurrence.py` (testy `tests/test_recurrence.py`). Volby při výkladu spec jsou v
`experiments/015_rare_stem_recurrence/notes.md`. Čísla níže jsou z `distance_profile.csv`,
`far_pairs_sections.csv`, `form_by_distance.csv`, `far_clusters.csv`, `label_stems.csv` a `summary.csv`
(hlavní sada: dekompozice (a), kmeny s 2–5 výskyty, `split`; ostatní sady tamtéž).

## Co se měřilo

Kmen = `stem` z `decompose_table` (a: afixy naučené na P textu; b: + poziční glyfy). Každá **dvojice výskytů
téhož vzácného kmene** má vzdálenost v pořadí stránek souboru (0 = táž stránka), atributy stránek (sekce =
`illustration`, jazyk, písař, kvatern, bifolio) a tvar (identický token / jiný tvar). Null podle spec: tokeny
kmene přiděleny náhodným stránkám s pravděpodobností úměrnou velikosti stránky, frekvence zachována; jeho
očekávání jsou spočítána **přesně** (každá dvojice = dva nezávislé tahy stránky), 200 vzorkovaných permutací
(`perm_tokens`) dává interval a kontrolu. Interval reálu = bootstrap přes kmeny (200×, 2,5–97,5 %).

Velikost vzorku (reál, dekompozice a): ZL split 480 kmenů s 2–5 výskyty (1 398 dvojic), 158 s 6–20 (10 102
dvojic), 41 popiskových (158 dvojic); ZL join 572 / 174 / 49; GC split 632 / 241 / 52; GC join 686 / 242 / 55.
Spec čekal ~3 000 vzácných kmenů – to platí pro celé tokeny, ne pro kmeny s naučenými afixy; intervaly jsou
proto širší, než spec předpokládal.

## 1. Rozdělení vzdáleností (`distance_profile.csv`, `distance_profile.png`)

Poměr pozorovaný/očekávaný podíl dvojic po binech (dekompozice a, kmeny 2–5):

| zdroj | 0 | 1 | 2–3 | 4–7 | 8–15 | 16–31 | 32+ | přebytek ≥ 16 (`x_far`) | ≥ 16 mezi ≥ 4 (`x_far_nl`) |
|---|---|---|---|---|---|---|---|---|---|
| **ZL split, reál** | 2,89 | 1,97 | 1,60 | 1,38 | 1,24 | 1,02 | 0,87 | −0,088 [−0,119; −0,057] | −0,053 [−0,081; −0,024] |
| **GC split, reál** | 3,09 | 1,73 | 1,44 | 1,45 | 1,34 | 1,09 | 0,85 | −0,097 [−0,124; −0,070] | −0,067 [−0,091; −0,043] |
| ZL join, reál | 3,17 | 2,25 | 1,64 | 1,33 | 1,20 | 1,04 | 0,87 | −0,089 | −0,047 |
| GC join, reál | 2,70 | 1,66 | 1,24 | 1,38 | 1,36 | 1,13 | 0,86 | −0,083 | −0,061 |
| permutace tokenů (ZL) | 1,01 | 1,03 | 0,99 | 1,01 | 1,00 | 1,00 | 1,00 | −0,001 [−0,022; 0,022] | −0,001 |
| permutace stránek globálně (ZL) | 2,89 | 0,98 | 0,98 | 1,00 | 0,98 | 0,97 | 0,99 | −0,014 | 0,000 |
| permutace stránek v sekci (ZL) | 2,89 | 1,33 | 1,37 | 1,29 | 1,19 | 1,05 | 0,91 | −0,062 | −0,038 |
| ts2020 cfg 21 (ZL) | 35,2 | 1,56 | 0,79 | 0,72 | 0,75 | 0,75 | 0,68 | −0,254 [−0,295; −0,204] | −0,008 [−0,031; 0,015] |
| šablona 005 s tématem (ZL) | 25,9 | 0,78 | 0,78 | 0,78 | 0,76 | 0,83 | 0,79 | −0,168 [−0,223; −0,131] | +0,004 [−0,026; 0,034] |
| šablona 005 bez tématu (ZL) | 25,3 | 0,72 | 0,68 | 0,74 | 0,75 | 0,80 | 0,81 | −0,156 | +0,010 |
| ac-tri 86 (ZL) | 52,0 | 4,38 | 0,58 | 0,52 | 0,50 | 0,50 | 0,51 | −0,400 | −0,001 |
| Plinius, próza (layout ZL) | 7,76 | 3,90 | 2,05 | 1,57 | 1,19 | 0,86 | 0,76 | −0,179 [−0,194; −0,165] | −0,083 [−0,096; −0,069] |
| Celsus, próza | 7,33 | 4,11 | 2,45 | 1,77 | 1,28 | 1,05 | 0,69 | −0,209 | −0,110 |
| Linnaeus po položkách | 37,3 | 4,67 | 1,38 | 0,96 | 0,70 | 0,66 | 0,56 | −0,343 | −0,046 |
| blábol (013, 60 stránek) | 13,1 | 4,54 | 0,83 | 0,49 | 0,45 | 0,38 | 0,30 | −0,329 | −0,070 |
| ZL ořez na 60 stránek (013) | 2,45 | 1,28 | 1,26 | 1,14 | 0,86 | 0,99 | 0,67 | −0,066 [−0,104; −0,011] | −0,026 [−0,068; 0,032] |

Kmeny 6–20 (kontrola) dávají stejný tvar s mírnějším vrcholem (ZL split 2,53 / 1,66 / 1,39 / 1,41 / 1,32 /
1,12 / 0,87; `x_far` −0,082 [−0,114; −0,055]; GC −0,067). Dekompozice (b): ZL `x_far` −0,073, GC −0,091,
tvar stejný.

Výklad: **dálkový přebytek podle spec je záporný ve všech sadách** – vzácné kmeny Voynichu se opakují blíž,
než dá náhoda, nikoli dál. Záporný je ale i u všech kontrol i nullů; rozlišuje **tvar**: reál má mírný vrchol
na téže stránce (2,9–3,2×) a **pozvolný pokles přes 1–15 stránek** (2,0 → 1,2), tj. lokalitu v měřítku
sousedních listů a kvaternu. Latinská próza má stejný tvar se strmějším začátkem (7–8× na stránce, 4× na
sousední). Kopírovací generátory a šablona mají **jen vrchol na stránce** (25–52×) a mimo ni plochý profil
(poměr ≈ 0,5–0,8 ve všech binech, `x_far_nl` ≈ 0): neumějí gradient 1–15 stránek, který reál sdílí s prózou.
Blábol a Linnaeus po položkách mají pokles jen do 1 stránky (položka/řádek), dál nic.

## 2. Dálkové dvojice (≥ 16 stránek) a sekce / jazyk / písař (`far_pairs_sections.csv`)

Podíl dálkových dvojic se stejným atributem vs očekávání párované podle vzdálenosti (dekompozice a, 2–5):

| zdroj | n dvojic | sekce: podíl [boot.] | očekáv. | přebytek | jazyk: přebytek | písař: přebytek |
|---|---|---|---|---|---|---|
| **ZL split, reál** | 1 010 | 0,236 [0,198; 0,277] | 0,155 | **+0,081** | +0,066 [+0,029; +0,107] | +0,068 [+0,025; +0,108] |
| **GC split, reál** | 1 249 | 0,207 [0,178; 0,235] | 0,149 | **+0,058** | +0,096 [+0,056; +0,133] | +0,063 [+0,036; +0,094] |
| ZL join / GC join, reál | 1 224 / 1 395 | 0,237 / 0,211 | 0,152 / 0,151 | +0,085 / +0,060 | +0,083 / +0,088 | +0,073 / +0,050 |
| permutace tokenů (ZL) | 1 132 | 0,139 [0,118; 0,161] | 0,140 | −0,001 | 0,00 | 0,00 |
| permutace stránek globálně (ZL) | 1 182 | 0,361 | 0,363 | −0,002 | – | – |
| permutace stránek v sekci (ZL) | 1 053 | 0,255 [0,238; 0,272] | 0,151 | +0,104 | – | – |
| ts2020 (ZL) | 1 259 | 0,158 [0,121; 0,199] | 0,143 | +0,014 | – | – |
| šablona s tématem (ZL / GC) | 818 / 1 170 | 0,165 / 0,166 | 0,143 / 0,138 | +0,022 / +0,028 | – | – |
| šablona bez tématu (ZL / GC) | 801 / 1 167 | 0,161 / 0,155 | 0,139 / 0,137 | +0,022 / +0,018 | – | – |
| ac-tri 86 (ZL) | 1 934 | 0,139 | 0,139 | −0,001 | – | – |
| Celsus (próza, sekce = kniha) | 3 332 | 0,176 [0,162; 0,193] | 0,108 | +0,069 | – | – |
| Plinius (próza, sekce = kniha) | 4 890 | 0,069 [0,061; 0,078] | 0,046 | +0,022 | – | – |
| Linnaeus po položkách | 1 835 | 0,162 [0,145; 0,181] | 0,154 | +0,009 | – | – |

Kmeny 6–20: sekce +0,044 [+0,021; +0,075] (ZL), +0,064 [+0,040; +0,099] (GC); dekompozice (b): +0,082 (ZL),
+0,055 (GC). Globální permutace stránek má po opravě (velikost stránky cestuje s obsahem) přebytek 0, tj.
null je v pořádku; permutace uvnitř sekce dává přebytek +0,10, o něco víc než reál.

Střední vzdálenosti (2–15 stránek): stejné **bifolio** reál 0,103 vs 0,077 (ZL; interval přebytku −0,009 až
+0,054) a 0,119 vs 0,072 (GC; +0,019 až +0,079), stejný kvatern 0,785 vs 0,737 (ZL, n.s.), 0,766 vs 0,736 (GC,
n.s.) – náznak vazby na fyzický list jen v GC a GC join (+0,010 až +0,065), v ZL ne.

Výklad: **dálkové dvojice vzácných kmenů sedí do téže sekce, téhož jazyka i téhož písaře nad náhodu** (ZL i GC,
split i join, obě dekompozice, obě frekvenční pásma). Šablona s tématem α = 5 dá z konstrukce jen +0,02–0,03,
ts2020 +0,01, ac-tri 0; reál je 2–4× nad šablonou a na úrovni Celsova *De medicina* (knihy = témata) a nad
Pliniem. Sekce, jazyk A/B a písař jsou zde ale zaměnitelné (všechny tři přebytky ~ +0,06–0,10) – test
nerozliší „obsah sekce“ od „slovník jazyka/písaře“ (006, 007: sekční informace kmenů ≈ jazyk).

## 3. Tvar opakování podle vzdálenosti (`form_by_distance.csv`, `form_by_distance.png`)

Podíl identických tokenů mezi dvojicemi téhož kmene minus očekávání z paradigmatu kmene (nezávislá volba
tvaru):

| zdroj | 0 | 1–3 | 4–15 | ≥ 16 |
|---|---|---|---|---|
| **ZL split, reál** (n = 34 / 105 / 249 / 1 010) | −0,016 [−0,124; +0,104] | +0,007 [−0,089; +0,096] | +0,013 [−0,043; +0,079] | −0,003 [−0,036; +0,036] |
| **GC split, reál** (47 / 122 / 345 / 1 249) | +0,026 [−0,081; +0,176] | +0,023 [−0,049; +0,100] | +0,006 [−0,042; +0,052] | −0,005 [−0,030; +0,026] |
| ts2020 (ZL) | +0,082 | −0,017 | −0,036 | −0,035 |
| ac-tri 86 (ZL) | +0,086 | +0,034 | −0,081 | −0,084 |
| šablona (ZL, s / bez tématu) | +0,008 / +0,005 | −0,003 / +0,002 | −0,005 / −0,003 | −0,002 / −0,002 |
| Celsus | +0,048 | +0,054 | +0,008 | −0,020 [−0,047; +0,004] |
| Plinius | +0,043 | +0,022 | +0,015 | −0,013 [−0,034; +0,009] |
| Linnaeus po položkách | +0,022 | +0,022 | −0,021 | −0,015 |
| blábol | +0,073 | +0,037 | −0,104 | −0,124 |

Výklad: **tvar opakování ve Voynichu nezávisí na vzdálenosti** – ani na dálku jiný tvar častěji, ani nablízko
přesná kopie častěji; všechny intervaly obsahují nulu (ZL, GC, join, 6–20 i dekompozice b totéž, viz CSV).
Kopírovací generátory a blábol mají na téže stránce přesné kopie nad paradigma (+0,08) a na dálku pod ním
(−0,04 až −0,12); latinská próza slabě totéž (+0,04–0,05 nablízko, −0,01 až −0,02 na dálku, jak čekáme u
skloňovaných jmen: blízké opakování téhož pádu). Voynich je v tomto bodě **„nezávislejší“ než latina i než
kopírování** – tvar je vybírán jako by nezávisle na místě, což odpovídá nálezu 002/011 (sufix rozptýlený,
varianty kmene nekoncentrované).

## 4. Shluky na dálku (`far_clusters.csv`)

Kmeny s ≥ 3 výskyty rozložené přes ≥ 3 kvaterny (resp. ≥ 3 sekce), pozorováno / očekáváno (poměr):

| zdroj | kvaterny 2–5 | sekce 2–5 | kvaterny 6–20 | sekce 6–20 |
|---|---|---|---|---|
| **ZL split, reál** | 140 / 168,0 (0,83; boot. 121–160) | 78 / 119,5 (0,65) | 148 / 156,8 (0,94) | 137 / 151,7 (0,90) |
| **GC split, reál** | 150 / 195,9 (0,77) | 94 / 142,2 (0,66) | 220 / 238,5 (0,92) | 197 / 229,6 (0,86) |
| permutace tokenů | 1,00 | 1,00 | 1,00 | 1,00 |
| permutace stránek globálně (ZL) | 0,97 | 0,96 | – | – |
| ts2020 (ZL) | 0,53 | 0,48 | 0,87 | 0,80 |
| šablona s / bez tématu (ZL) | 0,71 / 0,75 | 0,65 / 0,69 | 0,83 / 0,84 | 0,79 / 0,81 |
| ac-tri 86 (ZL) | 0,34 | 0,33 | – | – |
| Celsus / Plinius / Linnaeus | 0,71 / 0,72 / 0,50 | 0,48 / 0,62 / 0,42 | 0,85 / 0,91 / 0,67 | 0,74 / 0,87 / 0,56 |
| blábol | 0,28 | 0,28 | 0,54 | 0,62 |

Výklad: shluků přes ≥ 3 kvaterny je ve Voynichu **méně než náhodně** (0,77–0,83), ale víc než u ts2020 (0,53)
a ac-tri (0,34) a stejně jako u šablony (0,71–0,75) a latinské prózy (0,71–0,72). Metrika tedy jen opakuje
lokalitu z bodu 1 a nad náhodu nic nepřidává.

## 5. Popiskové kmeny (`label_stems.csv`, pás `l`)

Kmeny z popisků (Lz/Lf/Lc/Ls/Ln; ZL split 328 kmenů) s 2–5 výskyty v P textu: ZL split 41 (158 dvojic), ZL
join 49, GC split 52, GC join 55. Profil (ZL split, a): 4,52 / 2,49 / 1,80 / 0,48 / 0,82 / 1,11 / 0,94;
`x_far` −0,032 [−0,162; +0,055], `x_far_nl` **+0,036 [−0,030; +0,094]**; GC split −0,075 [−0,168; +0,023] /
−0,052; ZL join −0,050 / −0,004; GC join −0,074 / −0,047. Dálkové dvojice v téže sekci: ZL 0,187 vs 0,153
(interval přebytku −0,034 až +0,114), GC 0,222 vs 0,155 (−0,015 až +0,188). Shluky přes ≥ 3 kvaterny: ZL 20 /
20,9 (0,96), sekce 15 / 15,0 (1,00); GC 17 / 20,7 (0,82). Nejrozptýlenější popiskové kmeny ZL (5 stránek,
4–5 kvaternů): `olsh`, `kald`, `aiino`, `tald`, `teol` (tabulka v `label_stems.csv`).

Výklad: popiskové kmeny mají **plošší profil než ostatní vzácné kmeny** (menší lokální vrchol mimo stránku,
`x_far_nl` kolem nuly nebo kladný, shluky přes kvaterny na úrovni náhody, zatímco ostatní kmeny 0,8); to je
směr, jaký by měla „jména“, ale s 41–55 kmeny jsou intervaly tak široké, že se od ostatních kmenů (ZL
`x_far_nl` −0,053) staticky neoddělí. Sekční vazba dálkových dvojic je u nich stejná (+0,03 až +0,07, n.s.).

## Kontroly nullu a kritéria přijetí

* Permutace tokenů: průměrný přebytek −0,001 až +0,001 ve všech sadách (`distance_profile.csv`,
  `perm_tokens`), shoda vzorkované permutace s přesným očekáváním i u sekcí a shluků (poměr 1,00); test
  `tests/test_recurrence.py::test_permutation_matches_exact_expectation`. Žádný přebytek > 0,3 (max |x_far|
  0,40 je deficit ac-tri, ne přebytek).
* Kritérium 1 (profily pro ZL/GC × split/join × a/b × 2–5/6–20, nully, kontroly): splněno (4 829 běhů,
  `meta.json`); ts2020 jen pro ZL (generuje Eva); popiskový pás jen u reálu a generátorů v layoutu ZL/GC.
* Kritérium 2: (i) dálkový přebytek je **< 0** v ZL i GC (−0,09), menší deficit než próza (−0,18 až −0,21) i
  než ts2020 (−0,25) a šablona (−0,17); (ii) **ano**, dálkové dvojice sedí do téže sekce nad náhodu (+0,06 až
  +0,09 v ZL i GC, 3–4× víc než šablona s tématem, na úrovni Celsa), stejně i jazyk a písař; (iii) **ne**,
  tvar opakování nezávisí na vzdálenosti; (iv) popiskové kmeny mají plošší profil, nerozlišitelně od ostatních.

## Shrnutí

**Čísla:** přebytek dálkových dvojic (≥ 16 stránek) ZL −0,088 [−0,119; −0,057], GC −0,097 [−0,124; −0,070];
poměr pozorované/očekávané po binech ZL 2,9 / 2,0 / 1,6 / 1,4 / 1,2 / 1,0 / 0,87 (ts2020 35 / 1,6 / 0,8 / 0,7
/ 0,75 / 0,75 / 0,68; Plinius 7,8 / 3,9 / 2,1 / 1,6 / 1,2 / 0,86 / 0,76); dálkové dvojice v téže sekci 0,236
vs 0,155 (ZL), 0,207 vs 0,149 (GC), šablona s tématem +0,02, Celsus +0,07; identický tvar na dálku −0,003
[−0,036; +0,036] proti paradigmatu; shluky přes ≥ 3 kvaterny 0,83 (ZL) / 0,77 (GC) očekávání.

**Interpretace:** Vzácné kmeny Voynichu se vracejí s pozvolným poklesem přes 1–15 stránek, který mají i
latinské prózy a který kopírovací generátory ani šablona nedávají (ty mají jen vrchol na stránce a plochý
zbytek); dálková složka nad náhodu neexistuje, ale dálkové dvojice jsou vázané na sekci/jazyk/písaře víc, než
dá šablona s tématem, a stejně jako u tematicky členěné latinské prózy. Tvar opakování je nezávislý na
vzdálenosti – ani kopie nablízko, ani „jiný pád“ na dálku.

**Co výsledek neukazuje:** nerozliší sekci od jazyka a písaře (všechny tři přebytky stejné), neukazuje
„obsahové“ opakování na dálku (žádný bin ≥ 16 není nad náhodu), a popiskové kmeny jsou příliš málo početné
(41–55) na samostatný závěr; profil 1–15 stránek je slučitelný s řazením podle tématu i s postupnou změnou
slovníku psaní (drift), což tento experiment nerozliší.

**Další krok:** oddělit drift od tématu: profil vzdáleností uvnitř jedné sekce a jednoho písaře (herbál A)
proti permutaci stránek uvnitř sekce (zde: přebytek v sekci +0,10 vs reál +0,08, tj. pořadí stránek v sekci
nese ~ nic navíc), a párovat stránky podle ilustrace (stejná rostlina herbál/farma) místo sekce – to je test,
který šablona s tématem z konstrukce nemůže projít.
