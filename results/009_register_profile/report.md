# 009 – Registr: kde leží Voynich mezi narativem, katalogem a generátorem

Kód: `experiments/009_register_profile/run.py` (profil), `summarize.py` (tabulky), načítání a vysázení kontrol
`voytools/controls.py` (testy `tests/test_controls.py`). Deník rozhodnutí: `experiments/009_register_profile/notes.md`.
Soubory: `controls_used.csv`, `profile.csv` (metrika × zdroj: hodnota, interval, SD), `profile_z.csv`,
`distances.csv` + `ranking.csv`, `lists.csv`, `checks.csv`, `profile_heatmap.png`, `runs.csv` (2 090 řádků), `meta.json`.

## Co bylo spočteno

* **Zdroje.** Voynich: ZL split (hlavní), ZL join, GC split – bod + 100 subsamplů 80 % stránek.
  Kontroly (9, všechny vysázené do layoutu ZL split, 34 268 tokenů, 4 128 řádků, 207 stran) – bod + 100 subsamplů.
  Nully (30 seedů): ts2020 calB (ZL split, ZL join; Eva-only, pro GC chybí), šablona 005 plný model (ZL, GC),
  ac-tri 86 (ZL, GC), míchání v řádcích (ZL, GC).
* **Interval** reálu a kontrol: bod + 2 × kvantilové odchylky subsamplů (2,5–97,5 %, jako 007); SD = 2 × SD subsamplů.
  Nully: průměr a 2,5–97,5 % přes seedy.
* **Profil: 37 metrik v 7 skupinách** (`meta.json`, `GROUPS`): slovník (4), pořadí tokenů (2), pozice vs sekce
  (12 = 4 cíle × I(·;pozice), I(·;sekce), I(·;sekce|pozice)), paradigma (7), lokalita 006 (3), posloupnost
  MI(d) 007 (4), řádek (5).

### Kontroly (kritérium 1)

| kontrola | žánr | zdroj (soubor v klonu) | vysázení | sekce |
|---|---|---|---|---|
| `caesar_prose` | narativ, la | `…/latin_controls/caesar_la.txt` + `caesar_la2.txt` (Gutenberg; `.tess` kniha I–III jen 17,8 tis. slov) | próza, souvisle od začátku | 4 (knihy, sloučené) |
| `celsus_prose` | lékařská la | `data/voynich-units/latin/celsus.de_medicina.part.1–5.tess` | próza | 5 knih |
| `pliny_prose` | botanická la | `data/voynich-units/herbal/pliny_the_elder…part.4–6.tess` (knihy 16–30) | próza | 7 knih |
| `plinyP_prose` | dtto, 2. edice | Boxer `corpora/latin/legacy/pliny/pliny_lat0.csv` (Perseus) | próza | stejných 7 knih |
| `culpeper_prose` | herbář, en | `data/voynich-units/herbal/culpeper_en.txt` | próza | 7 (písmena abecedy) |
| `culpeper_items` | herbář, en | dtto | **položky = věty/klauze** (odstavců je méně než řádků) | 6 |
| `italian_prose` | próza, it (1895) | `data/voynich-units/italian/27608.txt` | próza | 7 bloků (oddíly nerozpoznatelné) |
| `linnaeus_prose` | katalog, la | `…/latin_controls/species1.txt` + `species2.txt` (*Species Plantarum*) | próza | 5 (řády, sloučené) |
| `linnaeus_items` | katalog, la | dtto | **položky = odstavce hesla** (jméno rodu, diagnóza, synonymum, Habitat) | 4 |

Nepoužito: Boxer `ellis_1854` (receptář, jen 10 565 slov), `italian/26166.txt` (21 tis. slov). Culpeper „po heslech“
nešel (≈ 2 000 odstavců na 4 128 řádků) – proto po větách. **Kritérium 1 splněno** (7 různých textů, dvě položkové varianty).
Kontroly podezřelých hodnot (`checks.csv`): žádná kontrola nemá `hapax_share` > 0,9 (max 0,669, Plinius) a žádná
mechanicky lámaná próza nemá I(·; pozice) vyšší než Voynich (max 0,003 bitu u kmene, Voynich 0,19).

## Výsledky

### Vzdálenost profilů (`ranking.csv`, RMS standardizovaných rozdílů, standardizace SD přes zdroje)

ZL split, všech 37 metrik: **line-shuffle 1,10 < ts2020 1,38 < Celsus 1,66 < šablona 1,77 < italština 1,78 <
Caesar 1,85 < Linnaeus položky 1,85 < ac86 1,86 < Plinius P 1,89 < Plinius 1,89 < Linnaeus próza 2,00 < Culpeper
próza 2,13 < Culpeper položky 2,23**. GC split: line 1,14 < Celsus < šablona < italština < ac86 < Caesar < Plinius ×2 <
Linnaeus položky < Linnaeus próza < Culpeper ×2 (Celsus 1,65 nejblíž ze skutečných textů). ZL join: ts2020 nejblíž,
pak Celsus.

Míchání v řádcích je nejblíž jen proto, že je to tentýž text (slovník, paradigma, lokalita a sekční informace
jsou z konstrukce totožné); smysluplné srovnání je se skutečnými texty a s generátory.

Nejmenší vzdálenost podle třídy zdroje (ZL split / GC split):

| metriky | próza | položky | katalog-próza | generátor | nejblíž |
|---|---|---|---|---|---|
| všechny | 1,66 / 1,65 | 1,85 / 1,80 | 2,00 / 1,91 | **1,38** / 1,68 | ts2020 (ZL), Celsus (GC, bez shuffle) |
| jen slovník | **0,48** / 0,64 | 1,17 / 1,36 | 1,19 / 1,39 | 0,59 / **0,61** | próza ≈ generátor |
| jen pořadí | **0,86** / **1,05** | 1,67 / 1,64 | 1,96 / 1,81 | 1,19 / 1,44 | Celsus |
| jen pozice | **1,97** / 1,92 | 2,12 / **1,89** | 2,24 / 1,96 | 2,04 / 2,19 | nic blízko (vše ≈ 2) |
| jen paradigma | 1,41 / 1,29 | 1,74 / 1,67 | 1,72 / 1,65 | **0,63** / **0,74** | generátor (ts2020 / šablona) |
| jen lokalita | 1,12 / **0,89** | 1,33 / 1,18 | **1,07** / **0,89** | 1,29 / 1,08 | nerozhoduje |
| jen posloupnost | **0,22** / **0,21** | 1,80 / 1,83 | 2,71 / 2,71 | 0,32 / 0,37 | próza ≈ generátor |
| jen řádek | 2,16 / 2,34 | 1,94 / 2,17 | 1,99 / **2,16** | **1,19** / 2,41 | ts2020 (ZL); položky před prózou |
| bez pozice | 1,38 / 1,40 | 1,70 / 1,76 | 1,86 / 1,88 | **0,92** / **1,37** | ts2020 |

Které skupiny co dělají: **posloupnost a pořadí** táhnou Voynich k próze (a stejně tak ke generátorům – obě
skupiny jsou u katalogu daleko, protože se tam opakují formule). **Paradigma** táhne ke generátorům. **Řádek**
je jediná skupina, kde je Voynich blíž položkovému textu než próze (ZL 1,94 vs 2,16; GC 2,17 vs 2,34) – a ještě
blíž ts2020 (ZL 1,19), které má pravidla pro začátky řádků. **Pozice** vzdaluje Voynich ode všeho (vyřazení
pozice sníží vzdálenost k próze z 1,66 na 1,38 a ke ts2020 z 1,38 na 0,92). Pořadí nejbližších skutečných textů
se mezi variantami „bez skupiny“ téměř nemění: vždy Celsus, italština, Caesar před Linnéem a Culpeperem.

### Seznamy (`lists.csv`; verdikt platí, je-li stejný v ZL split i GC split)

Pravidlo: (a) bod Voynichu leží v rozpětí intervalů skutečných textů; (b) mimo ně a mimo intervaly všech
generátorů (ts2020, šablona, ac86); (c) mimo skutečné texty, ale uvnitř intervalu některého generátoru.
Míchání v řádcích o (b)/(c) nerozhoduje (je to tentýž text).

**(b) Voynich-specifické – mimo všechny skutečné texty i všechny generátory (ZL i GC shodně 10 metrik; 3 smíšené):**

| metrika | ZL split [interval] | skutečné texty (min–max bodů) | generátory (ZL) |
|---|---|---|---|
| `prefix__I_pos` | 0,134 [0,124; 0,142] | 0,000–0,042 | ts2020 0,063 |
| `prefix__I_sec` / `I_sec_g_pos` | 0,085 / 0,104 | ≤ 0,010 / ≤ 0,011 | ts2020 0,012 / 0,011 |
| `suffix__I_sec` / `I_sec_g_pos` | 0,070 / 0,080 | ≤ 0,026 / ≤ 0,025 | ts2020 0,011 |
| `mean_suffixes_per_stem` | 5,25 [5,00; 5,46] | 1,71–2,36 | ts2020 3,78, šablona 4,47 |
| `H_suffix_given_stem` | 1,82 | 0,47–0,88 | ts2020 1,61, šablona 2,81 (nad) |
| `token_mi_excess_share` | 0,009 [0,006; 0,012] – **pod** | 0,022–0,171 | ts2020 0,005, ostatní 0 |
| `mi_stem_d1` | 0,051 [0,041; 0,062] – **pod** | 0,091–1,178 | ts2020 0,014, ostatní ≈ 0 |
| `jsd_last_mid` | 0,107 [0,095; 0,117] | 0,006–0,011 | ts2020 0,048 |
| `first_top5_excess`, `last_top5_excess` | −0,027 / 0,060 | GC (b), ZL (c: ts2020) – smíšené | |

Plus `frac_stems_multi_suffix` je (c) (šablona 0,908 ~ reál 0,912) a `hapax_share` je smíšené (ZL 0,676 leží
těsně v intervalu Plinia 0,669, GC 0,702 je nad všemi skutečnými texty i generátory).

**(c) mimo skutečné texty, uvnitř generátoru:** `frac_stems_multi_suffix` (šablona). Jinak žádná metrika – generátory
nevysvětlují žádnou vlastnost, kterou by skutečné texty neměly, kromě šíře paradigmatu v šabloně.

**(a) uvnitř rozpětí skutečných textů (23 metrik):** slovník (`n_types`, `top50_coverage`, Zipf), hrana
`edge_mi_excess` (0,186 – katalog má 0,30, latinská próza 0,02–0,07), všechny metriky kmene a délky v bloku
pozice (I(kmen; pozice) 0,19 leží mezi prózou ≈ 0 a položkovým Linnéem 0,25), `suffix__I_pos`, `mi_prefix_suffix`,
`cmi_excess`, `stem_dispersion`, `jsd_stem_sections`, lokalita (všechny 3), `mi_affix_d1`, plochy MI(d),
`jsd_first_mid`, `I_first_rest`. U `suffix/stem/length__I_pos`, `jsd_first_mid` a `I_first_rest` je Voynich
v rozpětí **jen díky položkovým kontrolám** (próza má ≈ 0; Linnaeus/Culpeper po položkách 0,16–0,26 u `jsd_first_mid`).

**Mimo všechny generátory (bez ohledu na skutečné texty), ZL i GC:** edge MI (z vůči ts2020 70), všechny poziční
a sekční informace prefixů/sufixů (z vůči ts2020 3,8 u I(sufix; pozice), jinak 13–25), `jsd_first_mid`, `jsd_last_mid`, `I_first_rest`,
`mi_stem_d1`, `mi_affix_d1`, `token_mi_excess_share`, `mean_suffixes_per_stem`, `H_suffix_given_stem`, hapax.
Uvnitř některého generátoru: `cmi_excess`, `mi_prefix_suffix`, `stem_dispersion`, `jsd_stem_sections`,
`loc_same_section`, plochy MI(d), Zipf (souhlasí s 004/SYNTHESIS).

## Odpovědi (kritérium 3)

1. **Je Voynich blíž položkovému textu než próze?** **Celkově ne.** Nejbližší skutečné texty jsou prózy (Celsus 1,66,
   italština, Caesar); Linnaeus po položkách je až za nimi (1,85; GC 1,80 vs 1,65). Položkovému textu se Voynich
   podobá jen v **řádkových konvencích**: první token řádku má jiné rozložení tříd než střed (`jsd_first_mid`
   0,19; próza ≤ 0,008, položky 0,16–0,26) a nese informaci o zbytku řádku (`I_first_rest` 0,057; próza ≤ 0,02,
   položky 0,11–0,19), a I(kmen; pozice) 0,19 (próza ≈ 0, Linnaeus položky 0,25). Posloupnost (MI kmenů na d = 1
   a za ním) je naopak prozaická až „chudší než próza“, ne katalogová (katalogy mají MI(d) 10–40× vyšší kvůli
   opakovaným formulím).
2. **Mimo všechny skutečné texty:** šíře paradigmatu (5,25 sufixu na kmen vs ≤ 2,4; H(sufix|kmen) 1,82 vs ≤ 0,88),
   síla pozičních a sekčních vazeb **afixů** (I(prefix; pozice) 0,134 vs ≤ 0,042; I(prefix; sekce | pozice) 0,10 vs
   ≤ 0,011), rozdíl posledního tokenu řádku od středu (`jsd_last_mid` 0,107 vs ≤ 0,011) – a dvě metriky **pod**
   všemi texty: MI sousedních tokenů (0,9 % H vs ≥ 2,2 %) a MI kmenů na d = 1 (0,051 vs ≥ 0,091).
3. **Mimo všechny nully:** hrana, všechna poziční a sekční informace afixů, řádkové konvence (první i poslední
   token, záhlaví), MI kmenů i afixů na d = 1, šíře paradigmatu a hapaxy. Všechny tyto vlastnosti ts2020 má
   ve stejném směru, ale menší (pozice prefixu 0,063 vs 0,134; `jsd_last_mid` 0,048 vs 0,107).

## Interpretace

Registrovým profilem stojí Voynich **mezi prózou a generátorem, ne u katalogu**: slovník, posloupnost a
lokalita odpovídají souvislému textu (nebo stejně dobře ts2020), paradigma odpovídá generátorům, a jedinou
„položkovou“ vlastností je chování začátku řádku. Co ho odděluje od všech skutečných textů i generátorů,
je **síla pozičních a sekčních vazeb afixů** a **šíře paradigmatu** – obojí řádově mimo latinu, angličtinu i
italštinu, ale ve směru, kterým se ubírá ts2020 (generátor s poziční pravidly je vždy blíž než prózy).

## Co výsledek neukazuje

* Kontroly jsou abecední písma; Voynich je tokenizován do glyfů Eva/v101 a afixy jsou učené stejným algoritmem
  na jiných jednotkách – část rozdílu v afixových metrikách (pozice, sekce, paradigma) může být rozdíl
  „afix glyfového písma vs koncovka slova“, ne registru. Hlavní závěry ale drží v ZL i GC.
* Sekce Voynichu jsou zaměněné s písařem a jazykem A/B (007); vysoká I(afix; sekce) u Voynichu je
  pravděpodobně z velké části rozdíl A/B, který kontroly z konstrukce nemají (mají jen kovariáty `hand`/`language`
  převzaté z layoutu, nezávislé na textu).
* Položkové vysázení je umělé (utínání položek na délku řádku ZL); skutečný katalog má řádky jiné délky.
  Položková Culpeperova varianta je po větách, ne po heslech.
* Vzdálenosti závisí na výběru metrik a jejich standardizaci (SD přes zdroje); pořadí skutečných textů je
  stabilní, ale absolutní velikosti nejsou testem.
* Žádná metrika tu neodlišuje „smysl“ od „nesmyslu“; jde o popis registru.

## Další krok

* 010 (řádek jako položka) rozvede řádkové metriky: hranice řádku, záhlaví, šablony, vertikální kopírování.
* Rozložit I(afix; sekce) na část vysvětlenou A/B a písařem (podmínit `hand`, jako 007), a to i pro ts2020 s
  dvěma startovními řádky (A/B).
* Kontrola, která má glyfové jednotky blíž Voynichu (např. Naibbe šifra z klonu R&T, nebo text v zkratkovém
  písmu), aby se oddělilo „afix písma“ od „koncovky jazyka“.
