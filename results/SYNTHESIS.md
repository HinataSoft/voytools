# SYNTHESIS – co která vlastnost textu vysvětluje a co ne

Průběžná tabulka: pro každou měřenou vlastnost, který nulový model ji reprodukuje a který ne.
Aktualizuje se po každém experimentu; čísla jsou z `results/NNN_*/`. **Stav po 001–010** (28. 10. 2026).

Nully: **ts2020** = port generátoru Timm & Schinner (004, config 21; kanonický kopírovací null; Eva-only →
jen ZL, split i join); **tmpl** = šablona s paradigmaty, tématem a lokalitou z konstrukce (005); **ac-tri** =
starší autokopírování (003); **slot**, **perm**, **pos-shuf** (míchání uvnitř poziční třídy), **line-shuf**.
Kladné kontroly (009): latinská próza (Caesar), lékařská (Celsus) a botanická (Plinius) latina, anglický herbář
(Culpeper), italská próza, Linnéův katalog – vysázené do layoutu ZL.
„Reprodukuje“ = |z| < 3 vůči SD nullu; u ts2020 pozor na velký rozptyl mezi seedy (uváděn i průměr).

## 1. Vlastnosti mimo dosah všech nullů (a zda je mají skutečné texty)

| Vlastnost | Reál ZL split | ts2020 | tmpl | skutečné texty (009) | Poznámka |
|---|---|---|---|---|---|
| Síla pozičních vazeb afixů: I(prefix; pozice) | 0,134 bit | 0,063 (47 %, z 17) | 0 | próza ≈ 0, položky ≤ 0,042 | **mimo všechny texty i nully**; ts2020 dává polovinu z pravidel pro začátky řádků |
| Sekční vazby afixů I(prefix/sufix; sekce \| pozice) | 0,104 / 0,080 | 0,013 / 0,013 | ≈ 0 | ≤ 0,011 / ≤ 0,025 | mimo vše; sekce ≈ A/B ≈ písař – viz A/B níže |
| Posun sufixů A/B **uvnitř písaře 3, po podmínění pozicí** | 0,047 (95 % neподmíněného), p 0,005 | 0 z konstrukce | 0 | – | dva registry/„jazyky“, ne dvě ruce, ne pozice (005, 008) |
| Šířka paradigmatu: sufixů na kmen; H(sufix \| kmen) | 5,25; 1,82 bit | 3,8; 1,65 | 4,4; 2,8 | ≤ 2,4; ≤ 0,88 | **mimo všechny texty i nully**: v jazycích jeden tvar převažuje, zde jsou varianty kmene četné a vyrovnané; není velikost slovníku (004) |
| Hapax (podíl typů) | 0,676 | 0,58 | 0,61 | 0,47–0,67 (Plinius 0,669) | nad všemi generátory; GC 0,70 nad vším |
| Okrajová MI (d = 1, uvnitř řádku) | 0,185 | 0,007 | 0 | próza 0,02–0,07; katalog 0,30 | uvnitř rozpětí textů; žádný generátor; R&T: půl dělá y→q |
| Konec řádku: `-m/-g` na posledním místě | 17,5× obohaceno | 7× | 1× | 1,1–2,4× | ortografická konvence konce řádku; `jsd_last_mid` 0,107 mimo všechny texty |
| Pokles sdílení kmenů se vzdáleností stránek | −0,0099 | −0,0028 (z −2,9 / −3,2, p 0,065) | −0,0009 | próza má také pokles | **na hraně**; kontrast „stejná sekce daleko vs jiná blízko“ nereprodukuje žádný null |
| BPE žlab @64 | 1,06 | 0,69 (z 5,7) | 0,76 | – | – |
| **Pod** všemi texty: pořadí tokenů; MI kmenů na d = 1 | 0,9 %; 0,051 | 0,44 %; 0,013 | 0; 0 | ≥ 2,2 %; ≥ 0,091 | Voynich má méně posloupnostní vazby než každá próza i katalog – ale víc než každý generátor |

## 2. Vlastnosti, které některý null reprodukuje (vyřazeny jako samostatné indicie)

| Vlastnost | Reál | Reprodukuje | Co ji vysvětlilo |
|---|---|---|---|
| `cmi_excess`, `mi_prefix_suffix`, produktivita, `jsd_stem_sections`, disperze kmenů | – | ts2020; tmpl (3–5 pravidel) | náhrady podobných glyfů v kopírování; paradigma/téma/lokalita jako parametry (004, 005) |
| Sekční informace kmenů nezávislá na písaři I(kmen; sekce \| písař) | 0,055 | ts2020 0,043 (80 %, z 1,6); tmpl | kopírování v souvislých blocích téhož písaře (008) |
| Tvar MI(d): vrchol při d = 1, pak nic | – | všechny nully | sdílená vlastnost, ne rozdíl; liší se jen výška vrcholu (008) |
| Polovina poziční informace tvaru slova | – | ts2020 (pravidla začátků řádků, stejná pozice nad) | 008 |
| Vertikální podobnost (token nad, 1,15×) | 0,005 přebytek | ts2020 1,14×; **pos-shuf 1,16×** | čistě poziční; Timmův mechanismus „stejná pozice nad“ není potřeba (010) |
| Záhlaví řádku, šablony řádků | 0,012; +0,007 | pos-shuf, line-shuf; próza | žádná položková struktura nad poziční profil (010) |
| Sekční složka sdílení nad vzdálenost; stopa popisků v textu | ≈ 0 | – | není co reprodukovat (006) |
| Zipf, délky, typy, H2, blízcí sousedé, Jaccard inventářů | – | ts2020 / všechny | cíle kalibrace / strukturální |

## 3. Registrový profil (009) a řádek (010)

* Celkově je Voynich nejblíž **próze** (Celsus) a ze všech zdrojů **ts2020**; katalog je dál. Slovník, posloupnost
  a lokalita = próza ≈ generátor; paradigma = generátor; jediná „položková“ vlastnost = chování začátku a konce
  řádku. Pozice vzdaluje Voynich ode všeho.
* Řádek je uzavřená jednotka: mezitokenová vazba (hrana, třídy) **končí na konci řádku** (poměr ≈ 0 proti
  0,6–1,3 u lámané prózy), ale uvnitř řádku nejsou záhlaví, šablony ani vertikální paralely nad to, co dá pozice.

## 4. Co z toho plyne (stav po deseti experimentech)

* **Variabilita textu je převážně formální:** pozice v řádku, registr A/B (uvnitř písaře), konec řádku a
  vnitroslovní obměna (mnoho vyrovnaných variant kmene). Posloupnostní vazby jsou slabší než v každém
  skutečném textu, položková struktura chybí, téma na dálku chybí. Jako **próza** se text nechová a jako
  **katalog** také ne.
* **Co žádný model nevysvětlil** a co by mohlo nést informaci: (a) šířka a vyrovnanost paradigmat, (b) síla
  pozičních a registrových vazeb afixů, (c) A/B jako registr uvnitř písaře, (d) hapaxy, (e) okrajová vazba a
  konec řádku. Všechno jsou vlastnosti *tvaru slov*, ne pořadí slov. Pokud text nese význam, nese ho spíš
  volba tvaru (afixů) než syntax – nebo je afixová variace řízena pozicí a registrem a význam v ní není.
* **Rozhodující otázka pro další kolo:** je vnitroslovní variace (sufixy daného kmene) vysvětlitelná pozicí
  v řádku + registrem + písařem? Pokud ano, „paradigma“ je poziční jev a zbývá jen hapax a A/B. Pokud ne,
  je to první kandidát na obsahovou volbu.
* Chybí dva nully, na které se otázka smyslu nejvíc odvolává: **lidský blábol** (Gaskell & Bowern, data na
  GitHubu) a **šifra skutečného textu** (Naibbe v klonu R&T) jako referenční body na registrovém profilu.

## 5. Otevřené body

011 rozklad vnitroslovní variace (H(sufix | kmen) vs pozice, registr, písař; totéž pro latinu a ts2020) ·
012 A/B jako registr: systematické mapování sufixů, sdílené kmeny, pozvolnost podle pořadí psaní, písař 3 ·
013 lidský blábol (Gaskell) a Naibbe šifra na registrovém profilu 009/010 · rozklad okrajové MI na páry ·
ts2020 s delší pamětí přes stránky · fyzické pořadí bifolií (Layfield & Davis 2026) · dvojice herbář–farmacie.
