# SYNTHESIS – co která vlastnost textu vysvětluje a co ne

Průběžná tabulka: pro každou měřenou vlastnost, který nulový model ji reprodukuje a který ne.
Aktualizuje se po každém experimentu; čísla jsou z `results/NNN_*/`. **Stav po 001–004** (27. 10. 2026).

Nully: **shuf** = míchání tokenů v řádcích; **ac-bi** = autokopírování s bigramovými editacemi (001,
config 573); **ac-tri** = s trigramovými editacemi (003, config 86); **slot** = slotová nezávislost
prefix/kmen/sufix (002); **ts2020** = port generátoru Timm & Schinner (004, config 21 „calB“ –
**kanonický kopírovací null od 004**); **tmpl** = strukturovaná šablona s paradigmaty (005, zatím –).
„Reprodukuje“ = |z| < 3 vůči SD generátoru přes 50 seedů, v obou sadách (ZL i GC, u ts2020 split i join).
Pozor: ts2020 má velký rozptyl mezi seedy, takže |z| < 3 někdy znamená jen „reál leží v rozptylu běhů“,
ne „průměr sedí“ – u takových řádků je uveden i průměr.

## 1. Vlastnosti, které zatím **žádný** null nereprodukuje

| Vlastnost | Reál (ZL split) | shuf | ac-bi | ac-tri | slot | ts2020 | Poznámka |
|---|---|---|---|---|---|---|---|
| Okrajová MI (poslední→první glyf) | 0,185 bit | 0 | 0 | 0 | – | 0,007 (z 40) | žádný null nemá mezislovní mechanismus; R&T: polovinu dělá pár y→q, celý token nepředpoví lépe než poslední glyf → spíš povrchový/ortografický jev; nutný rozklad (008) |
| Poziční JSD (první glyf řádku vs. střed) | 0,19 | – | 0,002 | 0,001 | – | 0,045 (z 22) | ts2020 má pravidla jen pro začátky řádků; zbytek = rozložení v řádku (Feaster) – spíš písařská/layoutová konvence než obsah |
| BPE žlab @64 | 1,06 | 0,92 | 0,60 | 0,76–0,88 | – | 0,69 (z 5,7) | i shuf dává 0,92 → část závisí na pořadí tokenů; ts2020 odpovídá R&T (0,65–0,86) |
| Sufixů na kmen (≥5 tok.) | 5,25 | = reál | 2,3 | 2,7 | 7,8 | 3,8 (z 10) | reál **mezi** kopírováním a nezávislou šablonou; nejde o velikost slovníku (004 krok 3) – **nejsilnější zbývající paradigmatický signál** |
| Hapax (podíl typů) | 0,676 | = reál | 0,53 | 0,56 | – | 0,58 (z 8) | žádný generátor nedá dost hapaxů; R&T totéž (59–60 %) |
| `jsd_suffix_AB` (frekvence sufixů A vs B) | 0,088 | = reál | 0,002 | – | 0,0003 | 0,004 (z 33) | **zaměněno s písařem** – kontrola `hand` v 005; generátor A/B nerozlišuje z konstrukce |

## 2. Vlastnosti, které ts2020 reprodukuje (vyřazeny jako samostatné indicie smyslu)

| Vlastnost | Reál | ts2020 (průměr, z) | Co ji vysvětlilo |
|---|---|---|---|
| `cmi_excess` I(P;S\|kmen) | 0,088 | 0,115 (−1,4) | náhrady podobných glyfů uvnitř slova; ac-bi/ac-tri dávaly dokonce 2–3× víc |
| `mi_prefix_suffix` (nepodmíněná) | 0,205 | 0,187 (0,7) | dtto – v 002 vypadala jako nejsilnější paradigmatický signál |
| sklon produktivity sufixů | 0,31 | 0,31 (0,2) | |
| `jsd_stem_sections` (kmeny podle sekce) | 0,23 | 0,24 (−0,4) | lokální kopírování + sekce jako souvislé bloky stačí |
| `stem_dispersion` | 0,46 | 0,33 (2,7) | **jen v rozptylu seedů**; průměr pořád pod reálem; ac-bi 0,005 |
| pořadí tokenů (excess/H) | 0,9 % | 0,44 % (2,5) | dtto, průměr poloviční |
| `H_suffix_given_stem` | 1,82 | 1,65 (1,9) | |
| podíl kmenů s >1 sufixem | 0,91 | 0,82 (3,0 / 3,5) | na hraně prahu (nekonzistentní split/join) |
| Zipf, délky, počet typů, podíl blízkých sousedů, H2 glyfů | – | trefeno kalibrací | cíle kalibrace (loss 0,05) |
| Jaccard inventáře sufixů mezi sekcemi | 0,96–0,98 | všechny nully | strukturální (uzavřený inventář) |

## 3. Co nully nedokážou trefit kalibrací (cíle)

| Cíl | Reál | ts2020 calB | ac-tri | Poznámka |
|---|---|---|---|---|
| H2 glyfů | 2,69 | 2,62 | 3,33 | **vyřešeno** v 004 (náhrady podobných glyfů, ne libovolné editace) |
| Hapax | 0,676 | 0,58 | 0,56 | otevřené; snad jiné parametry / delší šablona |
| loss (5 cílů) | – | 0,052 | 0,269 | kritérium 1 poprvé splněno |

## 4. Co z toho zatím plyne (opatrně)

* Po 004 se obraz změnil: většina „paradigmatických“ nálezů z 002 (vazba prefix–sufix, produktivita,
  tematické rozdělení kmenů) vzniká i kopírováním s náhradami podobných glyfů a pravidly pro začátky
  řádků. **Nejsou to indicie smyslu.**
* Zbývající nevysvětlené vlastnosti se dělí na dvě skupiny. (a) Takové, pro něž generátor nemá
  mechanismus a přidat ho je levné: okrajová MI (jeden parametr „začátek slova podle konce předchozího“),
  poziční profil (Feasterův drift), A/B rozdíl (dva startovní řádky). Tyhle samy o sobě smysl neprokážou;
  stanou se indicií, až se ukáže, že mají **strukturu závislou na obsahu** (na kmeni/třídě slova, ne jen
  na posledním glyfu). (b) Takové, které generátor měl vytvořit a nevytvořil: **šířka paradigmat**
  (5,3 vs 3,8 sufixu na kmen při stejné velikosti slovníku), hloubka BPE žlabu, hapaxy. To jsou
  kandidáti pro další práci.
* Varování k ts2020: tabulky podobnosti glyfů a pravidla „curve/line“ jsou ručně odvozené z rukopisu,
  generátor tedy zná tvary slov Voynichu; je to silný null, ne slabý.

## 5. Otevřené body

005 strukturovaná šablona (kolik ingrediencí potřebuje umělý systém na řádek 1), robustnost afixů,
A/B vs písař · 006 téma na dálku, popisky · 007 pozice vs obsah, MI(d) · **nové:** 008 struktura
paradigmat – jsou sady sufixů kmenů sdílené (flexivní třídy), nebo idiosynkratické (historie kopírování)?
Rozklad okrajové MI na páry a její závislost na třídě slova/sekci/A-B · později: A/B a poziční
mechanismus v ts2020, lidský blábol (Gaskell), dobové jazyky z Boxerova repa jako kladné kontroly.
