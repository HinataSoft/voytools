# SYNTHESIS – co která vlastnost textu vysvětluje a co ne

Průběžná tabulka: pro každou měřenou vlastnost, který nulový model ji reprodukuje (|z| < 3 v ZL i GC)
a který ne. Aktualizuje se po každém experimentu; čísla jsou z `results/NNN_*/`. Stav k experimentům 001–003
(27. 10. 2026). Zkratky nullů: **shuf** = míchání tokenů v řádcích; **ac-bi** = autokopírování s bigramovými
editacemi (001, config 573 / p_new 0,5); **ac-tri** = s trigramovými editacemi (003, config 86);
**slot** = slotová nezávislost prefix/kmen/sufix (002); **ts2020** = port Timm & Schinner (004, zatím –);
**tmpl** = strukturovaná šablona s paradigmaty (005, zatím –).

## 1. Vlastnosti, které zatím **žádný** null nereprodukuje (v ZL i GC)

| Vlastnost | Reál (ZL/GC, split) | shuf | ac-bi | ac-tri | slot | Poznámka |
|---|---|---|---|---|---|---|
| Okrajová MI (poslední→první glyf) | 0,185 / 0,226 bit | 0 | 0 | 0 | – | generátory nemají mezislovní mechanismus → z-skóre testuje konstrukci, ne text |
| Pořadí tokenů (excess / H) | 0,9 % / 1,3 % | 0 | 0 | 0 | – | slabé, ale nenulové; R&T: pod všemi prózami |
| Poziční JSD (první glyf řádku vs. střed) | 0,19 / 0,20 | – | 0,002 | 0,001 | – | generátory bez pozičního mechanismu |
| BPE žlab @64 | 1,06 / 1,13 | – | 0,60 | 0,76–0,88 | – | vázané editace přiblížily, ne do 3 SD; R&T ts2020: 0,65–0,86 |
| Sufixů na kmen (≥5 tok.) | 5,3 / 4,7 | = reál | 2,3 | 2,7 | 7,8 | reál **mezi** kopírováním a nezávislou šablonou |
| Podíl kmenů s >1 sufixem | 0,91 / 0,84 | = reál | 0,51 | 0,62 | 1,00 | dtto |
| `mi_prefix_suffix` (nepodmíněná) | 0,21 / 0,18 bit | = reál | 0,06–0,09 | – | 0,006 | nejsilnější „paradigmatický“ signál |
| `jsd_stem_sections` (kmeny podle sekce) | 0,23 / 0,27 | = reál | 0,49 | – | 0,13 | téma bez lokálního driftu – mezi nully |
| `stem_dispersion` (rozptyl variant kmene) | 0,46 / 0,48 | ≈ reál | 0,005–0,05 | – | 0,69 | varianty rozprostřené, mírně shluknuté (stránka/téma) |
| `jsd_suffix_AB` (frekvence sufixů A vs B) | 0,088 / 0,065 | = reál | 0,002 | – | 0,0003 | **zaměněno s písařem** – nutná kontrola `hand` (005) |

## 2. Vlastnosti, které nully reprodukují (nerozlišují)

| Vlastnost | Reál | Reprodukuje | Poznámka |
|---|---|---|---|
| Zipf, délky slov, počet typů | – | ac-bi, ac-tri, tmpl (očekáváno) | Gaskell & Bowern: i lidský blábol |
| Podíl blízkých sousedů (Timm) | 0,173 (poměr 1,18) | ac-tri 0,18 | cíl kalibrace |
| `cmi_excess` I(P;S\|kmen) nad náhodu | 0,09 | **ac-bi/ac-tri dávají 2–3× víc** (0,19–0,24) | samo o sobě paradigma neprokazuje; nad slot (0,035) ano |
| Jaccard inventáře sufixů mezi sekcemi | 0,96–0,98 | všechny nully | strukturální (uzavřený inventář 15 položek) |
| Sklon produktivity sufixů | 0,31–0,34 | slot (0,30) | reál se od šablony liší hladinou, ne sklonem |
| Gramatičnost tokenů (slot/trigram) | 0,99 / 1,0 | ac-tri z konstrukce | ac-bi jen 0,60–0,65 trigramově |

## 3. Vlastnosti, které zatím **nelze** trefit kalibrací generátoru

| Cíl | Reál | Nejlepší generátor | Proč |
|---|---|---|---|
| H2 glyfů (extended, scope all) | 2,69 | ac-tri 3,33; ac-bi 3,4–3,7 | editace skládají doložené trigramy v jiném poměru; R&T port dosahuje 2,6–2,7 → **004** |
| Hapax (podíl typů) | 0,676 | 0,56 | H2 a hapax jdou proti sobě (méně editací = méně typů) |
| H2 + hapax + share_near + n_types zároveň | – | loss 0,27 (ac-tri), 0,98 (ac-bi v mřížce) | kritérium 1 z 001/003 nesplněno |

## 4. Co z toho zatím plyne (opatrně)

* Kmeny se chovají jako lexikální jednotky s **omezeným** paradigmatem: méně variant než nezávislá
  šablona, víc než kopírování, se silnou nepodmíněnou vazbou prefix–sufix, s variantami rozptýlenými
  po celém textu a s rozdělením závislým na sekci. Platí ve dvou abecedách a třech politikách mezer.
* **Není rozhodnuto**, zda to je stopa významu, nebo jen *systému*: chybí null pro strukturovaný umělý
  jazyk (005) a kanonický kopírovací null (004). Část paradigmatických rozdílů může jít na vrub
  nekalibrované velikosti slovníku generátorů (2–3× víc kmenů s ≥ 5 výskyty).
* Okrajová MI, poziční profil a pořadí tokenů zůstávají mimo dosah všech generátorů, ale ty žádný
  poziční/mezislovní mechanismus nemají; interpretovatelné budou až proti ts2020 (má pravidla pro
  začátky řádků) a proti pozičnímu nullu (007).
* `cmi_excess` a Jaccard inventářů vyřazeny jako samostatné indicie.

## 5. Otevřené body k dalším experimentům

004 kanonický kopírovací null a kontrola velikosti slovníku · 005 strukturovaná šablona (kolik
„ručních“ pravidel potřebuje umělý systém, aby napodobil řádek 1) a robustnost afixů/parseru, A/B
vs písař · 006 téma na dálku, popisky · 007 pozice vs obsah, reziduální MI(d) · později: poziční
mechanismus v generátoru, lidský blábol (Gaskell data), dobové jazyky z Boxerova repa jako kladné
kontroly.
