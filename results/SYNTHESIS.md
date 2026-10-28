# SYNTHESIS – co která vlastnost textu vysvětluje a co ne

Průběžná tabulka; čísla z `results/NNN_*/`. **Stav po 001–013** (28. 10. 2026).

Nully a kontroly: **ts2020** (port Timm & Schinner, 004, config 21; Eva-only → jen ZL); **tmpl** (šablona s paradigmaty,
tématem, lokalitou, 005); **ac-tri** (003); **pos-shuf**, **line-shuf**, **perm** (permutace stránek/štítků);
**próza** = latinský narativ, lékařská a botanická latina, anglický herbář, italština (009); **katalog** = Linnaeus
po položkách (009); **blábol** = 38 pisatelů Gaskell & Bowern, 10 k tokenů, srovnáváno s Voynichem oříznutým (013);
**Naibbe** = šifra latiny (Greshko), referenční bod, ne hypotéza (013).

## 1. Vlastnosti, které má Voynich a nemá je žádný skutečný text, blábol, šifra ani generátor

| Vlastnost | Reál ZL split | nejbližší | Poznámka |
|---|---|---|---|
| **Poziční vazby afixů** I(prefix; pozice) | 0,134 bit | ts2020 0,063 (pravidla začátků řádků); próza ≤ 0,003; katalog 0,042; blábol 0,006; Naibbe 0 | jediná vlastnost, kterou nemá ani blábol, ani katalog, ani šifra; vysvětlí ale jen 2–5 % variace tvaru při daném kmeni (011) |
| **Šev slov**: okrajová MI (poslední→první glyf) a I(sufix_t; prefix_{t+1} \| kmen, pozice, …) | 0,185 bit; 0,038 bit | ts2020 0,007 / 0,004; próza 0,02–0,07 / 0,006–0,010; blábol 0,013 / 0,011; Naibbe 0,005 / 0,001 | jen d = 1, končí na hranici řádku; R&T: půl dělá y→q |
| **Konec řádku a uzavřený řádek**: `-m/-g` na konci 17,5×; vazba přes hranici řádku / uvnitř 0,02 | – | ts2020 7×; próza 1–2×, poměr 0,6–1,3; blábol 2×, 0,79; katalog po položkách 0,05–0,13 | řádek jako jednotka (Currier); položkový rys bez položkového obsahu (010) |
| **A/B jako dva slovníky s ostrou hranicí** (012): Jaccard kmenů A×B o 0,17–0,19 pod polovinami téhož jazyka i v herbáři; posun sufixů difúzní (víc tvarů A → B `-dy`), párování 1:1 vysvětlí < ½; skok na hranici 4,4× | – | ts2020/tmpl 0; blábol: pisatelé se liší 5× víc (ale každý jinak) | přežívá písaře (005), pozici (008); „jazyk vs. pisatelský okruh“ nerozlišeno (písař 3: 2 stránky A) |

## 2. Vlastnosti, které Voynich sdílí s lidským blábolem (a ne s jazykem ani generátory)

| Vlastnost | Reál | blábol | próza | ts2020 / tmpl | Poznámka |
|---|---|---|---|---|---|
| Šířka a vyrovnanost paradigmat: sufixů na kmen; H(sufix \| kmen) | 5,25 (oříznuto 4,9); 1,82 | 3,96; **1,89** | ≤ 2,4; ≤ 0,88 | 3,8; 1,65 / 4,4; 2,8 | Naibbe také (4,5; 2,0) – vzniká i z více zápisů téhož písmene |
| >90 % variace tvaru při daném kmeni nevysvětlí pozice, registr, písař, sekce (011) | 92 % zbytek | zbytek ~65 % (pisatel vysvětlí 30–35 %) | latina: zbytek ~98 % (variace je flexe, kovariáty ji neurčují) | ts2020 98 % | zbytek má jen slabou strukturu: šev (výše) a malá konzistence uvnitř stránky (menší než u kopírování) |
| Hapax | 0,676 | 0,85 | 0,47–0,67 | 0,58 / 0,61 | Naibbe opačně (0,41) |
| Slabé pořadí tokenů, MI kmenů na d = 1 | 0,9 %; 0,051 | ano (kmeny) | ≥ 2,2 %; ≥ 0,091 | 0,44 % / 0 | Naibbe také (0,012; 0,029) |

## 3. Vlastnosti vysvětlené nully (vyřazeny jako indicie)

`cmi_excess`, `mi_prefix_suffix`, produktivita, `jsd_stem_sections`, disperze kmenů (ts2020, tmpl) · sekční
informace kmenů v písaři (ts2020 80 %) · tvar MI(d) (všechny nully) · polovina poziční informace (ts2020) ·
vertikální podobnost (pos-shuf) · záhlaví a šablony řádků (pozice, próza) · sekční složka sdílení nad vzdálenost,
stopa popisků v textu (nic) · Zipf, délky, typy, H2, blízcí sousedé (kalibrace).

## 4. Obraz po třinácti experimentech

* Voynich se **nechová jako próza** (pořadí tokenů pod všemi texty, žádná posloupnostní vazba za sousedem, žádné
  téma na dálku), **ani jako katalog** (žádná záhlaví a šablony řádků), a **ani jako kopírovací generátor** (ten
  nedá šev, konec řádku, A/B, polovinu pozice).
* **Vnitroslovní variace je jako u lidského blábolu:** široké, vyrovnané sady tvarů při daném kmeni, mnoho hapaxů,
  volba tvaru nevysvětlená ničím měřitelným. Blábol ale **nemá** to, co Voynich odlišuje od všeho: vazby tvaru na
  místo v řádku, pravidla na švu slov, konvence konce řádku, uzavřený řádek, a sdílený (ne idiosynkratický) registr
  A/B. Tyto vlastnosti nemá ani šifra skutečného textu typu Naibbe.
* Jinak řečeno: rukopis vypadá jako **systém s pevnými pravidly vázanými na místo (řádek, pozice, šev) a volnou,
  blábolu podobnou obměnou tvarů uvnitř slov**, se dvěma ostře oddělenými slovníky A a B. Nic z toho zatím není
  stopa obsahu; nic z toho ale také obsah nevylučuje – jen říká, že pokud tam je, nenese ho pořadí slov ani
  sekční rozložení kmenů v míře, kterou by dala próza.
* Zbývají tři otázky, které lze ještě zodpovědět daty: (a) jsou poziční a švová pravidla **malá sada konvencí**
  (několik glyfových pravidel, které jdou přidat generátoru), nebo difúzní? (b) opakují se **vzácné kmeny** na
  vzdálených místech nad náhodu (jediný typ signálu, který by kopírování ani blábol nedaly a obsah ano)?
  (c) jsou kmeny jen v B **glyfové varianty** kmenů A (změna zápisu), nebo nová slova?

## 5. Otevřené body

014 inventář místních pravidel (rozklad pozice a švu na glyfy/páry; kolik pravidel dá ts2020 reál) · 015 opakování
vzácných kmenů na dálku · 016 A/B kmeny: varianty nebo nová slova · blábol jednoho pisatele na dlouhém textu
(není korpus) · fyzické pořadí bifolií · průběžná zpráva (001–013).
