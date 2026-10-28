# 010 – Řádek jako položka: nezávislost řádků, šablony řádků, vertikální kopírování

Kód: `experiments/010_line_as_item/run.py`, `summarize.py`; nové funkce `voytools.stats.line_boundary_mi` a
`voytools.stats.vertical_similarity` (testy `tests/test_line_stats.py`); kontroly z 009 (`voytools/controls.py`).
Deník: `experiments/010_line_as_item/notes.md`. Tabulky (dlouhý formát: metrika × zdroj, hodnota, interval, SD):
`boundary_vs_inside.csv`, `line_header.csv`, `line_templates.csv` (+ `signatures_top20.csv`), `same_length.csv`,
`vertical_copy.csv`, `line_end.csv`, ke každé `*_z.csv` (z reálu vůči nullům stejné sady a kontrolám),
`summary.csv`, `runs.csv`, `meta.json`; grafy `boundary.png`, `vertical_copy.png`, `templates.png`.

**Zdroje.** Reál ZL/GC × split/join: bod + 100 subsamplů 80 % stránek (interval = bod + 2 × kvantilové odchylky,
jako 007/009). Nully po 30 seedech: ts2020 calB (ZL split, join), šablona 005 (ZL, GC split), ac-tri 86 (ZL, GC),
míchání v řádku (`line`), míchání v poziční třídě 007 (`posclass`), permutace celých řádků stejné délky
(`samelen`). Kontroly z 009 ve stejném layoutu: Caesar, Plinius, Celsus, italština, Culpeper (próza / po větách),
Linnaeus (próza / po položkách), bod + 100 subsamplů. Kritérium 1 splněno (všech 6 metrik pro všechny zdroje).
Varianta afixů (b) z 008 neprovedena – 008 nemá výsledky.

Hodnoty níže: ZL split [95% interval]; GC split v závorce, kde se liší.

## 1. Hranice řádku vs vnitřek (`boundary_vs_inside.csv`)

Korigovaná MI na d = 1 (bity; korekce permutací uvnitř každého druhu dvojic):

| jednotka | uvnitř řádku | přes hranici řádku | (b)/(a) | próza (b)/(a) | položky (b)/(a) |
|---|---|---|---|---|---|
| hrana (poslední→první glyf) | 0,195 [0,181; 0,207] | −0,002 [−0,016; 0,010] | −0,01 [−0,08; 0,05] (GC 0,04) | 0,62–1,00 | 0,05 (Culpeper), 0,12 (Linnaeus) |
| afixová třída | 0,078 | −0,002 | −0,03 [−0,22; 0,16] (GC 0,14) | 0,78–1,27 | 0,05–0,13 |
| kmen (top 200) | 0,127 | 0,028 | 0,22 [−0,06; 0,52] (GC −0,12) | 0,35–0,66 | 0,08–0,19 |

Přes hranici odstavce (546 dvojic) je hranová MI 0,025 [−0,050; 0,096] – neodlišitelná od nuly i od hranice řádku.
Generátory (ts2020, šablona, ac86) mají už **uvnitř** řádku hranovou MI ≈ 0 (ts2020 0,008), takže jejich poměr
nic neříká; přes hranici mají všechny ≈ 0. Míchání v řádku dává přes hranici afixů 0,015 [0,004; 0,026] – sousední
řádky v reálu sdílejí třídy (sekce, A/B, písař), přesto reál přes hranici nemá nic (−0,002). Míchání v poziční
třídě dává 0,002 [−0,006; 0,010], tj. stejně jako reál.

## 2. Záhlaví řádku (`line_header.csv`; null = permutace prvních tokenů mezi řádky stejné délky)

| I(třída 1. tokenu; …) | ZL split | p (bod) | posclass | ts2020 | `line` | próza | položky |
|---|---|---|---|---|---|---|---|
| každého dalšího tokenu (`head_rest`) | 0,012 [0,007; 0,016] | 0,005 | 0,010 | 0,007 | 0,018 | 0,002–0,006 (Linnaeus próza 0,023) | 0,045–0,058 |
| 2. tokenu | 0,041 [0,022; 0,058] | 0,005 | 0,022 | 0,026 | 0,019 | 0,027–0,103 (Linnaeus próza 0,25) | 0,45–0,63 |
| 3. tokenu | 0,020 [0,003; 0,038] | 0,005 | 0,013 | 0,007 | 0,021 | 0,008–0,031 | 0,21–0,35 |
| posledního | 0,003 [−0,008; 0,013] | 0,26 | 0,000 | 0,007 | 0,020 | −0,007–0,013 | 0,027–0,044 |

## 3. Šablony řádků (`line_templates.csv`, `signatures_top20.csv`)

* Zisk informace o třídě z (délka řádku, pozice): 0,088 [0,079; 0,102] bitu (GC 0,091) – posclass 0,084, ts2020 0,061,
  ac86/šablona 0,04, próza 0,037–0,040, položky 0,15–0,19. Zisk z předchozí třídy v řádku: 0,088 [0,079; 0,097] –
  posclass 0,036, ts2020 0,035, próza 0,034–0,130, položky 0,24–0,27.
* Podpisy řádků 4–8 tokenů: 932 různých na 1 000 řádků; 20 nejčastějších pokrývá 0,122 [0,108; 0,136] řádků,
  po „sloupcovém“ zamíchání (tokeny se stejnou délkou řádku a pozicí mezi řádky) 0,115 → přebytek 0,007 (GC 0,009);
  ts2020 0,105 vs 0,099 (0,006), míchání v řádku 0,120 vs 0,114, próza ~0 (0,059 vs 0,058), položky 0,088 vs 0,064
  a 0,126 vs 0,072 (přebytek 0,024–0,054). Nejčastější podpisy jsou kombinace nejčastějších tříd (`signatures_top20.csv`).
* zlib: poměr 0,477; zisk po sloupcovém zamíchání 0,012 [0,009; 0,014] (GC 0,009) – posclass 0,007, `line` 0,007,
  ts2020 0,004, šablona 0, próza 0,002–0,021 (Linnaeus próza 0,041), položky 0,044–0,050.

## 4. Řádky stejné délky (`same_length.csv`)

Přebytek Jaccardovy podobnosti tříd u stejné délky minus u různé délky (po odečtení globálního zamíchání):
0,001 [−0,016; 0,016] (GC −0,004); všechny zdroje v [−0,007; 0,012]. **Žádný efekt** – nikde.

## 5. Vertikální kopírování (`vertical_copy.csv`, `vertical_copy.png`)

| zdroj | shoda ≤ 1 editace s tokenem nad | očekávání (náhodná pozice nad) | přebytek | poměr nad / jinde nad |
|---|---|---|---|---|
| ZL split | 0,046 [0,042; 0,049] | 0,041 | 0,005 [0,003; 0,008], p = 0,005 | 1,15 [1,09; 1,21] |
| GC split / ZL join / GC join | 0,056 / 0,042 / 0,048 | | 0,007 / 0,006 / 0,006 | 1,15 / 1,20 / 1,15 |
| ts2020 ZL split | 0,053 | 0,049 | 0,004 [0,002; 0,007] | 1,14 [1,05; 1,23] |
| **posclass** ZL / GC | 0,045 / 0,054 | | **0,006 [0,005; 0,007] / 0,007** | **1,16 / 1,15** |
| samelen, ac86, šablona, `line` | | | 0,002 / 0 / 0 / 0 | 1,07 / 1,00 / 1,00 / 1,00 |
| próza | 0,009–0,031 | | 0,000–0,001 | 0,97–1,11 |
| položky (Linnaeus, Culpeper) | 0,037 / 0,054 | | 0,016 | 2,13 / 1,54 |

Šikmo nad (i ± 1) je reál také mírně nad očekáváním (0,044 / 0,045 vs 0,041).

## 6. Konce řádků (`line_end.csv`)

Podíl tokenů končících `-m`/`-g` (Eva): poslední v řádku 0,172 [0,153; 0,191], jinde 0,009 (obohacení 17,5×; GC
nejvíc obohacený koncový glyf 18,6×). ts2020 0,253 vs 0,035 (7×), šablona 0,010 vs 0,011, ac86 0,005 vs 0,005.
Poslední token **není kratší**: rozdíl délky −0,02 [−0,08; 0,04] glyfu (GC −0,01; ZL join +0,13 [0,06; 0,20] – delší).
Kontroly: obohacení koncového písmene 1,1–2,4×.

## Odpovědi (kritérium 2; ZL i GC shodně, není-li řečeno jinak)

1. **Je řádek uzavřenější jednotka než u prózy a než u generátorů?** **Ano vůči próze:** závislost sousedních
   tokenů (hrana, afixová třída) se na hranici řádku zastaví úplně (poměr ≈ 0 vs 0,6–1,3 u mechanicky lámané
   prózy), podobně jako u položkového vysázení (0,05–0,13). **Vůči generátorům nelze rozhodnout poměrem** –
   ty nemají závislost ani uvnitř řádku; reál se od nich liší tím, že silnou vnitřřádkovou hranu (0,195 vs ts2020
   0,008) má a za hranicí řádku ji ztrácí. Hodnota přes hranici (≈ 0) je stejná jako po zamíchání v poziční třídě,
   tj. tokeny na krajích řádku jsou dané pozicí, ne sousedem.
2. **Nese první token informaci o zbytku řádku?** **Ne nad rámec pozice.** `head_rest` 0,012 je nad permutací
   prvních tokenů (p = 0,005), ale stejnou hodnotu dává míchání v poziční třídě (0,010) a víc míchání v řádku (0,018 –
   homogenita řádku). Vazba na 2. token (0,041) je jen vazba sousedů a leží v rozpětí prózy (0,027–0,103).
   Katalog po položkách má 4–15× víc.
3. **Existují šablony řádků nad náhodu a nad ts2020?** **Ne v katalogovém smyslu.** Pokrytí 20 nejčastějšími
   podpisy přesahuje sloupcové zamíchání jen o 0,007 (ts2020 0,006, položky 0,024–0,054); zlib zisk 0,012 je nad
   ts2020 (0,004) a nad míchanými nully (0,007), ale v rozpětí prózy a 4× pod položkami. Řádková „šablona“ je
   poziční profil tříd (zisk z pozice 0,088 = posclass 0,084), ne opakované vzory.
4. **Je podobnost k tokenu nad zvýšená, a reprodukuje ji ts2020?** **Ano, mírně** (1,15×, přebytek 0,005, p = 0,005),
   a ts2020 ji reprodukuje (1,14×) – ale **stejný přebytek dává i míchání v poziční třídě** (1,16×), které žádné
   kopírování nemá. Zvýšená podobnost k tokenu nad je tedy vysvětlena tím, že tokeny na stejném indexu sousedních
   řádků sdílejí poziční třídu (relativní pozice, začátek/konec řádku). Mechanismus „kopie ze stejné pozice v řádku
   nad“ tím **není potvrzen jako zdroj** a není ani potřeba; reál nemá nic navíc oproti poziční nule.

Podezřelé hodnoty (kritérium 3): poměr (b)/(a) > 1,5 nemá žádný reál ani kontrola (max 1,27, Plinius afixy);
velké poměry u šablony a ac86 (až 140) jsou dělení téměř nulou (vnitřní MI ≈ 0), ne nález. Podíl „nad ≤ 1 editace“
je u všech zdrojů ≤ 0,056 (< 0,5); kompresní poměr ≥ 0,46 (> 0,2).

## Interpretace

Řádek Voynichu je uzavřený v jednom smyslu: **mezitokenová vazba (hrana, třídy) končí na konci řádku**, zatímco u
lámané prózy pokračuje – v tom se Voynich chová jako text sázený po položkách. Uvnitř řádku ale nejsou vidět
položkové vlastnosti (záhlaví, opakované šablony, vertikální paralely) nad to, co vysvětlí poziční třída tokenu;
vertikální podobnost, kterou má ts2020 z pravidla „stejná pozice“, má reál ve stejné míře, ale tutéž hodnotu dává
i zamíchání v poziční třídě.

## Co výsledek neukazuje

* Neříká, *proč* vazba končí na konci řádku: je to stejně slučitelné s položkovým zápisem jako s řádkovým
  generováním (každý řádek nový start), nebo s tím, že hranová MI je převážně poziční/ortografická (009, R&T).
* Položkové kontroly jsou uměle vysázené (utínání na délku řádku ZL), proto jejich „záhlaví“ a „šablony“ jsou
  horní mez katalogu, ne odhad skutečného katalogu.
* Hranice odstavce má jen ~550 dvojic – intervaly jsou široké, rozdíl hranice odstavce vs řádku neprůkazný.
* Afixové třídy jsou z naučených afixů (bez pozičních glyfů z 008); jiná definice tříd může posunout čísla.
* Nully `posclass` a `line` vycházejí z reálu (zachovávají slovník), generátory ne; ts2020 je Eva-only (GC nemá).

## Další krok

* Zopakovat vertikální test s podmíněním poziční třídou (očekávání = shoda s tokeny řádku nad ve *stejné poziční
  třídě*, ne s libovolnou pozicí) – rozhodne, zda po odečtení pozice zbývá vůbec něco pro kopírování.
* Hranovou MI uvnitř řádku rozložit podle pozice (první pár řádku vs vnitřní páry) – je vazba rovnoměrná, nebo
  ji nesou páry u začátku řádku?
* Až bude 008: zopakovat 1–3 s variantou afixů (b) (poziční glyfy mimo kmen).
