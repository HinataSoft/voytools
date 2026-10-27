# 010 – Řádek jako položka: nezávislost řádků, šablony řádků, vertikální kopírování

## Otázka
007 ukázalo, že za bezprostředním sousedem není žádná posloupnostní závislost a že pozice v řádku
určuje tvar slova víc než cokoli jiného; Currier (1976) řádek označil za funkční jednotku. Otázka:
**chová se řádek jako samostatná položka (heslo, recept, popisek) s vnitřní šablonou?** A konkrétní
test kopírovací hypotézy: Timm & Schinner předpokládají, že zdrojem kopie je slovo na *stejné pozici
v řádku nad* – **je podobnost k tokenu přímo nad zvýšená nad podobnost k ostatním tokenům
předchozího řádku?**

## Hypotézy a predikce
- **Položka:** závislost (hrana, afixová třída, kmen) přes hranici řádku je výrazně slabší než uvnitř řádku
  na d = 1; první token řádku se chová jako „záhlaví“ (nese informaci o složení zbytku řádku); posloupnosti
  afixových tříd v řádku se opakují jako šablony (málo častých vzorů, nízká entropie na pozici);
  řádky stejné délky si jsou strukturně podobnější. Próza (kontroly 009, mechanicky lámané): žádný rozdíl
  uvnitř/přes hranici, žádné šablony. Katalog (Linnaeus po položkách): podobné Voynichu.
- **Vertikální kopírování (ts2020, Timm):** podobnost k tokenu nad (stejný index v předchozím řádku)
  > podobnost k ostatním tokenům předchozího řádku; reál: pokud ne, mechanismus „stejná pozice“
  je vyvrácen jako hlavní zdroj; pokud ano, je to silný argument pro kopírování.
- **Šablona 005:** žádná řádková struktura (nemá řádková pravidla) → slouží jako „bez položek“ null.

## Data a filtry
ZL3b a GC2a, `P`, `clean`, split (hlavní) a join; odstavce a řádky podle tabulky (`par_index`,
`line_in_par`, `token_index`). Dekompozice naučenými afixy (+ varianta (b) z 008 s pozičními glyfy,
pokud 008 ukáže, že záleží). Afixová třída = (prefix, sufix) → pro šablony řádků zredukovat na
~20 nejčastějších dvojic + `other`.

## Metriky
1. **Hranice vs vnitřek:** korigovaná MI na d = 1 pro hranu (poslední→první glyf), afixovou třídu a kmen
   (top 200) zvlášť pro dvojice (a) uvnitř řádku, (b) přes hranici řádku uvnitř odstavce, (c) přes hranici
   odstavce. Poměr (b)/(a) a (c)/(a) s intervaly (100 subsamplů stránek).
2. **Záhlaví řádku:** I(třída prvního tokenu; multiset tříd zbytku řádku) a I(třída prvního; třída druhého,
   třetího, posledního); null: permutace prvních tokenů mezi řádky stejné délky (200×).
3. **Šablony řádků:** řádek jako posloupnost afixových tříd; (i) entropie třídy na pozici 1…k (k = délka)
   a podmíněná entropie H(třída_i | třída_{i−1}) vs H(třída_i | pozice); (ii) počet různých „podpisů“
   řádku (posloupnost tříd) na 1 000 řádků a pokrytí 20 nejčastějšími podpisy, zvlášť pro délky 4–8;
   (iii) kompresní poměr (zlib) posloupnosti tříd po řádcích vs po řádcích zamíchaných uvnitř (pos-shuf)
   a vs tokeny zamíchané mezi řádky stejné délky.
4. **Řádky stejné délky:** průměrná Jaccardova podobnost tříd mezi řádky stejné délky vs různé délky,
   po odečtení očekávání z marginál.
5. **Vertikální kopírování:** pro každý token s indexem i v řádku l (l > 0 v odstavci): editační vzdálenost
   k tokenu (l−1, i) vs minimum přes ostatní tokeny řádku l−1 vs tokeny řádku l vlevo; podíl shod ≤ 1
   editace „nad“ proti „jinde v řádku nad“ (a očekávání po permutaci indexů uvnitř řádku l−1, 200×).
   Totéž pro (l−1, i±1). Reportovat pro reál, ts2020 (má pravidlo p_same_position 14/28), ac-tri 86,
   šablonu.
6. **Konce řádků:** podíl tokenů končících na `-m`/`-g`/`-am` na posledním místě řádku vs jinde, a zda
   je poslední token řádku kratší (Steckley 2024) – deskriptivně, pro všechny zdroje.

## Nulové modely a kontroly
ts2020 calB (30 seedů), šablona 005 (30), ac-tri 86 (převzít/30), míchání v řádku, míchání uvnitř
poziční třídy (007), permutace řádků stejné délky; kladné kontroly z 009 (latinská próza mechanicky
lámaná; Linnaeus/Culpeper po položkách), jsou-li 009 výstupy v `main` – jinak Plinius a Linnaeus
načíst stejně jako v 009 (sdílet kód přes `voytools/controls.py`, pokud existuje).

## Výstupy (`results/010_line_as_item/`)
`boundary_vs_inside.csv`, `line_header.csv`, `line_templates.csv` (+ `signatures_top20.csv`),
`same_length.csv`, `vertical_copy.csv`, `line_end.csv`, `runs.csv`/`meta.json`, grafy
(`boundary.png`, `vertical_copy.png`, `templates.png`), `report.md`, `notes.md`.

## Kritéria přijetí
1. Všech 6 metrik pro reál (ZL/GC, split/join) a všechny nully/kontroly, s intervaly.
2. Report odpoví ano/ne (ZL, GC): (i) je řádek uzavřenější jednotka než u prózy a než u generátorů?
   (ii) nese první token řádku informaci o zbytku? (iii) existují šablony řádků nad náhodu a nad ts2020?
   (iv) je podobnost „k tokenu nad“ zvýšená – a reprodukuje ji ts2020 tak, jak ji má reál?
3. Podezřelé: poměr (b)/(a) > 1,5; podíl „nad ≤ 1 editace“ > 0,5; kompresní poměr < 0,2.
4. Testy procházejí; Python 3.10 / pandas 2.3; bez scipy; nové funkce (např. `stats.line_boundary_mi`,
   `stats.vertical_similarity`) s docstringy a testy.

## Rozpočet
Metriky jsou lehké (řádkové tabulky); 100 subsamplů × 4 sady + 90 vzorků nullů + kontroly ≈ 1 h.
Jedna session, jedna větev, PR.
