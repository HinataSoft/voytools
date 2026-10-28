# 016 – Kmeny jen v B: glyfové varianty kmenů A, nebo nová slova?

## Otázka
012: A a B sdílejí méně kmenů než poloviny téhož jazyka, posun sufixů je difúzní (víc tvarů → `-dy`), přechod ostrý.
Výhrada 012: pokud B mění i glyfy uvnitř „kmene“ (`ed`, `qo`, `ch`/`sh`), parser rozdělí stejné slovo jinak a sdílení
klesne uměle. Otázka: **jsou kmeny vyskytující se jen v B (a jen v A) systematické glyfové varianty kmenů druhého
jazyka (změna zápisu), nebo nové kmeny bez protějšku (nový slovník)?**

## Hypotézy a predikce
- **Změna zápisu:** většina B-only kmenů (vážená tokeny) má protějšek v A ve vzdálenosti ≤ 1–2 editace, přičemž editace
  tvoří malou sadu opakovaných operací (vložení `e`/`ed`, `o`→`qo`, `ch`↔`sh`, `ol`→`ol`+`d`…); po aplikaci naučené
  substituce Jaccard A×B dosáhne polovin téhož jazyka; totéž platí pro A-only vůči B.
- **Nový slovník:** protějšky nad náhodu nejsou (poměr k náhodným párům kmenů ≈ 1), nebo jsou, ale editace jsou
  různorodé (žádná malá sada); Jaccard po substituci zůstává pod polovinami.

## Data a filtry
ZL3b a GC2a, `P`, clean, split a join; dekompozice (a), (b) a navíc **bez dekompozice** (celé tokeny) – aby parser
nebyl ve hře. Jazyk A/B stránek; výřezy celkem a herbář. Vyrovnané velikosti jako 012 (200 losů).

## Metriky
1. **Protějšky:** pro každý B-only kmen (≥ 2 výskyty) nejbližší A kmen podle Levenshteina na glyfech (preset
   `extended`; u GC v101 znaky); podíl B-only tokenů s protějškem ve vzdálenosti ≤ 1 a ≤ 2; očekávání: totéž pro
   B-only kmeny proti **polovině B** (vnitřní kontrola: jak často má kmen náhodou blízký jiný kmen) a proti náhodně
   zamíchaným glyfům kmenů A (null tvaru). Symetricky A-only vs B.
2. **Inventář editací:** pro protějšky ≤ 2 spočítat operace (vložení/smazání/substituce konkrétních glyfů, s pozicí
   ve slově: začátek/vnitřek/konec); seřadit; pokrytí 50/80/95 % tokenů; srovnat s inventářem editací mezi
   polovinami téhož jazyka (šum parseru) a s R&T/Currierovými rozdíly (`ed`, `qo`, `chol`).
3. **Substituce a sdílení:** naučit mapování A→B z top-k editací (k = 1, 3, 5, 10), aplikovat na kmeny A, spočítat
   Jaccard/pokrytí A'×B vs poloviny; kolik editací je třeba, aby sdílení dosáhlo úrovně polovin (pokud vůbec).
4. **Celé tokeny:** totéž (1–3) na celých tokenech bez parseru; a varianta na BPE jednotkách (64 sloučení).
5. **Sémantická kontrola:** protějšky by měly zachovat **sekční rozložení**: kmen A a jeho B varianta se vyskytují
   v týchž typech sekcí/ilustrací nad náhodu? (srovnat JSD sekčních profilů dvojic protějšek vs náhodné dvojice).

## Nulové modely
Poloviny téhož jazyka (200×), ts2020 (A/B podle layoutu; 30 seedů) – protějšky a editace mezi jeho „A“ a „B“ jsou
šum, šablona 005 (30), kontroly: latina – kniha vs kniha (Plinius) jako „stejný jazyk, jiné téma“ a Caesar vs Celsus
jako „stejný jazyk, jiný registr“ (co dělá registr v jazyce s flexí).

## Výstupy (`results/016_ab_stem_variants/`)
`counterparts.csv` (kmen, protějšek, vzdálenost, operace, tokeny), `counterpart_rates.csv`, `edit_inventory.csv`,
`substitution_sharing.csv`, `whole_tokens.csv`, `section_check.csv`, grafy, `runs.csv`/`meta.json`, `report.md`,
`notes.md`. Funkce do `voytools/registers.py` (rozšíření; s testy).

## Kritéria přijetí
1. Všechny míry pro ZL i GC, split i join, tři reprezentace (afixy, celé tokeny, BPE), s nully a kontrolami.
2. Report odpoví: (i) podíl B-only (A-only) tokenů s protějškem ≤ 1 / ≤ 2 nad vnitřní kontrolu; (ii) kolik editací nese
   80 %; jsou to Currierovy znaky? (iii) kolik editací vrátí sdílení na úroveň polovin; (iv) zachovávají protějšky
   sekční profil? – ano/ne, ZL/GC; věta: změna zápisu, nový slovník, nebo obojí (s podíly).
3. Podezřelé: podíl protějšků ≤ 1 u polovin téhož jazyka > 0,9 (parser příliš jemný – uvést), inventář s jedinou
   editací > 80 % (zkontrolovat ligatury).
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
Levenshtein mezi ~2 k × 2 k kmeny na glyfech: sekundy–minuty v čistém Pythonu (omezit kandidáty délkou ±2);
200 losů × 3 reprezentace × 2 soubory ≈ 30–45 min; nully ≈ 20 min. ≈ 1–1,5 h.
