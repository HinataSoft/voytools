# 019 – Tvary jen v B: jména, jiná morfologie, nebo slovník druhého systému?

## Otázka
016: pětina tokenů B (0,21–0,23; poloviny téhož jazyka 0,05–0,06) jsou tvary, které A nemá a které nejsou jeho
glyfovými variantami; jedinou systematickou operací je vložení `d` (`-ey → -edy`, v101 `8`; 28–46 % protějšků), zbytek
je dlouhý ohon. Latina téhož jazyka s jiným tématem (Plinius, jiné knihy) i s jiným registrem (Caesar vs Celsus) dává
stejný obraz. Otázka: **z čeho se ta pětina skládá – z `-edy` variant sdílených kmenů (jiný zápis/paradigma), z nových
tvarů sdílených kmenů (jiná morfologie), nebo z nových kmenů – a chovají se nové kmeny jako jména/položky (úzké
paradigma, shluky na stránkách, popisky), jako registrová slova (široké paradigma, rovnoměrně), nebo stejně jako
ostatní slovník B?** A je `-ey → -edy` deterministické pravidlo kontextu (změna zápisu), nebo volba tvaru (011)?

## Hypotézy a predikce
- **Jména / položky (obsah):** nové kmeny tvoří velkou část pětiny; mají úzká paradigmata (1–2 tvary), jsou shlukované
  (málo stránek na výskyt, nízká disperze), sekčně soustředěné, častěji v popiscích a na krajích řádků; glyfová
  stavba se neliší od zbytku (nejsou „cizí“). Vzor Plinius-téma (jména rostlin a zvířat v jiných knihách).
- **Registr / morfologie:** pětinu nesou nové tvary sdílených kmenů a `-edy` varianty; nové kmeny mají široká
  paradigmata jako sdílené, jsou rozprostřené; vzor Caesar/Celsus.
- **Slovník druhého systému:** nové kmeny se při vyrovnané frekvenci neliší od sdílených v ničem měřitelném.
- **Pravidlo `-edy`:** deterministické – P(`-edy` | kontext uvnitř slova) blízko 0 nebo 1 pro většinu kontextů v B
  a kmeny nemají oba tvary; volba tvaru – P v pásmu 0,3–0,7, kmeny mají oba tvary i na téže stránce.

## Data a filtry
ZL3b a GC2a, `P`, clean, split (hlavní), join; celé tokeny (nosná reprezentace 016) a dekompozice (a) pro kmeny.
Vyrovnané velikosti jako 012/016 (celé stránky, 200 losů); protějšky a operace z `voytools.registers` (016). Třídy
tvarů B (symetricky A): **sdílené** (≥ 1 v A), **B-only** (≥ 2 v B, 0 v A) rozdělené na (i) `d`-varianty (protějšek
≤ 1 operací `ins d` / `ins e(e)` / `ch↔sh` / `q`), (ii) jiný protějšek ≤ 1, (iii) bez protějšku ≤ 1; hapaxy B zvlášť
(iv). Kmen tvaru: (a); **nový kmen** = kmen bez jediného tokenu v A.

## Metriky
1. **Skladba pětiny:** podíl tokenů B-only v (i)–(iv); podíl B-only tokenů s novým kmenem vs novým tvarem sdíleného
   kmene; totéž pro poloviny B (šum) a pro A-only.
2. **Profil tříd** (každá třída vs sdílené tvary **vyrovnané frekvencí**, 200 losů; a vs táž třída z polovin B):
   paradigma kmene (tvarů na kmen, H(sufix | kmen)); disperze (stránek na výskyt, Griesův DP, podíl výskytů na
   nejčastější stránce); sekční entropie; podíl v popiscích (`L` loci) a na krajích řádků (`is_first`, `is_last`,
   `is_par_first_line`); souvýskyt dvou B-only tvarů v témž řádku vs náhoda; glyfová stavba: log-věrohodnost bigramů
   pod modelem sdílených tvarů B (jsou „cizí“?); šev: rozdělení prvního glyfu následujícího tokenu.
3. **Pravidlo `-edy`:** všechny tokeny A i B končící `-ey` / `-edy` (v101 zřejmě `-c9` / `-c89`; ověřit podle 016
   a `glyphs.py`, uvést):
   P(`-edy` | glyf před `-ey`, prefix, pozice v řádku, poslední glyf předchozího tokenu) v B a v A; entropie volby
   při daném kontextu vs entropie bez kontextu; podíl kmenů B s oběma tvary, a z nich podíl s oběma na téže stránce;
   srovnat s jinou dvojicí sufixů téhož paradigmatu (`-y`/`-dy`, `-ol`/`-or`) jako měřítkem běžné variace 011.
4. **Latinské vzory:** Plinius téma a Caesar/Celsus registr (016, 200 losů): tytéž třídy a profily (bod 2) –
   referenční hodnoty pro „jména“ a „registr“; ts2020 a šablona (A/B podle layoutu, 30 seedů) jako šum.
5. **Seznamy:** 30 nejčastějších nových kmenů B s počtem tvarů, stránkami, sekcemi, popisky; totéž A.

## Nulové modely a kontroly
Poloviny B (200×), vyrovnání frekvencí (200×), ts2020 config 21 (30), šablona 005 (30), Plinius téma, Caesar/Celsus
registr; u bodu 3 permutace kontextů (200×).

## Výstupy (`results/019_b_only_forms/`)
`composition.csv`, `classes.csv` (tvar, třída, počty, protějšek, operace, kmen nový/sdílený), `class_profiles.csv`
(zdroj × třída × metrika × hodnota × vyrovnaná základna × CI), `edy_rule.csv`, `edy_contexts.csv`, `latin_reference.csv`,
`top_new_stems.csv`, grafy (`composition.png`, `class_profiles.png`, `edy_rule.png`), `runs.csv`/`meta.json`,
`report.md`, `notes.md`. Nový modul `voytools/bforms.py` (používá `registers.py`; testy `tests/test_bforms.py`).

## Kritéria přijetí
1. Vše pro ZL i GC, split i join (join alespoň body 1 a 3), oba směry (B-only hlavní, A-only kratší).
2. Report odpoví: (i) skladba pětiny v podílech (d-varianty / nové tvary sdílených kmenů / nové kmeny / hapaxy);
   (ii) chovají se nové kmeny jako Plinius-téma (jména), jako Caesar/Celsus (registr), nebo jako sdílené (druhý
   systém) – po metrikách; (iii) je `-edy` deterministické, nebo volba tvaru. Věta s podíly.
3. Podezřelé: třída z polovin B s profilem odlišným od sdílených nad interval (chyba vyrovnání), bigramová
   věrohodnost nových kmenů = sdílených na 3 desetinná místa (chyba modelu), P(`-edy`) bez rozptylu přes kontexty.
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
Tabulkové operace; 200 losů × 2 soubory × 2 politiky ≈ 20 min; latina 200 losů ≈ 10 min; generátory 60 vzorků ≈ 20 min.
≈ 1–1,5 h.
