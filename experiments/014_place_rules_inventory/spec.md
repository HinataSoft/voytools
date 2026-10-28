# 014 – Inventář místních pravidel: kolik konvencí dělá pozici a šev slov

## Otázka
Po 013 zůstávají jako Voynich-specifické jen vlastnosti vázané na místo: I(prefix; pozice) 0,13 bitu (ts2020 dá
polovinu, próza, katalog, blábol ani šifra nic), okrajová MI na švu slov 0,19 bitu a vazba sufix → prefix dalšího
slova 0,038 bitu, konec řádku (`-m/-g` 17×) a uzavřený řádek. Otázka: **jsou to malé sady konvencí (několik glyfových
pravidel, která lze přidat generátoru), nebo difúzní jevy rozprostřené přes celý inventář? A jsou stejné napříč
sekcemi, jazyky A/B a písaři (konvence), nebo se mění (něco jiného)?**

## Hypotézy a predikce
- **Konvence:** 80 % každé z těchto informací nese < 10 glyfů/párů/tříd (např. začátek řádku: `d-`, `y-`, `s-`, gallows
  `p-/f-` v prvním řádku odstavce; konec: `-m`, `-g`, `-am`; šev: `y→q`, `r/s→a`, `-n→ch/sh` vyhýbání); pravidla jsou
  stabilní napříč A/B, písaři a sekcemi; ts2020 doplněný o tyto tabulky reál dožene.
- **Difúzní jev:** informace je rozprostřená, 80 % vyžaduje desítky položek, liší se mezi A/B nebo sekcemi; doplněný
  generátor nedožene.

## Data a filtry
ZL3b a GC2a, `P`, clean, split a join; dekompozice (a) a (b). Poziční třídy jako 007 (22) + jemné `rel_pos`.

## Kroky
1. **Rozklad I(prefix; pozice) a I(sufix; pozice)** na příspěvky hodnot (pointwise: p(x,y)·log(p(x,y)/(p(x)p(y))) sečtené
   po afixu a po poziční třídě); seřadit; kolik afixů / tříd pokryje 50 %, 80 %, 95 %. Totéž pro první glyf tokenu
   (místo afixu) a pro poslední glyf. Oddělit: (i) první řádek odstavce vs ostatní, (ii) první/poslední token řádku vs
   vnitřní (`rel_pos`) – kolik z poziční informace je „kraje řádku“ a kolik „drift uvnitř“ (Feaster).
2. **Rozklad okrajové MI** (poslední→první glyf uvnitř řádku) na páry; kolik párů dává 50/80/95 %; totéž po odečtení
   y→q; a vazba sufix → prefix dalšího tokenu (011) na páry (sufix, prefix). Test, zda je šev **jen glyfový**:
   I(poslední glyf; první glyf dalšího) vs I(sufix; prefix dalšího | poslední glyf, první glyf) – zbývá něco nad
   glyfový pár? (jazyk: shoda tvarů by dala zbytek; konvence: nic).
3. **Stabilita:** všechny rozklady zvlášť pro A, B, každého písaře (≥ 3 000 tokenů), každou sekci; JSD mezi
   příslušnými tabulkami (např. rozdělení prefixů na začátku řádku v A vs B) vs JSD mezi náhodnými polovinami
   téže skupiny (200×). Vyjádřit jako „podíl pravidel sdílených“.
4. **Doplněný generátor** (jediná větev; nový modul `voytools/placerules.py`, `nulls.py` neměnit): post-procesor,
   který vezme výstup ts2020 (nebo šablony 005) a aplikuje naučené tabulky: (R1) rozdělení prefixů/prvních glyfů
   pro první token řádku a pro první řádek odstavce (přepis prefixu podle tabulky při zachování kmene a sufixu),
   (R2) sufixy/poslední glyfy posledního tokenu řádku, (R3) švová tabulka – pro každý token vylosovat první glyf/
   prefix podmíněně na posledním glyfu předchozího tokenu (přepis prefixu). Varianty: R1, R1+R2, R1+R2+R3.
   Pro každou změřit celý profil (009 metriky + 011 rozklad + 010 řádkové metriky), 30 seedů, a zjistit, které
   Voynich-specifické vlastnosti doplněný generátor dožene (|z| < 3, split i join) a co se tím **pokazí**
   (cíle kalibrace, hapax, paradigma).
5. **Blábol a katalog pro srovnání** (z 013/009): stejný rozklad kraje řádku u Linnéa po položkách a u blábolu ve
   vlastních řádcích – jak vypadá „konvence“ u položkového textu.

## Výstupy (`results/014_place_rules_inventory/`)
`position_contributions.csv`, `seam_contributions.csv`, `coverage.csv` (kolik položek na 50/80/95 %), `stability.csv`,
`augmented_profile.csv` + `augmented_lists.csv`, grafy (`position_pareto.png`, `seam_pareto.png`, `augmented.png`),
`runs.csv`/`meta.json`, `report.md`, `notes.md`.

## Kritéria přijetí
1. Rozklady a pokrytí pro ZL i GC, split i join; stabilita pro A/B, písaře, sekce.
2. Report odpoví: (i) kolik položek nese 80 % poziční informace prefixu, sufixu, prvního/posledního glyfu; (ii) kolik
   párů nese 80 % švu a zda nad glyfovým párem zbývá něco; (iii) jsou pravidla sdílená napříč A/B, písaři, sekcemi
   (podíl, JSD vs poloviny); (iv) které vlastnosti oddílu 1 syntézy doplněný ts2020 dožene a co se pokazí.
3. Podezřelé: pokrytí 80 % jednou položkou (chyba výpočtu), JSD mezi skupinami < JSD polovin.
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
Rozklady: tabulky, minuty. Doplněný generátor 3 varianty × 30 seedů × profil ≈ 1 h. Celkem ≈ 1,5–2 h.
