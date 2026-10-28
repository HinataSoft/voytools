# 017 – Dá lokální mechanismus pozici, šev a soudržnost slova zároveň?

## Otázka
014: ts2020 + naučené tabulky R1 (začátek řádku) / R2 (konec řádku) / R3 (prefix | poslední glyf předchozího tokenu)
dožene poziční vazbu prefixu i sufixu, vazbu sufix → prefix dalšího tokenu a 2/3 okrajové MI (0,126 vs 0,186), ale
**pokazí** vazbu prefix–sufix (`mi_prefix_suffix` 0,094 vs 0,205), sekční vazbu prefixu, záhlaví řádku a přidá
~2 000 typů; zbytek švu sedí v tokenech bez prefixu (R3 jim první glyf nemění) a obohacení `-m/-g` na konci řádku
zůstává 6–8× místo 17× (generátor píše `-m` i uvnitř řádku). R3 byl přepis **nezávislý** na zbytku slova – nejjednodušší
možná realizace. Otázka: **stačí přepis podmíněný i na zbytku slova (sufix, první glyf kmene), aby kopírovací generátor
dožene pozici, šev, konec řádku a vazbu prefix–sufix najednou bez nafouknutí slovníku – nebo některá kombinace
vlastností lokálním přepisem nejde a je to ona, co Voynich odlišuje?**

## Hypotézy a predikce
- **Konvence (lokální gramatika stačí):** varianta R3′ (prefix | poslední glyf předchozího, sufix) vrátí
  `mi_prefix_suffix` na reál (|z| < 3) a šev drží; R3″ (první glyf kmene u tokenů bez prefixu | poslední glyf
  předchozího) dožene zbylou třetinu švu; R2′ (odstranění koncových glyfů uvnitř řádku) dožene 17×; počet typů a hapax
  zůstanou v intervalu ts2020 (± 10 %). Pak oddíl 1 syntézy (kromě A/B) = konvence reprodukovatelné lokálním přepisem.
- **Podpis (lokální přepis nestačí):** každá varianta, která dožene šev, pokazí něco jiného (prefix–sufix, sekční
  vazbu prefixu, záhlaví řádku, typy, H(sufix | kmen)); kombinace *pozice + šev + prefix–sufix + šířka paradigmatu*
  zůstane mimo dosah; report ji pojmenuje jako nejmenší nedosažitelnou kombinaci.

## Data a filtry
ZL3b (ts2020 je Eva-only; GC nelze – uvést), `P`, clean, split i join (přepis se učí a hodnotí zvlášť pro obě politiky);
dekompozice (a). Tabulky se učí z reálu (jako 014); u podmíněných tabulek **backoff**: buňka s < 5 pozorováními →
méně podmíněná tabulka (R3′ → R3 → nepodmíněné); podíl tokenů přepsaných z backoffu hlásit (`table_support.csv`).

## Kroky
1. **Varianty post-procesoru** (rozšíření `voytools/placerules.py`; R1–R3 z 014 nechat beze změny, přidat funkce):
   - **R3′a**: prefix | (poslední glyf předchozího tokenu, sufix tokenu); kmen a sufix zachovány.
   - **R3′b**: prefix | (poslední glyf předchozího, první glyf kmene).
   - **R3′c**: prefix | (poslední glyf předchozího, sufix, poziční třída 007 hrubá: první/vnitřní/poslední token,
     první řádek odstavce ano/ne) – nejpodmíněnější člen rodiny, horní mez.
   - **R3″**: u tokenů bez prefixu přepis prvního glyfu kmene | poslední glyf předchozího (tabulka z tokenů bez prefixu).
   - **R2′**: R2 + u tokenů uvnitř řádku přepis posledního glyfu z tabulky „uvnitř řádku“, pokud je glyf z množiny
     koncových (`m`, `g`; naučit: glyfy s poměrem konec/vnitřek > 5).
   - **Licence typů** (přepínač): přepis se přijme jen tehdy, když výsledný token existuje ve slovníku reálu
     (≥ 2 výskyty); jinak losovat znovu (max 5×), pak nechat původní. Měřit, kolik přepisů licence zamítne.
   Sady k vyhodnocení (každá 30 seedů): ts2020; +R1+R2 (014, znovu spočítat touž cestou); +R1+R2+R3 (014);
   +R1+R2′+R3′a; +R1+R2′+R3′b; +R1+R2′+R3′c; +R1+R2′+R3′a+R3″; totéž poslední s licencí typů. Celkem 8 sad.
2. **Profil** každé sady: metriky 009 (37) + řádkové 010 + rozklad 011 (vysvětlený podíl H(sufix | kmen)) + rozklad
   švu 014 (páry, zbytek nad glyfovým párem) + `mi_prefix_suffix`, H(sufix | kmen), sufixů na kmen, typy, hapax,
   sekční vazba prefixu, záhlaví řádku; z-skóre proti reálu podle SD 30 seedů **a** podle bootstrapu reálu (100×;
   přísnější z obou uvést – 014 upozorňuje, že SD seedů je u některých metrik příliš malá).
3. **Seznamy:** pro každou vlastnost oddílu 1 syntézy (bez A/B) a pro cíle kalibrace 004 (loss, H2, typy, hapax,
   paradigma): dosaženo / pokaženo / beze změny (|z| < 3 split i join) – tabulka vlastnost × sada. Z ní **nejmenší
   kombinace** vlastností, kterou žádná sada nedá najednou (pokud existuje).
4. **Totéž na šabloně 005** (s tématem; 30 seedů) pro sady +R1+R2′+R3′a a +R1+R2′+R3′a+R3″: šablona má vazbu
   prefix–sufix z konstrukce (paradigmata) – zachová ji podmíněný přepis? Závěr nesmí záviset na základním generátoru.
5. **Cena přepisu:** podíl přepsaných tokenů, podíl z backoffu, podíl zamítnutý licencí, počet nových typů proti
   základnímu generátoru; graf `rewrite_cost.png`.

## Nulové modely a kontroly
Reál ZL split/join s bootstrapem (100×); ts2020 config 21 (30); šablona 005 (30). Kontrola správnosti: R3′ s tabulkou
naučenou z **ts2020 samotného** (ne z reálu) nesmí změnit profil ts2020 (|z| < 1 u všech metrik) – test, že přepis
neinjektuje nic mimo tabulky.

## Výstupy (`results/017_conditional_seam_rewrite/`)
`augmented_profile.csv` (sada × metrika × průměr × SD × z_seed × z_boot), `augmented_lists.csv` (vlastnost × sada:
dosaženo/pokaženo/beze změny), `minimal_unreachable.csv`, `table_support.csv`, `rewrite_cost.csv`,
`template_check.csv`, `selfcheck.csv`, grafy (`augmented.png`, `rewrite_cost.png`), `runs.csv`/`meta.json`,
`report.md`, `notes.md`. Funkce do `voytools/placerules.py` (R1–R3 beze změny; testy `tests/test_placerules.py`).

## Kritéria přijetí
1. Všech 8 sad × 30 seedů, split i join; z-skóre oběma způsoby; samokontrola (ts2020 tabulky na ts2020) v mezích.
2. Report odpoví: (i) vrátí R3′ vazbu prefix–sufix při zachování švu? (ii) dožene R3″ zbylou třetinu švu a R2′ 17×?
   (iii) kolik to stojí typů/hapaxu a kolik tokenů jde z backoffu; (iv) existuje kombinace vlastností, kterou žádná
   sada nedá – která, a platí to i na šabloně? Věta: konvence (lokální přepis stačí) / podpis (nestačí; kombinace).
3. Podezřelé: sada s |z| < 3 u všech metrik včetně typů (zkontrolovat, zda přepis nekopíruje reál přes příliš jemné
   tabulky – podíl buněk s < 5 pozorováními, podíl z backoffu), samokontrola mimo meze (chyba přepisu).
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
014: 3 sady × 30 seedů ≈ 1 h. Zde 8 + 2 sady ≈ 2–2,5 h; bootstrap reálu minuty. Běhy po sadách přes `run_seeds`,
commit po každé sadě.
