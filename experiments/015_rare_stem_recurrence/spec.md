# 015 – Opakování vzácných kmenů na dálku

## Otázka
Sdílení slovníku mezi stránkami řídí blízkost (006), ale měřeno kosinem frekvenčních vektorů, kde dominují časté kmeny.
Obsah (jména rostlin, hvězd, surovin) by se projevil jinde: **vzácný kmen se vrátí na vzdáleném místě, kde se vrací
i jeho „věc“**, zatímco kopírování dá jen lokální opakování a blábol žádné nad náhodu. Otázka: **opakují se vzácné kmeny
(2–5 výskytů) na vzdálených stránkách nad náhodu, a vrací se jako přesný token, nebo jako jiný tvar téhož kmene?**

## Hypotézy a predikce
- **Obsah:** rozdělení vzdáleností mezi výskyty vzácného kmene má navíc k lokálnímu vrcholu i dálkovou složku nad
  permutací (víc dvojic výskytů ≥ 16 stránek od sebe, než dá permutace tokenů kmene mezi stránkami při zachování
  frekvence); dálkové dvojice sedí častěji do téže sekce / téhož typu ilustrace než náhodné dvojice stránek; opakování
  na dálku jsou častěji **jiný tvar téhož kmene** (jméno v jiném pádu) než u blízkých (kopie).
- **Kopírování (ts2020):** opakování jen lokální (stejná stránka, sousední), dálkové ne nad náhodu; tvar spíš přesná kopie.
- **Blábol / šablona bez tématu:** žádná dálková složka nad náhodu; šablona s tématem α = 5 (005) ji má z konstrukce –
  porovnat velikost.

## Data a filtry
ZL3b a GC2a, `P`, clean, split (hlavní), join; dekompozice (a) a (b); kmeny s 2–5 výskyty (hlavní), 6–20 (kontrola);
kmeny s ≥ 1 výskytem v popiscích (`L`) zvlášť (jsou mezi nimi jména?). Vzdálenost = pořadí stránek v souboru;
sekce = `illustration`; navíc „stejný bifolio/kvatern“ z `$B`/`$Q`.

## Metriky
1. **Rozdělení vzdáleností** dvojic výskytů téhož vzácného kmene (v stránkách; 0 = táž stránka), biny 0 / 1 / 2–3 /
   4–7 / 8–15 / 16–31 / 32+; podíl dvojic v binu vs očekávání z permutace: tokeny kmene přiděleny náhodným stránkám
   s pravděpodobností úměrnou velikosti stránky (200×); **dálkový přebytek** = podíl(≥ 16) − očekávání, s intervalem.
2. **Dálkové dvojice a sekce:** podíl dálkových dvojic (≥ 16) se stejnou sekcí / stejným jazykem / stejným písařem vs
   očekávání z náhodných dvojic stránek stejné vzdálenosti (vzdálenostně párovaná permutace, 200×).
3. **Tvar opakování:** pro dvojice výskytů téhož kmene podíl identických tokenů (stejný prefix i sufix) vs jiný tvar,
   jako funkce vzdálenosti (0, 1–3, 4–15, ≥ 16); srovnat s očekáváním z paradigmatu kmene (nezávislá volba tvaru).
4. **Shluky na dálku:** počet trojic a více výskytů téhož vzácného kmene rozložených přes ≥ 3 různé kvaterny vs
   permutace.
5. **Popisky:** kmeny vyskytující se v popiscích (Lz, Lf, Lc, Ls, Ln) a zároveň 2–5× v textu: vzdálenostní profil
   zvlášť (jsou to kandidáti na jména).

## Nulové modely a kontroly
Permutace stránek pro tokeny kmene (výše), permutace **celých stránek** uvnitř sekce a globálně (006), ts2020 config 21
(30 seedů), šablona 005 s tématem a bez tématu (30 + 30), ac-tri 86 (30); kontroly: latinská próza (Plinius, Celsus) a
Linnaeus po položkách ve layoutu ZL (009), blábol (013; krátký – jen biny do 8 stránek).

## Výstupy (`results/015_rare_stem_recurrence/`)
`distance_profile.csv` (zdroj × bin × podíl × očekávání × přebytek × interval), `far_pairs_sections.csv`,
`form_by_distance.csv`, `far_clusters.csv`, `label_stems.csv`, grafy (`distance_profile.png`, `form_by_distance.png`),
`runs.csv`/`meta.json`, `report.md`, `notes.md`. Funkce do `voytools/recurrence.py` s testy.

## Kritéria přijetí
1. Profily pro reál (ZL/GC, split/join, obě dekompozice, 2–5 i 6–20), všechny nully a kontroly.
2. Report odpoví: (i) je dálkový přebytek > 0 (ZL i GC) a jak velký proti próze a proti ts2020/šabloně s tématem?
   (ii) sedí dálkové dvojice do téže sekce nad náhodu? (iii) jsou dálková opakování jiným tvarem častěji než blízká?
   (iv) popisková jména: mají dálkový profil?
3. Podezřelé: dálkový přebytek > 0,3 (chyba očekávání), permutace s nenulovým průměrným přebytkem (chyba nullu).
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
Tabulkové operace nad ~3 000 vzácnými kmeny; 200 permutací × 8 sad ≈ 15 min; generátory 120 vzorků ≈ 30 min. ≈ 1 h.
