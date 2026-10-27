# 008 – Doplnění kanonického nullu ts2020 do 005–007 a oprava „pozičního kmene“

## Otázka
Experimenty 005, 006 a 007 běžely bez `voytools.selfcitation` (004 bylo sloučeno později) a jejich
kopírovací null je jen provizorní `autocopy_like` config 86. Které z jejich nálezů (sekční informace
kmenů nezávislá na písaři, strmý pokles sdílení se vzdáleností, nulová reziduální MI za d = 2, posun
sufixů A/B uvnitř písaře) **přežijí nejsilnější null ts2020**? A je „pozice nese o kmeni víc než sekce“
(007) artefakt parseru, který nechává poziční glyfy v kmeni?

## Nulové modely
`selfcitation_like(real, seed=…, **config21)` s parametry z `results/004_selfcitation_port/best_configs.csv`
(řádek `ts2020_calB`: p_copy 0,1, replace_weights 20/50/30, p_same_position 14), 50 seedů
(006 a 007 stačí 30). ts2020 je Eva-only → jen ZL (split a join); u GC uvést, že null chybí.
Dále šablona 005 (`ParadigmTemplate`, plný model rank 1 s α = 5, λ = 0,3; z `results/005_…/best_configs.csv`)
30 seedů, ZL i GC – do 006 a 007 (sekční informace, pokles se vzdáleností), kde chyběla.

## Kroky
1. **007 znovu s ts2020 a šablonou:** `experiments/007_position_vs_content/run.py` s přidanými zdroji
   (`Context.corpus`; viz jeho `notes.md`) – rozklad informace (4 cílové proměnné) a MI(d) s korekcemi
   `pc` a `line`. Reportovat, kolik procent reálné hodnoty ts2020 dává u I(Y; pozice), I(Y; sekce | pozice),
   I(kmen; sekce | písař) a zda má ts2020 na d = 1 hranu (má pravidla jen pro začátky řádků a volbu
   zdroje „stejná pozice v řádku nad“).
2. **Parser:** zopakovat rozklad 007 (bod + 100 subsamplů, ZL i GC, split) se třemi variantami kmene:
   (a) výchozí naučené afixy; (b) afixy doplněné o poziční glyfy – prefixy `p`, `f`, `t`, `k`, `cp`, `cf`,
   `ct`, `ck` (gallows a jejich ligatury) a sufixy `m`, `g`, `am`, `om` (konce řádků); (c) kmen = prostřední
   BPE jednotka(y) po odstranění první a poslední jednotky z `learn_bpe(merges=64)`. Pro každou variantu
   I(kmen; pozice) vs I(kmen; sekce) vs I(kmen; sekce | písař). Závěr: je poziční převaha u kmene stabilní?
3. **006 znovu s ts2020 a šablonou:** `experiments/006_topic_at_distance/run.py pages --sources …`
   (model `spec+size`, K = 30, kontrast, křivka). ts2020 kopíruje jen uvnitř stránky: očekáváme plochou
   mezistránkovou křivku – pak reálný pokles 3–4× strmější než ac-tri **není** kopírováním vysvětlitelný
   vůbec. Šablona s λ = 0,3 a tématem dává sekční složku z konstrukce – porovnat velikost s reálem.
4. **A/B a pozice:** `jsd_suffix_AB` z 005 (celkově a u písaře 3) spočítat **podmíněně na poziční třídě**
   (22 tříd z 007; vážený průměr JSD uvnitř tříd) a zvlášť pro vnitřní tokeny řádku. Pokud posun A/B
   po podmínění pozicí zůstane ≥ ½ původního, není poziční. Null: permutace štítku jazyka mezi
   stránkami uvnitř písaře (200×).
5. **Tabulka pro syntézu:** pro každý řádek `results/SYNTHESIS.md` oddílu 1 doplnit sloupec ts2020
   (hodnota, z, reprodukuje ano/ne) a tmpl, kde chyběl.

## Výstupy (`results/008_ts2020_recheck/`)
`decomp_ts2020.csv`, `mi_distance_ts2020.csv`, `parser_variants.csv`, `pages_ts2020.csv` (OLS + kontrast
+ křivka), `ab_by_position.csv`, `synthesis_update.csv` (řádek syntézy × null × hodnota × z × verdikt),
`report.md`, `notes.md`. Žádné nové moduly v balíku; pokud je třeba, jen přidat zdroje do existujících
`run.py` (bez změny jejich výstupů pro dosavadní zdroje).

## Kritéria přijetí
1. Všech pět kroků má tabulku; `synthesis_update.csv` pokrývá všech 9 řádků oddílu 1 syntézy.
2. Report odpoví ano/ne (ZL; GC kde lze): přežije ts2020 (a) sekční informace kmenů nezávislá na písaři,
   (b) strmý pokles sdílení, (c) „žádná MI za d = 2“ jako *rozdíl* od nullu (nebo je to u nullu stejné),
   (d) A/B uvnitř písaře po podmínění pozicí; (e) je poziční převaha kmene artefakt parseru?
3. Reprodukce dřívějších čísel: reálné hodnoty se shodují s 005/006/007 (uvést max. rozdíl).
4. Testy procházejí; Python 3.10 / pandas 2.3; bez scipy.

## Rozpočet
ts2020 vzorek 5–15 s; 50 + 30 + 30 vzorků s metrikami ≈ 30 min na 4 procesech; parser varianty
3 × (1 + 100) × 2 soubory ≈ 20 min; A/B ≈ 10 min. Celkem ≈ 1–1,5 h, jedna session, jedna větev, PR.
