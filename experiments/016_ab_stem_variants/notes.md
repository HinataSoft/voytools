# 016 – deník session

## Rozhodnutí a interpretace spec

* **Kód.** Rozšíření `voytools/registers.py` (žádná změna existujících signatur): `levenshtein`, `edit_script`
  (jedno optimální zarovnání; sousední vložení / smazání sloučena do jedné operace, např. `ins:ed@end`),
  `nearest_counterparts` (index symetrických smazání – SymSpell: dvě posloupnosti ve vzdálenosti ≤ d vždy sdílejí
  řetězec vzniklý ≤ d smazáními na každé straně; přesné pro d ≤ 2, bez průchodu všech dvojic), `counterpart_rates`,
  `edit_scripts`, `edit_inventory`, `coverage_levels`, `substitution_sharing`, `shuffle_units`, `profile_test`.
  Testy `tests/test_registers.py` (Levenshtein proti hrubé síle, úplnost hledání do 2 editací proti všem dvojicím,
  inventář a substituce na ručním příkladu). Drobná oprava: `transfer_matrix` kopíruje pole z `crosstab`
  (pandas 3 vrací read-only pohled; na pandas 2.3 beze změny).
* **Jednotky.** Kmeny / tokeny jako n-tice jednotek: ZL glyfy presetu `extended` (`ch`, `sh`, `cth`…, `iin`, `in`, `ee`
  = jedna jednotka), GC znaky v101, latina písmena, reprezentace `bpe` = jednotky BPE se 64 sloučeními naučenými
  na celém daném korpusu (reál: celý text P; nully: každý seed na sobě). Jedna editace = jedna jednotka.
* **Reprezentace:** `a` kmen dekompozice (a), `b` kmen (b) (+ poziční glyfy, jen ZL/GC), `tok` celý token,
  `bpe` celý token v BPE jednotkách. Prázdné kmeny vynechány.
* **X-only kmen** = ≥ 2 tokeny v X a 0 tokenů v Y (ve stejném losu); kandidáti = kmeny Y s ≥ 2 tokeny. Inventáře
  pro Jaccard také ≥ 2 (012 používal ≥ 3; hodnoty Jaccard proto nejsou s 012 přímo srovnatelné).
* **Směry.** `x_only` = X-only kmeny proti Y, `y_only` = Y-only proti X. U A×B je X = A, Y = B, tedy `y_only`
  = **B-only vs A** (hlavní otázka) a `x_only` = A-only vs B. Operace se vždy zapisují ve směru A → B
  (X → Y): u B-only kmene je to skript z A protějšku do B kmene. Pozice `begin`/`mid`/`end` vzhledem ke
  zdrojovému (A) tvaru.
* **Stejná velikost** jako 012: N tokenů na stranu = polovina menšího jazyka ve výřezu (celé stránky v náhodném
  pořadí, poslední oříznuta), 200 losů; poloviny A|A, B|B = dva disjunktní losy téhož jazyka. Výřezy celkem a herbář.
* **Null tvaru (`SHUF`).** Kandidáti = inventář A (≥ 2) s glyfy promíchanými napříč kmeny při zachování délek
  (zachová frekvence jednotek a délky, zničí skladbu); dotazy = B-only kmeny definované proti *reálnému* A.
* **Inventář editací** se počítá na skriptech dvojic ≤ 2 editace, vážených tokeny X-only kmene; `explained(k)`
  = podíl tokenů dvojic, jejichž celý skript leží v top-k operacích (operace řazeny podle tokenové váhy);
  pokrytí 50/80/95 % = nejmenší k. V `edit_inventory.csv` je inventář **sdružený přes 20 losů** (AB i poloviny
  stejným postupem); v `runs.csv` jsou `ops_50/80/95` z každého losu zvlášť.
* **Substituce** (metrika 3): v každém losu se naučí pořadí operací (in-sample) a X-only kmeny, jejichž skript
  leží v top-k, se **přepíší na svůj protějšek** (jen kmeny s protějškem ≤ 2; sdílené kmeny se nemění). Jaccard
  přepsaného inventáře proti druhé straně; `k = 0` = bez přepisu, `all` = všechny dvojice ≤ 2 sloučeny (horní
  mez). Stejný postup u polovin – ty tím rostou také, proto se A×B srovnává s polovinami při témže k i s
  polovinami při k = 0 (úroveň, kterou má „změna zápisu" dohnat). Plošná aplikace pravidla na všechny kmeny
  (např. každé `o-` → `qo-`) by vyžadovala kontext u vložení uvnitř slova; licencovaná varianta je k hypotéze
  změny zápisu velkorysá (horní mez).
* **Currierovy operace** (jen Eva): vložení/smazání složené jen z `e`, `ee`, `d`, `q` (tj. `ed`, `eed`, `d`, `qo`
  bez `o`…) nebo záměna `ch` ↔ `sh`; `currier_token_share` v `edit_coverage.csv`. U GC (v101) se značky
  neklasifikují automaticky, inventář se čte ručně.
* **Sekční kontrola** (metrika 5) na plném textu (A plné × B plné), reprezentace `a` a `tok`: JSD profilů dvojic
  (protějšek ≤ 1 / ≤ 2) vs 500 permutací přiřazení; referenční hladina = **týž kmen v A a v B** (sdílené kmeny
  ≥ 2 na obou stranách). Profil `sections` (všechny kódy ilustrací) je ovládán tím, že A a B mají jiné sekce
  (A: H, P; B: B, S, H, T), takže ani sdílené kmeny nic neukážou; proto přidán profil `herbal_binary`
  (podíl tokenů kmene v herbáři vs mimo něj – jediná sekce s hmotou v obou jazycích).
* **Nully:** ts2020 config 21 (ZL split/join) a šablona 005 (ZL/GC split), 30 seedů, A/B podle layoutu, 1 los
  na seed; generované korpusy cachované v `scratch/016_cache/` (gitignored).
* **Latinské kontroly** (`voytools.controls`, korpusy `lrozanova/voynich-units`, layout ZL split jako v 009):
  Plinius knihy 16, 18, 21 = „A", 23, 25, 28, 30 = „B" (stejný jazyk, jiné téma: stromy a rostliny vs léky);
  Caesar = „A", Celsus = „B" (dva kontrolní korpusy spojené, stránky s prefixem). Reprezentace `a`, `tok`, `bpe`;
  výřez jen celkem; N = polovina menší strany (Plinius ≈ 5,4 k, Caesar/Celsus 17 k tokenů).

## Průběh

* Jeden los ZL split ≈ 1,2 s (4 reprezentace × 2 výřezy × 2 směry × 4 srovnání); 200 losů ≈ 4–8 min na sadu.
  Všechny sady, nully i kontroly běžely paralelně na pozadí (`nohup`), výstupy v `results/016_ab_stem_variants/runs/`.
* Prostředí: Python 3.13 / pandas 3.0.5 (cloud); kód drží 3.10 / pandas 2.3 (žádná nová syntaxe, `copy=True`
  u `to_numpy`).
