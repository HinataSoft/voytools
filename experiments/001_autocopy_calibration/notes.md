# 001 – deník session

## Rozhodnutí a výklad spec

* **Generátor:** `autocopy_like(real, allowed_from=<všechny reálné tokeny>, seed=..., **config)`.
  Seed tokeny = prvních 100 tokenů (výchozí), ale bigramový model legálních editací se odhaduje
  ze všech reálných slov (jak doporučuje docstring `autocopy_generate` a `examples/03`); s modelem
  ze 100 seed tokenů by editace byly téměř nemožné. `p_new` tedy čerpá jen ze 100 seed tokenů.
* **Kódování běhu:** `seed = config_id*100 + rep`; tentýž `seed` jde do generátoru i do
  míchání uvnitř metrik. Fáze A (rep 0–2) a fáze B (rep 0–49) sdílejí prvních 5 vzorků téže
  konfigurace (deterministicky stejné).
* **Fáze A** počítá cílové metriky + edge/token MI (spec dovoluje vynechat BPE a paradigmata).
  Paralelně 4 procesy (`phaseA/part0..3`), sloučeno `merge_runs` do `runs.csv`.
* **Fáze B:** 10 nejlepších konfigurací z fáze A × 50 seedů, všechny metriky, pro ZL split
  (hlavní) a jako kontroly ZL join, GC split, GC join (stejných 10 konfigurací, vybraných na ZL).
  `best_configs.csv` = pořadí podle loss fáze B na ZL split.
* **Loss složky:** práh 0,05 platí pro nevážený čtverec relativní odchylky ((x−t)/t)² (↔ 22 %).
* **H2 pro GC:** v101 jsou glyfy už jednotlivé znaky; preset `extended` (Eva kompozity) by
  některé v101 znaky chybně slučoval (140 533 vs 140 940 symbolů), proto GC používá preset `None`.
  `autocopy_like` na GC substituuje mezi libovolnými glyfy seedu (prázdná `edit_table`), omezeno
  bigramovým modelem.
* **BPE:** `bpe_dependence_curve(lines, (0,32,64,128))` s výchozím `preset=None` (znaky), jak
  odpovídá referenční hodnotě 1,05 ve spec (ověřeno: 1,058).
* **Paradigmata:** afixové seznamy se učí na reálném textu (pro každý soubor/politiku zvlášť) a
  stejné seznamy se použijí pro generované vzorky a míchání.
* **Bootstrap reálného textu:** 200× převzorkování stránek s vracením (duplicitní stránky dostanou
  odlišné `folio`). Interval = x ± 1,96·SD bootstrapu (normální), protože n_types a hapax jsou
  při převzorkování s vracením vychýlené dolů (duplikované stránky), percentilový interval by byl
  posunutý.
* **Rozanova & Temerev tab. 7:** hodnoty nejsou v repu a spec nedovoluje stahovat; v `profile.csv`
  chybí (uvedeno v reportu).
* **Push:** session dostala repo jako lokální bundle (`clone: from /home/user/.seed.bundle`) bez `origin`;
  origin doplněn dodatečně (`https://github.com/HinataSoft/voytools`) a větev pushnuta na konci.

## Průběh

* Pilot (seed 0, výchozí konfigurace): share_near 0,59, n_types 11 069, H2 3,58. Rohové
  konfigurace mřížky: H2 vždy 3,5–3,65 → H2 bude pravděpodobně neslučitelné s ostatními cíli.
* Fáze A (4 paralelní procesy, ~1 h): min loss 0,964; H2 v celé mřížce 3,52–3,69; nejlepší
  konfigurace všechny na okraji mřížky (window 32, p_new 0,15).
* Diagnostika mimo mřížku (`diag`): bez editací H2 ≈ reál → přebytek H2 tvoří editace; p_new 0,5
  dá loss ~0,20 s jedinou netrefenou složkou H2. Doplněn `diagfull` (50 seedů, všechny metriky,
  všechny 4 sady), aby šlo oddělit vliv nekalibrace na necílové metriky: BPE žlab se na ZL srovná,
  paradigmata a MI/poziční metriky ne.
* Čekací smyčka `pgrep -f "run.py grid"` našla sama sebe a doběhla až na časový limit (bez vlivu
  na výpočet); příště čekat na PID.
