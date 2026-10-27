# 002 – deník session

## Rozhodnutí a výklad spec

* **Kombinace:** spec píše „8 kombinací (2 soubory × 3 politiky, + 2 dekompozice u ZL)“, to je
  ve skutečnosti 9: ZL × {split, join, strict} × {learned, reference} + GC × {split, join, strict} ×
  learned. Spočítáno všech 9.
* **Autokopírování:** podle zadání vlastníka dva nulové modely z prvních dvou řádků
  `results/001_autocopy_calibration/best_configs.csv` – `autocopy_573` (rank 1: window 32,
  p_copy_exact 0,5, p_new 0,15, p_two_edits 0,5, recency 0, edit 0,4/0,3/0,3) a `autocopy_536`
  (rank 2: totéž s p_copy_exact 0,35 a edit 0,6/0,2/0,2; nejlepší na ZL join a GC v 001).
  Volání stejné jako v 001: `autocopy_like(real, allowed_from=<všechny reálné tokeny>, seed, **cfg)`.
  Na GC (preset `None`) je `edit_table` prázdná, substituce jen podle bigramového modelu.
* **Afixové seznamy** se učí jednou na reálném textu dané kombinace (nebo referenční seznamy) a
  stejné seznamy se použijí pro bootstrap i všechny nulové vzorky (fér srovnání, viz docstring
  `decompose_table`).
* **strict:** tokeny obsahující `,` jsou vyřazeny (ZL 2 104, GC 1 921 tokenů; uloženo v
  `meta.json` každé části jako `n_dropped_comma`).
* **Slotová nezávislost** (`voytools.morphology.slot_independent_sample`): pro každý token se
  prefix, kmen a sufix losují nezávisle s vracením z reálných tokenových marginál, slepí se a
  vzorek projde **stejnou** `decompose_table` se stejnými seznamy jako reál (greedy re-parse může
  slepeninu rozdělit jinak než los). Rozložení řádků, stránek a metadata zůstávají.
  `nulls.py` nebyl měněn (pracuje na něm jiná session), proto je sampler v `morphology.py`.
* **Stem-conditional shuffle** = `cmi_shuffled` uvnitř `paradigm_metrics` (20 permutací prefixů
  uvnitř kmene); v `effects.csv` jako null `stem_cond_shuffle`, rozdíl = `cmi_excess`.
* **Inventáře:** inventář skupiny = hodnoty (sufix včetně prázdného / kmen) s ≥ 3 výskyty v té
  skupině. Sekce = všechny hodnoty `illustration` v P textu (H, S, B, P, T, C, A), průměr přes všech
  21 dvojic (v bootstrapu někdy méně, když vypadne malá sekce A); jazyk = A vs B (tokeny bez jazyka
  vynechány). Navíc `jsd_*` (Jensen–Shannon frekvenčních rozdělení, bez prahu), protože sufixů je
  jen 15 a Jaccard sufixů je proto téměř vždy ≈ 1 (strukturálně, ne empiricky).
* **Kontrola frekvence u `frac_stems_multi_suffix`:** navíc `frac_multi_suffix_5to9` – podíl
  kmenů s > 1 sufixem jen mezi kmeny s 5–9 tokeny.
* **`productivity_slope`:** OLS ln(# různých sufixů) ~ ln(frekvence), kmeny ≥ 5; uložen i
  průsečík a R².
* **`stem_dispersion`:** medián přes kmeny (≥ 5) mediánové mezery mezi po sobě jdoucími výskyty
  / (N/k). Pozor: při náhodných nezávislých pozicích je mezera přibližně geometrická a poměr
  ≈ ln 2 ≈ 0,69, ne 1; rovnoměrné rozmístění dá 1.
* **Intervaly v `effects.csv`:** rozdíl = bodový odhad reálu − průměr nulového modelu;
  95% interval = rozdíl ± 1,96·√(SD_boot² + SD_null²); `z_null` = rozdíl / SD_null. Bodový odhad
  (ne průměr bootstrapu), protože typové metriky (počty kmenů atd.) jsou při převzorkování stránek
  s vracením vychýlené (duplicitní stránky), stejně jako v 001.
* **Seedy:** `seed` = číslo opakování uvnitř zdroje (real 0; real_boot 0–199; slot_indep,
  autocopy_* 0–49; shuffle_lines 0–4). Tentýž seed jde do vzorku i do permutací v metrikách.

## Průběh

* Pilot ZL strict learned (seedy 0–1): reál cmi_excess 0,078, slot 0,022, autocopy 0,18–0,24;
  stem_dispersion reál 0,46, slot 0,68, autocopy 0,006–0,03.
* Běh: 3 bloky po 3 kombinacích paralelně (`run.py run --combo ...`), commit + push po bloku.
