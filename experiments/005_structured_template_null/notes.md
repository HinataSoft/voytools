# 005 – deník session

## Balík
* `voytools/templates.py`: `ParadigmTemplate.learn(dfd, **params)`, `.sample(template_df, seed, pool_seed=None)`,
  `.with_params(**params)` (sdílí naučené statistiky, gramatiky a cache kmenů; pro mřížku), `.stem_pool(seed, n)`.
  `nulls.py` beze změny. Testy `tests/test_templates.py`.
* `morphology.decompose` / `decompose_table`: nový parametr `order='prefix_first'|'suffix_first'`
  (výchozí beze změny; `suffix_first` = zrcadlový greedy parser). Test `test_decompose_suffix_first_order`.

## Rozhodnutí a výklad spec
* **Umělé kmeny**: `WordGrammar.sample` (slot + trigram) naučená na reálných kmenech. Pořadí (= Zipfův rank)
  podle prvního výskytu: nejdřív procházka gramatikou naučenou na *tokenech* kmenů (pravděpodobné krátké řetězce
  dostanou vysoké ranky), po vyčerpání (2 000 neúspěchů za sebou) gramatika naučená na *typech* kmenů (dlouhý
  chvost). Krátké umělé kmeny nutně splývají s reálnými (`share_stems_real` v `runs.csv`).
* **Velikost paradigmatu k**: z reálných kmenů se stejným `floor(log2 count)`; count = *realizovaný* počet
  tokenů kmene ve vzorku (ne očekávaný – očekávaný systematicky podhodnocoval).
* **Téma**: θ_kmen ~ Dirichlet(α·K·base), base = podíl tokenů sekce; α = ∞ bez tématu.
* **Lokalita stránky**: slovník stránky = 20 tahů z rozdělení sekce stránky (`page_vocab=20`, pevné, mimo mřížku);
  s pravděpodobností λ se kmen tokenu bere odtud.
* **Další „ruční pravidla“** (mimo mřížku, jen v ablacích/adicích): `prefix_dependence='suffix'`
  (P(prefix|sufix) z reálu), `'stem_suffix'` (sada prefixů kmene vážená P(prefix|sufix)),
  `suffix_by='language'` (frekvence sufixů zvlášť pro A a B).
* **Společná náhoda**: `pool_seed = rep`, tj. konfigurace se stejným rep sdílejí inventář umělých kmenů.
* **Kalibrace** na ZL split i GC split zvlášť (každý soubor svou mřížku a fázi B), join jen kontrola nejlepší
  konfigurace (50 seedů).
* **Kolize jmen**: metrika 002 `n_stems` (počet kmenů po dekompozici) přepsala v řádcích fáze A a části fáze B
  parametr `n_stems`. Parametr je proto uložen jako `inventory_n_stems`; u starších řádků se dopočítá z
  `config_id` (`load_dir`). Metriky nejsou dotčené.
* **Ablace – výchozí „plný“ model**: cílové metriky téma α a lokalitu λ prakticky neurčují (loss je přes α, λ plochá).
  Plný model = parametry kmenů/paradigmatu/prefixu z fáze B (rank 1) + (α, λ) z mřížky fáze A (stejné ostatní
  parametry) s nejmenší relativní chybou na `jsd_stem_sections` a `stem_dispersion` (krok „stage 2“).
  Pokud plný model nemá vazbu prefixů, ablace `no_prefix_link` by byla totožná – místo ní běží `add_prefix_stem`.
* **Robustnost**: podvzorky 80 % stránek bez vracení; SD celého textu ≈ SD podvzorku · √(0,8/0,2) = 2× (delete-d
  jackknife). Bodový odhad = celý text. BPE varianta: 64 sloučení (`learn_bpe`, preset souboru), prefix = první
  jednotka tokenu, pokud je mezi 12 nejčastějšími počátečními jednotkami a zbude ≥ 1 jednotka; sufix = poslední
  jednotka mezi 14 nejčastějšími koncovými. Slotový null se pro každou variantu tvoří z její dekompozice;
  autocopy_536 je jeden vzorek na seed vyhodnocený všemi variantami.
* **A/B × písař**: skupina = jazyk + ruka (`A1`, `B2`, …). Interval: podvzorky stránek (jako výše); základní
  linie: permutace štítků mezi stránkami dvojice (200×) – JSD malých skupin je vychýlená nahoru, rozhoduje
  přebytek nad permutací. Sekce uvnitř písaře jen s ≥ 500 tokeny.
