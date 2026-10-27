# 005 – Null „strukturovaný umělý jazyk“ a robustnost experimentu 002

## Otázka
Experiment 002 našel, že kmeny mají užší paradigma než nezávislá šablona a širší než kopírování,
silnou nepodmíněnou vazbu prefix–sufix, rozptýlené výskyty a tematicky závislé rozdělení po sekcích.
Otázka: **odliší tyto metriky reálný text i od systému, který paradigmata a tematickou volbu kmenů
má *z konstrukce*, ale žádný význam nenese** (systematická glosolálie, vymyšlený jazyk s tabulkami)?
Pokud ne, jsou nálezy 002 důkazem *systému*, ne *smyslu*, a směr 1 je vyčerpaný. Druhá část:
ověřit, že 002 nezávisí na volbě afixových seznamů a parseru.

## Hypotézy a predikce
- Šablonový generátor s vestavěnými paradigmaty reprodukuje `mean_suffixes_per_stem`,
  `frac_stems_multi_suffix` a (s tematickým parametrem) `jsd_stem_sections`.
- Nereprodukuje bez dalšího: nepodmíněnou `mi_prefix_suffix` (0,18–0,21 bitu), `cmi_excess`
  (0,09), disperzi kmenů 0,46 (šablona bez tématu dává 0,69), posun sufixových frekvencí A/B.
  Každou z těchto čtyř lze do generátoru přidat parametrem – a tím se zjistí, **kolik nezávislých
  „ručních“ pravidel** by umělý systém potřeboval, aby reál napodobil. To je výstup pro syntézu.
- Robustnost: závěry 002 (i)–(iii) platí při 8 i 20 sufixech, délce ≤ 3, BPE jednotkách jako afixech
  a parseru s preferencí sufixu; jinak uvést, které ne.

## Změny balíku (jediná větev `exp/005-structured-template`)
- Nový modul `voytools/templates.py` (nesahat na `nulls.py`, paralelně se mění v 004):
  `ParadigmTemplate.learn(dfd, **params)` z dekomponované tabulky reálu a
  `ParadigmTemplate.sample(template_df, seed) -> DataFrame` (přes `nulls.as_table`, layout reálu).
  Mechanismus: inventář prefixů a sufixů s frekvencemi z reálu; **umělý** inventář kmenů o velikosti
  `n_stems` s Zipfovým rozdělením (exponent `zipf_s`) – kmeny jsou náhodné gramatické řetězce
  (použít `wordgrammar.WordGrammar.sample` nebo trigramovou procházku), ne reálné kmeny; každému
  kmeni přidělena sada sufixů velikosti k ~ rozdělení naučené z reálu (podle frekvence kmene) –
  parametr `paradigm='fixed'|'free'` (free = každý token losuje sufix nezávisle = slotová
  nezávislost); prefix losován globálně, nebo `prefix_dependence='stem'|'suffix'|'none'`
  (zapne vazbu prefix–sufix); `topic_concentration` α: rozdělení kmenů po sekcích z Dirichletovy
  směsi (α→∞ = bez tématu, malé α = silné téma); `page_locality` λ: s pravděpodobností λ se kmen
  losuje z rozdělení stránky (vytváří disperzi < 0,69 bez kopírování). Bez kopírování tokenů.
- Testy: determinismus, layout, gramatičnost kmenů, `paradigm='free'` dává `cmi_excess` ≈ slotový null z 002.

## Data a filtry
Jako 002: ZL3b a GC2a, `P`, `clean`, politiky `split` (hlavní) a `join`; dekompozice naučená
(`decompose_table` výchozí) – stejné seznamy pro reál i vzorky.

## Kroky
1. **Kalibrace šablony** na cílové metriky: `n_types`, `hapax_share`, `mean_len`, `n_stems_k`,
   `mean_suffixes_per_stem`, `frac_stems_multi_suffix` (váhy 1,1,1,2,2,2; loss jako v 001).
   Mřížka: `n_stems ∈ {2k, 4k, 8k}`, `zipf_s ∈ {0.9, 1.1, 1.3}`, `paradigm ∈ {fixed}`,
   `prefix_dependence ∈ {none, stem}`, `topic_concentration ∈ {∞, 5, 1, 0.2}`,
   `page_locality ∈ {0, 0.3, 0.6}` → 216 konfigurací × 3 seedy; fáze B: 5 nejlepších × 50 seedů,
   plné metriky 002 (`metrics` z `experiments/002_paradigm_coherence/run.py` – importovat nebo
   zkopírovat beze změn).
2. **Ablace** na nejlepší konfiguraci (50 seedů každá): bez tématu, bez lokality, bez vazby
   prefixů, `paradigm='free'`. Pro každou metriku 002 tabulka „kterou ingredienci potřebuje“.
3. **Robustnost 002** (reál, 200 subsamplů po stránkách *bez* vracení – 80 % stránek, místo
   bootstrapu s vracením, který typové metriky vychyluje; viz 002 report): afixové seznamy
   `n_suffixes ∈ {8, 14, 20}`, `max_len ∈ {3, 4}`, afixy = BPE jednotky z `learn_bpe(merges=64)`
   (prefix/sufix = první/poslední jednotka), parser s preferencí sufixu (`decompose` zrcadlově:
   nejdřív nejdelší sufix, pak prefix – přidat parametr `order='prefix_first'|'suffix_first'`
   do `morphology.decompose`/`decompose_table`, výchozí beze změny, test). Reportovat (i)–(iii)
   z 002 pro každou variantu proti slotovému nullu a `autocopy_536`.
4. **A/B versus písař:** `jsd_suffix_AB` rozložit: mezi písaři uvnitř jazyka B (ruce 2/3/4/5 podle
   `hand`), mezi sekcemi uvnitř písaře, mezi jazyky uvnitř písaře (pokud existuje). Porovnat s
   šablonou s tématem i bez něj.

## Výstupy (`results/005_structured_template_null/`)
`runs.csv`, `meta.json`, `grid_summary.csv`, `best_configs.csv`, `ablation.csv`, `robustness.csv`,
`ab_by_hand.csv`, `template_vs_real.png`, `report.md`, `notes.md`.

## Kritéria přijetí
1. Kalibrace: existuje konfigurace s loss < 0,15 (šablona má dost volnosti; pokud ne, uvést proč).
2. Tabulka „metrika × ingredience“: pro každou metriku 002 (6 hlavních + disperze + JSD) uvést,
   zda ji nejlepší šablona reprodukuje (|z| < 3, ZL i GC) a která ablace ji rozbije.
3. Robustnost: závěry (i)–(iii) z 002 označeny jako „stabilní“ / „závislé na variantě X“.
4. `ab_by_hand.csv` s intervaly; věta, zda posun sufixů A/B přežije kontrolu písaře.
5. Testy procházejí; Python 3.10 / pandas 2.3; bez scipy. Podezřelé: Jaccard = 1 (uvést jako
   strukturální jako v 002), `cmi_excess` > 1, disperze > 10.

## Rozpočet
Vzorek šablony ≈ 1–2 s, metriky 002 ≈ 2 s → mřížka 648 vzorků ≈ 40 min na 4 procesech, fáze B
a ablace ≈ 450 vzorků ≈ 30 min, robustnost 200 × ~10 variant × 2 soubory ≈ 1 h. Celkem 2–3 h,
jedna session, jedna větev, PR na konci.
