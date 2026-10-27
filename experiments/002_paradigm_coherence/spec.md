# 002 – Paradigma, nebo šablona? Koherence kmenů a afixů

## Otázka
Chovají se „kmeny“ voynichských tokenů jako lexikální jednotky s paradigmatem (tentýž kmen se vrací
s různými afixy, kombinace prefix–sufix je při daném kmeni omezená, inventář afixů je napříč
sekcemi stálý, zatímco inventář kmenů se mění), nebo jde jen o slotovou šablonu, kde se afixy
a kmeny kombinují nezávisle?

## Hypotézy a predikce
- **Paradigma (jazyk):** `cmi_excess` (I(prefix; sufix | kmen) nad náhodu) > 0 a výrazně nad
  nulovými modely; vysoký podíl kmenů s více sufixy (`frac_stems_multi_suffix`) i po kontrole
  frekvence kmene; Jaccard inventářů sufixů mezi sekcemi/jazyky A–B **vyšší** než Jaccard kmenů;
  produktivita sufixů (počet různých sufixů vs. frekvence kmene) roste pomaleji než u nezávislé
  šablony (afixy nejsou u každého kmene libovolné).
- **Šablona / generátor:** slotově nezávislý vzorek (viz nulové modely) dá `cmi_excess ≈ 0`
  a produktivitu odpovídající náhodnému párování; autokopírovací generátor dá vysokou
  „paradigmatičnost“ jen jako vedlejší efekt lokálního kopírování – odlišit tím, že u něj
  varianty téhož kmene leží blízko sebe (střední vzdálenost mezi výskyty kmene v tokenech),
  zatímco u jazyka jsou rozprostřené po celém textu.
- **Artefakt:** efekt zmizí v GC (v101) nebo závisí na volbě afixových seznamů (naučené vs.
  referenční) nebo politice mezer.

## Data a filtry
- `data/ZL3b-n.txt` a `data/GC2a-n.txt`; `generic_type="P"`, `clean=True`; politiky mezer
  `split`, `join`, `strict` (u `strict` tokeny s `,` uvnitř vyřadit a uvést počet).
- Dekompozice: (a) `decompose_table(df)` s naučenými afixy (výchozí parametry), (b)
  `decompose_table(df, prefixes="reference", suffixes="reference")`; u GC jen (a).
- Podskupiny pro inventáře: `illustration` (H, B, S, P, …) a `language` (A, B).

## Metriky
Z `paradigm_metrics(dfd, min_stem_count=5)`: `n_stems_k`, `mean_suffixes_per_stem`,
`frac_stems_multi_suffix`, `H_suffix_given_stem`, `mi_prefix_suffix`, `cmi_excess`.
Navíc (implementovat v `run.py`, případně přidat do `voytools.morphology` s testem):
- `jaccard_suffix_sections` vs. `jaccard_stem_sections`: průměrný Jaccard inventářů (sufixů /
  kmenů s ≥ 3 výskyty) přes všechny dvojice sekcí; totéž pro A vs. B;
- `productivity_slope`: sklon log(počet různých sufixů kmene) ~ log(frekvence kmene) přes kmeny
  s ≥ 5 výskyty (OLS, numpy);
- `stem_dispersion`: medián přes kmeny (≥ 5 výskytů) mediánové vzdálenosti (v tokenech) mezi po
  sobě jdoucími výskyty téhož kmene, dělený očekáváním při rovnoměrném rozložení (≈ N/k).

## Nulové modely (stejný kód, stejné filtry)
1. **Slotová nezávislost:** pro každý token vylosovat prefix, kmen a sufix nezávisle z reálných
   marginálních rozdělení (zachovat počty tokenů a délku řádků; vzorky s pravděpodobností
   úměrnou frekvenci). 50 seedů.
2. **Stem-conditional shuffle:** sufixy permutovat mezi tokeny téhož kmene (to je
   `cmi_shuffled` uvnitř `paradigm_metrics`); uvést explicitně.
3. **Autokopírování:** `autocopy_like(real, **best)` s parametry z
   `results/001_autocopy_calibration/best_configs.csv` (první řádek); pokud 001 ještě neběželo,
   použít výchozí parametry a v reportu to označit. 50 seedů.
4. **Míchání v řádcích** (`shuffle_within_lines`): 5 seedů – jen pro doložení, že metriky na
   úrovni typů jsou invariantní a `stem_dispersion` ne.

## Nejistota
Reálný text: bootstrap po stránkách, 200 opakování (jeden běh = jeden bootstrap vzorek, seed =
číslo opakování). Nulové modely: 50 seedů každý. V `runs.csv` sloupce `source`, `file`,
`space_policy`, `decomposition`, pak metriky; `seed` kóduje opakování.

## Výstupy (`results/002_paradigm_coherence/`)
- `runs.csv`, `meta.json`, `summary.csv` (metrika × zdroj × soubor × politika: průměr, 2,5 %, 97,5 %),
- `effects.csv` – rozdíl reál − null pro každou metriku a null, s intervaly a z-skóre,
- `inventories.csv` – Jaccard sufixů a kmenů pro každou dvojici sekcí (reál),
- `paradigm.png` – reál vs. nulové modely (cmi_excess, frac_multi_suffix, jaccard, dispersion),
- `report.md`.

## Kritéria přijetí
1. Všech 8 kombinací (2 soubory × 3 politiky, + 2 dekompozice u ZL) má kompletní tabulku.
2. Report odpoví na tři otázky zvlášť: (i) je `cmi_excess` reálu nad slotovou nezávislostí i nad
   autokopírováním? (ii) jsou afixy mezi sekcemi stálejší než kmeny? (iii) jsou varianty kmene
   rozprostřené (jazyk), nebo shluknuté (kopírování)? U každé: ZL ano/ne, GC ano/ne.
3. Podezřelé hodnoty: `cmi_excess` > 1 bit, Jaccard = 1 nebo 0, `stem_dispersion` ≫ 10.
4. Nové funkce v balíku mají docstring a test; `python -m pytest -q` prochází.

## Rozpočet
Jedna dekompozice + metriky ≈ 1–3 s. Reál 200 bootstrapů × 8 kombinací ≈ 1 h; nulové modely
(50+50+5) × 8 ≈ 1 h. Rozdělit na 2–3 session (`seed_range`), výsledky do `partK/`, sloučit.
