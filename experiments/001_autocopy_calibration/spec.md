# 001 – Kalibrace autokopírovacího generátoru a co zůstane nevysvětleno

## Otázka
Jak dobře dokáže jednoduchý generátor „kopíruj a obměň“ (`voytools.nulls.autocopy_like`, ve
stylu Timm & Schinner 2020) napodobit ZL3b v těch vlastnostech, pro které byl navržen – a které
vlastnosti textu **zůstanou nevysvětlené** i po nejlepší kalibraci? Výsledek je zároveň nulový
model pro všechny další experimenty (parametry z `results/001_autocopy_calibration/best_configs.csv`).

## Hypotézy a predikce
- Generátor lze nakalibrovat tak, aby trefil **cílové vlastnosti**: podíl blízkých sousedů
  (`share_near`, reálně 0,173, poměr k míchání 1,18), podíl hapaxů (0,676), počet typů (~6 800),
  průměrnou délku tokenu (4,99 znaků) a podmíněnou entropii glyfů H2 (2,69 bitu, preset
  `extended`, scope `all`). Podle Rozanova & Temerev (2026) to jejich implementace dokázala.
- **Necílové vlastnosti**, u kterých očekáváme, že generátor selže (tj. kde by se měl projevit
  obsah/gramatika): vazba okrajových glyfů (`edge_glyph_mi` excess 0,185 bitu; generátor ≈ 0),
  podíl pořadí tokenů (`adjacent_token_mi` excess share 0,9 %; generátor ≈ 0), hloubka BPE žlabu
  při 64 sloučeních (reálně 1,05 bitu, generátor výrazně míň), poziční profil (JSD prvního glyfu
  řádku vs. uprostřed řádku), a paradigmatické metriky (`paradigm_metrics`: `cmi_excess`,
  `frac_stems_multi_suffix`).
- Pokud generátor po kalibraci trefí i necílové vlastnosti, je to důležitý negativní výsledek:
  tyto metriky pak nerozlišují smysl od kopírování a v dalších experimentech je nepoužijeme
  jako indicii.

## Data a filtry
- `data/ZL3b-n.txt`, tokeny `filter_tokens(tables.tokens, generic_type="P", clean=True)`,
  politika mezer `split` (hlavní) a `join` (kontrola).
- Kontrola abecedy: totéž na `data/GC2a-n.txt` (v101; `autocopy_like` tam používá substituce
  mezi libovolnými glyfy seedu – uvést v reportu).
- Generátor dostane vždy stejný layout jako reálný text (`autocopy_like(real, ...)`), seed
  tokeny = prvních 100 tokenů reálného textu (výchozí chování).

## Metriky (všechny z `voytools`)
Cílové: `neighbor_similarity(window=5)` → `share_near`, `ratio`; `vocabulary_summary` →
`n_types`, `hapax_share`, `mean_len`; `glyph_entropy_report(preset="extended", scope="all")` → `H2`.
Necílové: `edge_glyph_mi(shuffles=20)` → `excess_bits`; `adjacent_token_mi(cap=2000, shuffles=20)`
→ `excess_share`; `bpe_dependence_curve(lines, merge_checkpoints=(0,32,64,128))` → `gap` při 64;
`positional_profile(which="first")` → JSD `line_first` vs `mid`; `paradigm_metrics(decompose_table(df))`
→ `cmi_excess`, `frac_stems_multi_suffix`, `mean_suffixes_per_stem`.
Každou metriku spočítat identickým kódem pro reálný text (jednou, seed 0) i pro každý vzorek.

## Parametrická mřížka
`window ∈ {4, 8, 16, 32}`, `p_copy_exact ∈ {0.1, 0.2, 0.35, 0.5}`, `p_new ∈ {0.02, 0.05, 0.15}`,
`p_two_edits ∈ {0.1, 0.25, 0.5}`, `recency ∈ {0, 0.3}`, `edit_type_probs ∈ {(0.6,0.2,0.2), (0.4,0.3,0.3)}`
→ 576 konfigurací. Fáze A: 3 seedy na konfiguraci (1 728 vzorků). Fáze B: 10 nejlepších
konfigurací × 50 seedů. Identifikace běhu v `runs.csv`: `seed = config_id * 100 + rep`
(`config_id` 0..575 v pořadí `itertools.product` v uvedeném pořadí os); sloupce s parametry
konfigurace do každého řádku.

## Ztrátová funkce
`loss = Σ_f w_f · ((x_f − t_f) / t_f)²` přes cílové vlastnosti f, váhy: `share_near` 2, `hapax_share` 2,
`H2` 2, `n_types` 1, `mean_len` 1; `x_f` = průměr přes seedy konfigurace, `t_f` = reálná hodnota.
Vyhodnocovat loss i zvlášť po složkách (žádná složka > 0,05, tj. relativní odchylka > 22 %, je
„netrefeno“).

## Nejistota
Fáze B dává pro každou metriku průměr a 2,5/97,5% kvantil přes 50 seedů; pro reálný text
bootstrap po stránkách (200 opakování) u metrik, kde to dává smysl (hapax, n_types, H2, share_near).
Rozdíl reál − generátor vyjádřit v jednotkách SD generátoru (z-skóre) i absolutně.

## Výstupy (`results/001_autocopy_calibration/`)
- `runs.csv`, `meta.json` (fáze A i B; fázi B lze dát do podsložky `phaseB/` a sloučit `merge_runs`),
- `grid_summary.csv` – konfigurace × průměrné metriky × loss, seřazeno,
- `best_configs.csv` – 10 nejlepších s parametry (tento soubor čtou další experimenty),
- `profile.csv` – pro reálný text, nejlepší konfiguraci, čisté míchání v řádcích a (pokud
  existuje) referenční hodnoty z Rozanova & Temerev tab. 7: všechny metriky s intervaly,
- `profile.png` – cílové vs. necílové metriky (reál vs. generátor, s intervaly),
- `report.md`.

## Kritéria přijetí
1. Existuje konfigurace s loss < 0,15 a žádnou složkou > 0,05; jinak report vysvětlí, která cílová
   vlastnost je neslučitelná s ostatními (to je samo o sobě výsledek).
2. `profile.csv` má řádek pro každou metriku a každý zdroj; report explicitně vyjmenuje
   necílové metriky, u nichž je reálný text > 3 SD od nejlepšího generátoru (v ZL i GC), a ty,
   kde generátor reálný text reprodukuje.
3. Testy balíku procházejí; `run.py` jde spustit znovu a dopočítá jen chybějící seedy.
4. Podezřelé hodnoty (bug): `share_near` generátoru < 0,05 nebo > 0,9; `H2` < 1,5 nebo > 4;
   záporný `hapax_share`.

## Rozpočet
Jeden vzorek + všechny metriky ≈ 3–6 s (BPE a paradigmata jsou nejdražší; v fázi A je lze
vynechat a počítat jen cílové + edge/token MI). Fáze A ≈ 1 728 × 3 s ≈ 1,5 h na jednu VM →
rozdělit na 4 session (`--seeds 0:57600/4:k` při kódování `config_id*100+rep` – run.py musí
rozsah interpretovat přes config_id, viz `_template/run.py: parse_seeds`). Fáze B ≈ 500 × 6 s ≈ 50 min.
