# 003 – Autokopírování s editacemi vázanými gramatikou slova

## Otázka
Experiment 001 ukázal, že generátor „kopíruj a obměň“ trefí hapaxy, počet typů i podíl blízkých
sousedů, ale ne podmíněnou entropii glyfů (H2 3,4–3,7 proti 2,69): libovolné editace (bigramový
model) vyrábějí posloupnosti glyfů, které rukopis nezná. Timm & Schinner (2020) obměňují slova
podle tabulky povolených skupin glyfů, tedy v mezích slovní šablony, a reimplementace u Rozanové
a Temereva (2026) H2 2,6–2,7 dosáhla. Otázka: **stačí vázat editace gramatikou slova, aby
generátor trefil všechny cílové vlastnosti zároveň – a které necílové vlastnosti tím začne
reprodukovat?** To druhé je vlastní vědecký výsledek: cokoli, co kopírování + šablona vysvětlí,
přestává být indicií smyslu.

## Hypotézy a predikce
- **Vázané editace** (trigram, slot) sníží H2 na úroveň rukopisu při zachování hapaxů/typů;
  bigram (001) zůstane referencí.
- Necílové vlastnosti, které gramatika slova pravděpodobně „dožene“: BPE žlab při 64 sloučeních,
  podíl kmenů s více sufixy, počet sufixů na kmen, `cmi_excess`.
- Necílové vlastnosti, které bez pozičního a mezislovního mechanismu zůstanou nevysvětlené:
  okrajová MI, podíl pořadí tokenů, poziční JSD. (Pokud je generátor začne reprodukovat, je to
  překvapení a report to zdůrazní.)

## Změny balíku (větev `feat/autocopy-strict-edits`, pak merge do `exp/003-...`)
V `voytools.nulls.autocopy_generate` (a `autocopy_like`) přidat parametry s výchozími hodnotami,
které **nemění stávající chování**:
- `edit_model='bigram'|'trigram'|'slot'`:
  - `trigram`: substituce/inzerce/delece se přijme jen, když všechny trigramy glyfů výsledného
    tokenu (s okrajovými značkami `^`, `$`) jsou doloženy v reálném slovníku; po `max_tries`
    neúspěších se token zkopíruje beze změny (a počítá se do `p_copy_exact` – uvést v reportu
    skutečný podíl nezměněných kopií).
  - `slot`: z reálných tokenů se naučí slotová gramatika – pořadí slotů podle průměrné relativní
    pozice glyfu v tokenu (alternativně pevné pořadí podle Zattery 2022, pokud je v balíku
    snadno zapsatelné); token je gramatický, když jeho glyfy jdou v neklesajícím pořadí slotů,
    žádný slot se neopakuje mimo doložené výjimky (`ee`, `ii`) a všechny bigramy jsou doloženy.
    Editace se přijme jen do gramatického tokenu.
- `new_token_source='seed'|'grammar'`: `seed` = dosavadní chování (losování z `seed_tokens`);
  `grammar` = nové slovo se vygeneruje náhodnou procházkou po doložených trigramech od `^` po
  `$` s délkou vzorkovanou z reálného rozdělení délek a s podmínkou gramatičnosti (`slot`).
- Funkce `grammaticality(df_or_tokens, model)` vracející podíl tokenů, které gramatika připouští
  (na reálném textu má být ≈ 1,0 pro model naučený na něm samém; uvést pro ZL i GC).
- Testy: výchozí volání dává stejný výstup jako dřív (regresní test na pevném seedu);
  `trigram` nevyrobí nedoložený trigram; `grammar` sampler vrací jen gramatické tokeny.

## Data a filtry
Jako 001: `data/ZL3b-n.txt` a `data/GC2a-n.txt`, `generic_type="P"`, `clean=True`, politiky
`split` (hlavní) a `join`. Seed tokeny = prvních 100 tokenů; gramatika a trigramy se učí
z celého reálného textu (uvést).

## Cílové metriky a ztráta
Shodné s 001 (`share_near`, `hapax_share`, `H2` extended/scope all – u GC preset `None`,
`n_types`, `mean_len`; váhy 2,2,2,1,1; loss = Σ w·((x−t)/t)²). **Omezení:** `p_new` při
`new_token_source='seed'` nejvýš 0,3 – vyšší hodnoty jen převzorkovávají reálný slovník a nejsou
generativním vysvětlením. Bod `p_new=0,5` z 001 slouží jen jako horní mez (uvést v tabulce).

## Mřížka
`edit_model ∈ {trigram, slot}` × `new_token_source ∈ {seed, grammar}` × `p_new ∈ {0.05, 0.15, 0.3}`
× `window ∈ {16, 32, 64}` × `p_copy_exact ∈ {0.35, 0.5, 0.65}` × `p_two_edits ∈ {0.25, 0.5}`
× `recency ∈ {0, 0.3}` → 432 konfigurací; `edit_type_probs` pevně (0.6, 0.2, 0.2).
Fáze A: 3 seedy na konfiguraci; fáze B: 10 nejlepších × 50 seedů (ZL/GC × split/join) + míchání.
`seed = config_id*100 + rep` jako v 001, sloupce s parametry v každém řádku.
Bigramovou referenci nepočítat znovu – převzít z `results/001_autocopy_calibration/`.

## Necílové metriky
Jako 001: `edge_glyph_mi`, `adjacent_token_mi` (excess share), `bpe_dependence_curve` (`gap`
při 64), `positional_profile` JSD `line_first` vs `mid`, `paradigm_metrics` (`cmi_excess`,
`frac_stems_multi_suffix`, `mean_suffixes_per_stem`). Navíc `grammaticality` reálného textu
i vzorků (kontrola, že `slot`/`trigram` opravdu drží gramatiku) a skutečný podíl nezměněných kopií.

## Výstupy (`results/003_autocopy_strict/`)
`runs.csv`, `meta.json`, `grid_summary.csv`, `best_configs.csv` (sloupce `edit_model`,
`new_token_source`, všechny parametry, `loss`, složky ztráty, `provisional` 0/1),
`profile.csv` + `profile.png` (zdroje: reál, bigram z 001, nejlepší trigram, nejlepší slot,
míchání v řádcích, referenční hodnoty z `docs/reference_values.md`), `report.md`, `notes.md`.

## Kritéria přijetí
1. Konfigurace s loss < 0,15 a každou složkou < 0,05 **při `p_new ≤ 0,3`**; jinak vysvětlit,
   která cílová vlastnost je s ostatními neslučitelná i při vázaných editacích.
2. `profile.csv` kompletní pro všechny zdroje a 4 sady (ZL/GC × split/join); report má tři
   seznamy necílových metrik: (a) reprodukované vázaným generátorem (|z| < 3 v ZL i GC),
   (b) nereprodukované v obou, (c) nekonzistentní mezi ZL a GC.
3. `best_configs.csv` první řádek = doporučený nulový model pro další experimenty; pokud
   kritérium 1 splněno, `provisional=0`, jinak `provisional=1` s poznámkou.
4. Regresní test výchozího chování prochází; `python -m pytest -q` prochází; Python 3.10 /
   pandas 2.3 kompatibilita zachována (žádné scipy).
5. Podezřelé hodnoty: `grammaticality` reálného textu < 0,95 (pak je gramatika špatně naučená);
   `share_near` vzorku > 0,9 (mimo `window` 16 to signalizuje chybu); H2 < 2,0.

## Rozpočet
Vázané editace jsou pomalejší (opakované pokusy): odhad 4–8 s na vzorek s plnými metrikami.
Fáze A: ve fázi A počítat jen cílové metriky + edge/token MI (≈ 2–3 s) → 432 × 3 × 3 s ≈ 1 h na
4 procesech ≈ 20 min. Fáze B ≈ 500 × 8 s ≈ 1,2 h na 4 procesech ≈ 20 min. Celkem se změnami
balíku 2–3 h v jedné session; dělit mezi session jen pokud fáze A překročí hodinu.
