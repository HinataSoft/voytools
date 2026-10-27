# 004 – Port generátoru Timm & Schinner (self-citation) jako kanonický kopírovací null

## Otázka
Náš generátor `autocopy_generate` netrefuje H2 ani hapaxy ani s editacemi vázanými gramatikou
(003: loss 0,27, H2 3,3 proti 2,69). Rozanová a Temerev (2026) mají v repu `voynich-units` (MIT)
Python port referenční Java implementace Timma a Schinnera (2020), který na stejném textu dává
h2 2,6–2,7, hapax 59–60 %, D64 0,65–0,86, podíl pořadí tokenů 0,5–0,6 %, edge MI 0,006–0,009
(`docs/reference_values.md`). Otázka: **reprodukuje tento mechanismus naše cílové vlastnosti, a co
z necílových (paradigmata, BPE žlab, disperze kmenů) zůstane nevysvětleno i jemu?** Výsledek
nahradí provizorní null z 001/003.

## Hypotézy a predikce
- Port reprodukuje publikované hodnoty (jinak je chyba v portu, ne v mechanismu).
- Na našich cílových metrikách trefí H2 a hapax (to naše editace neuměly); `share_near` a `n_types`
  mohou zůstat mimo.
- Necílové: edge MI, pořadí tokenů, poziční JSD – generátor má pravidla pro začátky řádků
  a odstavců (pravidlo VIII, gallows) a volbu zdroje „stejná pozice v řádku“, takže poziční JSD
  může částečně reprodukovat; edge MI nejspíš ne. Paradigmatické metriky a disperze kmenů:
  predikce z 002/003 je, že zůstanou nereprodukované (disperze kmenů ≪ reál).

## Změny balíku (jediná větev `exp/004-selfcitation-port`)
- Nový modul `voytools/selfcitation.py`: port třídy `SelfCitation` a datových tabulek
  (`LIGATURES`, podobnostní tabulky glyfů, `can_before`/`can_after`, `G`) ze souboru
  `analysis/reproduce_selfcitation_control.py` repa <https://github.com/lrozanova/voynich-units>
  (naklonovat `git clone --depth 1` do dočasného adresáře; GitHub je v cloudu dostupný).
  V hlavičce modulu uvést původ a licenci (MIT, Rozanova & Temerev 2026; mechanismus Timm &
  Schinner 2020, Java reference github.com/TorstenTimm/SelfCitationTextgenerator). Převzít
  generátor beze změny logiky; vypustit jen části driveru (entropie, BPE, útok), které už máme.
- API: `selfcitation_generate(template_lines, *, seed, p_copy=0.0, replace_weights=(30,50,20),
  p_add_remove=20, p_combine_split=30, p_same_position=28, p_reuse_last=10, p_line_final=0.35,
  p_suggest=40, max_repeat=3) -> list[list[str]]`, kde `template_lines` je seznam řádků s počtem
  tokenů a příznakem `is_par_first_line` (z reálné tabulky); a `selfcitation_like(template_df,
  **params) -> DataFrame` přes `nulls.as_table` (stejná metadata a layout jako reál, jako
  `autocopy_like`). Výstupní tokeny jsou v základní Eva; ligatury generátoru rozepsat zpět na Eva
  znaky, aby `collapse_eva` fungoval stejně jako na reálu.
- Testy: determinismus podle seedu; výstup má přesně layout šablony; na 2 000 tokenech jsou
  všechny glyfy z Eva abecedy; regresní test, že `autocopy_generate` je nedotčen.
- Nic v `nulls.py` neměnit kromě případného re-exportu; `morphology.py`, `stats.py` neměnit.

## Data a filtry
`data/ZL3b-n.txt`, `generic_type="P"`, `clean=True`, politika `split` (hlavní) a `join`.
Generátor pracuje jen v Eva (jeho tabulky jsou Eva-specifické), takže vzorek pro GC2a nelze vyrobit;
v profilu uvést reálné hodnoty GC jen pro srovnání a robustnost zde nahradit dvojicí split/join.

## Kroky
1. **Validace portu (bez kalibrace):** výchozí parametry R&T, 5 seedů, šablona = ZL split.
   Spočítat náš profil (`glyph_entropy_report` extended/scope all, `vocabulary_summary`,
   `bpe_dependence_curve` gap@64, `adjacent_token_mi` excess share, `edge_glyph_mi`,
   `neighbor_similarity`). Zároveň tytéž metriky na jejich třech hotových vzorcích
   `data/controls/selfcitation/selfcite_seed{1,2,3}.txt` z klonu (prostý text, řádek = řádek; načíst
   jako tokeny a obalit `as_table` se syntetickými metadaty). Rozdíl „jejich vzorky vs. náš port“
   odhalí chyby portu; rozdíl „jejich vzorky vs. jejich tabulka 7“ odhalí rozdíly metrik.
2. **Kalibrace A (jejich ztráta):** mřížka `p_copy ∈ {0, 0.1, 0.2, 0.3}` × `replace_weights ∈
   {(30,50,20), (50,35,15), (20,50,30), (60,30,10)}`, 5 seedů, ztráta R&T (délka tokenu /0,25,
   typy /500, hapax /0,03). **Kalibrace B (naše ztráta z 001):** stejná mřížka + `p_same_position ∈
   {14, 28, 42}`, 3 seedy, loss = Σ w·((x−t)/t)² se složkami share_near, hapax, H2, n_types,
   mean_len (váhy 2,2,2,1,1). Fáze B: nejlepší bod každé kalibrace × 50 seedů, všechny metriky
   včetně `paradigm_metrics`, `stem_dispersion`, `inventory_jaccard`, `productivity_slope`
   (funkce z `voytools.morphology`, dekompozice `decompose_table` s naučenými afixy z reálu –
   stejné seznamy pro reál i vzorky).
3. **Velikost slovníku jako konfound:** pro paradigmatické metriky navíc (a) reál podvzorkovaný
   po stránkách tak, aby `n_stems_k` odpovídalo generátoru (20 opakování), (b) generátor na
   šabloně prodloužené tak, aby `n_stems_k` odpovídalo reálu. Uvést, jak moc se rozdíly z 002
   změní.
4. **Profil:** `profile.csv` ve formátu 003 (zdroje: reál, `ts2020_default`, `ts2020_calA`,
   `ts2020_calB`, `trigram_86` z 003, `bigram_573` z 001, míchání, `ref_RT2026_*`), z-skóre,
   tři seznamy (a)/(b)/(c) jako v 003 – u GC jen poznámka, že generátor je Eva-only; místo „ZL i GC“
   zde platí „split i join“.

## Výstupy (`results/004_selfcitation_port/`)
`validation.csv` (port vs. jejich vzorky vs. tabulka 7), `runs.csv` + `meta.json` (fáze),
`grid_A.csv`, `grid_B.csv`, `best_configs.csv` (řádky `ts2020_calA`, `ts2020_calB`, s `provisional`
0/1 a poznámkou, které cíle míjí), `vocab_matching.csv`, `profile.csv`, `profile.png`, `report.md`,
`notes.md`.

## Kritéria přijetí
1. Validace: náš port na výchozích parametrech leží v rozmezí jejich vzorků (h2 ±0,1, hapax ±0,05,
   D64 ±0,1, edge MI ±0,01, podíl pořadí ±0,3 p. b.); jinak hledat chybu portu, dokud nesedí, nebo
   doložit, že stejný rozdíl mají i jejich hotové vzorky (pak je to rozdíl metrik, ne portu).
2. `best_configs.csv` první řádek = doporučený null (`ts2020_calB`, pokud loss < 0,15 a složky < 0,05,
   jinak lepší z A/B s `provisional=1` a výčtem minutých cílů).
3. Tři seznamy necílových metrik (reprodukované / nereprodukované / nekonzistentní split vs. join)
   plus odstavec o vlivu velikosti slovníku na paradigmatické rozdíly.
4. `python -m pytest -q` prochází; Python 3.10 / pandas 2.3 kompatibilita; bez scipy.
5. Podezřelé: validace mimo tolerance u více než dvou metrik současně; H2 < 2,0; `share_near` > 0,9.

## Rozpočet
Generátor je pomalejší než náš (Python port, rich-get-richer): odhad 5–15 s na vzorek 34 k tokenů.
Kalibrace A 16×5 + B 48×3 ≈ 224 vzorků ≈ 30–60 min na 4 procesech; fáze B 100 vzorků s plnými
metrikami ≈ 30 min; celkem 2–3 h v jedné session. Jedna větev, PR na konci.
