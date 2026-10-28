# CLAUDE.md – Voynich text analysis (voytools)

Project owner: Václav ("Vasa"). He writes in **Czech**; reply in Czech. Code, docstrings,
commit messages and file names stay in English. Be concise and cite file paths.

## What this project is

Statistical analysis of the Voynich manuscript text (Beinecke MS 408) using the IVTFF
transliterations in `data/`. The research question is **whether the text carries meaningful
information**, approached through properties that meaningless or mechanically generated text does
not reproduce. It is *not* a decipherment project: we never propose readings, sound values,
plaintext languages or "solutions".

Working hypotheses and exclusions decided by the owner:

* Excluded: mathematical/algorithmic ciphers, number sequences, and the fringe "decipherment"
  literature (Cheshire, Pahlavi, Book of Dunstan, Jama, Zenodo-style "99.8 % decoded" claims).
  Do not cite or build on them.
* Allowed: some *obfuscation* of a natural-language source – redundant parts of "words", units
  that are not single letters, spaces that are not all word boundaries.
* Hardest competing hypotheses, always kept on the table: text generated without meaning
  (Timm & Schinner auto-copy, Rugg grille), human gibberish/glossolalia (Gaskell & Bowern 2022).

Three research directions (details in the owner's notes; summary here):

1. **Morphology as paradigm, not template** – do stems recur with different affixes, are
   prefix–suffix combinations constrained given the stem, are affix inventories stable across
   sections while stem inventories change? Templates and copy-and-modify generators give slot
   structure without paradigms.
2. **Position vs content** – which variance in token form is explained by position in line/
   paragraph (Currier, Feaster, Stafford) and which by section/topic; what remains after
   positional effects are removed (residual dependencies at distance 1, 2, 3 …).
3. **Topic at a distance and links to illustrations** – shared stems on non-adjacent pages that
   share content (same plant in herbal and pharma sections, zodiac vs astronomical pages); label
   stems (`ok-`/`ot-` labels) reappearing in running text.

## Status after experiments 001–007 (read `results/SYNTHESIS.md` for the table)

* Canonical copying null = `voytools.selfcitation.selfcitation_like` with config 21 from
  `results/004_selfcitation_port/best_configs.csv` (Timm & Schinner port, calibrated: loss 0.05).
  The older `autocopy_like` config 86 (003) stays only for continuity.
* Paradigm-level findings (002) are reproduced by copying (004) and by a structured template (005):
  they are evidence of a *system*, not of meaning. Text structure is dominated by line/paragraph
  position and by proximity; no sequential dependence beyond the immediate neighbour (007); stem
  sharing between pages follows distance, not section (006).
* Working frame for the next round: if the text carries meaning, the profile (hapax-rich, low token
  order, line as unit, strong positional conventions) points to a **catalogue / register / recipe**
  kind of text (names and formulaic entries), not continuous prose. Tests target that frame
  (008 recheck with ts2020, 009 register profile vs real genres, 010 line as item).

## Status after experiments 008–013 (see `results/SYNTHESIS.md`)

* Voynich-specific after all controls (prose, catalogue, human gibberish, Naibbe cipher, ts2020, template):
  place-bound properties only – affix–position binding, word-seam coupling, line-end convention, closed line –
  plus A/B as two sharply separated vocabularies. Within-word form variation (wide, even paradigms, hapax-rich,
  >90 % unexplained by position/register) is shared with human gibberish. No token-order, topic-at-distance or
  item/catalogue structure. Treat earlier "paradigm" signals as system, not meaning.
* Next round (014–016) asks whether place rules are a small set of conventions, whether rare stems recur at a
  distance above chance (the one signal copying and gibberish cannot give), and whether B-only stems are glyph
  variants of A stems.

## Status after experiments 014–016 (see `results/SYNTHESIS.md`)

* Place rules are a **small shared system**: ~10 glyph rules (8 prefixes / 4 seam pairs carry 80 %), same direction
  across A/B, scribes and sections, different weights in A vs B; nothing beyond the glyph pair at the seam.
  ts2020 + learned tables (`voytools/placerules.py`, R1–R3) reproduces position and 2/3 of the seam but breaks
  prefix–suffix coupling and adds types: copying + local rewrite does not give position, seam and word cohesion
  at once (open: 017).
* Rare stems recur with a **gradual decay over 1–15 pages, like Latin prose**; generators and the template give
  only a same-page peak and a flat remainder; far pairs are section/language/scribe-bound as in topic-structured
  Latin; no far excess above chance; form of repetition is independent of distance. Drift vs topic is undecided
  (open: 018, 020).
* B-only forms are **1/5 of B tokens** and are not glyph variants of A (fewer close counterparts than halves of one
  language); the only systematic operation is `-ey → -edy` (v101 `8`). Latin with another topic or register shows
  the same pattern; generators give nothing. Content vs register vs scribal circle undecided (open: 019).
* Working frame: the manuscript behaves like a fixed place-bound system with gibberish-like form variation, two
  vocabularies and prose-like vocabulary locality. The three properties no generator or gibberish gave (locality
  over pages, section-bound far pairs, new vocabulary in B) are each compatible with drift/register; the next
  round separates drift from topic blocks and names from morphology.

## Conventions that every analysis follows

* **Primary file `data/ZL3b-n.txt` (Eva), control `data/GC2a-n.txt` (v101).** A result counts
  only if it holds in both; what disappears when the alphabet changes is a transliteration
  artefact.
* **Units are not assumed.** Eva characters are glyph transcriptions, not letters (Rozanova &
  Temerev 2026): use `preset='extended'` composites or BPE units for glyph-level work, and run
  token-level work under all three space policies (`split`, `join`, `strict`); uncertain
  spaces (`,`) behave like word-internal junctures, never pool them silently with `.`.
* **Paragraph text (`generic_type == 'P'`) and labels (`'L'`) are analysed separately.**
  Circular/radial text (`C`, `R`) is excluded unless the question is about it.
* **Covariates:** Currier language (`language` A/B), Davis scribe (`hand`), section
  (`illustration`), paragraph-first line and line-first/last position. Scribe and language are
  confounded; report both.
* **Null models always run alongside the real text**, in the same code path: within-line and
  within-page shuffles (`voytools.nulls`), the auto-copy generator (`autocopy_generate`, after
  calibration in experiment 001), and slot-independence samplers where relevant. Type-level
  metrics (vocabulary, affix tables, paradigms) are invariant under token shuffles – say so and
  use generator-based nulls for them.
* **Uncertainty:** bootstrap over pages or quires, permutation tests, many seeds. Report effect
  sizes with intervals, not just p-values. Every random procedure takes a `seed`.
* **Reproducibility:** results written by scripts in `experiments/NNN_name/run.py`, using
  `voytools.experiment.run_seeds` (checkpointed, resumable) and committed to
  `results/NNN_name/`. Numbers in reports must come from files in `results/`.

Known baselines (full tables in `VERIFICATION.md`): ZL3b has 227 pages, 5 385 loci, 39 026
tokens (`split`), 36 278 (`join`), 34 268 clean paragraph tokens; glyph H(X_{t+1}|X_t) = 2.69
bits (extended preset, scope all), adjacent-token MI excess ≈ 0.9 % of H(X), edge-glyph MI ≈
0.20 bits, near-neighbour share 0.173 vs 0.146 shuffled. Do not re-derive these in reports;
cite them.

## Repository layout

```
voytools/        the package (parser, tables, glyphs, stats, morphology, nulls, experiment, cli, plots)
tests/           pytest suite (python -m pytest -q); real-file tests use data/
examples/        runnable examples (01_summary, 02_compare_sections, 03_morphology_and_nulls)
data/            IVTFF transliterations (read-only; never edit; see data/README.md)
experiments/     one folder per experiment: spec.md (the task), run.py (the code), notes
results/         committed outputs per experiment: runs.csv, summary.csv, report.md, *.png
docs/            cloud-setup.md (environment, workflow), further notes
scripts/         setup_env.sh (cloud environment setup)
README.md        package manual (Czech); VERIFICATION.md measured numbers
```

## Working rules for sessions (cloud or local)

* Start by reading `experiments/README.md` and the `spec.md` of the experiment you were given.
  If a spec is ambiguous, state the interpretation you chose in `report.md`; don't stall.
* Keep the package API stable; add functions rather than changing signatures. New public
  functions get a docstring and a test. `python -m pytest -q` must pass before every push
  (real-file tests take ~30 s).
* Dependencies: numpy, pandas only (matplotlib optional, imported lazily). **No scipy.**
  The owner runs Python 3.10 / pandas 2.3 locally; keep code compatible (no 3.11+ syntax,
  no pandas-3-only APIs).
* Never modify files in `data/`. Never delete or overwrite another session's `results/`
  folder; parallel sessions use disjoint seed ranges (`seed_range`) and separate subfolders,
  merged with `merge_runs`.
* Commit early and often on a branch named `exp/NNN-short-name` (or `feat/...` for package
  changes): push after each finished chunk of seeds and at least every ~15 minutes while a long
  run is in progress, so a paused or reclaimed VM loses nothing. Keep each results file small
  (CSV/MD/PNG, < 5 MB); large intermediates are gitignored.
* Long computations: use `run_seeds` chunks that finish within the command timeout
  (default 10 min foreground; see `docs/cloud-setup.md` for raising it), or run in the
  background with `nohup` and poll; never rely on a single multi-hour command.
* `report.md` ends with: numbers (from `results/`), interpretation in two or three sentences,
  what the result does *not* show, and suggested next step. Czech is fine for reports.
* When done, write `experiments/NNN_*/result.md` (template `experiments/_template/result.md`:
  a `Stav:` line and a 2–4 sentence `Souhrn:`). **Do not edit `experiments/README.md`** or any
  other experiment's folder – status lives only in `result.md`, so parallel PRs never conflict.
  Shared files (`voytools/`, `tests/`) change only by addition. `python scripts/status.py` prints
  the overview.
* Do not fetch papers or web pages unless the spec asks; the literature is summarised in the
  owner's `papers/README.md` outside the repo and in the spec files.

## Key references (for context, not to be re-read every session)

Currier 1976; D'Imperio 1978; Reddy & Knight 2011; Montemurro & Zanette 2013; Timm & Schinner
2020; Bowern & Lindemann 2021; Lindemann & Bowern 2020 (entropy); Gaskell & Bowern 2022
(gibberish); Zattera 2022 (slot alphabet); Feaster 2022; Stolfi 2000 (crust–mantle–core);
Rozanova & Temerev 2026, arXiv:2608.17096 (units; code at github.com/lrozanova/voynich-units);
Greshko 2025 (Naibbe cipher, Cryptologia).
