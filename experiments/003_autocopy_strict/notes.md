# 003 – deník session

## Změny balíku (větev `feat/autocopy-strict-edits`, sloučeno do `exp/003-autocopy-strict`)

* Nový modul `voytools/wordgrammar.py`: `learn_word_grammar` / `WordGrammar` (doložené bigramy a trigramy
  s okrajovými značkami, slotová gramatika, sampler `sample`), `grammaticality(df_or_tokens, model)`.
* `autocopy_generate` / `autocopy_like`: nové parametry `edit_model`, `new_token_source`, `word_grammar`
  (předem naučená gramatika, aby se neučila pro každý vzorek) a `info` (skutečné počty; `autocopy_like` je
  ukládá do `attrs['autocopy_info']`). Výchozí hodnoty nemění výstup: regresní test
  `tests/test_wordgrammar.py::test_default_autocopy_output_is_unchanged` a `..._on_real_files_...` porovnávají
  otisky (sha256) výstupu zaznamenané se starým kódem (syntetický seed, ZL i GC `autocopy_like`).
* Testy prošly i v prostředí Python 3.10.20 / pandas 2.3.3 (uv venv), tedy u vlastníka.

## Výklad spec

* **Trigram:** editace se přijme, jen když jsou doložené všechny trigramy výsledného tokenu (`^ g1 g2`, …,
  `gn-1 gn $`; jednoglyfový token `^ g $`). Implementačně: kandidátní glyf substituce/inzerce se losuje
  (s bigramovými vahami jako dřív) jen mezi kandidáty, kteří dají přípustný token; když pro zvolenou pozici a typ
  editace žádný není (nebo delece není přípustná), je to neúspěšný pokus; po `max_tries` (10) neúspěších zůstane
  token nezměněn. Skutečný podíl nezměněných tokenů se zapisuje (`p_unchanged`, `p_edit_failed`).
* **Slot:** zvolena naučená gramatika (ne pevné pořadí podle Zattery – v balíku by šlo jen o ručně opsanou
  tabulku pro Eva, pro v101 nepoužitelnou). Výchozí pořadí slotů = průměrná relativní pozice glyfu v tokenu,
  jeden slot na glyf; to na reálném textu připouští jen 36 % (ZL) / 49 % (GC) tokenů, protože např. `o`, `d`, `k`
  stojí na začátku i uprostřed slova. Proto se greedy přidávají další sloty (glyf, místo vložení), dokud
  gramatika nepřipouští 99 % reálných tokenů (nebo zisk kroku < 0,02 %). Opakovatelné glyfy (`e`, `i`, `h` na ZL
  basic; `c` na GC) = glyf následovaný sám sebou v ≥ 5 % případů; jejich běh zabírá jeden slot do nejdelší
  doložené délky. Plus všechny bigramy doložené.
* **Grammar sampler:** náhodná procházka po doložených trigramech od `^` do `$`, délka z reálného rozdělení
  (token-vážené délky v glyfech), v každém kroku jen glyfy, které udrží slotovou gramatičnost; poslední glyf musí
  mít doložené `… $`. Slepá ulička = nový pokus (až 100), pak záložní seed token (počítá se, `n_grammar_fallback`).
* Gramatika i trigramy se učí z **celého** reálného textu dané sady (P, clean), preset generátoru (`basic` ZL,
  `None` GC); seed tokeny = prvních 100 tokenů.
* Okrajové značky jsou interní řetězce, ne `^`/`$` (v101 má `$` jako glyf). Původní bigramový kód používá
  `^`/`$` – na GC to může vzácně splést glyf `$` s koncem slova; neměněno kvůli regresi výchozího chování.
* **Fáze B:** „10 nejlepších“ = 5 nejlepších konfigurací fáze A pro každý `edit_model`, aby profil měl
  nejlepší trigram i nejlepší slot s 50 seedy. Stejných 10 konfigurací na všech 4 sadách. Nejlepší
  trigram/slot/doporučená konfigurace se volí na ZL split a drží se pevně i pro ostatní sady.
* **Převzato z 001** (nepočítáno znovu): bigramový generátor (nejlepší konfigurace fáze B dané sady), bod
  p_new = 0,5 (jen horní mez), míchání v řádcích, page bootstrap reálného textu. Bodové hodnoty reálného textu
  přepočteny a shodují se s 001 na 1e-15 (`real/*/check_vs_001.json`). Gramatičnost bigramových vzorků 001
  dopočtena na týchž vzorcích (stejné seedy, stejné volání, 5 seedů; `bigram001/`).
