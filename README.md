# voytools

Malý, testovaný Python balíček pro statistickou analýzu přepisů Voynichova rukopisu ve formátu
**IVTFF 2.0** (Zandbergen). Čte přepisové soubory (ZL3b, GC2a, RF1b, IT2a, interlineární LSI),
převádí je na tabulky `pandas` (tokeny, loci, stránky) a nabízí základní statistiky: frekvence,
entropie glyfů, vzájemnou informaci sousedních tokenů, polohové profily, podobnost sousedů,
afixy a paradigmata, BPE a nulové modely (zamíchání, generátor autokopírování).

Závislosti jsou pouze `numpy` a `pandas` (žádné `scipy`); `matplotlib` je volitelný. Balíček nepoužívá
síť a nikde neobsahuje pevné absolutní cesty. Všechny náhodné procedury mají `seed`.

## Obsah

1. [Instalace](#instalace)
2. [Rychlý start](#rychlý-start)
3. [Příkazová řádka](#příkazová-řádka)
4. [Tabulka tokenů: slovník sloupců](#tabulka-tokenů-slovník-sloupců)
5. [Tabulky `loci` a `pages`](#tabulky-loci-a-pages)
6. [Pravidla pro mezery a další volby tokenizace](#pravidla-pro-mezery-a-další-volby-tokenizace)
7. [Abecedy a presety glyfů](#abecedy-a-presety-glyfů)
8. [Filtrování, „clean“, interlineární soubory, odstavce](#filtrování-clean-interlineární-soubory-odstavce)
9. [Přehled funkcí](#přehled-funkcí)
10. [Příklady](#příklady)
11. [Testy a ověření](#testy-a-ověření)
12. [Experimenty, výsledky a cloudové session](#experimenty-výsledky-a-cloudové-session)
13. [Známá omezení](#známá-omezení)
14. [Zdroje](#zdroje)

## Instalace

Potřebujete Python 3.10 nebo novější. V adresáři balíčku (např. `D:\Work\Voynich\analysis\`):

```
pip install -e .            # numpy + pandas
pip install -e .[plot]      # navíc matplotlib (grafy)
pip install -e .[dev]       # navíc pytest
pip install -e .[parquet]   # navíc pyarrow (export do Parquet)
```

Bez instalace lze balíček použít i ze zdrojů: příklady v `examples/` si cestu k balíčku doplní samy
a `python -m voytools ...` funguje z kořenového adresáře balíčku.

Přepisové soubory se předpokládají v sousedním adresáři `D:\Work\Voynich\transliterations\`.
Příklady a testy je hledají na několika místech (viz [Příklady](#příklady)); cestu lze vždy zadat
parametrem nebo proměnnou prostředí `VOYTOOLS_DATA` (adresář se soubory `ZL3b-n.txt` atd.).

## Rychlý start

```python
from voytools import load_tables, filter_tokens
from voytools.stats import vocabulary_summary, glyph_entropy_report, adjacent_token_mi
from voytools.morphology import decompose_table, paradigm_metrics

tables = load_tables(r"D:\Work\Voynich\transliterations\ZL3b-n.txt")    # tokens, loci, pages
text = filter_tokens(tables.tokens, generic_type="P")                    # odstavcový text, čisté tokeny
print(vocabulary_summary(text, by="language"))                           # slovník podle Currierova jazyka
print(glyph_entropy_report(text, "extended", scope="all"))               # H(X), H(X_t+1|X_t), H3
print(adjacent_token_mi(text, shuffles=100, seed=0))                     # MI sousedních tokenů
dec = decompose_table(text)                                              # prefix + kmen + sufix
print(dec.attrs["prefixes"], dec.attrs["suffixes"])                      # naučené seznamy afixů
print(paradigm_metrics(dec).summary)                                     # metriky paradigmat
```

`load_tables` vrací objekt `Tables` s atributy `tokens`, `loci`, `pages`, `header`, `warnings`
(seznam upozornění parseru; chybně utvořený identifikátor locu vyvolá `IVTFFError`)
a `transcribers`. Tabulky lze uložit: `tables.export("vystup", fmt="csv")`.

## Příkazová řádka

Po instalaci jako `voytools <příkaz> SOUBOR [volby]`, jinak `python -m voytools <příkaz> SOUBOR [volby]`.
Výstup jsou zarovnané textové tabulky v čisté ASCII (znaky z privátní oblasti Unicode, kterými se
reprezentují vysoké ASCII znaky `@nnn;`, se vypisují zpět jako `@nnn;`, takže výstup nedělá potíže ani
ve starších konzolích Windows); volba `--csv` vypíše CSV. Chyby (chybějící soubor, špatný formát)
končí na stderr a návratovým kódem 2.

| Příkaz | Co dělá |
|---|---|
| `summary SOUBOR` | stránky, loci podle typu, tokeny, typy, hapaxy, délky, entropie glyfů, rozpis podle typu/jazyka/ilustrace/ruky |
| `export SOUBOR VÝSTUP [--tables]` | uloží tabulku tokenů do CSV (nebo `.parquet`); `--tables` přidá `VÝSTUP_loci.csv` a `VÝSTUP_pages.csv` |
| `freq SOUBOR [--top N] [--by SLOUPEC]` | nejčastější tokeny (případně po skupinách) |
| `entropy SOUBOR [--preset basic\|extended\|none\|auto] [--scope line\|all] [--by SLOUPEC] [--boundary SYM]` | `H(X)`, `H(X_{t+1}\|X_t)`, `H3`; bez `--preset` oba presety pro Evu |
| `positions SOUBOR [--which first\|last\|both] [--top N] [--permutations N] [--seed S]` | počáteční/koncové glyfy tokenů podle polohy v řádku, Jensen-Shannonova divergence tříd, volitelně permutační p-hodnoty |
| `affixes SOUBOR [--side prefix\|suffix\|both] [--max-len N] [--min-count N] [--top N] [--by SLOUPEC]` | produktivní prefixy a sufixy |
| `mi SOUBOR [--shuffles N] [--cap N] [--window N] [--no-neighbor] [--seed S]` | MI sousedních tokenů a krajních glyfů s korekcí zamícháním, podobnost sousedů |

Společné volby (výběr tokenů a tokenizace):

| Volba | Význam |
|---|---|
| `--generic-type T` | typ locu: `P` (odstavce), `L` (popisky), `C`, `R`, seznam `P,L` nebo `all` (statistické příkazy mají výchozí `P`, `summary` a `export` `all`) |
| `--language A\|B`, `--illustration H`, `--hand 1` | výběr podle proměnných stránky `$L`, `$I`, `$H` |
| `--transcriber X` | ID přepisovatele v interlineárních souborech (výchozí: nejúplnější, s upozorněním na stderr) |
| `--clean` / `--no-clean` | vynechá tokeny s nečitelnými znaky, alternativními čteními a vysokými ASCII znaky (výchozí: zapnuto pro statistiky a `summary`, vypnuto pro `export`) |
| `--space-policy split\|join\|strict` | zacházení s nejistou mezerou `,` (viz níže) |
| `--alternatives first\|mark` | zacházení s alternativami `[a:b]` |
| `--infer-paragraphs` | dopočítá odstavce v souborech bez značek `<%>` a `<$>` (např. RF1b) |
| `--csv` | tisk ve formátu CSV |

Příklady:

```
python -m voytools summary D:\Work\Voynich\transliterations\ZL3b-n.txt
python -m voytools entropy ZL3b-n.txt --scope all
python -m voytools mi ZL3b-n.txt --shuffles 100
python -m voytools positions ZL3b-n.txt --which first --permutations 200 --language B
python -m voytools affixes GC2a-n.txt --side suffix --min-count 50 --top 20
python -m voytools freq LSI_ivtff_0d.txt --transcriber H --top 20
python -m voytools export ZL3b-n.txt tokens.csv --tables
```

Na ZL3b i GC2a trvá `summary` asi 2 s a `mi --shuffles 100` (včetně podobnosti sousedů) asi 4 s.

## Tabulka tokenů: slovník sloupců

`tables.tokens` má jeden řádek na token v pořadí souboru a přesně tyto sloupce (41):

| Sloupec | Typ | Význam |
|---|---|---|
| `idx` | int | pořadové číslo tokenu v souboru (od 0); po filtrování zůstává původní |
| `folio` | str | **název stránky** podle IVTFF (`f1r`, `f67r1`, `fRos`), nikoli číslo folia |
| `quire` | str / NA | `$Q` – kvintern (A–T) |
| `page_in_quire` | str / NA | `$P` – stránka v kvinternu |
| `folio_in_quire` | str / NA | `$F` – folio v kvinternu |
| `bifolio` | str / NA | `$B` – bifolium v kvinternu |
| `illustration` | str / NA | `$I` – typ ilustrace: A astronomická, B biologická, C kosmologická, H herbář, P farmaceutická, S jen hvězdy na okraji, T jen text, Z zodiak |
| `language` | str / NA | `$L` – Currierův jazyk (A, B) |
| `hand` | str / NA | `$H` – písař (1–5) |
| `currier_hand` | str / NA | `$C` – Currierova ruka |
| `extraneous` | str / NA | `$X` – cizí písmo na stránce (C, M, O, S, V) |
| `locus_id` | str | identifikátor locu bez locatoru, typu a přepisovatele (`f1r.3`; u čísel jako `44a` včetně písmene). V interlineárních souborech **není jedinečný napříč přepisovateli** |
| `locus_num` | int | číslo locu |
| `locator` | str | locator: `@` nezávislá poloha (první locus stránky), `+` pod předchozím, `*` začátek řádku u levého okraje, `=` na stejném řádku jako předchozí, `&` totéž na kruhové linii, `~` totéž s horším zarovnáním, `/` vpravo nahoře, `!` text, který v rukopisu není |
| `locus_type` | str | úplný typ locu (`P0`, `Pt`, `L0`, `Lz`, `Cc`, `Ri`, ...) |
| `generic_type` | str | obecný typ: `P` lineární text v odstavcích, `L` popisky, `C` text po obvodu kruhu, `R` text po poloměru |
| `subtype` | str | podtyp (`0`, `t`, `z`, ...) |
| `transcriber` | str / NA | ID přepisovatele (jen interlineární soubory) |
| `par_index` | Int64 / NA | pořadí odstavce na stránce (od 0); NA mimo odstavce a pro loci jiného typu než `P` |
| `line_in_par` | Int64 / NA | pořadí řádku v odstavci (od 0); NA jako `par_index` |
| `is_par_first_line` | bool | první řádek odstavce |
| `is_par_last_line` | bool | poslední řádek odstavce (značka `<$>` nebo konec před dalším odstavcem) |
| `line_index` | int | pořadí locu na stránce (od 0; u interlineárních souborů zvlášť pro každého přepisovatele) |
| `token_index` | int | pořadí tokenu v lociu (od 0) |
| `n_tokens` | int | počet tokenů v lociu |
| `pos_from_end` | int | vzdálenost od konce locu (poslední token = 0) |
| `rel_pos` | float | `token_index / (n_tokens - 1)`; 0 pro jednoslovný locus |
| `is_first`, `is_last` | bool | první / poslední token locu |
| `token` | str | normalizovaný token: bez komentářů a značek, s prvním čtením alternativ (viz `--alternatives`), ligatury `{..}` bez závorek |
| `token_raw` | str | token tak, jak je v souboru (včetně `[a:b]`, `{..}`) |
| `len_chars` | int | délka ve znacích (nejistá mezera `,` v politice `strict` se nepočítá) |
| `len_glyphs` | int | délka v glyfech podle presetu tabulky (`attrs['glyph_preset']`) |
| `has_unreadable` | bool | token obsahuje `?` |
| `has_alternative` | bool | token obsahuje alternativní čtení `[a:b]` |
| `has_ligature` | bool | token obsahuje ligaturu `{..}` |
| `has_high_ascii` | bool | token obsahuje vysoký ASCII znak `@nnn;` (128–255) |
| `after_uncertain_space`, `before_uncertain_space` | bool | token sousedí s nejistou mezerou `,` (jen politika `split`) |
| `after_drawing_break`, `before_drawing_break` | bool | token sousedí s přerušením kresbou `<->` nebo `<~>` |

Konvence pro chybějící hodnoty: proměnná stránky, která není uvedena (nebo má hodnotu `@`), je `NA`
(`None` / `NaN` podle verze pandas; v seskupených statistikách se objeví jako skupina `NA`).
Sloupce `par_index` a `line_in_par` mají typ `Int64`. Tabulky mají `df.attrs['alphabet']` a
`df.attrs['glyph_preset']`; z nich si glyfové funkce odvozují `preset='auto'`.

Po `decompose_table` přibudou sloupce `prefix`, `stem`, `suffix` a v `attrs` použité seznamy afixů.

## Tabulky `loci` a `pages`

`tables.loci` (jeden řádek na locus) má stejné sloupce stránky a locu jako tokeny a navíc
`num_suffix` (písmeno za číslem, `44a`), `n_chars`, `paragraph_start`, `paragraph_end`, `comments`
(volné komentáře `<!...>`), `inline_tags` (textové značky `<@L=B>`), `raw_text`, `text` (text po odstranění
komentářů a značek), `line_no` (číslo řádku v souboru) a `continued` (locus se skládá z řádků spojených `/`).

`tables.pages` (jeden řádek na stránku): `folio`, `page_index`, `line_no`, proměnné stránky
(`quire` ... `extraneous`), `n_loci`, `n_loci_P`, `n_tokens`, `n_comment_lines`.

Textové značky `<@H=2>` se použijí na řádek, na kterém stojí, a na všechny další loci stránky.

## Pravidla pro mezery a další volby tokenizace

Mezery mezi tokeny jsou `.` (jistá) a `,` (nejistá). Politika se volí parametrem
`space_policy` (`build_tables`, `load_tables`) nebo `--space-policy`:

| Politika | Chování | Počet tokenů v ZL3b (všechny loci, bez filtru) |
|---|---|---|
| `split` (výchozí) | `.` i `,` dělí; sousední tokeny dostanou příznaky `after_/before_uncertain_space` | 39 026 |
| `join` | `,` se odstraní a obě části se spojí v jeden token | 36 278 |
| `strict` | dělí jen `.`; `,` zůstává v tokenu jako literální znak (délky a glyfy ho ignorují) | 36 278 (2 449 tokenů s čárkou) |

Přerušení kresbou `<->` a `<~>` dělí vždy (a nastaví příznaky `*_drawing_break`). Více oddělovačů za sebou
(po odstranění vycpávek ve starých souborech) se zredukuje na nejsilnější z nich.

Další volby:

* `alternatives='first'` (výchozí) nahradí `[a:b:c]` prvním čtením; `'mark'` ponechá `[a:b:c]` doslova
  (token je pak zvláštní typ, který se snadno odfiltruje; při počítání glyfů je celá alternativa **jedním**
  pseudo-glyfem). Příznak `has_alternative` se nastaví vždy.
* Vysoké ASCII znaky `@nnn;` (nnn = 128–255) se v tokenu ukládají jako jediný znak z privátní oblasti Unicode
  `chr(0xE000 + nnn)`; `voytools.render_high_ascii` je převede zpět na `@nnn;`, stejně to dělá CSV export a CLI.
* Ligatury `{..}` se v `token` zapisují bez závorek, v `token_raw` se závorkami; příznak `has_ligature`.
* Starší soubory (IVTFF < 2.0 s abecedou Eva, např. `LSI_ivtff_0d.txt`, verze 1.5) používají `!` a `%` jako
  vycpávky a mezery pro zarovnání; `parse_ivtff(..., strip_fillers='auto')` je odstraní. V abecedě `v101` (GC2a)
  jsou `!` a `%` písmena a zůstávají.
* Nečitelné `?` a `???` zůstávají v tokenu (příznak `has_unreadable`).

Počty tokenů závisejí na politice a na filtrování: například `join` dává 36 278 tokenů ve všech loci a 32 612
v odstavcích (`P`), `split` po vyčištění 37 978 ve všech loci a 33 205 v `P0`. Podrobnosti jsou ve `VERIFICATION.md`.

## Abecedy a presety glyfů

Eva zapisuje některé skupiny tahů jako více písmen. Statistiky na úrovni glyfů proto umějí tyto skupiny
slučovat (`voytools.glyphs`):

| Preset | Sloučené skupiny (nejdelší shoda zleva) |
|---|---|
| `basic` | `cth ckh cph cfh ch sh` |
| `extended` | `basic` + `iin in ee` (`iin` se zkouší před `in`) |
| `none` / `None` | identita: každý znak je glyf |

`preset='auto'` (výchozí ve statistikách) zvolí `basic` pro abecedy `Eva-` a `EvaT` a **identitu pro ostatní
abecedy, zejména `v101`** (GC2a), kde `c`, `h`, `s` nejsou části složených glyfů. Výsledek odpovídá
postupné náhradě `str.replace` v pořadí `cth, ckh, cph, cfh, ch, sh, iin, in, ee` (ověřeno testem na všech
tokenech ZL3b). Další funkce: `collapse_eva(token, preset)`, `glyph_sequence(df, preset, separator=None,
scope='line'|'all')`, `to_symbol_string`, `to_symbol_strings` (každý glyf jako jeden znak, např. pro externí nástroje).

## Filtrování, „clean“, interlineární soubory, odstavce

`filter_tokens(df, generic_type, language, illustration, hand, *, transcriber, quire, locus_type,
subtype, folio, clean=True)` vybírá podle metadat; každý selektor přijme hodnotu, seznam nebo řetězec
oddělený čárkami (`'P,L'`), `'all'`/`None` znamená bez omezení. Řádky s chybějící hodnotou ve vybraném sloupci
nikdy nevyhovují. Výsledek je kopie s novým indexem `0..n-1` (původní číslo zůstává v `idx`).

`clean=True` (výchozí) vynechá tokeny s příznaky `has_unreadable`, `has_alternative` nebo `has_high_ascii`;
`clean=False` ponechá vše; seznam (`['unreadable', 'ligature']`) vybere příznaky explicitně
(k dispozici i `'alternative'`, `'high_ascii'`).

**Interlineární soubory** (LSI, 18 přepisovatelů) mají zvlášť číslování řádků pro každého přepisovatele.
Před statistikami vyberte jednoho: `load_tables(path, transcriber="H")`. CLI bez `--transcriber`
vybere nejúplnějšího (u LSI je to `H`, stejný počet tokenů jako IT2a, tj. 37 919) a oznámí to na stderr;
`export` ponechává všechny.

**Odstavce.** `par_index` počítá odstavce na stránce podle značek `<%>` a `<$>`. Soubory bez značek
(RF1b) lze zpracovat s `infer_paragraphs=True` (`--infer-paragraphs`): locus typu `P`, který leží mimo
označený odstavec, začíná nový implicitní odstavec, má-li locator `@` nebo `*`.

## Přehled funkcí

Všechny funkce, které berou tabulku tokenů, přijímají volitelně `by=` (sloupec nebo seznam sloupců) a vrací
`DataFrame` se sloupci skupin vpředu. Hodnoty chybějící ve sloupci skupiny tvoří skupinu `NA`.

**`voytools.ivtff`** – parser: `parse_ivtff(path)`, `parse_ivtff_text(text)`, `tokenize_locus`, modelové třídy
`Header`, `Page`, `Locus`, `Token`, `Transliteration`, `render_high_ascii`. Neznámé konstrukce jen vyvolají
upozornění v `warnings` (uložených je nejvýše 20 na kategorii, počty v `warning_counts`).

**`voytools.table`** – `build_tables`, `load_tables`, `filter_tokens`, `export(df, path, fmt='auto')`
(Parquet bez `pyarrow` se zapíše jako CSV s varováním `RuntimeWarning`).

**`voytools.stats`**

| Funkce | Výsledek |
|---|---|
| `token_frequency(df, by)` | `token, count, freq, rank` |
| `vocabulary_summary(df, by)` | `n_tokens, n_types, n_hapax, hapax_share` (hapaxy / typy), `top50_coverage`, `ttr`, `mean_len`, `mean_len_glyphs` |
| `zipf_fit(freq)`, `zipf_slope(freq)` | regrese log četnosti na log pořadí |
| `length_distribution(df, unit)` | rozdělení délek ve znacích nebo glyfech |
| `entropy`, `perplexity`, `conditional_entropy(seq, order)` | plug-in odhady v bitech |
| `glyph_entropy_report(df, preset, by, boundary, scope)` | `H1 = H(X)`, `H2 = H(X_{t+1}\|X_t)`, `H3`, `perplexity = 2^H2`, `mutual_info = H1 - H2` pro text bez mezer; `scope='line'` (výchozí) nepřekračuje konce řádků, `'all'` bere celý text jako jeden proud (tak to dělají Rozanova a Temerev) |
| `adjacent_token_mi(df, cap=2000, shuffles=100, seed)` | MI sousedních tokenů uvnitř řádku; typy nad `cap` se sloučí do `<other>`; `excess_bits` = pozorovaná MI minus průměr z `shuffles` zamíchání tokenů uvnitř řádků, `excess_share` = `excess_bits / H_cap` (podíl entropie tokenu vysvětlený levým sousedem) |
| `edge_glyph_mi(df, shuffles, seed, preset)` | MI posledního glyfu tokenu a prvního glyfu následujícího, stejná korekce |
| `positional_profile(df, which='first'\|'last', unit, preset, permutations, seed)` | rozdělení počátečních/koncových glyfů ve třídách `line_first, mid, line_last` a `parfirst_*` (řádky začínající odstavec), matice JSD, volitelně permutační p-hodnoty; jednoslovné řádky se ignorují |
| `neighbor_similarity(df, window=5, max_edit=1, shuffles, seed)` | podíl tokenů s „blízkým“ (editační vzdálenost ≤ `max_edit`, včetně shody) tokenem mezi `window` předchozími; základ je zamíchání uvnitř stránek (Timm a Schinner) |
| `vocabulary_overlap(df, col, a, b)`, `labels_vs_text(df)` | průnik slovníků dvou skupin, popisky proti odstavcovému textu |
| `keyness(df, col, unit='token'\|'stem')`, `top_keyness` | Dunningovo G² pro charakteristické tokeny/kmeny skupin |

**`voytools.morphology`**

| Funkce | Výsledek |
|---|---|
| `affix_inventory(df, side, max_len, min_count, preset, by)` | `affix, n_glyphs, count, share, n_types, n_remainders, n_free` (zbytky, které jsou samy tokeny), `exact`, `ext_share` (podíl pokračování nejčastějším delším afixem; blízko 1 = afix není samostatná jednotka, např. `q` → `qo`) |
| `learn_affixes(df, side, top_n, ...)` | seznam produktivních afixů podle `n_free`, závislé afixy vynechány |
| `decompose(token, prefixes, suffixes)` | chamtivé dělení nejdelší shodou: nejdřív prefix, pak sufix; `decompose("qokeedy", ["qo","o"], ["edy","dy","y"])` = `('qo', 'ke', 'edy')` |
| `decompose_table(df, prefixes=None\|'reference'\|seznam, suffixes=...)` | kopie tabulky se sloupci `prefix`, `stem`, `suffix`; seznamy v `attrs`; `learn_on=` naučí seznamy na jiné tabulce |
| `paradigm_table(dfd, counts, min_stem_count, top_stems, top_suffixes, axis)` | kontingenční tabulka kmen × sufix (prázdný afix je `-`) |
| `paradigm_metrics(dfd, min_stem_count=5, shuffles, seed)` | počet kmenů, sufixů na kmen, `H(sufix\|kmen)`, `I(prefix; sufix)` a `I(prefix; sufix \| kmen)` včetně zamíchaných základů a `cmi_excess` |
| `learn_bpe`, `apply_bpe`, `bpe_dependence_curve` | BPE uvnitř tokenů (remíza: větší dvojice) a „závislostní mezera“ `H(U) - H(U_{t+1}\|U_t)` pro několik počtů sloučení |

Výchozí seznamy `'reference'` (`REFERENCE_PREFIXES`, `REFERENCE_SUFFIXES`) vznikly naučením na odstavcovém textu ZL3b:

* prefixy: `che cho qo ch ol sh d l o s t y`
* sufixy: `aiin ain chy dy al am ar ol or d l r s y`

**`voytools.nulls`** – `shuffle_within_lines(df, seed)`, `shuffle_within_pages(df, seed)` (vlastnosti tokenu
cestují s tokenem, polohové sloupce zůstávají), `autocopy_generate(seed_tokens, n_tokens, window, ...)`,
`as_table(lines, template=None)`, `autocopy_like(template, ...)`.
Parametry `edit_model='bigram'|'trigram'|'slot'` (editace se přijme jen do tokenu, který připouští trigramový
nebo slotový model slova) a `new_token_source='seed'|'grammar'` (nové slovo z gramatiky místo ze seed slovníku);
výchozí hodnoty dávají přesně původní výstup. `info={}` vrátí skutečné počty (např. `p_unchanged`, podíl
nezměněných kopií), `autocopy_like` je ukládá do `attrs['autocopy_info']`.

**`voytools.wordgrammar`** – `learn_word_grammar(tokens, preset)` (bigramy, trigramy a slotová gramatika slova
naučené z reálných tokenů; `WordGrammar.accepts(glyphs, model)`, `WordGrammar.sample(rng)`),
`grammaticality(df_or_tokens, model='slot', grammar=None, reference=None, weight='token')` – podíl tokenů,
které model připouští.
Zamíchání zachovává slovník, četnosti, délky řádků a vše, co závisí jen na slovech samotných (afixy, paradigmata,
vnitroslovní statistiky): **ty se na zamíchané kopii nezmění**, takže slouží jako nulový model pouze pro efekty
pořadí (MI sousedů, MI krajních glyfů, podobnost sousedů, polohové profily a entropie glyfů textu bez mezer,
která na hranicích slov závisí na pořadí slov).

**`voytools.plots`** (volitelně, matplotlib): `plot_zipf`, `plot_positional`, `plot_bpe_curve`; funkce vracejí
`Figure` a při zadání `path` uloží obrázek.

## Příklady

Ve složce `examples/` (každý skript přijme volitelný argument s cestou k souboru):

| Skript | Obsah |
|---|---|
| `01_summary.py [soubor] [--plot zipf.png]` | souhrn souboru, slovník, Zipf, délky, entropie glyfů |
| `02_compare_sections.py [soubor]` | Currierův jazyk A proti B: slovník, entropie, klíčová slova a kmeny, afixy, popisky proti textu |
| `03_morphology_and_nulls.py [soubor] [--seed N] [--shuffles N]` | rozklad na prefix-kmen-sufix, paradigmata, zamíchání a autokopírování se stejnými seznamy afixů, statistiky závislé na pořadí, BPE |

Bez argumentu skript hledá `ZL3b-n.txt` v tomto pořadí: `data/` v repu, `<examples>/../../transliterations/` (rozložení
`D:\Work\Voynich\analysis\examples` + `D:\Work\Voynich\transliterations`), `<examples>/../transliterations/`,
`./transliterations/`, `../transliterations/` a adresář z proměnné `VOYTOOLS_DATA`; když nic nenajde, vypíše
všechny zkoušené cesty. Příklad 3 trvá asi 10 s.

## Testy a ověření

```
pip install -e .[dev]
python -m pytest
```

Testy používají syntetický minimální soubor IVTFF (`tests/conftest.py`) a „smoke“ testy na skutečných
souborech (ZL3b, GC2a, RF1b, IT2a, LSI); pokud soubory nejsou nalezeny (hledá se `VOYTOOLS_DATA`, `data/` v repu,
`../transliterations`, ...), přeskočí se. Čísla naměřená na skutečných souborech a porovnání s publikovanými
hodnotami jsou ve `VERIFICATION.md`; soubor lze znovu vytvořit příkazem
`python tools/make_verification.py [ADRESÁŘ_SE_SOUBORY]`.

## Experimenty, výsledky a cloudové session

Repo je připravené pro práci v Claude Code (lokálně i v cloudových session nad GitHubem):

* `CLAUDE.md` – kontext projektu, výzkumné směry, konvence a pravidla pro každou session.
* `data/` – kopie transliterací (ZL3b, GC2a, RF1b, IT2a, LSI), aby session nepotřebovaly síť.
* `experiments/NNN_name/spec.md` + `run.py` – zadání a kód experimentu; šablona v `experiments/_template/`,
  protokol v `experiments/README.md`.
* `results/NNN_name/` – commitované výstupy (`runs.csv`, `summary.csv`, `report.md`, grafy).
* `voytools.experiment` – `run_seeds` (navazovatelný běh po seedech, zápis po každém běhu),
  `seed_range` (rozdělení seedů mezi paralelní session), `merge_runs`, `summarize_runs`.
* `docs/cloud-setup.md` – založení repa, nastavení cloudového prostředí (`scripts/setup_env.sh`,
  timeouty), spouštění experimentů a cesta výsledků zpět.

## Známá omezení

* **Verze pandas.** Balíček byl vyvíjen a testován s Pythonem 3.13, numpy 2.x a pandas 3.0. Kód se vyhýbá API
  zrušeným v pandas 3 i API dostupným jen v pandas 3 a nespoléhá na dtype `object` / `str`, ale s pandas 2.3
  (a Pythonem 3.10) nebyl spuštěn; spusťte na svém stroji `python -m pytest`.
* **Odhady entropie** jsou prosté (plug-in), bez korekce na malý vzorek. `H3` na ~140 000 glyfech je již
  podhodnocené. Korekci obstarávají základy z mícháním (`excess_bits`), ne analytická oprava. Odhady MI
  sousedních tokenů závisejí na `cap` (výchozí 2000 typů).
* **Řádky.** Loci s locatorem `=`, `&`, `~` leží na stejném fyzickém řádku jako předchozí locus (např.
  nadpis `Pt`), ale počítají se jako samostatné řádky (`line_index`), takže polohové statistiky a MI uvnitř
  řádku je berou odděleně. Přerušení kresbou `<->` dělí token, ne řádek.
* **Apostrof** `'` (součást některých ligatur, např. `{c'y}`) zůstává v tokenu a počítá se jako samostatný glyf.
* **Výchozí rozsah entropie glyfů** je `scope='line'` (bez bigramů přes konce řádků); publikované hodnoty
  (Rozanova a Temerev) odpovídají `scope='all'`. Rozdíl u ZL3b je asi 0,05 bitu.
* **Jazyk `NA`.** Stránky bez proměnné `$L` tvoří ve statistikách podle jazyka skupinu `NA`
  (v ZL3b 531 čistých odstavcových tokenů).
* **Generátor autokopírování** je záměrně jednoduchý základ (Timm a Schinner 2020 jen v hrubých rysech):
  kopie z posledních `window` tokenů s 0-2 úpravami glyfů omezenými doloženými dvojicemi glyfů. Nereprodukuje
  publikovaný algoritmus; jeho entropie glyfů je vyšší a podíl blízkých sousedů mnohem vyšší než ve skutečném
  textu. Pro smysluplný výstup předávejte `allowed_from=` se všemi slovy skutečného textu.
* **Parquet** vyžaduje `pyarrow`; jinak se zapíše CSV.
* **Windows.** CLI vypisuje jen ASCII (vysoké znaky jako `@nnn;`); CSV je v UTF-8. V Excelu importujte
  s kódováním UTF-8 (tokeny samotné jsou ASCII).
* **Staré formáty.** Z formátů před IVTFF 2.0 byl vyzkoušen jen `LSI_ivtff_0d.txt` (1.5); `locus_id` nese
  písmeno čísel typu `44a`, ale není jedinečný napříč přepisovateli.

## Zdroje

* Currier, P. H. (1976). *Some important new statistical findings.* In: D'Imperio, M. E. (ed.), *New Research on
  the Voynich Manuscript: Proceedings of a Seminar* (rozdělení na jazyky A a B, polohové efekty v řádku).
* Reddy, S., Knight, K. (2011). *What we know about the Voynich manuscript.* Proceedings of the 5th ACL-HLT
  Workshop on Language Technology for Cultural Heritage, Social Sciences, and Humanities.
* Timm, T., Schinner, A. (2020). *A possible generating algorithm of the Voynich manuscript.* Cryptologia 44(1),
  1-19, doi:10.1080/01611194.2019.1596999.
* Rozanova, M., Temerev, A. (2026). Entropie glyfů, MI sousedních tokenů a BPE uvnitř tokenů; použito podle
  referenčních skriptů z balíku `voynich-units` (`analysis/reproduce_*.py`). Úplný bibliografický údaj
  nebyl k dispozici.
* Zattera, M. (2022). *A new transliteration alphabet brings new evidence of word structure and multiple
  writing systems in the Voynich manuscript.* Cryptologia.
* Zandbergen, R. *IVTFF – Intermediate Voynich MS Transliteration File Format 2.0* (`docs/IVTFF_format_2.0.txt`).
