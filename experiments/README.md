# experiments/ – protokol

Každý experiment má vlastní složku `NNN_short_name/` s těmito soubory:

| Soubor | Obsah |
|---|---|
| `spec.md` | zadání: otázka, hypotézy a co by ukázal smysluplný vs. generovaný text, data a filtry, metriky, nulové modely, počet seedů/bootstrapů, výstupy, kritéria přijetí |
| `run.py` | kód; používá `voytools` a `voytools.experiment.run_seeds`; parametry přes argparse (`--seeds`, `--out`, `--chunk`), výchozí hodnoty = hodnoty ze spec |
| `notes.md` | (volitelné) deník session: rozhodnutí, odchylky od spec, slepé uličky |
| `result.md` | **stav a souhrn** – píše runner na konci: řádek `Stav:` (`hotovo` / `hotovo s výhradou` / `nedokončeno`), datum, větev, odkaz na report; pak 2–4 věty: hlavní čísla, odpověď na otázku, výhrady. Jediné místo, kde se stav eviduje (tento README se **nemění**, aby PR nekolidovaly). Šablona v `_template/result.md`. |

Výsledky jdou do `results/NNN_short_name/` (stejné číslo a jméno):

| Soubor | Obsah |
|---|---|
| `runs.csv` | jeden řádek na seed/běh (zapisuje `run_seeds`, průběžně, navazovatelně) |
| `summary.csv` | agregace (`summarize_runs`) nebo vlastní tabulka efektů s intervaly |
| `report.md` | čísla (jen z `results/`), interpretace ve 2–3 větách, co výsledek **neukazuje**, další krok |
| `*.png` | grafy (matplotlib), malé |
| `meta.json` | parametry, časy, počet seedů (zapisuje `run_seeds`) |

Souhrn stavu napříč experimenty: `results/SYNTHESIS.md` (udržuje Cowork po každém kole z `result.md` a reportů).
Přehled stavů na jeden pohled: `python scripts/status.py` (vypíše tabulku z `result.md`; nic neukládá).

Pravidla:

1. **Nulové modely ve stejném běhu.** Každá metrika se počítá pro reálný text i pro každý nulový
   model ve stejné funkci a stejném filtru; výsledky v jedné tabulce (sloupec `source`).
2. **ZL i GC.** Pokud spec neříká jinak, všechno dvakrát (`data/ZL3b-n.txt`, `data/GC2a-n.txt`).
3. **Tři politiky mezer** tam, kde jde o tokeny (`split`, `join`, `strict`); u glyfů preset
   `extended` a případně BPE jednotky.
4. **Seedy a bloky.** `run_seeds` zapisuje po každém seedu; bloky volte tak, aby jeden příkaz
   skončil do 10 minut. Paralelní session = disjunktní `seed_range`, podsložky `partK/`, nakonec
   `merge_runs`.
5. **Commit a push** po každém dokončeném bloku a nejméně každých ~15 minut během dlouhého běhu.
   Větev `exp/NNN-short-name`.
6. **Report píše ten, kdo počítal.** Bez čísel z `results/` žádný závěr. Interpretace opatrná:
   co by vysvětlil generátor, co pozice, co je artefakt abecedy.
7. **Nikdy** nenavrhovat čtení, zvukové hodnoty ani „řešení“; nepracovat s okrajovou literaturou
   (viz `CLAUDE.md`).
8. **Stav jen do `result.md` vlastního experimentu.** Runner needituje tento README ani cizí složky;
   sdílené soubory (`voytools/`, `tests/`) mění jen přidáváním (nové funkce, nové testy).

Šablona: `_template/` (`spec.md`, `run.py`, `result.md`). Zkopírujte ji a doplňte.

## Seznam experimentů

Stav a souhrn každého experimentu: `NNN_*/result.md`; čísla: `results/NNN_*/report.md`.

| # | Název | Otázka |
|---|---|---|
| 001 | autocopy_calibration | kalibrace autokopírovacího generátoru na ZL; které vlastnosti textu zůstanou nevysvětlené |
| 002 | paradigm_coherence | chovají se kmeny a afixy jako paradigma, nebo jako šablona? |
| 003 | autocopy_strict | editace vázané gramatikou slova (trigram/slot): trefí generátor H2, a které necílové vlastnosti tím začne reprodukovat? |
| 004 | selfcitation_port | port generátoru Timm & Schinner (R&T 2026, MIT) jako kanonický kopírovací null; validace proti tab. 7; velikost slovníku jako konfound |
| 005 | structured_template_null | null „strukturovaný umělý jazyk“ s paradigmaty a tématem z konstrukce; kolik ingrediencí potřebuje k napodobení 002; robustnost afixů/parseru; A/B vs písař |
| 006 | topic_at_distance | sdílení kmenů mezi stránkami: sekce vs vzdálenost vs písař; popisky (ok-/ot-) a jejich stopa v textu |
| 007 | position_vs_content | rozklad informace o afixech/kmenech mezi pozici, sekci a písaře; reziduální MI(d) po odečtení pozice |
| 008 | ts2020_recheck | doplnit ts2020 a šablonu do 005–007; je „poziční kmen“ artefakt parseru; A/B po podmínění pozicí |
| 009 | register_profile | registrový profil Voynichu proti latinské próze, lékařské/botanické latině, herbáři a katalogu (voynich-units) a nullům – kterému žánru je nejblíž? |
| 010 | line_as_item | řádek jako položka: hranice vs vnitřek, záhlaví řádku, šablony řádků, vertikální kopírování (test Timmova mechanismu) |
| 011 | suffix_variation | kolik z H(sufix \| kmen) vysvětlí pozice, registr A/B, písař, sekce; má zbytek strukturu (shoda, konzistence)? latina a ts2020 jako srovnání |
| 012 | ab_register | A/B: sdílení kmenů, mapování sufixů (řídké = zápis, difúzní = obsah), pozvolnost podle pořadí psaní, písař 3, A/B vs sekce |
| 013 | gibberish_and_cipher_controls | lidský blábol (Gaskell & Bowern) a Naibbe šifra latiny na profilu 009/010/011: které Voynich-specifické vlastnosti mají |
| 014 | place_rules_inventory | rozklad poziční informace a švu na glyfy/páry; stabilita napříč A/B, písaři, sekcemi; ts2020 doplněný o naučené místní tabulky – co dožene |
| 015 | rare_stem_recurrence | opakování vzácných kmenů (2–5×) na vzdálených stránkách nad náhodu; vazba dálkových dvojic na sekci; tvar opakování; popisková jména |
| 016 | ab_stem_variants | kmeny jen v B/A: glyfové varianty (změna zápisu) nebo nová slova; inventář editací; sdílení po substituci; celé tokeny a BPE bez parseru |
| 017 | conditional_seam_rewrite | ts2020 + přepis podmíněný sufixem a předchozím glyfem (R3′) a přepis prvního glyfu kmene: dá lokální mechanismus pozici, šev a vazbu prefix–sufix zároveň? kde přesně selže |
| 018 | locality_drift_vs_blocks | lokalita slovníku: hladký pokles (drift; ts2020 s pamětí přes stránky) vs bloky (segmentace, hranice sekcí); fyzické pořadí (kvatern, bifolio, anomálie sousedů) |
| 019 | b_only_forms | tvary jen v B rozdělené na `-edy` varianty a zbytek: jména/položky (úzké paradigma, popisky, hapax, pozice) nebo jiná morfologie; pravidlo `-ey → -edy` deterministické? |
| 020 | paragraph_as_unit | odstavec a stránka jako jednotky obsahu: sdílení kmenů mezi řádky podle vzdálenosti, skok na hranici odstavce/stránky proti nullům a próze s vlastními odstavci |
