# experiments/ – protokol

Každý experiment má vlastní složku `NNN_short_name/` s těmito soubory:

| Soubor | Obsah |
|---|---|
| `spec.md` | zadání: otázka, hypotézy a co by ukázal smysluplný vs. generovaný text, data a filtry, metriky, nulové modely, počet seedů/bootstrapů, výstupy, kritéria přijetí |
| `run.py` | kód; používá `voytools` a `voytools.experiment.run_seeds`; parametry přes argparse (`--seeds`, `--out`, `--chunk`), výchozí hodnoty = hodnoty ze spec |
| `notes.md` | (volitelné) deník session: rozhodnutí, odchylky od spec, slepé uličky |

Výsledky jdou do `results/NNN_short_name/` (stejné číslo a jméno):

| Soubor | Obsah |
|---|---|
| `runs.csv` | jeden řádek na seed/běh (zapisuje `run_seeds`, průběžně, navazovatelně) |
| `summary.csv` | agregace (`summarize_runs`) nebo vlastní tabulka efektů s intervaly |
| `report.md` | čísla (jen z `results/`), interpretace ve 2–3 větách, co výsledek **neukazuje**, další krok |
| `*.png` | grafy (matplotlib), malé |
| `meta.json` | parametry, časy, počet seedů (zapisuje `run_seeds`) |

Souhrn stavu napříč experimenty: `results/SYNTHESIS.md` (udržuje Cowork po každém experimentu).

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

Šablona: `_template/` (`spec.md` + `run.py`). Zkopírujte ji a doplňte.

## Seznam experimentů

| # | Název | Stav | Otázka |
|---|---|---|---|
| 001 | autocopy_calibration | hotovo (krit. 1 nesplněno: H2 nedosažitelné bigramovými editacemi; provizorní null p_new 0,5, loss 0,20) | kalibrace autokopírovacího generátoru na ZL; které vlastnosti textu zůstanou nevysvětlené |
| 002 | paradigm_coherence | hotovo (kmeny: paradigma užší než šablona, širší než kopírování; varianty rozptýlené; cmi_excess sám nerozlišuje) | chovají se kmeny a afixy jako paradigma, nebo jako šablona? |
| 003 | autocopy_strict | hotovo (krit. 1 nesplněno: loss 0,27, H2 3,3; žádná necílová vlastnost nereprodukována) | editace vázané gramatikou slova (trigram/slot): trefí generátor H2, a které necílové vlastnosti tím začne reprodukovat? |
| 004 | selfcitation_port | hotovo | port generátoru Timm & Schinner (R&T 2026, MIT) jako kanonický kopírovací null; validace proti tab. 7; velikost slovníku jako konfound |
| 005 | structured_template_null | hotovo | null „strukturovaný umělý jazyk“ s paradigmaty a tématem z konstrukce; kolik ingrediencí potřebuje k napodobení 002; robustnost afixů/parseru; A/B vs písař |
| 006 | topic_at_distance | hotovo | sdílení kmenů mezi stránkami: sekce vs vzdálenost vs písař; popisky (ok-/ot-) a jejich stopa v textu |
| 007 | position_vs_content | hotovo | rozklad informace o afixech/kmenech mezi pozici, sekci a písaře; reziduální MI(d) po odečtení pozice |
| 008 | ts2020_recheck | hotovo (ts2020 reprodukuje I(kmen; sekce \| písař) a „nic za d = 2“; pokles sdílení na hraně; A/B uvnitř písaře není poziční; poziční kmen není artefakt parseru) | kanonický null ts2020 a šablona 005 v 006/007; „poziční kmen“ třemi parsery; A/B podmíněně na pozici; doplněk syntézy |
| 009 | register_profile | hotovo (7 kontrol + 2 položkové varianty; Voynich nejblíž próze/ts2020, ne katalogu; Voynich-specifické: poziční/sekční vazby afixů, šíře paradigmatu, konec řádku) | registrový profil Voynichu vs skutečné texty (narativ, lékařství, botanika, herbář, katalog) a nully |
