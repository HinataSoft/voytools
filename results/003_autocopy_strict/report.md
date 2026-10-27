# 003 – Autokopírování s editacemi vázanými gramatikou slova

Kód: `experiments/003_autocopy_strict/run.py`. Změny balíku: `voytools/wordgrammar.py` a parametry
`edit_model` / `new_token_source` v `voytools/nulls.py` (větev `feat/autocopy-strict-edits`). Deník a výklad
spec: `experiments/003_autocopy_strict/notes.md`. Všechna čísla pocházejí ze souborů v této složce
(`grid_summary.csv`, `phaseB_configs.csv`, `best_configs.csv`, `profile.csv`, `diag_summary.csv`,
`real/*/point.csv`, `summary_info.json`). Bigramová reference je převzatá z `results/001_autocopy_calibration/`.

## Co bylo spočítáno

| Krok | Rozsah | Soubor |
|---|---|---|
| Reálný text: bodové hodnoty + gramatičnost | ZL3b, GC2a × `split`, `join`; P text, clean | `real/<sada>/point.csv` (shoda s 001 na 1e-15: `check_vs_001.json`) |
| Fáze A: 432 konfigurací × 3 seedy, cílové metriky + edge/token MI + gramatičnost | ZL split | `runs.csv` (z `phaseA/part0..3`), `grid_summary.csv` |
| Fáze B: 10 konfigurací (5 nejlepších z A pro `trigram` a 5 pro `slot`) × 50 seedů, všechny metriky | všechny 4 sady | `phaseB/<sada>/runs.csv`, `phaseB_configs.csv`, `best_configs.csv` |
| Gramatičnost vzorků bigramového generátoru 001 (stejné seedy, tytéž vzorky, 5 seedů) | všechny 4 sady | `bigram001/<sada>/runs.csv` |
| **Mimo spec:** rozklad přebytku H2 (5 variant kolem nejlepší trigram/slot konfigurace × 3 seedy) | ZL split | `diag/ZL_split/`, `diag_summary.csv` |
| Převzato z 001 | bigram (nejlepší konfigurace fáze B dané sady), bod p_new 0,5 (horní mez), míchání v řádcích, page bootstrap | `profile.csv` (zdroje `bigram_001`, `bigram_001_pnew05_upper`, `shuffle_lines_001`) |

Gramatika slova i trigramy se učí z **celého** reálného textu dané sady. Seed tokeny jsou prvních 100 tokenů.
Výklad `slot`: naučená gramatika, ne pevné pořadí podle Zattery. Výchozí pořadí slotů dává průměrná relativní
pozice glyfu. Jeden slot na glyf ale připouští jen 36 % (ZL) / 49 % (GC) reálných tokenů, a tak se greedy
přidávají další sloty pro glyfy, které stojí na dvou místech slova, dokud gramatika nepřipouští 99 % tokenů.
`ee`, `ii` = běh opakovatelného glyfu v jednom slotu. Podrobnosti v `notes.md`.

## Čísla

**Gramatika reálného textu** (`real/*/point.csv`, token-vážené; v závorce podíl typů):

| | ZL split | ZL join | GC split | GC join |
|---|---|---|---|---|
| slot | 0,990 (0,952) | 0,990 (0,961) | 0,991 (0,960) | 0,990 (0,963) |
| trigram / bigram | 1 / 1 | 1 / 1 | 1 / 1 | 1 / 1 |
| počet slotů | 63 | 70 | 113 | 122 |

Kritérium 5 (gramatičnost reálu ≥ 0,95) je splněno. Vázané generátory gramatiku opravdu drží
(`summary_info.json`): `trigram` má ve všech 1 296 bězích fáze A trigramovou gramatičnost 1,0, `slot` má
slotovou gramatičnost ≥ 0,987 (zbytek tvoří kopie 1 % negramatických reálných seed tokenů). Sampler
`grammar` nepotřeboval žádnou zálohu (0 fallbacků). Bigramový generátor 001 má slotovou gramatičnost
0,88–0,91 a trigramovou 0,60–0,65 (`profile.csv`, `bigram_001`).

**Kalibrace (kritérium 1 – nesplněno).** Nejlepší konfigurace (fáze B, ZL split, 50 seedů) je
`config_id 86`: `trigram`, `seed`, p_new 0,3, window 32, p_copy_exact 0,35, p_two_edits 0,5, recency 0, s
**loss 0,269**. Složky: share_near 0,004, hapax 0,031, **H2 0,057**, **n_types 0,086**, mean_len 0,001.
Skutečný podíl nezměněných tokenů je 0,380 (nominálně 0,35; 4,8 % pokusů o editaci selže). Nejlepší
`slot` (`config_id 316`, p_copy_exact 0,5, window 64) má loss 0,293 (H2 0,075, hapax 0,065). Pro srovnání:
bigram 001 v mřížce 0,978, bigram 001 s p_new 0,5 (jen horní mez, mimo povolené p_new) 0,202.

| zdroj (ZL split) | loss | share_near | hapax | H2 | n_types | mean_len | p_unchanged |
|---|---|---|---|---|---|---|---|
| reál | – | 0,173 | 0,676 | 2,689 | 6 796 | 4,99 | – |
| trigram, config 86 | 0,269 | 0,183 | 0,557 | 3,329 | 8 785 | 4,86 | 0,380 |
| slot, config 316 | 0,293 | 0,160 | 0,504 | 3,424 | 6 998 | 4,73 | 0,510 |
| bigram 001, config 573 | 0,978 | 0,256 | 0,534 | 3,609 | 9 753 | 4,90 | 0,516 |
| bigram 001, p_new 0,5 (**horní mez**) | 0,202 | 0,161 | 0,565 | 3,365 | 7 482 | 4,75 | – |

Na ostatních sadách je config 86 také nejlepší nebo téměř nejlepší (`phaseB_configs.csv`): loss ZL join
0,250, GC split 0,223 (nejlepší tam config 96, 0,222), GC join 0,204 (jako jediná sada se všemi složkami
< 0,05).

Proč to nejde (`grid_summary.csv`, `diag_summary.csv`):

1. **Žádná z 432 konfigurací nemá loss < 0,15 ani všechny složky < 0,05.** Nejmenší maximum složek má
   config 96 (0,052). Práh 0,05 pro H2 odpovídá relativní odchylce 22 %, tedy H2 ≤ 3,29. Tento práh
   překoná 63 konfigurací, ale skutečné H2 je i tak nejméně 3,17 bitu (trigram), resp. 3,31 (slot), proti
   2,69 reálně.
2. **H2 a hapax jdou proti sobě.** Všech 63 konfigurací s H2 pod prahem minou hapax (složka ≥ 0,072,
   hapax ≤ 0,58). Nejvyšší hapax v celé mřížce je 0,640, reálně 0,676. Nižší H2 vyžaduje méně editací (vyšší
   p_copy_exact), ale editace jsou hlavní zdroj nových typů a hapaxů.
3. **Rozklad mimo spec** (`diag_summary.csv`, ZL split, kolem config 86 / 316): kopie bez editací mají H2
   2,667 / 2,675, ale jen 87 typů. Samotná slova z gramatiky (p_new 1, bez kopírování) mají H2 2,99 a hapax
   0,696 (≈ reál), share_near ale jen 0,060. Samotné editace seed slov bez kopírování dají H2 3,25 / 3,35.
   Přebytek H2 tedy nevzniká tím, že by editace porušovaly gramatiku (tu teď drží). Vzniká tím, že editovaná
   slova používají doložené trigramy v jiném poměru než reálný text. Generátor, který by trefil H2 i hapax,
   by musel brát nová slova z gramatiky ve velkém podílu. Pak ale ztratí podíl blízkých sousedů.
4. `new_token_source='grammar'` při p_new ≤ 0,3 nepomáhá (nejlepší loss 0,402 trigram, 0,467 slot).
   Gramatická slova přidávají typy a zhoršují n_types. `recency 0,3` je výrazně horší (nejlepší loss 5,9,
   share_near 0,5–0,8).

**Necílové metriky** (`profile.csv`; z = (reál − průměr) / SD generátoru přes 50 seedů; konfigurace
vybrané na ZL split a stejné na všech sadách):

| metrika | reál ZL / GC split | trigram 86 ZL / GC split | slot 316 ZL / GC split | bigram 001 ZL / GC split | z trigram (4 sady) | z slot (4 sady) |
|---|---|---|---|---|---|---|
| edge_mi_excess | 0,185 / 0,226 | 0,000 / 0,000 | 0,000 / 0,000 | 0,000 / 0,000 | 174 až 263 | 175 až 248 |
| token_mi_excess_share | 0,0089 / 0,0129 | 0,0001 / 0,0003 | −0,0003 / 0,0001 | −0,0003 / 0,0000 | 5,5 až 12,1 | 6,9 až 14,0 |
| pos_jsd_first_mid | 0,186 / 0,195 | 0,001 / 0,003 | 0,001 / 0,003 | 0,002 / 0,003 | 325 až 425 | 363 až 416 |
| bpe_gap64 | 1,058 / 1,126 | 0,759 / 0,812 | 0,879 / 0,853 | 0,601 / 0,586 | 8,1 až 28,0 | 3,7 až 21,4 |
| cmi_excess | 0,088 / 0,089 | 0,209 / 0,226 | 0,230 / 0,239 | 0,229 / 0,191 | −25,5 až −20,2 | −25,1 až −20,4 |
| frac_stems_multi_suffix | 0,912 / 0,839 | 0,620 / 0,694 | 0,550 / 0,622 | 0,507 / 0,698 | 9,9 až 19,2 | 13,8 až 24,2 |
| mean_suffixes_per_stem | 5,25 / 4,72 | 2,71 / 3,25 | 2,48 / 2,95 | 2,28 / 3,15 | 27,6 až 64,6 | 32,6 až 69,0 |

**Kritérium 2 – tři seznamy** (vázaný generátor = trigram 86 i slot 316; „v ZL i GC“ = ve všech čtyřech
sadách):

* (a) **reprodukované (|z| < 3 v ZL i GC):** žádná.
* (b) **nereprodukované v obou:** všech sedm – `edge_mi_excess`, `token_mi_excess_share`, `pos_jsd_first_mid`,
  `bpe_gap64`, `cmi_excess`, `frac_stems_multi_suffix`, `mean_suffixes_per_stem`. Platí pro trigram i slot.
* (c) **nekonzistentní mezi ZL a GC:** žádná. Nejblíže k reprodukci je `bpe_gap64` u slotu na ZL join
  (z 3,7). Na GC zůstává z 19–21.

Predikce ze spec: u okrajové MI, pořadí tokenů a poziční JSD se nereprodukce potvrdila, překvapení nenastalo.
U BPE žlabu, podílu kmenů s více sufixy, sufixů na kmen a `cmi_excess` se predikce „gramatika je dožene“
**nepotvrdila**. Vázání editací posunulo BPE žlab k reálu (0,60 → 0,76 / 0,88 na ZL split), ale ne do 3 SD.
Paradigmatické metriky se oproti bigramu změnily jen málo, a to ve všech sadách.

**Kritérium 3.** `best_configs.csv`, řádek 1 = config 86 (trigram, seed, p_new 0,3, window 32,
p_copy_exact 0,35, p_two_edits 0,5, recency 0) s **`provisional=1`**, protože kritérium 1 není splněno.
Poslední řádek je bod 001 s p_new 0,5, označený jen jako horní mez.

**Kritérium 4.** `python -m pytest -q`: 182 passed, 1 skipped (Python 3.13 / pandas 3.0). Totéž platí pro
Python 3.10.20 / pandas 2.3.3. Regresní test výchozího chování
(`tests/test_wordgrammar.py::test_default_autocopy_output_is_unchanged`, `..._on_real_files_is_unchanged`)
prochází. Žádné scipy.

**Kritérium 5 (podezřelé hodnoty).** Gramatičnost reálu ≥ 0,990 ve všech sadách. Fáze A ani B nemají žádný
běh se share_near > 0,9 nebo H2 < 2,0 (`summary_info.json`). Maximum share_near ve fázi A je 0,82 (recency
0,3, window 16).

**Literatura** (`profile.csv`, zdroje `ref_RT2026_*`, jen ZL split). Self-citation generátor Rozanové
a Temereva má h2 2,6–2,7, hapax 59–60 %, D64 0,65–0,86 a edge MI 0,006–0,009. Náš trigram 86 má H2
3,33, hapax 0,557 a BPE žlab 0,759 (v jejich rozmezí) a edge MI 0,000. Jejich nízké h2 tedy naše vázání
editací nereprodukuje. Rozdíl je nejspíš v algoritmu (jejich generátor nemáme), ne v gramatice slova.

## Interpretace

Vázat editace gramatikou slova **nestačí** k tomu, aby generátor trefil všechny cílové vlastnosti najednou.
Pomůže hodně: loss 0,98 → 0,27 v povoleném rozsahu p_new. H2 ale zůstává o 0,6 bitu nad reálem a
podíl hapaxů pod ním. Editovaná slova jsou gramatická, ale skládají doložené trigramy v jiném poměru než
reálný text. Pro výzkumnou otázku je podstatné, že gramatika slova **nevysvětlila žádnou** z necílových
vlastností. Paradigmatická struktura (víc sufixů na kmen, nižší `cmi_excess`), BPE žlab, vazba na hranici
slov a poziční profil zůstávají ve všech čtyřech sadách mimo dosah vázaného generátoru, stejně jako u
bigramového.

## Co výsledek neukazuje

* Neukazuje, že text nese význam. Ukazuje jen, že **tento** generátor (kopírování + slotová/trigramová
  gramatika slova bez pozičního a mezislovního mechanismu) paradigmata a BPE žlab nereprodukuje. Generátor
  nemá žádný mechanismus pro edge MI, pořadí tokenů a poziční JSD, takže velká z tam testují spíš jeho
  konstrukci.
* Kalibrace není úplná: H2 a hapax jsou mimo cíl. Část rozdílu v BPE žlabu a paradigmatech proto může
  pocházet z nekalibrace. Proti tomu mluví bod 001 s p_new 0,5: je kalibrovaný jinak a paradigmata u něj
  vycházejí podobně (`bigram_001_pnew05_upper`).
* Slotová gramatika je naučená z dat a dost volná (63–122 slotů). Pevné Zatterovo pořadí nebylo
  testováno. Volnost gramatiky se projevuje tím, že slotové editace H2 snižují méně než trigramové.
* Gramatika se učí z celého textu, který se pak měří. Generátor tak „zná“ všechny doložené trigramy. To je
  pro něj výhoda, ne handicap.
* z-skóre jsou vůči SD generátoru. Page bootstrap reálného textu (z 001) existuje jen pro cílové metriky.
* Rozklad H2 (`diag`) je mimo spec, jen 3 seedy a jen ZL split.

## Doporučený další krok

1. Jako nulový model pro další experimenty brát **config 86** (`best_configs.csv`, provisional) a vždy
   uvádět, že H2 a hapax netrefuje. Paradigmatické výsledky experimentu 002 srovnávat s ním i s bigramem 001.
2. Zjistit, proč editovaná gramatická slova mají vyšší H2. Například vážit editace podle četnosti trigramů
   v reálném textu (ne podle bigramových vah) nebo editovat celé afixové skupiny místo glyfů.
3. Přidat poziční mechanismus (řádkové prefixy) a mezislovní vazbu (konec slova → začátek dalšího). Teprve
   pak brát edge MI a poziční JSD jako indicie.
