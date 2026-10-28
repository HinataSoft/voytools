# 020 – deník session

Větev `exp/020-paragraph-as-unit`, 28. 10. 2026. Kód: `run.py`, `summarize.py`, modul `voytools/paragraphs.py`
(testy `tests/test_paragraphs.py`), rozšíření `voytools/controls.py` (`tess_paragraphs`, `merge_short_items`,
`layout_paragraphs`, `paragraph_control`; testy přidány do `tests/test_controls.py`).

## Volby při výkladu spec

* **Jednotka a sdílení.** Řádek = lokus `P` (clean). Sdílení dvojice řádků = počet *různých* kmenů společných oběma
  řádkům; očekávání = Σ_s (1−(1−p_s)^{n_i})(1−(1−p_s)^{n_j}) (tahy s vracením z frekvencí kmenů), kde p_s jsou
  frekvence kmenů **na dané stránce** (hlavní; struktura nad slovníkem stránky) nebo globální (vedlejší). Poměr
  s = Σ pozorované / Σ očekávané přes dvojice skupiny (ne průměr poměrů – dvojice s očekáváním 0 by jinak padaly).
  Globální očekávání má u permutace řádků nenulové Δ (≈ +0,05: stránky s krátkými odstavci mají jiný poměr
  stejných/cizích dvojic i jiné sdílení) – proto je hlavní očekávání stránkové, u něhož permutace dává Δ ≈ 0.
* **Kmeny.** `decompose_table` s naučenými afixy (dekompozice a), kmeny s ≥ 2 výskyty v korpusu (hapax nelze
  sdílet a nejsou ani v očekávání). Varianty: `st` (vše), `st50` (bez 50 nejčastějších), `tok` (celé tokeny ≥ 2),
  `ng` (kmeny s odstraněnými počátečními gallows `k t p f cth ckh cph cfh`, v101 `g j f u k h` + velká) – kontrola
  artefaktu prvního řádku. Generátory jsou dekomponovány seznamy afixů reálného textu (jako 015), kontroly vlastními.
* **Odstavce.** Index odstavce přepočítán jako *souvislý běh* řádků se stejným `par_index` v pořadí na stránce;
  v ZL se tím f84v rozpadne na 5 odstavců místo 4 (řádky f84v.42–51, jednoslovné, nesou `par_index` 0 za
  odstavcem 3). Jinak 740 → 741 odstavců.
* **Null (a) přelosování hranic** = permutace pořadí délek odstavců na stránce (počet a délky zachovány). Pozor:
  stránka s odstavci stejné délky je pod tímto nullem beze změny a i jinak se část dvojic nemění – null (a) tedy
  **nedává Δ = 0**, ale dolní mez „kolik Δ zbývá, když hranice posunu“; úplný null bez struktury je (b) permutace
  řádků (Δ ≈ 0 z konstrukce; kritérium 3).
* **Metrika 3 (první řádek).** Dvojice jen dopředu (dřívější řádek i, pozdější j). `first_own` = první řádek
  odstavce vs další řádky vlastního odstavce ve vzdálenosti d, `first_next` = první řádek vs řádky *následujícího*
  odstavce v témže d (d = 1…6; při d ≤ 3 jsou takové dvojice jen u krátkých odstavců), `mid_own` / `mid_next`
  totéž pro prostřední řádky (lip ≥ 1, ne poslední). „Ohlašuje“ = (first_own − first_next) > (mid_own − mid_next);
  pouhý kraj sníží first_own i first_next stejně. Navíc Δ bez prvních řádků (`stepnf`, dvojice bez
  prvního řádku odstavce) – odečítá artefakt generátoru, který první řádky kopíruje z jiného zásobníku.
* **Metrika 4 (hranice stránky).** Dvojice (a-tý poslední řádek stránky i, b-tý první řádek stránky i+1), d = a+b−1
  ≤ 3, jen po sobě jdoucí stránky `P` v pořadí souboru (`rank` z celé tabulky; mezilehlá stránka bez `P` textu dvojici
  ruší); očekávání globální, srovnání s dvojicemi uvnitř stránky v témž d (globální očekávání), poměr přes/uvnitř.
  Podskupiny: stránka i+1 téhož jazyka **a** sekce / jiného; `across_open` = poslední odstavec stránky bez značky
  `<$>` – v ZL jen 2 stránky (f68r1, f84v), v GC podle `paragraph_end`; u kontrol s vlastními odstavci =
  odstavec pokračuje na další stránce (`par_continues`). Blábol s vlastními stránkami (`gibberish_native`) jen zde.
* **Metrika 5.** Mezi dvojicemi výskytů společného kmene podíl identických tokenů minus očekávání z paradigmatu kmene
  (Σ_f q_f² přes korpus), stejné / přes hranici, d ≤ 3 a sdruženě (d = 0 v tabulce).
* **Bootstrap** reálu = převzorkování stránek s vracením (váhy stránek), 200×.

## Kontroly s vlastními odstavci (`paragraph_control`)

Řádky a stránky jako ZL split P (4 128 řádků, 207 stránek), odstavce z pramene: konec odstavce **ukončí řádek**
(zbytek šablony řádku se nevysází, řádek je kratší – jako poslední řádek odstavce v rukopise), další odstavec začíná
novým řádkem. Index odstavce na stránce začíná vždy 0 (jako v tabulce ZL); zda odstavec pokračuje přes stránku, nese
`par_continues`.

| kontrola | zdroj odstavců | tokenů | odstavců (medián řádků) | stránek s ≥ 2 odst. | stránek s otevřeným koncem |
|---|---|---|---|---|---|
| `celsus_chap` | kapitoly Celsa (řádky `.tess`), knihy 1–5 | 33 926 | 275 (12) | 52 | 203 |
| `pliny_chap` | kapitoly Plinia, knihy 16–30 | 33 202 | 461 (8) | 148 | 197 |
| `celsus_sent40` | věty seskupené na ≥ 40 tokenů uvnitř kapitoly (slabá kontrola) | 31 562 | 831 (5) | 201 | 184 |
| `pliny_sent40` | totéž | 31 563 | 850 (5) | 203 | 193 |
| `linnaeus_par` | položka = odstavec (položky < 3 tokeny sloučeny) | 23 507 | 2 776 (1) | 207 | 114 |
| `pliny_prose` (009) | próza v layoutu ZL **s hranicemi ZL** (nic neznamenají) | 34 268 | 741 (4) | 159 | – |
| `gibberish_pooled` (013) | blábol v ořezu 60 stránek s hranicemi ZL | 10 024 | 235 (4) | 44 | – |
| `gibberish_native` (013) | vlastní řádky a stránky pisatelů, 1 odstavec = stránka | 10 032 | 88 | 0 | – |

Reál ZL split: 741 odstavců (medián 4 řádky), 159 stránek s ≥ 2 odstavci. Kapitoly Celsa jsou 3× delší než
odstavce Voynichu (medián 12 řádků) a dávají jen 52 stránek s ≥ 2 odstavci – proto i varianta `sent40`.
Přepisy blábolu (Gaskell & Bowern) **nenesou odstavce**: 38 souborů má jen značky stránek (`__________`) a 85
prázdných řádků uvnitř stránek celkem (zjevně mezery, ne členění) – blábol je tu jen jako negativní kontrola s
hranicemi ZL a pro hranici stránky s vlastními stránkami.

## Průběh

* Jeden vzorek: `line_table` + `pair_content` ≈ 1,3 s (ZL split, 59 519 dvojic), všechny metriky 4 variant ≈ 1 s;
  layoutové nully ≈ 1 s, ts2020 ≈ 12 s, šablona ≈ 45 s na vzorek. Čtyři paralelní procesy (ZL/GC × split/join a
  po nich kontroly).
* První běh ukázal **Δ u ts2020 (0,23–0,25) nad reálem (0,146)** – port kopíruje první řádky odstavců z jiného
  zásobníku (`choose_source`, režim `paragraph_initial`) a přidává gallows, tj. skok z konstrukce přes první řádek.
  Proto přidána varianta Δ bez prvních řádků (`stepnf`) a `ng`.
