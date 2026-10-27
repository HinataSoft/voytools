# 007 – deník session

## Stav vstupů
* `main` při startu = 7d6a7a9; 004 (selfcitation) ani 005 (šablona) nejsou sloučeny → nully
  `selfcitation_like` a šablona 005 **vynechány**; běží autocopy config 86 z 003 (seedy 8600–8629 =
  stejné vzorky jako 003 phase B), míchání v řádcích a nové míchání v poziční třídě.
* `voytools/nulls.py`: přidána jen `shuffle_within_position_class` (+ 2 testy v `tests/test_nulls.py`).
* Prostředí session: Python 3.13; pandas připnut na 2.3.3 (jako u vlastníka), numpy 2.x; bez scipy.

## Interpretace spec (rozhodnutí)
* **Poziční třída** = místo v řádku (first/mid/last; `single` pro jednotokenový řádek) × `rel_pos`
  v 5 binech × `line_in_par` 0/1/2+ (2+ i pro chybějící). `parfirst` z `positional_profile` je
  totožný s `line_in_par == 0`, proto není zvlášť. V ZL split 22 obsazených tříd.
* **Odstavec** = stránka + `par_index`; řádek bez `par_index` (14 tokenů v ZL) je vlastní jednotka.
* **Cílové proměnné**: prefix / sufix z naučených seznamů reálného textu (12 + 14, `decompose_table`,
  preset basic; nully rozkládány týmiž seznamy); kmen = top 200 kmenů + `<other>` (kategorie z daného
  korpusu); délka = `len_glyphs` (basic) zastropená na 10.
* **Korekce mícháním** (metrika 1): pro I(Y;X|Z) je bias = průměr plug-in hodnoty po permutaci Y
  uvnitř vrstev Z (pro nepodmíněné globálně); 50 permutací pro bodový odhad, 10 pro subsample a nully.
* **Permutační null** (metrika 1, bodový odhad, 200 opakování): (a) Y permutováno uvnitř stránky
  (dle spec) – testuje jen poziční členy, sekce/písař/jazyk jsou vlastnosti stránky, takže pro ně
  je tento null z konstrukce totožný s reálem (p = 1 je očekávané, ne výsledek); (b) proto navíc
  permutace záznamů (sekce, písař, jazyk) mezi stránkami – test stránkových členů.
* Navíc proti spec: I(Y; jazyk) a I(Y; písař | sekce) (CLAUDE.md: písař a jazyk jsou zaměněné, uvádět obojí).
* **Jednotky MI(d)**: (a) poslední glyf tokenu t → první glyf tokenu t+d (preset extended / v101 znaky);
  (b) totéž pro BPE jednotky (64 sloučení naučených na reálném textu, aplikováno i na nully):
  poslední jednotka → první jednotka; (c) afixová třída = dvojice prefix|sufix s ≥ 20 výskyty, jinak
  `<other>` (stejná proměnná na obou koncích); (d) kmen s ≥ 20 výskyty + `<other>`.
  Vzdálenost d v tokenech uvnitř odstavce (přes řádky ano, přes odstavce ne).
* **Korekce MI(d)**: odečten průměr ze 100 (bod) / 20 (subsample, nully) míchání: `pc` = poziční
  třída uvnitř odstavce (spec), `line` = uvnitř řádku (spec); citlivostní `pcpage` = stejné třídy
  uvnitř stránky; `pcfine` = jemné třídy (rel_pos 10 binů × line_in_par 0..5/6+) uvnitř stránky
  (kovariát driftu, metrika 3).
* **Zamrznutí**: v ZL split zůstane při `pc` na místě 41 % tokenů (řádky 0 a 1 odstavce: ~70 %),
  protože v jednom odstavci je v dané třídě často jen jeden token. Korekce `pc` je proto
  konzervativní (stlačuje korigovanou MI k nule) – proto citlivostní `pcpage` (12 % zamrzlých).
* **Intervaly reálu**: 200 subsamplů 80 % stránek bez vracení; odchylky od průměru subsamplů škálovány
  √(m/(n−m)) = 2 a přičteny k bodovému odhadu z celého textu (prostý percentil 80% subsamplů by
  interval podhodnotil zhruba dvakrát). Nully: 2,5–97,5 % přes 30 seedů.
* **Polovina poklesu**: první d, kde korigovaná křivka klesne na ≤ 50 % hodnoty při d = 1
  (lineární interpolace); 10,5 = nedosaženo do d = 10; NaN, je-li hodnota při d = 1 ≤ 0.
* **Drift** (metrika 3): podíly 12 nejčastějších glyfů (+ other) podle `rel_pos` (5 binů, všechny
  tokeny), `rel_pos_mid` (jen vnitřní tokeny řádku – bez efektu prvního/posledního slova) a
  `line_in_par` (0..4, 5+); I(glyf; bin) korigovaná permutací binů mezi tokeny.

## Průběh
