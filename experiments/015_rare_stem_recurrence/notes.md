# 015 – poznámky ze session

## Interpretace spec (kde byla volnost)

* **Vzácné kmeny:** kmen = `stem` z `decompose_table` (a: naučené afixy na reálném P textu dané kombinace;
  b: + poziční glyfy 008/011). Prázdné kmeny vynechány. S naučenými afixy má ZL split jen **480 kmenů s 2–5
  výskyty** (1 331 tokenů, 1 398 dvojic) a 158 s 6–20; spec počítal s ~3 000 – to by odpovídalo celým
  tokenům, ne kmenům. Rozpočet je proto mnohem menší, intervaly širší.
* **Vzdálenost:** pořadí stránek v souboru přes všechny stránky (i ty jen s popisky), stránky = všechny
  s ≥ 1 čistým P tokenem (žádný práh 40 tokenů jako v 006; váhy velikostí to řeší).
* **Permutační null (tokeny kmene → náhodné stránky, váha ∝ velikost stránky, s vracením):** očekávání
  jsou **přesná** (`recurrence.expected_profile`, `same_by_distance`, `far_clusters`), protože každá
  dvojice výskytů jsou dva nezávislé tahy stránek. 200 vzorkovaných permutací (`perm_tokens`) dává interval
  a kontrolu nullu (průměrný přebytek musí být ≈ 0).
* **Dálkový přebytek** podle spec = podíl(≥ 16) − očekávání (`x_far`). Protože reál má silný lokální vrchol,
  je tento přebytek nutně záporný; proto navíc `x_far_nl` = podíl(≥ 16 | vzdálenost ≥ 4) − očekávání, tj.
  tvar rozdělení mimo lokální vrchol, a poměr pozorované/očekávané po binech (`ratio`).
* **Interval reálu:** bootstrap přes kmeny (200×; dvojice jsou vnořené v kmeni), ne přes stránky.
* **Dálkové dvojice a sekce:** očekávání párované podle vzdálenosti = vážený podíl stejných atributů mezi
  dvojicemi stránek přesně v té vzdálenosti, průměr přes dálkové dvojice (přesně); interval z bootstrapu
  kmenů a z permutace tokenů. Navíc rozsah `mid` (2–15 stránek) pro kvatern/bifolio.
* **Tvar opakování:** identický tvar = identický token (při stejném kmeni ⇔ stejný prefix i sufix);
  očekávání z paradigmatu kmene = Σ n_f(n_f−1)/(k(k−1)) na kmen, průměr přes dvojice třídy.
* **Shluky:** kmeny s ≥ 3 výskyty rozložené přes ≥ 3 kvaterny (a navíc ≥ 3 sekce); očekávání přesně
  (inkluze–exkluze přes ≤ 2 kategorie).
* **Popisky:** kmeny z `L` tokenů podtypů z/f/c/s/n (čisté), dekompozice stejnými seznamy; pás `l` =
  tyto kmeny s 2–5 výskyty v P textu (ZL split: 41 kmenů z 328 popiskových).
* **Generátory:** ts2020 config 21 (seedy 2100+, 30) jen Eva; šablona 005 `best` dané kombinace (seedy
  100000+, pool_seed = rep) s tématem a s `topic_concentration = inf` (varianta `no_topic` z 005; lokalita
  stránek λ ponechána); ac-tri 86 (seedy 8600+, 30). Generované vzorky se dekomponují seznamy reálu;
  popiskový pás = reálné popiskové kmeny ve vzorku.
* **Kontroly:** 009 `pliny_prose`, `celsus_prose`, `linnaeus_items` v layoutu ZL split (vlastní naučené
  afixy, jen dekompozice a); 013 `gibberish_pooled` a `voynich_crop` (60 stránek, vzdálenosti jen v ořezaném
  pořadí – biny 16+ téměř prázdné).
* **Permutace stránek (006):** `perm_global` / `perm_section` (100×): obsah stránky na jinou pozici,
  metadata zůstávají s pozicí.

## Prostředí

* Cloud: Python 3.13, pandas 3.0.5, numpy 2.5.3 (vlastník 3.10 / 2.3). Test `test_registers.py::
  test_transfer_and_mapping_fit` padal už na `main` (pandas 3 vrací read-only pole z `to_numpy`); opraveno
  `copy=True` v `voytools/registers.py` (kompatibilní s pandas 2.3).
* Externí kontroly: `git clone --depth 1` `lrozanova/voynich-units` do `../lrozanova/voynich-units` a
  `danielgaskell/voynich` do `../ext/gaskell` (cesty jako v 009/013).
* Časy: reál/permutace < 0,2 s; bootstrap 0,8 s; ts2020 3 s; šablona 18 s; ac-tri 2,4 s na vzorek.
