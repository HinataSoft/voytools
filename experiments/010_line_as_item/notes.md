# 010 – deník session

* Větev `exp/010-line-as-item` vychází z `exp/009-register-profile` (sdílí `voytools/controls.py` a
  načítání kontrol z `experiments/009_register_profile/run.py`); 009 v době běhu nebylo v `main`.
* **Nové funkce v balíku:** `stats.line_boundary_mi` (MI na d = 1 zvlášť uvnitř řádku / přes hranici řádku
  v odstavci / přes hranici odstavce na stejné stránce; korekce permutací pravých štítků uvnitř každého druhu
  dvojic, takže bias se rovná pro různé počty dvojic) a `stats.vertical_similarity` (shoda ≤ k editací s
  tokenem přímo nad, šikmo nad, jinde v řádku nad, vlevo ve vlastním řádku; očekávání = přesná střední
  hodnota permutace řádku nad; permutační p).
* **Afixová třída** = 20 nejčastějších dvojic (prefix, sufix) + `other`, afixy naučené na celém korpusu
  (`decompose_table`), nully Voynichu s reálnými seznamy. Kmen = 200 nejčastějších + `other`.
* **Záhlaví (metrika 2):** „multiset tříd zbytku řádku“ je pro plug-in MI příliš řídký; použito
  I(třída prvního; třída každého dalšího tokenu řádku) (párová verze, `head_rest`) a zvlášť 2., 3. a
  poslední token. Null: permutace prvních tříd mezi řádky stejné délky (bod 200×, subsample 20×).
* **Šablony (metrika 3):** (i) H(třída), H(třída | délka řádku, pozice) pro řádky 4–8 tokenů,
  H(třída_i | třída_{i−1}) uvnitř řádku; zisky `gain_pos`, `gain_prev`. (ii) podpisy řádků 4–8 tokenů:
  různé podpisy na 1 000 řádků, pokrytí 20 nejčastějšími (součet přes délky), totéž po „sloupcovém“
  zamíchání (tokeny se stejnou délkou řádku a pozicí permutované mezi řádky). (iii) zlib (úroveň 9)
  posloupnosti tříd s oddělovačem řádků; zisk po zamíchání v řádku a po sloupcovém zamíchání.
* **Stejná délka (metrika 4):** 20 000 náhodných dvojic řádků, Jaccard množin tříd; očekávání = totéž po
  globálním zamíchání tokenů (délky řádků zachovány); `jac_samelen_effect` = přebytek stejné délky minus
  přebytek různé délky.
* **Konce řádků (metrika 6):** `-m`/`-g` = poslední glyf `m` nebo `g` (zahrnuje `-am`), jen Eva; pro
  všechny zdroje délka posledního tokenu (řádky ≥ 2 tokeny) vs ostatní a maximální obohacení koncového glyfu
  (glyfy s ≥ 50 výskyty).
* **Varianta (b) z 008** (poziční glyfy v afixech) neprovedena – 008 zatím nemá výsledky.
* Nully: ts2020 calB (ZL split, join), šablona 005 (ZL, GC split), ac-tri 86 (ZL, GC split), míchání v řádku,
  míchání v poziční třídě (007), permutace celých řádků stejné délky (`samelen`), vše 30 seedů.
  Kontroly 009: Caesar, Plinius, Celsus, italština, Culpeper (próza a po větách), Linnaeus (próza a po
  položkách), bod + 100 subsamplů.
