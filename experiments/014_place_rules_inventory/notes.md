# 014 – deník session

## Rozhodnutí a interpretace spec

* **Rozklad MI na položky** (`voytools.placerules`): pointwise členy p(x,y)·log₂ p(x,y)/(p(x)p(y)) sečtené po
  hodnotě afixu/glyfu (= p(x)·KL(p(pos|x)‖p(pos)), nikdy záporné) a po poziční třídě; pokrytí 50/80/95 % se čte
  z položek seřazených sestupně. U párů jsou členy i záporné (vyhýbání, např. `-n → ch`); kladné členy dosáhnou
  každé úrovně dřív, než je potřeba záporný, počet záporných je v `n_negative`. Pokrytí počítáme z hrubého
  (plug-in) rozkladu; korekce na bias (permutace Y, 50×) je uvedena vedle (`mi_corrected`, `position_chain.csv`)
  a činí 0,005–0,012 bitu, tj. pořadí položek nemění.
* **„Kraj řádku“ vs „drift uvnitř“**: řetězový rozklad I(Y; pos) = I(Y; místo v řádku) + I(Y; řádek v odstavci |
  místo) + I(Y; rel_pos | místo, řádek v odstavci); pos je třída 007 (místo × 5 binů rel_pos × 0/1/2+), takže
  součet je přesný.
* **Dekompozice (a) a (b)** jako v 011: (a) naučené afixy, (b) + poziční glyfy (`p f t k cph cfh cth ckh` /
  `m g am om`). Glyfové cíle (první/poslední glyf) na seznamech nezávisí (variant `-`).
* **Šev bez y→q**: člen nejsilnějšího páru (v Eva y→q, v GC `9→4`-ekvivalent, viz `removed_pair`) je odečten
  z čitatele i z celku; pokrytí zbytku je na zbylé informaci.
* **Test „jen glyfový pár“**: I(sufix_t; prefix_{t+1} | poslední glyf_t, první glyf_{t+1}) s korekcí permutací
  tokenů uvnitř řádků (jako `edge_glyph_mi`).
* **Stabilita**: JSD mezi skupinami vs JSD mezi náhodnými polovinami (řádky jako jednotky) téže skupiny, 200×;
  aby se velikosti vzorků shodovaly, „mezi skupinami“ se počítá též na polovinách (půl A vs půl B).
  `ratio` = JSD mezi / JSD uvnitř; 1 = skupiny se neliší víc než poloviny. Podmíněné tabulky (šev) se srovnávají
  jako vážený průměr JSD řádků tabulky přes kontexty s ≥ 20 páry v obou skupinách. „Podíl pravidel sdílených“ =
  z pointwise členů, které pokryjí 80 % Σ|příspěvek| sdružené tabulky, podíl těch se stejným znaménkem
  (obohacení/vyhýbání) v obou skupinách. Písaři s ≥ 3 000 tokeny (1, 2, 3), sekce s ≥ 2 000 (B, H, P, S, T);
  jazyk A/B.
* **Post-procesor R1–R3** (`apply_place_rules`): tabulky učené na reálu s afixy (b) (jinak by `-m/-g`, `p-/f-`
  nebyly afixy a R2/R1 by je neuměly zapsat). Přepis afixu je *minimální vazba*: token si afix nechá
  s pravděpodobností min(1, q(a)/p(a)), jinak losuje z přebytku (q − p)₊ – mění se jen tolik tokenů, kolik
  je nutné, aby rozdělení v kontextu bylo cílové. Kmen a druhý afix zůstávají. R1: kontexty první token
  řádku × první řádek odstavce (`first|pf`, `first|other`, `inner|pf`); R2: sufix posledního tokenu řádku
  (≥ 2 tokeny); R3: prefix každého vnitřního tokenu podmíněně na posledním glyfu předchozího tokenu
  (kontexty s ≥ 20 páry v reálu; R3 přepisuje po R1). Poslední glyf tokenu se přepisem prefixu nemění
  (kmen má ≥ 1 glyf), takže R3 lze aplikovat vektorově z původního textu. R3 mění jen prefix – první glyf
  tokenu bez prefixu (kmen začíná) nechává, proto doplněný generátor okrajovou MI dohnat zcela nemůže
  (viz report).
* **Profil doplněného generátoru**: stejný kód jako 009 (`M9.profile`), 010 (`M10.line_profile`),
  011 (`M11.measure`, varianty a i b), 30 seedů (stejné seedy ts2020 jako v 009–011, takže referenční ts2020
  je převzat z `results/009..011/runs.csv` a nepočítá se znovu). Navíc kalibrační cíle 004 (H2 glyfů, délka,
  typy, hapax) v `calibration.csv`.
* **Krok 5 (katalog, blábol)**: klony `lrozanova/voynich-units` a `danielgaskell/voynich` v `../lrozanova`,
  `../ext/gaskell` jako v 009/013 (nejsou v repu). Použit Linnaeus po položkách a blábol ve vlastních řádcích.
* Kritérium „dožene“: |z| < 3 v split i join, z = (Voynich − průměr nullu)/SD nullu (konvence 009).
  „Pokazí“: ts2020 mělo |z| < 3 a varianta ≥ 3.

## Odchylky a poznámky

* `tests/test_registers.py::test_transfer_and_mapping_fit` padal už na `main` v tomto prostředí
  (pandas 3.0: `to_numpy()` vrací read-only pole, in-place `/=` selže). Opraveno `copy=True` v
  `voytools/registers.py` (chování beze změny). Lokálně (pandas 2.3) se to neprojevilo.
* Přepsaný token = prefix + kmen + sufix; opětovná dekompozice stejnými seznamy může dát jiný řez
  (`ka`+`m` → `k`+`am`). Profilové metriky dekomponují novým řezem, jako u každého jiného korpusu.
* Prostředí: Python 3.13 / pandas 3.0 / numpy 2.5 (cloud); kód drží 3.10 / pandas 2.3 kompatibilitu.
