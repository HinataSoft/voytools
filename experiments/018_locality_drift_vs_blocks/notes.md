# 018 – poznámky ze session

## Interpretace spec (kde byla volnost)

* **Stránky:** všechny stránky s ≥ 1 čistým `P` tokenem (ZL split 207, GC split 207), pořadí a vzdálenost = pořadí
  v souboru přes všechny stránky (jako 015). Metadata `quire`, `bifolio` (= `quire/bifolio` jako v 015), `illustration`
  (sekce), `language`, `hand` z tabulky stránek `R.page_table`.
* **Míra (a), hlavní:** pro každou dvojici stránek počet dvojic výskytů téhož vzácného kmene (2–5 výskytů, dekompozice
  (a) = naučené afixy reálu) `N_ij` a jeho očekávání při permutaci tokenů na stránky s vahou ∝ velikosti (015)
  `W_ij = 2 w_i w_j P`. Podobnost množiny dvojic = sdružený poměr Σ N / Σ W (přesně poměr pozorované/očekávané z 015;
  pro d = 1 vychází 1,97 jako v 015). Míry (b) kosinus log(1+f) bez 30 nejčastějších kmenů (006) a (c) Jaccard
  množin kmenů s ≥ 2 výskyty v korpusu (hapaxy nelze sdílet) bez 30 nejčastějších mají `W = 1`. Kosinus a Jaccard
  nesou závislost na velikosti stránky (006); permutace stránek ji sdílí, generátory v layoutu reálu také.
* **Hladký model vzdálenosti f(d):** sdružený poměr po vzdálenostech 1…40 sloučený do binů (1, 2, 3, 4, 5, 6–7, 8–9,
  10–11, 12–15, 16–19, 20–23, 24–31, 32–40, 41+), izotonicky nerostoucí (PAVA vážené Σ W). Biny místo jednotlivých d:
  při polovině vzácných kmenů připadá na d = 1 jen několik dvojic a izotonický fit šumu je na krátkých vzdálenostech
  systematicky nadhodnocený (první blok PAVA = běžící maximum). Reziduum `R = N − f(d) W`.
* **Segmentace (DP):** prostý součet reziduí v blocích (doslovné čtení spec) je degenerovaný – optimum je jeden velký
  blok plus samostatné stránky (singleton nic nestojí, odříznutí záporných dvojic se vyplatí vždy). Použito proto
  **jednostranné deviance skóre bloku** sčítané přes bloky: Poisson pro (a) (`O log(O/E) − (O − E)` pro O > E, jinak 0;
  profilová věrohodnost bloku s vlastní mírou sdílení) a Gauss pro (b), (c) (`(ΣR)²/(2σ²·n_dvojic)` pro ΣR > 0).
  Blok s poměrem ≈ 1 přispívá ≈ 0, takže DP utrácí hranice jen tam, kde je evidence o vyšší podobnosti uvnitř; stejné
  kritérium jako kernel change-point (Harchaoui & Cappé 2007; Arlot, Celisse & Harchaoui 2019). K = 1…30.
* **Zisk bloků:** podle spec sdružené reziduum uvnitř bloků Σ R / Σ W (pro (a) v jednotkách poměru k náhodě: 0,10 =
  o 10 % víc dvojic uvnitř bloků, než dá vzdálenost). **Křížové ověření:** kmeny (sloupce matice stránky × kmeny)
  náhodně rozděleny na poloviny; bloky nalezeny na polovině 1 (reziduum proti f celého textu), zisk měřen na polovině 2
  **proti binovanému sdruženému profilu poloviny 2 bez izotonie** (reziduum má v každém binu nulový součet, takže
  bloky nezávislé na polovině 2 mají očekávaný zisk 0 při libovolné struktuře vzdálenosti). Sdílený model f mezi
  polovinami i izotonický fit poloviny dávaly u permutace stránek zisk systematicky **záporný** (−0,1 při K = 30):
  na krátkých vzdálenostech leží skoro všechny dvojice uvnitř bloků a výběr bloků podle poloviny 1 se promítá do f.
  Reál: 50 rozdělení = 50 řádků (`seed`), interval = kvantily přes rozdělení; nully a generátory: 5 rozdělení na vzorek
  zprůměrováno, interval přes vzorky. Navíc `fit_K` = zisk na všech kmenech bez křížového ověření (vždy ≥ 0, přeučený).
* **K_ref:** počet souvislých úseků sekce v pořadí stránek (ZL P text: 24 úseků, 23 hranic; sekce 7 kódů, ale H a B
  se střídají a T leží uvnitř), u herbáře A počet úseků kvaternů (9). Body změny z plné segmentace (všechny kmeny);
  shoda ± 1 stránka; náhodná základna = K − 1 náhodných hranic (200×).
* **Výřezy:** `all` = celý text; `within` = jen dvojice v téže sekci (maska na W; f i bloky z nich); `herbalA` =
  stránky H & jazyk A (ZL 95, GC 95), vzdálenosti v původním pořadí, váhy a vzácné kmeny přepočítané uvnitř výřezu.
  U kontrol jen `all` a `within` (sekce = kniha/řád/pisatel).
* **Fyzické pořadí:** (a)(b) reziduum R (pro (a) dělené f·W, tj. v jednotkách poměru) dvojic na témž bifoliu / v témž
  kvaternu proti dvojicím ostatním **na stejné vzdálenosti** (vážený průměr po d), d ≤ 15. (c) sousedé = 3 stránky
  téhož jazyka s nejvyšším standardizovaným reziduem `(O − fE)/√(fE)` (a) resp. reziduem (b, c); anomálie = všichni
  tři > 8 stránek. (d) matice kvatern × kvatern = sdružené reziduum; pořadí kvaternů spektrální seriací (Fiedlerův
  vektor Laplaciánu) a |Spearman| s pořadím čtení; podíl kvaternů, jejichž nejpodobnější kvatern je soused v pořadí čtení.
  Permutační interval = zdroj `perm_global` (metadata zůstávají s pozicí, obsah i velikost stránky se stěhují).
* **ts2020-mem:** `selfcitation.py` – v režimu „local“ (zdroj z předchozích řádků téže stránky) s pravděpodobností
  `p_prev_page` zdrojový řádek z dřívější stránky, počet stránek zpět geometrický se střední hodnotou `page_memory`
  (oříznuto na první stránku). Při `p_prev_page = 0` se generátor náhodných čísel nevolá → bitově shodný výstup
  (test `test_memory_defaults_are_bit_identical`). Kalibrace: mřížka 4 × 4, 5 seedů, ztráta = Σ (poměr_bin − reál)² přes
  biny 1 / 2–3 / 4–7 / 8–15 profilu 015 (dekompozice a, kmeny 2–5); vítěz 30 seedů + ztráta 004 (loss B proti
  `results/004_selfcitation_port/real/ZL_split/point.csv`) a profil 015 (biny, přebytek v sekci ≥ 16, tvar opakování).
* **Zdroje a seedy:** real 50 rozdělení; perm_global / perm_section 200; ts2020 cfg 21 seedy 2100+ (30); ts2020-mem
  2300+ (30); šablona 005 s tématem `best` dané kombinace (ZL_split, GC_split; seedy 100000+, pool_seed = rep, 30).
  Kontroly 009 (Plinius, Celsus, Linnaeus po položkách; layout ZL split, vlastní afixy, sekce = kniha / řád), 013
  blábol (60 stránek, body 1–2); real + permutace.

## Prostředí

* Cloud: Python 3.13, pandas 3.0.5, numpy 2.5.3 (vlastník 3.10 / 2.3); kód bez novější syntaxe i API.
* Externí kontroly: `lrozanova/voynich-units` → `../lrozanova/voynich-units`, `danielgaskell/voynich` → `../ext/gaskell`
  (cesty jako v 009/013).
* Časy: reál (1 rozdělení) 1,8 s; permutace stránek (5 rozdělení) 3 s; ts2020 vzorek 3 s + metriky; šablona 18 s.
* Běhy: 4–6 procesů najednou s `OMP_NUM_THREADS=1` (bez toho se vlákna BLAS přetahovala, 8 s na řádek místo 3 s).
  Jeden řetězec úloh přežil ukončení a počítal Celsa/Linnaea souběžně s frontou – duplicitní seedy v `runs.csv`
  (identické hodnoty, ověřeno) byly odstraněny (`drop_duplicates('seed')`); žádný poškozený řádek.
* Nejistota zisku (a): po polovinách kmenů je na d = 1 jen ~7 dvojic, proto biny v f(d) a hodnocení proti
  binovanému profilu poloviny 2; intervaly přes rozdělení jsou široké, `se` průměru 0,02–0,06.
* Kosinus/Jaccard: zisk bloků roste s K stejně u všech zdrojů v témž layoutu (velikost stránek, jazyk podle pozice) –
  pro bloky nepoužitelné, pro fyzické pořadí (bifolio) ano, protože permutace dávají přesně 0.
* `neighbours_list.csv`: u surového poměru (a) mají malé stránky poměry 0 / 20 / 40 a mnoho remíz (`f1r f1v f2r`),
  seznam anomálií proto podle surového kosinu.
