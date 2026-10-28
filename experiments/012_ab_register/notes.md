# 012 – deník session

## Rozhodnutí a interpretace spec

* **Kód.** `voytools/registers.py` (nové, ne v `morphology.py`): `draw_tokens`, `split_pages`, `stem_sharing`,
  `transfer_matrix`, `max_weight_matching` (Hungarian, O(n³), bez scipy; test proti hrubé síle),
  `mapping_fit`, `jsd`/`jsd_contributions`, `page_scores`, `spearman`, `order_drift`, `boundary_jumps`.
  Testy `tests/test_registers.py` (mj. syntetická „změna zápisu" x→y, y→z, z→x: párování ji najde a JSD po
  párování ≈ 0).
* **Stejná velikost.** Sdílení kmenů i mapování silně závisí na velikosti vzorku, proto všechna srovnání
  (A×B, poloviny A, poloviny B, permutace štítků, nully) běží na **N tokenech na stranu**, N = polovina
  menšího jazyka v daném výřezu (celkem ZL 5 497, herbář 1 676, písař 3 366). Celé stránky v náhodném pořadí,
  poslední stránka oříznuta na přesně N. 200 losů. Plné A × plné B (nestejné velikosti) a souvislé poloviny
  podle pořadí stránek jsou jednorázové body (`ABfull`, `AAcontig`, `BBcontig`).
* **Písař 3** má v A jen 2 stránky (733 tokenů, sekce S); polovina = 366 tokenů → mapování na velmi málo
  kmenech (uvedeno v `n_stems`). Je to nejslabší výřez.
* **Permutace štítku jazyka uvnitř písaře:** písaři 1 (jen A), 2, 4, 5 (jen B) mají jediný jazyk, takže
  permutace mění jen stránky písaře 3; null je proto skoro totožný s reálem a jeho informační hodnota je
  malá (uvedeno v reportu). Spec ho žádá, ponechán.
* **Mapování:** M[a, b] = Σ_kmen p(a | kmen, A) p(b | kmen, B) přes kmeny s ≥ 5 tokeny na obou stranách
  (neváženě, podle spec); diagonální hmota u dvou polovin téhož jazyka *není* 1 (volba sufixu má i v rámci
  jazyka 1,8 bitu, 011), proto je srovnávacím bodem vždy null „poloviny". Párování = maximální váhové
  párování (Hungarian), zbytková JSD = JSD rozdělení sufixů A po přeznačení podle párování vs B.
* **Přesuny** (`transfers_top.csv`): mimodiagonální buňky s největším přebytkem průměrné normalizované
  matice A→B nad průměrem polovin (50 losů na každou stranu, stejné N).
* **Pozvolnost:** skóre stránky = Σ p_str(s) log2(p_B(s)/p_A(s)), p_A, p_B bez dané stránky (add-0,5),
  stránky ≥ 20 tokenů. Pořadí = pořadí stránek v souboru (kvaternion, folio); fyzické pořadí bifolií
  může být jiné (Layfield & Davis 2026). Pořadí je silně spjaté se sekcí, takže „drift" uvnitř jazyka
  nelze oddělit od střídání sekcí.
* **Sekce:** JSD kmenů a sufixů mezi skupinami sekce × jazyk s ≥ 1 500 tokeny (B_B lázně, S_B hvězdy,
  H_B, T_B text, H_A, P_A farmacie), náhodné vzorky 1 500 tokenů (bez stránkové struktury), 100×;
  „same_group" = dva disjunktní vzorky téže skupiny (šumová hladina).
* **Nully:** ts2020 (ZL split/join) a šablona 005 (ZL/GC split), 30 seedů, generátory 009; vzorky se cachují
  v `scratch/012_cache/` (gitignored). A/B = štítky layoutu (generátor jazyk nezná).

## Průběh

* První spuštění: bodový řádek (seed −1) měl jiné sloupce než losy → `run_seeds` odmítl; bod se teď píše do
  `runs/point__<sada>/`, losy do `runs/real__<sada>/`.
