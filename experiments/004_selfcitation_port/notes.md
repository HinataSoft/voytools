# 004 – deník session

## Port

* Zdroj: `analysis/reproduce_selfcitation_control.py` z github.com/lrozanova/voynich-units, commit `956a7c4`
  (MIT). Řádky 51–969 (tabulky, `tokenize`, `can_before`/`can_after`, `G`, `SelfCitation`) převzaty doslova do
  `voytools/selfcitation.py`; vypuštěn jen driver (entropie, BPE, pořadí tokenů, Levenshtein, útok, kalibrace).
  Přidáno API `selfcitation_generate`, `selfcitation_template`, `selfcitation_like`.
* „Ligatury rozepsat zpět na Eva“: generátor vrací `G.s`, tedy zřetězení ligaturních jednotek – už je to prostý
  Eva řetězec (písmena a–z). Žádný převod nebyl potřeba; test `test_output_is_basic_eva_on_2000_tokens` ověřuje,
  že `collapse_eva(…, 'extended')` dává zpět týž řetězec.
* Náhodnost: generátor používá `random.Random(seed)` (jako originál), ne numpy – jinak by nešla ověřit shoda
  token po tokenu.

## Validace – zásadní zjištění

* Jejich tři hotové vzorky (`data/controls/selfcitation/selfcite_seed{1,2,3}.txt`) **nejsou** vyrobené s Java
  výchozími parametry (`replace_weights` 30/50/20), ale s výsledkem jejich kalibrace: `p_copy = 0`,
  `replace_weights = (0, 0, 100)` (vždy trojnásobná náhrada). Zjištěno prohledáním jejich kalibrační mřížky
  (16 bodů) jejich vlastní třídou na jejich šabloně (ZL3b z jejich bundlu, jejich čištění tokenů, šablona
  prodloužená cyklením na 37 500 tokenů): bod (0, 0/0/100) dává seed 1 shodně ve všech 4 653 řádcích, ostatní body
  v žádném.
* Náš port s týmiž parametry na téže šabloně reprodukuje **všechny tři** vzorky token po tokenu
  (`results/004_selfcitation_port/validation/identity.csv`). Šablona je v `run.py` (`rt_template`) přepsaná z jejich
  kódu, jejich kód se v `run.py` nespouští.
* Výklad spec „výchozí parametry R&T“: počítám obojí – `port_rt_params` (jejich publikované parametry, tj. to, co
  odpovídá tabulce 7) a `port_java_default` (výchozí hodnoty API = Java). Kritérium 1 hodnotím na obou.
* Jejich vzorky jsou měřeny na prvních 34 268 tokenech (= velikost naší ZL split šablony), stránky a příznaky
  prvních řádků odstavců z jejich šablony, ostatní metadata syntetická (`table_from_lines`). První verze měla
  příznak prvního řádku odstavce syntetický (první řádek každé 25řádkové stránky) a `pos_jsd_first_mid` vyšel
  0,080 místo 0,047 – opraveno.
* `near_ratio` (naše definice: okno 5 tokenů, Levenshtein ≤ 1, míchání v řádcích) dává na jejich vzorcích 1,00,
  tabulka 7 uvádí „edit-1 poměr“ 1,06–1,23 (sousední tokeny, jiná normalizace). Je to rozdíl definic metrik,
  ne portu (jejich vlastní vzorky dávají totéž co port).

## Mřížka

* Spec mřížka A (p_copy × 4 váhy, 5 seedů) a B (× p_same_position, 3 seedy) se počítá jako jedna sada běhů
  (`seed = config_id*100 + rep`): body s p_same_position 28 mají 5 opakování, ostatní 3. Kalibrace A čte
  `rep < 5` u p_same_position 28, kalibrace B `rep < 3` u všech.
* Mimo spec přidána váha (0, 0, 100) – R&T bod. V tabulkách `in_spec_grid = 0`, do výběru calA/calB nevstupuje.
* Cíle obou ztrát = reálný ZL split (34 268 tokenů), generátor na stejné šabloně (stejný počet tokenů).
  R&T kalibrovali na 32 747 tokenech; délka tokenu je v obou případech ve znacích Eva.
