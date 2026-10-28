# 012 – A a B jako registr: změna zápisu, nebo změna obsahu?

Kód: `experiments/012_ab_register/run.py`, `summarize.py`; funkce `voytools/registers.py` (testy
`tests/test_registers.py`). Rozhodnutí a odchylky: `experiments/012_ab_register/notes.md`.

## Co se počítalo

* ZL3b a GC2a, `P`, clean, split (hlavní) a join; dekompozice (a) naučené afixy a (b) + poziční glyfy.
* Výřezy: **celkem**, **herbář** (sekce H: A 95 stránek / 7 779 tokenů, B 32 / 3 353 – nejčistší kontrast),
  **písař 3** (A jen 2 stránky f58r/f58v, 733 tokenů – slabý výřez).
* Všechna srovnání A × B na **stejné velikosti** N = polovina menšího jazyka (celkem 5 497, herbář 1 676,
  písař 3 366 tokenů na stranu), 200 losů celých stránek. Nully: dvě poloviny téhož jazyka (A|A, B|B; náhodné
  stránky a souvislé poloviny podle pořadí), permutace štítku jazyka uvnitř písaře, ts2020 (30) a šablona 005
  (30) s A/B podle layoutu. `AB−poloviny` = A×B minus průměr (A|A, B|B) v témže losu.
* Hodnoty: průměr [2,5 %; 97,5 %] přes losy (reál) nebo seedy (nully).

## (i) Sdílí A a B kmeny jako poloviny téhož jazyka? (`stem_sharing.csv`)

Dekompozice (a), kmeny s ≥ 3 výskyty, split.

| sada / výřez | Jaccard A×B | A\|A | B\|B | A×B − poloviny | souvislé poloviny A / B | pokrytí A→B / B→A (tokeny) | totéž u polovin |
|---|---|---|---|---|---|---|---|
| ZL celkem | 0,422 | 0,616 | 0,601 | **−0,187 [−0,244; −0,124]** | 0,493 / 0,487 | 0,826 / 0,855 | 0,864–0,883 |
| GC celkem | 0,383 | 0,577 | 0,533 | **−0,172 [−0,213; −0,119]** | 0,462 / 0,481 | 0,782 / 0,823 | 0,824–0,844 |
| ZL herbář | 0,492 | 0,626 | 0,697 | **−0,169 [−0,239; −0,080]** | 0,631 / 0,692 | 0,774 / 0,795 | 0,812–0,834 |
| GC herbář | 0,444 | 0,557 | 0,573 | **−0,121 [−0,194; −0,034]** | 0,548 / 0,604 | 0,731 / 0,720 | 0,755–0,764 |
| ZL písař 3 | 0,457 | 0,315 | 0,506 | +0,047 [−0,193; 0,334] | – | 0,577 / 0,594 | – |

* Nully: ts2020 A×B − poloviny −0,017 [−0,047; 0,019], šablona −0,005 [−0,038; 0,028] (ZL celkem) – generátor
  bez jazyka dá podle očekávání nulu.
* Join a dekompozice (b) stejně (−0,17 až −0,19 celkem; herbář −0,13 až −0,17).
* **Odpověď (i): ne.** A a B sdílejí výrazně méně kmenů (inventáře) než dvě poloviny téhož jazyka – i v herbáři,
  tj. při stejné sekci, a i proti souvislým polovinám (které zahrnují posun sekcí). Rozdíl je v řidších kmenech:
  pokrytí tokenů klesá jen o 3–5 p. b. (časté kmeny jsou společné). Uvnitř písaře 3 nelze rozhodnout (A = 2 stránky,
  interval přes nulu).

## (ii) Je mapování sufixů řídké? Kolik JSD zbude po nejlepším párování? (`mapping_fit.csv`, `transfers_top.csv`)

Sufix, celkem, (a), split; kmeny s ≥ 5 tokeny v A i B (ZL průměrně 64 kmenů na los).

| sada | diag. hmota A×B / poloviny | entropie řádků A×B / poloviny | JSD sufixů A×B | po párování | poloviny |
|---|---|---|---|---|---|
| ZL | 0,366 / 0,43–0,45 | 2,36 / 2,34–2,38 | 0,092 | **0,060** | 0,003–0,004 |
| GC | 0,375 / 0,41–0,43 | 2,51 / 2,34–2,41 | 0,069 | **0,036** | 0,003–0,004 |
| ZL herbář | 0,368 / 0,46–0,51 | 1,99 / 2,03–2,21 | 0,109 | 0,064 | 0,006 |
| ZL prefix | 0,251 / 0,32–0,34 | 2,91 / 2,73–2,81 | 0,070 | 0,058 | 0,002–0,004 |

* ts2020 (A/B podle layoutu): JSD sufixů A×B 0,008, poloviny 0,004–0,005 – identita + šum.
* Největší přesuny nad úrovní polovin (ZL; GC stejné páry ve v101): **A `-y` → B `-dy`** (přebytek 0,057 hmoty;
  GC `9 → 89`), **`-ol` → `-dy`** (0,041), **`-or` → `-dy`** (0,028); prefixy `d-` → ∅, ∅ → `ch-`, `cho-` → `o-`,
  `o-` → `qo-`.
* **Odpověď (ii): difúzní, ne řídké.** Řádky mapování nejsou ostřejší než u polovin téhož jazyka (entropie
  2,36 vs 2,34–2,38), diagonála klesá jen o 7 p. b. a nejlepší párování 1 : 1 sníží JSD sufixů jen o 35 % (ZL)
  až 48 % (GC); zbytek 0,036–0,060 bitu je 10–15× nad úrovní polovin. Hlavní posun je „více tvarů A → jeden
  tvar B" (`-y`, `-ol`, `-or` → `-dy`), což jednoznačná náhrada znaků není.

## (iii) Skokový, nebo pozvolný přechod? (`page_scores.csv`, `gradient.csv`, `gradient.png`)

Skóre stránky na ose A–B (sufixy, bez dané stránky), stránky ≥ 20 tokenů (ZL 114 A, 83 B).

| | ZL | GC | ts2020 ZL (30 seedů) | šablona ZL |
|---|---|---|---|---|
| oddělení A/B (rozdíl průměrů / SD) | **4,15** | **3,71** | 0,27 [−0,12; 0,71] | 0,08 |
| skok na hranici A/B / změna uvnitř | **4,41** | **4,51** | 1,14 [0,86; 1,44] | 1,06 |
| chybně zařazené stránky | 7 z 197 | 8 | 86 | 93 |
| drift uvnitř A (Spearman, p) | 0,26 (0,025) | 0,07 (0,45) | 0,25 [0,00; 0,54] | −0,02 |
| drift uvnitř B | −0,07 (0,48) | −0,40 (0,005) | 0,11 | 0,04 |

* Prefixy: drift uvnitř A i B silný v obou sadách (ρ 0,41–0,55, p 0,005) – pořadí stránek je ale spjaté se
  sekcemi (herbář → lázně → farmacie → hvězdy), takže „drift" nelze oddělit od střídání sekcí.
* Chybně zařazené stránky ZL: f25r, f51r, f53v (A, herbář, malé stránky), **f58r, f58v (A, písař 3)**, f86v4 (B, C),
  f88v (A, farmacie). Obě A stránky písaře 3 leží na straně B (+0,09 / +0,17) – A písaře 3 je „B-ovitější" než A
  písaře 1.
* **Odpověď (iii): skokový.** Na hranicích A/B skáče skóre 4,4× víc než mezi sousedními stránkami téhož jazyka
  (nully ~1,1); oddělení 3,7–4,2 SD. Drift uvnitř jazyka je nekonzistentní mezi ZL a GC (A jen ZL, B jen GC)
  a ve velikosti, kterou dává i ts2020 – průkazný není.

## (iv) Je A/B větší rozdíl než rozdíl sekcí uvnitř B? (`sections_within.csv`)

JSD na 1 500 tokenech na stranu (náhodné tokeny, 100×); šumová hladina = dva vzorky téže skupiny.

| | sufixy ZL | sufixy GC | kmeny ZL | kmeny GC |
|---|---|---|---|---|
| A×B uvnitř herbáře (H_A × H_B) | **0,108** | **0,097** | 0,244 | 0,303 |
| sekce uvnitř B (lázně, hvězdy, herbář B, text) | 0,007–0,036 | 0,013–0,030 | 0,157–0,197 | 0,200–0,242 |
| sekce uvnitř A (herbář × farmacie) | 0,029 | 0,013 | 0,248 | 0,282 |
| šum (táž skupina) | 0,003–0,004 | 0,003 | 0,099–0,143 | 0,134–0,189 |

* **Odpověď (iv): ano u sufixů, jen mírně u kmenů.** Rozdíl A/B v sufixech je 3–15× větší než rozdíl mezi
  sekcemi uvnitř B (ZL i GC). U kmenů je A/B v herbáři o něco větší než rozdíly sekcí uvnitř B
  (0,244 vs ≤ 0,197; GC 0,303 vs ≤ 0,242), ale stejně velký jako rozdíl herbáře a farmacie uvnitř A.

## Co se mění (`contributors.csv`)

Největší příspěvky k JSD A/B (ZL, celkem): sufixy `-dy` (42 % JSD; A 4,7 % → B 20 %), ∅, `-ain`, `-ol`, `-or`;
prefixy `cho-` (29 %), `l-`, `d-`, `qo-` (A 9,5 % → B 17,6 %); glyfové bigramy `ed` (13 %; A 0,04 % → B 2,9 %),
`cho`, `dy`, `eed`. To jsou známé Currierovy znaky (`qo`, `ed`, `chol`) – v ZL i GC (v101 `c8`, `1o`, `89`).

## Kontroly podezřelých hodnot

* Jaccard A×B > Jaccard polovin: ne (ve všech sadách a výřezech kromě písaře 3, kde interval obsahuje nulu).
* Diagonální hmota > 0,95: ne (0,37); posun existuje, v souladu s 005/008.
* Permutace štítku uvnitř písaře dává skoro stejná čísla jako reál (např. Jaccard 0,437 vs 0,422): jazyk mají
  oba jen u písaře 3, permutace tedy mění jen jeho stránky – tento null je slabý, uvádím ho pro úplnost.

## Interpretace

A a B se neliší jen zápisem tvarů při stejném slovníku: sdílejí méně kmenů než dvě poloviny téhož jazyka i uvnitř
herbáře, posun sufixů je difúzní (víc tvarů A přechází na `-dy`) a jednoznačné přeznačení vysvětlí méně než
polovinu rozdílu. Přechod je ostrý (stránky jsou A nebo B), což spíš odpovídá dvěma konvencím/pisatelským
okruhům než pozvolnému driftu; generátor bez jazyka (ts2020, šablona) nic z toho nedává.

## Co výsledek neukazuje

* Neodděluje jazyk od písaře: písař 3 je jediný kontrast a jeho A má 2 stránky (výsledky uvnitř písaře 3
  jsou nerozhodné). „Změna obsahu" a „jiný pisatelský okruh se svým slovníkem" se tu nerozliší.
* Pořadí stránek v souboru není pořadí psaní; drift je zmatený se sekcemi.
* Nižší sdílení kmenů může zčásti vzniknout parserem: pokud A/B mění i glyfy uvnitř „kmene" (`ed`), parser
  rozdělí stejné slovo na jiné kmeny. Varianta (b) to nemění, ale jiný parser by mohl.
* Mapování je párování přes kmeny (Σ p·p); nezachytí kontextové náhrady (sufix závislý na kmeni).

## Další krok

* 013: zda lidský blábol ukazuje podobný „registr" mezi pisateli (I(afix; pisatel)) – kolik A/B může udělat
  prostě jiný člověk.
* Rozložit sdílení kmenů podle frekvence (které kmeny jsou jen v A / jen v B) a ověřit, zda jde o glyfové
  varianty téhož kmene (`ed`-varianty) – test parserového artefaktu.
