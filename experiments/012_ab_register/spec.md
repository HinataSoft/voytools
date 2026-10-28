# 012 – A a B jako registr: změna zápisu, nebo změna obsahu?

## Otázka
Posun sufixů mezi Currierovými jazyky A a B přežívá kontrolu písaře (písař 3 píše obojí, 005) i podmínění
pozicí (008). Žádný null ho nemá. Otázka: **je A→B systematická náhrada tvarů při stejném slovníku kmenů
(změna konvence/zápisu), nebo změna slovníku (jiné kmeny – jiný obsah nebo jiný jazyk)? Je přechod pozvolný
podle pořadí psaní, nebo skokový?**

## Hypotézy a predikce
- **Změna zápisu:** A a B sdílejí kmeny stejně jako dvě poloviny téhož jazyka; posun sufixů lze popsat
  řídkým mapováním (každý sufix A přechází převážně na jeden sufix B); u sdílených kmenů se celé sady sufixů
  posouvají konzistentně; v pořadí psaní může být přechod skokový (rozhodnutí) nebo pozvolný (drift).
- **Změna obsahu/jazyka:** kmeny sdílené výrazně méně než mezi polovinami téhož jazyka; mapování sufixů
  difúzní; sekce uvnitř A (herbář A) a B (herbář B, lázně, hvězdy) mají vlastní slovník nad rámec A/B.
- **Generátor:** ts2020 nemá A/B; jeho „A“ a „B“ stránky (podle layoutu) musí dávat identitní mapování
  a sdílení jako dvě poloviny – to je null pro všechny míry.

## Data a filtry
ZL3b a GC2a, `P`, `clean`, split (hlavní), join (kontrola). Stránky s `language` A / B; uvnitř písaře 3 zvlášť.
Dekompozice naučenými afixy (a) a s pozičními glyfy (b). Pořadí psaní = pořadí stránek v souboru (kvatern,
folio) – uvést, že fyzické pořadí bifolií může být jiné (Layfield & Davis 2026).

## Metriky
1. **Sdílení kmenů:** Jaccard a pokrytí (podíl tokenů A, jejichž kmen se vyskytuje v B, a naopak) pro
   kmeny s ≥ 3 výskyty; srovnání s dělením téhož jazyka na dvě poloviny stejné velikosti (náhodné stránky,
   200×; a souvislé poloviny podle pořadí – kontrola driftu). Celkově, uvnitř písaře 3, a uvnitř herbáře
   (herbář A vs herbář B: stejná sekce, různý jazyk – nejčistší kontrast).
2. **Mapování sufixů:** pro kmeny s ≥ 5 tokeny v A i v B matice M[s_A, s_B] = Σ_kmen p(s_A | kmen, A)·p(s_B | kmen, B)
   (přenos v rámci kmene), normalizovaná po řádcích; entropie řádků (řídkost), diagonální hmota (podíl „beze
   změny“), nejlepší párování (maximální váhové párování bez scipy: Hungarian nebo greedy s uvedením) a
   zbytková JSD sufixů A/B po aplikaci párování. Null: poloviny téhož jazyka (identita + šum) a ts2020.
   Totéž pro prefixy. Seznam 10 největších přesunů (s_A → s_B, váha).
3. **Pozvolnost:** pro každou stránku skóre S = Σ_sufix p_stránka(sufix)·log(p_B(sufix)/p_A(sufix)) (projekce na
   osu A–B, s vynecháním stránky z odhadu p_A/p_B); graf S podle pořadí stránek, obarveno jazykem, písařem,
   sekcí; test driftu uvnitř A (Spearman S × pořadí, permutace stránek 200×) a uvnitř B; velikost skoku na
   hranicích A/B vs rozptyl uvnitř; „přechodné“ stránky (S blízko 0) – seznam.
4. **Co se mění:** příspěvky sufixů k JSD A/B (top 10) a zda jde o glyfové substituce (např. `-ol` ↔ `-edy`,
   `ch-` ↔ `qo-`), tj. zda mapování odpovídá známým A/B rozdílům (Currier: `qo`, `ed`, `chol`).
5. **Sekce uvnitř jazyka:** JSD kmenů a sufixů mezi sekcemi uvnitř B (lázně vs hvězdy vs herbář B) proti
   JSD A/B uvnitř herbáře – je A/B větší rozdíl než rozdíl sekcí?

## Nulové modely
Poloviny téhož jazyka (náhodné a souvislé, 200×), permutace štítku jazyka mezi stránkami uvnitř písaře
(200×), ts2020 (30 seedů, A/B podle layoutu), šablona 005 (30; má téma po sekcích, ne jazyk).

## Výstupy (`results/012_ab_register/`)
`stem_sharing.csv`, `suffix_mapping_ZL.csv` / `_GC.csv` (matice), `mapping_fit.csv`, `transfers_top.csv`,
`page_scores.csv`, `gradient.png`, `contributors.csv`, `sections_within.csv`, `runs.csv`/`meta.json`,
`report.md`, `notes.md`.

## Kritéria přijetí
1. Všechny míry pro ZL i GC, celkově i uvnitř písaře 3 a uvnitř herbáře, s nully.
2. Report odpoví: (i) sdílí A a B kmeny jako poloviny téhož jazyka? (ii) je mapování sufixů řídké
   (změna zápisu) nebo difúzní? kolik JSD zbude po nejlepším párování? (iii) je přechod skokový nebo pozvolný;
   existuje drift uvnitř A/B? (iv) je A/B větší rozdíl než rozdíl sekcí uvnitř B? – ano/ne, ZL/GC.
3. Podezřelé: Jaccard kmenů A/B > Jaccard polovin (chyba), diagonální hmota mapování > 0,95 u reálu
   (žádný posun – rozpor s 005).
4. Testy, Python 3.10 / pandas 2.3, bez scipy; nové funkce do `voytools/registers.py` (ne do `morphology.py`).

## Rozpočet
Tabulkové operace; 200 permutací × 2 soubory × 3 výřezy ≈ 20 min; nully ≈ 20 min. Jedna session, jedna větev, PR.
