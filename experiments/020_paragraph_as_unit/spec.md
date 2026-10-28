# 020 – Odstavec a stránka jako jednotky obsahu

## Otázka
Řádek je uzavřená jednotka pro vazby mezi tokeny (010), slovník má lokalitu v měřítku stránek a listů (006, 015).
Mezi tím leží **odstavec** (ZL: 740 odstavců na 207 stránkách `P`, medián 4 řádky; 159 stránek má ≥ 2 odstavce) a
**stránka**. V próze je odstavec/kapitola tematická jednotka: sdílení slovníku mezi řádky klesá skokem na hranici
odstavce nad to, co dá vzdálenost řádků, a text pokračuje přes hranici stránky. Kopírovací generátor bere zdroje
z předchozích řádků stránky bez ohledu na odstavce (žádný skok), šablona má téma na stránku (žádný skok, ale skok na
hranici stránky), katalog po položkách má skok na každé položce. Otázka: **klesá sdílení kmenů mezi řádky na hranici
odstavce nad to, co dá vzdálenost; „ohlašuje“ první řádek slovník odstavce; a pokračuje slovník přes hranici stránky
jako v próze, nebo se láme jako u generátorů?**

## Hypotézy a predikce
- **Odstavec = jednotka obsahu:** při stejné vzdálenosti řádků je sdílení uvnitř odstavce vyšší než přes hranici
  (skok Δ(d) > 0 nad přelosování hranic a nad ts2020/šablonu); uvnitř odstavce je sdílení ploché (pytel slov),
  ne klesající; přes hranici stránky sdílení při d = 1 ≈ uvnitř stránky (pokračování) – jako Plinius/Celsus s vlastními
  odstavci.
- **Bez jednotky (kopírování / drift):** sdílení je hladká funkce vzdálenosti řádků, hranice odstavce nic nemění
  (Δ ≈ 0, jako ts2020), přes hranici stránky skok dolů (jako ts2020 a šablona).
- **Artefakt:** Δ > 0 i u ts2020 (gallows v prvním řádku odstavce mění kmeny přes parser) – odečíst, uvést.

## Data a filtry
ZL3b a GC2a, `P`, clean, split (hlavní) a join; dekompozice (a); kmeny ≥ 2 výskyty v korpusu; varianta **bez 50
nejčastějších kmenů** (aby skok nedělala funkční slova) a varianta celé tokeny. Řádky = lokusy `P`; `par_index`,
`line_in_par`, `is_par_first_line`, `line_index`, `page_index` z tabulky. Jen stránky s ≥ 2 odstavci pro body 1–2.

## Metriky
1. **Sdílení řádků podle vzdálenosti a hranice:** pro dvojice řádků téže stránky ve vzdálenosti d = 1…8 sdílení
   = počet společných kmenů / očekávání při délkách řádků (očekávání z frekvencí kmenů **na dané stránce**, aby šlo
   o strukturu nad slovník stránky; navíc globální očekávání pro srovnání); zvlášť pro dvojice v témž odstavci a
   přes hranici; **skok** Δ(d) = s_same(d) − s_diff(d), d = 1…4, s intervalem z nullů: (a) přelosování hranic odstavců
   na stránce při zachování počtu a délek odstavců (200×), (b) permutace řádků na stránce (200×).
2. **Uvnitř odstavce:** s(d) pro d = 1…6 jen v témž odstavci – pokles s d (kopírování) nebo plocho (pytel)? Sklon
   s intervalem vs ts2020 a próza.
3. **První řádek:** sdílení prvního řádku odstavce se zbytkem odstavce vs sdílení prostředního řádku se zbytkem
   (při stejných vzdálenostech); a sdílení prvního řádku s **následujícím** odstavcem vs s vlastním (ohlašuje, nebo
   je to jen kraj?). Zvlášť s gallows-prefixy odstraněnými z kmenů (kontrola artefaktu).
4. **Hranice stránky:** poslední k řádků stránky i vs první k řádků stránky i + 1 (k = 1…3; vzdálenost řádků přes
   hranici) vs dvojice ve stejné vzdálenosti uvnitř stránky; poměr přes/uvnitř; zvlášť když stránka i + 1 je téhož
   jazyka/sekce a když ne; zvlášť když poslední odstavec stránky nemá značený konec (`<$>` chybí – pokračuje?).
5. **Tvar opakování** (015 v měřítku řádků): u společných kmenů podíl identických tokenů uvnitř odstavce vs přes
   hranici (vyrovnaná vzdálenost), proti paradigmatu kmene.

## Nulové modely a kontroly
(a) přelosování hranic, (b) permutace řádků na stránce (200×); ts2020 config 21 (30; layout včetně odstavců, mechanismus
je nezná – měří artefakt); šablona 005 s tématem (30); **próza s vlastními odstavci**: Plinius a Celsus vysázené do
layoutu ZL (řádky a stránky jako ZL), ale hranice odstavců z pramene (kapitoly; pokud zdroj kapitoly nenese, použít
hranice vět po ~ 40 tokenech a uvést jako slabou kontrolu) – rozšíření `voytools/controls.py` o volbu „vlastní
odstavce“ (malé, s testem); Linnaeus po položkách (položka = odstavec; pozitivní kontrola skoku); blábol (013), mají-li
přepisy odstavce (ověřit, uvést).

## Výstupy (`results/020_paragraph_as_unit/`)
`line_pair_sharing.csv` (zdroj × d × same/diff × hodnota × null × CI), `boundary_step.csv`, `within_par_decay.csv`,
`first_line.csv`, `page_boundary.csv`, `form_by_boundary.csv`, grafy (`boundary_step.png`, `page_boundary.png`),
`runs.csv`/`meta.json`, `report.md`, `notes.md`. Nový modul `voytools/paragraphs.py` (+ testy); malé rozšíření
`controls.py`.

## Kritéria přijetí
1. Vše pro ZL i GC, split i join, s/bez top-50 kmenů; všechny nully a kontroly (próza s vlastními odstavci povinně).
2. Report odpoví: (i) Δ(d) > 0 nad přelosování a nad ts2020 (ZL i GC)? jak velký proti próze a katalogu?
   (ii) uvnitř odstavce plocho, nebo pokles? (iii) ohlašuje první řádek? (iv) hranice stránky: pokračování, nebo
   zlom – a jak u prózy a u generátorů? Věta: odstavec/stránka jako jednotka obsahu – ano / ne / artefakt.
3. Podezřelé: Δ > 0 u permutace řádků (chyba očekávání), Δ u ts2020 ≥ Δ reálu (artefakt – uvést a odečíst),
   Linnaeus bez skoku (chyba kontroly).
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
Tabulkové operace nad 4 000 řádky; nully 200× ≈ 10 min; generátory 60 vzorků ≈ 20 min; kontroly minuty. ≈ 1 h.
