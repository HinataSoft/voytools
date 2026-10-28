# 018 – Lokalita slovníku: drift, nebo bloky? A fyzické pořadí stránek

## Otázka
015: vzácné kmeny se vracejí s pozvolným poklesem přes 1–15 stránek (poměr k náhodě 2,0 → 1,2), stejně jako latinská
próza; kopírovací generátory a šablona dají jen vrchol na stránce a plochý zbytek (ts2020 bere zdroje jen z téže
stránky). Permutace stránek uvnitř sekce dá plochých 1,3 – reál má lokalitu i **uvnitř** sekcí. Dálkové dvojice sedí
do téže sekce nad náhodu jako u tematicky členěné latiny. Dvě čtení: **drift** (slovník psaní se mění postupně
podél pořadí psaní; kopírování s pamětí přes stránky by to dalo z konstrukce) a **bloky** (stránky seskupené podle
věci: podobnost vysoká uvnitř bloku a skok na hranici bez ohledu na vzdálenost). Otázka: **je lokalita hladká funkce
vzdálenosti, nebo má bloky a skoky; dá ji kopírování s pamětí přes stránky; a sleduje podobnost pořadí čtení, nebo
fyzické jednotky (kvatern, bifolio)?**

## Hypotézy a predikce
- **Drift / paměť:** ts2020 s pamětí přes stránky (kalibrovaný na biny 1–15) reprodukuje celý profil 015 včetně
  lokality uvnitř sekcí; zisk segmentace na bloky nad hladký model vzdálenosti (křížově ověřený) = ts2020-mem;
  body změny nesedí na hranice sekcí lépe než u ts2020-mem.
- **Bloky (téma):** křížově ověřený zisk bloků nad vzdálenost výrazně nad ts2020-mem (jako Plinius s knihami);
  body změny sedí na hranice sekcí / jazyka / kvaternů nad náhodu; ts2020-mem kalibrovaný na 1–15 nedá přebytek
  dálkových dvojic v téže sekci (+0,08) – ten zůstane jako blokový zbytek.
- **Fyzické pořadí:** podobnost sleduje pořadí čtení (stejný bifolio bez přebytku při stejné vzdálenosti); anomálie
  sousedů (stránky, jejichž nejbližší jsou daleko) ne nad permutaci – nebo nad ni, a pak se seznam porovná se známými
  hypotézami o přeskládání (herbář B uprostřed herbáře A, f57, farmaceutické bifolio).

## Data a filtry
ZL3b a GC2a, `P`, clean, split (hlavní), join; dekompozice (a). Stránky v pořadí souboru; `quire`, `bifolio`,
`page_in_quire`, `illustration`, `language`, `hand` z tabulky. Výřezy: celý rukopis, uvnitř sekcí, **herbář A** (jedna
sekce, jeden jazyk, převážně jeden písař – hlavní výřez pro drift vs bloky).

## Metriky
1. **Matice podobnosti stránek** (tři míry): (a) sdílení vzácných kmenů 2–5 výskytů, pozorované/očekávané při
   velikostech stránek (015); (b) kosinus všech kmenů (006); (c) Jaccard kmenů ≥ 2. Hlavní = (a).
2. **Hladký model vzdálenosti:** f(d) = průměr podobnosti po vzdálenosti (izotonicky klesající; d = 1…40, dál
   konstanta), odhadnutý na celém rukopisu a zvlášť uvnitř sekcí; reziduum r_ij = s_ij − f(d_ij).
3. **Bloky:** optimální segmentace posloupnosti stránek na K souvislých bloků maximalizující průměrné reziduum uvnitř
   bloků (dynamické programování, K = 2…30); **křížové ověření**: vzácné kmeny náhodně rozdělit na dvě poloviny,
   bloky najít na polovině 1, zisk (průměrné reziduum uvnitř bloků − 0) měřit na polovině 2; 50 rozdělení.
   Křivka zisku G(K) pro reál vs nully; přebytek nad ts2020-mem při K = počet sekcí, 2×, 3×. Body změny (K = počet
   sekcí a K = 2× počet): shoda s hranicemi sekcí, jazyka, písaře, kvaternů (± 1 stránka; přesnost/úplnost) vs
   náhodné body.
4. **ts2020 s pamětí přes stránky** (`voytools/selfcitation.py`: nové volitelné parametry `p_prev_page` – pravděpodobnost,
   že zdrojový řádek je z předchozí stránky místo z aktuální, a `page_memory` – střední vzdálenost geometrického
   rozdělení stránek zpět; výchozí 0 → výstup **bitově shodný** s dosavadním, stávající testy beze změny; nový test
   shody). Kalibrace: mřížka p_prev_page ∈ {0,1; 0,2; 0,3; 0,5} × page_memory ∈ {2, 4, 8, 16}, 5 seedů, nejmenší
   čtvercová odchylka poměrů binů 1 / 2–3 / 4–7 / 8–15 od reálu (015); u vítěze 30 seedů; zkontrolovat, že loss 004
   zůstane < 0,1 (jinak uvést). Pro vítěze: celý profil 015 (biny, přebytek ≥ 16 v téže sekci, tvar opakování) a body
   2–3.
5. **Fyzické pořadí:** (a) dvojice na témž bifoliu (stejný `quire` + `bifolio`) vs dvojice se stejnou vzdáleností
   v pořadí čtení jinde: přebytek podobnosti; (b) stejný kvatern vs jiný při stejné vzdálenosti (hranice kvaternu
   jako skok?); (c) **anomálie sousedů**: pro každou stránku 3 nejpodobnější stránky téhož jazyka (reziduum);
   podíl se vzdáleností ≤ 2 / ≤ 5 / > 16 vs permutace stránek (200×); seznam stránek, jejichž všichni tři sousedé jsou
   > 8 stránek daleko, se sekcí/jazykem/kvaternem; (d) matice kvatern × kvatern (průměrné reziduum) a z ní pořadí
   kvaternů podle největší podobnosti sousedů – shoduje se s pořadím čtení? (deskriptivně, s permutačním intervalem).

## Nulové modely a kontroly
Permutace stránek globálně a uvnitř sekce (200×); ts2020 config 21 (30) a ts2020-mem (30); šablona 005 s tématem
(30); kontroly: **Plinius ve layoutu ZL s hranicemi knih** (knihy = známé bloky; pozitivní kontrola bodů změny –
016 knihy rozlišuje, hranice knih zaznamenat v tabulce stránek), Celsus, Linnaeus po položkách; blábol krátký (jen
body 1–2, do 8 stránek).

## Výstupy (`results/018_locality_drift_vs_blocks/`)
`similarity_profile.csv`, `segmentation_gain.csv` (zdroj × K × zisk × CI), `changepoints.csv`, `crossval_blocks.csv`,
`memory_calibration.csv`, `physical_order.csv`, `neighbours.csv`, `quire_matrix.csv`, grafy (`distance_profile.png`,
`segmentation_gain.png`, `similarity_matrix.png` – reziduální matice stránek s hranicemi sekcí a kvaternů,
`quire_matrix.png`), `runs.csv`/`meta.json`, `report.md`, `notes.md`. Nový modul `voytools/locality.py` (+ testy);
rozšíření `selfcitation.py` (viz výše).

## Kritéria přijetí
1. Vše pro ZL i GC (ts2020-mem jen ZL), split i join (join alespoň body 1–3), tři výřezy.
2. Report odpoví: (i) reprodukuje ts2020-mem profil 1–15 a zároveň přebytek v téže sekci a lokalitu uvnitř sekcí?
   (ii) je křížově ověřený zisk bloků nad ts2020-mem (a kde stojí Plinius s knihami)? (iii) sedí body změny na
   hranice sekcí/kvaternů nad náhodu? (iv) sleduje podobnost pořadí čtení, nebo bifolia/kvaterny; které stránky
   jsou anomálie? Věta: drift / bloky / obojí, s čísly.
3. Podezřelé: zisk bloků u permutace stránek > 0 nad interval (chyba křížového ověření), ts2020-mem s p_prev_page = 0
   neshodný s ts2020 (chyba rozšíření), Plinius bez bodů změny na hranicích knih (chyba hranic).
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
Matice 207 × 207 a DP: minuty. Mřížka ts2020-mem 16 × 5 + 30 seedů ≈ 40 min; ostatní nully ≈ 20 min; kontroly minuty.
≈ 1,5–2 h.
