# 007 – Pozice versus obsah: report

Spec: `experiments/007_position_vs_content/spec.md`; kód `experiments/007_position_vs_content/run.py`
(+ `summarize.py`); rozhodnutí a odchylky `experiments/007_position_vs_content/notes.md`.
Všechna čísla jsou z `info_decomposition.csv`, `mi_distance.csv`, `mi_summary.csv` a `drift.csv`
v této složce (index běhů `runs.csv`, parametry `meta.json`, jednotlivé běhy v `decomp/`, `mi/`, `drift/`).

**Co běželo.** ZL3b a GC2a, `P` text `clean`, politiky `split` a `join` (4 sady). Reál: bodový odhad
z celého textu + 200 subsamplů po stránkách (80 % bez vracení). Intervaly reálu jsou 95%, odchylky
subsamplů škálované √(m/(n−m)) = 2. Nully po 30 seedech: `autocopy_like` config 86 z 003 (seedy
8600–8629, stejné vzorky jako 003 phase B), míchání tokenů v řádku a nové
`nulls.shuffle_within_position_class` (míchání uvnitř poziční třídy v odstavci). Intervaly nullů jsou
2,5–97,5 % přes seedy. **Vynechané nully:** `selfcitation_like` (004) a šablona 005 nejsou v `main`
(stav 7d6a7a9), tedy chybí oba nully, které spec podmiňuje sloučením.

**Poziční třída** = místo v řádku (první / vnitřní / poslední / jediný token) × `rel_pos` v 5 binech ×
`line_in_par` 0 / 1 / 2+ (22 obsazených tříd). Sekce = `illustration` (7), písař = `hand` (5),
jazyk = `language`. Všechny MI jsou plug-in v bitech s korekcí mícháním: odečten průměr hodnoty po
permutaci Y uvnitř vrstev podmínky. „Korigovaná“ hodnota je tedy přebytek nad tím, co dává
konečný vzorek.

## Kritérium 3 – podezřelé hodnoty

Žádné. I(Y; pozice) > H(Y): 0 případů; korigovaná MI < −0,01 bitu: 0 případů v rozkladu a 0 v
`mi_distance.csv` (`meta.json`, `flags`). Nejzápornější korigované hodnoty MI(d) leží kolem −0,003
(bod) a dolní meze intervalů kolem −0,01.

## 1. Rozklad informace (kritérium 1)

Tabulky: korigovaná MI v bitech, reál [95% interval] / autocopy 86 (průměr). `I_sec_g_pos` =
I(Y; sekce | pozice) atd. Podíl z H(Y) je ve sloupci `share` v `info_decomposition.csv`
a v grafu `decomposition.png`.

Permutační testy bodového odhadu (200 permutací): poziční členy proti Y permutovanému uvnitř
stránky, stránkové členy proti permutaci záznamů (sekce, písař, jazyk) mezi stránkami.
Všech 128 členů má p ≤ 0,03, 125 z nich p = 0,005 (minimum). Permutace uvnitř stránky stránkové
členy z konstrukce nemění; proto druhý null (`notes.md`).

**ZL_split** – korigovaná MI v bitech, reál [95% interval] / autocopy 86 (průměr 30 seedů)

| Y | H(Y) | I_pos | I_sec | I_hand | I_lang | I_sec_g_pos | I_pos_g_sec | I_sec_g_hand | I_hand_g_sec |
|---|---|---|---|---|---|---|---|---|---|
| prefix | 3,372 | 0,134 [0,125; 0,142] / -0,000 | 0,085 [0,072; 0,101] / 0,004 | 0,073 [0,059; 0,088] / 0,003 | 0,061 [0,050; 0,073] / 0,001 | 0,104 [0,089; 0,121] / 0,004 | 0,152 [0,141; 0,163] / -0,000 | 0,031 [0,023; 0,041] / 0,005 | 0,019 [0,009; 0,026] / 0,004 |
| suffix | 3,277 | 0,048 [0,042; 0,055] / 0,000 | 0,070 [0,057; 0,084] / 0,004 | 0,088 [0,073; 0,105] / 0,003 | 0,076 [0,062; 0,092] / 0,002 | 0,080 [0,066; 0,094] / 0,003 | 0,058 [0,050; 0,066] / -0,000 | 0,020 [0,014; 0,028] / 0,005 | 0,037 [0,026; 0,048] / 0,004 |
| stem | 5,542 | 0,191 [0,179; 0,203] / 0,001 | 0,119 [0,101; 0,136] / 0,025 | 0,104 [0,087; 0,122] / 0,017 | 0,080 [0,064; 0,093] / 0,009 | 0,117 [0,098; 0,132] / 0,012 | 0,192 [0,181; 0,205] / -0,000 | 0,055 [0,042; 0,072] / 0,025 | 0,039 [0,028; 0,049] / 0,017 |
| length | 2,746 | 0,024 [0,020; 0,027] / 0,000 | 0,022 [0,016; 0,029] / 0,004 | 0,021 [0,014; 0,027] / 0,003 | 0,016 [0,010; 0,021] / 0,001 | 0,027 [0,020; 0,035] / 0,004 | 0,029 [0,025; 0,033] / 0,000 | 0,006 [0,003; 0,009] / 0,005 | 0,005 [0,002; 0,008] / 0,004 |

**GC_split** – korigovaná MI v bitech, reál [95% interval] / autocopy 86 (průměr 30 seedů)

| Y | H(Y) | I_pos | I_sec | I_hand | I_lang | I_sec_g_pos | I_pos_g_sec | I_sec_g_hand | I_hand_g_sec |
|---|---|---|---|---|---|---|---|---|---|
| prefix | 3,307 | 0,097 [0,090; 0,105] / 0,000 | 0,072 [0,060; 0,084] / 0,004 | 0,064 [0,053; 0,077] / 0,002 | 0,054 [0,045; 0,066] / 0,001 | 0,086 [0,073; 0,099] / 0,004 | 0,111 [0,103; 0,118] / 0,000 | 0,026 [0,020; 0,033] / 0,005 | 0,019 [0,010; 0,025] / 0,003 |
| suffix | 3,246 | 0,032 [0,028; 0,037] / -0,000 | 0,049 [0,039; 0,059] / 0,004 | 0,067 [0,055; 0,080] / 0,003 | 0,055 [0,044; 0,067] / 0,001 | 0,056 [0,044; 0,066] / 0,003 | 0,040 [0,035; 0,045] / -0,000 | 0,016 [0,010; 0,020] / 0,005 | 0,034 [0,024; 0,043] / 0,004 |
| stem | 5,485 | 0,176 [0,165; 0,186] / 0,001 | 0,112 [0,098; 0,128] / 0,025 | 0,099 [0,083; 0,115] / 0,017 | 0,071 [0,056; 0,084] / 0,008 | 0,107 [0,092; 0,122] / 0,012 | 0,172 [0,161; 0,183] / 0,001 | 0,052 [0,039; 0,065] / 0,024 | 0,038 [0,028; 0,047] / 0,016 |
| length | 2,592 | 0,031 [0,027; 0,035] / 0,000 | 0,010 [0,007; 0,013] / 0,003 | 0,015 [0,011; 0,020] / 0,002 | 0,010 [0,006; 0,015] / 0,001 | 0,016 [0,011; 0,021] / 0,002 | 0,036 [0,031; 0,041] / -0,000 | 0,004 [0,002; 0,006] / 0,004 | 0,009 [0,003; 0,012] / 0,003 |

**ZL_join** – korigovaná MI v bitech, reál [95% interval] / autocopy 86 (průměr 30 seedů)

| Y | H(Y) | I_pos | I_sec | I_hand | I_lang | I_sec_g_pos | I_pos_g_sec | I_sec_g_hand | I_hand_g_sec |
|---|---|---|---|---|---|---|---|---|---|
| prefix | 3,413 | 0,147 [0,135; 0,156] / 0,000 | 0,089 [0,076; 0,104] / 0,005 | 0,076 [0,062; 0,092] / 0,003 | 0,063 [0,052; 0,075] / 0,002 | 0,109 [0,093; 0,127] / 0,004 | 0,166 [0,153; 0,178] / 0,000 | 0,033 [0,024; 0,042] / 0,006 | 0,020 [0,009; 0,027] / 0,005 |
| suffix | 3,337 | 0,051 [0,044; 0,059] / 0,000 | 0,075 [0,062; 0,090] / 0,004 | 0,093 [0,078; 0,110] / 0,003 | 0,080 [0,066; 0,098] / 0,002 | 0,086 [0,072; 0,101] / 0,004 | 0,061 [0,053; 0,069] / 0,000 | 0,021 [0,015; 0,030] / 0,006 | 0,039 [0,028; 0,050] / 0,004 |
| stem | 5,516 | 0,196 [0,186; 0,209] / 0,000 | 0,122 [0,101; 0,139] / 0,028 | 0,109 [0,092; 0,129] / 0,019 | 0,084 [0,067; 0,100] / 0,009 | 0,116 [0,098; 0,133] / 0,012 | 0,194 [0,181; 0,209] / -0,000 | 0,055 [0,042; 0,071] / 0,026 | 0,041 [0,031; 0,051] / 0,017 |
| length | 2,773 | 0,026 [0,021; 0,029] / 0,000 | 0,025 [0,018; 0,033] / 0,004 | 0,024 [0,017; 0,031] / 0,003 | 0,018 [0,012; 0,024] / 0,001 | 0,031 [0,023; 0,040] / 0,004 | 0,032 [0,027; 0,036] / 0,000 | 0,007 [0,004; 0,010] / 0,005 | 0,006 [0,002; 0,009] / 0,004 |

**GC_join** – korigovaná MI v bitech, reál [95% interval] / autocopy 86 (průměr 30 seedů)

| Y | H(Y) | I_pos | I_sec | I_hand | I_lang | I_sec_g_pos | I_pos_g_sec | I_sec_g_hand | I_hand_g_sec |
|---|---|---|---|---|---|---|---|---|---|
| prefix | 3,333 | 0,129 [0,120; 0,137] / 0,000 | 0,076 [0,063; 0,089] / 0,004 | 0,069 [0,055; 0,082] / 0,003 | 0,059 [0,048; 0,071] / 0,001 | 0,093 [0,078; 0,108] / 0,003 | 0,145 [0,136; 0,156] / 0,000 | 0,027 [0,020; 0,036] / 0,005 | 0,020 [0,011; 0,026] / 0,004 |
| suffix | 3,293 | 0,044 [0,039; 0,050] / 0,000 | 0,054 [0,045; 0,066] / 0,004 | 0,074 [0,061; 0,089] / 0,003 | 0,064 [0,051; 0,077] / 0,001 | 0,065 [0,052; 0,077] / 0,003 | 0,054 [0,048; 0,062] / -0,000 | 0,017 [0,011; 0,023] / 0,005 | 0,037 [0,026; 0,048] / 0,004 |
| stem | 5,503 | 0,171 [0,161; 0,182] / 0,000 | 0,120 [0,102; 0,135] / 0,025 | 0,105 [0,089; 0,121] / 0,017 | 0,076 [0,059; 0,092] / 0,008 | 0,110 [0,095; 0,122] / 0,011 | 0,164 [0,153; 0,176] / -0,001 | 0,056 [0,043; 0,069] / 0,025 | 0,040 [0,030; 0,051] / 0,017 |
| length | 2,598 | 0,031 [0,027; 0,035] / 0,000 | 0,014 [0,009; 0,019] / 0,004 | 0,017 [0,013; 0,023] / 0,003 | 0,013 [0,009; 0,018] / 0,001 | 0,019 [0,014; 0,025] / 0,003 | 0,036 [0,030; 0,041] / -0,000 | 0,004 [0,002; 0,007] / 0,005 | 0,008 [0,003; 0,011] / 0,003 |

**Věty ke kritériu 1** (platí ve všech čtyřech sadách, není-li řečeno jinak):

* **Prefix:** pozice > sekce ≥ písař (ZL split 0,134 vs 0,085 vs 0,073 bitu; GC split 0,097 vs 0,072
  vs 0,064). Sekce zůstává i po podmínění pozicí (0,086–0,109) a pozice po podmínění sekcí
  (0,111–0,166). Kopírovací null nereprodukuje nic: poziční informace ≈ 0 (|·| ≤ 0,0002), sekční
  0,004 bitu (≈ 5 % reálu).
* **Sufix:** písař > sekce > pozice (ZL split 0,088 / 0,070 / 0,048; GC split 0,067 / 0,049 / 0,032).
  Sufix nese víc informace o písaři/jazyce než o místě v řádku a I(sufix; písař | sekce) = 0,034–0,039
  je dvojnásobná proti prefixu (0,019–0,020) a srovnatelná s kmenem (0,038–0,041). Kopírovací null: pozice ≈ 0, sekce 0,004 (≈ 6–8 % reálu).
* **Kmen (top 200):** pozice > sekce > písař (ZL split 0,191 / 0,119 / 0,104; GC split 0,176 / 0,112 /
  0,099). I(kmen; sekce | písař) = 0,052–0,056, tedy sekce nese o kmeni informaci i v rámci písaře.
  Kopírovací null dává sekční informaci 0,025 a I(·; sekce | písař) 0,024–0,026, tedy asi 20 %
  resp. 45 % reálu, čistě z lokálního driftu kopírování. Poziční informaci kmene nedává (≈ 0).
* **Délka tokenu:** v ZL pozice ≈ sekce ≈ písař (0,024 / 0,022 / 0,021; join 0,026 / 0,025 / 0,024),
  v GC pozice > písař > sekce (0,031 / 0,015 / 0,010). Rozdíl ZL/GC je zčásti věcí abecedy (délka
  se počítá v glyfech dané transliterace). Kopírovací null: pozice ≈ 0, sekce 0,003–0,004.
* **Hypotézy ze spec:** „kmen: sekce > pozice“ **neplatí** v žádné sadě (pozice nese o kmeni o 40–60 %
  víc). „Afixy: pozice ≥ sekce“ platí pro **prefix**, pro **sufix ne**. Hypotéza o kopírovacím nullu
  (kmen ani afix nenesou sekční informaci nad rámec sousedství) se potvrzuje pro null, ne pro reál:
  reálná sekční informace kmenů je 4,4–4,8× a afixů 12–21× větší než to, co vyrobí lokální kopírování.
  Šablonový null 005 neběžel.

## 2. Reziduální MI na vzdálenost d = 1–10 (kritérium 2)

Graf `mi_distance.png` (korekce `pc`). Tabulka: korekce `pc` = odečten průměr ze 100 (bod) / 20
(subsample, nully) míchání uvnitř poziční třídy v odstavci; `line` = míchání v řádku. „max d“ =
největší d, do kterého má korigovaná MI dolní mez intervalu > 0 při všech menších d. Plocha = součet
d = 2…10. Plná tabulka, včetně citlivostních korekcí `pcpage` (stejné třídy, celá stránka) a
`pcfine` (jemné biny), je v `mi_summary.csv` a po jednotlivých d v `mi_distance.csv`.

| sada | U | reál pc: MI d=1 | d=2 [int.] | d=3 dolní mez | plocha 2–10 [int.] | max d s dolní mezí > 0 | reál line: plocha 2–10 | autocopy pc: d=1 / plocha | míchání v řádku pc: plocha | poz. třída pc: plocha |
|---|---|---|---|---|---|---|---|---|---|---|
| ZL_split | edge | 0,143 | -0,000 [-0,002; 0,003] | -0,004 | -0,002 [-0,008; 0,005] | 1 | -0,003 [-0,011; 0,004] | 0,000 / 0,001 [-0,003; 0,004] | 0,004 | -0,000 |
| ZL_split | bpe | 0,194 | 0,004 [-0,002; 0,011] | -0,003 | 0,009 [-0,009; 0,031] | 1 | -0,008 [-0,027; 0,015] | 0,001 / 0,004 [-0,004; 0,015] | 0,025 | 0,001 |
| ZL_split | affix | 0,123 | 0,010 [-0,000; 0,019] | -0,006 | 0,043 [0,012; 0,075] | 1 | -0,006 [-0,032; 0,025] | 0,003 / 0,038 [0,026; 0,049] | 0,032 | 0,002 |
| ZL_split | stem | 0,051 | 0,007 [0,001; 0,014] | -0,007 | 0,031 [0,013; 0,050] | 2 | -0,002 [-0,021; 0,018] | 0,002 / 0,029 [0,018; 0,044] | 0,028 | 0,002 |
| GC_split | edge | 0,180 | -0,002 [-0,006; 0,002] | -0,004 | -0,005 [-0,014; 0,006] | 1 | -0,007 [-0,016; 0,004] | 0,000 / 0,003 [-0,002; 0,008] | 0,004 | -0,001 |
| GC_split | bpe | 0,222 | 0,001 [-0,006; 0,011] | -0,007 | 0,008 [-0,013; 0,033] | 1 | -0,000 [-0,023; 0,022] | 0,001 / 0,007 [-0,005; 0,019] | 0,020 | 0,002 |
| GC_split | affix | 0,115 | 0,005 [-0,004; 0,014] | -0,006 | 0,015 [-0,013; 0,044] | 1 | -0,011 [-0,036; 0,020] | 0,003 / 0,029 [0,015; 0,047] | 0,018 | -0,000 |
| GC_split | stem | 0,049 | 0,006 [0,000; 0,013] | -0,005 | 0,023 [0,006; 0,042] | 2 | -0,001 [-0,019; 0,021] | 0,002 / 0,020 [0,007; 0,031] | 0,018 | 0,000 |
| ZL_join | edge | 0,113 | 0,000 [-0,002; 0,003] | -0,005 | 0,003 [-0,004; 0,011] | 1 | 0,001 [-0,006; 0,008] | -0,000 / 0,002 [-0,003; 0,006] | 0,003 | -0,001 |
| ZL_join | bpe | 0,157 | 0,005 [-0,003; 0,013] | -0,005 | 0,016 [-0,005; 0,041] | 1 | -0,003 [-0,025; 0,018] | 0,001 / 0,006 [-0,009; 0,016] | 0,026 | -0,003 |
| ZL_join | affix | 0,106 | 0,009 [-0,001; 0,019] | -0,003 | 0,043 [0,008; 0,073] | 1 | -0,008 [-0,043; 0,023] | 0,003 / 0,035 [0,021; 0,053] | 0,031 | 0,002 |
| ZL_join | stem | 0,040 | 0,007 [0,001; 0,013] | -0,005 | 0,025 [0,000; 0,046] | 2 | -0,011 [-0,031; 0,012] | 0,002 / 0,029 [0,009; 0,045] | 0,029 | -0,002 |
| GC_join | edge | 0,141 | -0,003 [-0,006; 0,001] | -0,005 | 0,002 [-0,007; 0,010] | 1 | -0,001 [-0,010; 0,009] | 0,000 / 0,002 [-0,002; 0,008] | 0,004 | 0,000 |
| GC_join | bpe | 0,176 | 0,004 [-0,002; 0,012] | -0,010 | 0,020 [-0,003; 0,043] | 1 | 0,008 [-0,016; 0,028] | 0,001 / 0,008 [-0,001; 0,021] | 0,022 | -0,000 |
| GC_join | affix | 0,100 | 0,011 [0,000; 0,020] | -0,002 | 0,027 [0,001; 0,056] | 2 | -0,013 [-0,047; 0,016] | 0,002 / 0,032 [0,015; 0,047] | 0,019 | -0,000 |
| GC_join | stem | 0,044 | 0,004 [-0,004; 0,010] | -0,007 | 0,023 [0,004; 0,045] | 1 | 0,001 [-0,024; 0,022] | 0,001 / 0,022 [0,004; 0,035] | 0,022 | 0,000 |

**Polovina poklesu.** U reálu klesne korigovaná křivka pod polovinu hodnoty při d = 1 už mezi
d = 1 a 2 ve všech 16 kombinacích. Interpolovaná hodnota 1,49–1,61 (intervaly 1,45–1,73) je jen
artefakt lineární interpolace: při d = 2 zbývá 0–19 % hodnoty z d = 1. U autocopy není polovina
poklesu dobře definovaná, protože jeho hodnota při d = 1 je jen 0,000–0,003 bitu a křivka je plochá
(medián 1,8–4,8, intervaly 1,1 až > 10).

**Odpověď kritéria 2: ne.** Po odečtení pozičních efektů reál **nedrží** závislost dál než na d = 2,
a to v ZL ani v GC. Při d = 3 je dolní mez korigované MI záporná ve všech 16 kombinacích (U × sada)
pro obě korekce ze spec (`pc` i `line`). Kladná dolní mez při d = 2 se objeví jen u kmenů (s `pc`: ZL split,
ZL join, GC split; s `line`: ZL split), u afixových tříd (GC join, `pc`) a jednou u BPE (GC join,
`line`); velikost je 0,004–0,011 bitu. Hrana poslední→první glyf má závislost **jen** při d = 1
(0,11–0,21 bitu), BPE jednotky s jedinou výjimkou také (0,16–0,26 bitu).
Plocha d = 2…10 je s korekcí `pc` kladná u kmenů (0,023–0,031 bitu, všechny 4 sady) a u afixů
(3 sady), ale **stejně velká je u autocopy (0,020–0,038) a u textu míchaného v řádku (0,018–0,032)**.
Je to tedy plošný příspěvek toho, že tokeny jednoho řádku/okolí si jsou podobné, ne pomalu klesající
sekvenční závislost. S korekcí `line` je plocha reálu ve všech 16 kombinacích nulová (interval
obsahuje 0). Predikce spec („reál klesá pomaleji než kopírovací null“) se obrací: reál má ostrý vrchol
při d = 1 a pak nic, kopírovací null má nízké ploché plato (okno 32 tokenů) bez vrcholu.

**Citlivost** (`mi_summary.csv`). S korekcí na celé stránce (`pcpage`) drží reál dolní mez > 0
do d = 2–5. Stejně nebo víc ji ale drží i míchání v řádku (d = 7–9) a autocopy (d = 10) a plocha reálu
(kmen 0,062–0,081, afix 0,083–0,116) je srovnatelná s mícháním v řádku a menší než u autocopy. Tato korekce tedy přičítá k „posloupnosti“
i shlukování na úrovni odstavce/řádku. Korekce `pc` zamrazí 41–43 % tokenů (v řádcích 0 a 1 odstavce
~70 %, `notes.md`), takže je spíš konzervativní; `pcpage` zamrazí 12–13 %. Závěr „nic za d = 2, co
by nevysvětlilo složení řádku“ platí pro obě.

## 3. Drift glyfů (Feaster) jako kovariát

`drift.csv`: podíly 12 nejčastějších glyfů (+ other) podle `rel_pos`, `rel_pos` jen u vnitřních tokenů
(`rel_pos_mid`) a `line_in_par`, plus korigovaná I(glyf; bin):

| sada | rel_pos | rel_pos_mid | line_in_par | autocopy (všechny osy) |
|---|---|---|---|---|
| ZL split | 0,0080 [0,0070; 0,0092] | 0,0037 [0,0029; 0,0047] | 0,0032 [0,0023; 0,0044] | ≤ 0,0001 |
| ZL join | 0,0084 [0,0073; 0,0094] | 0,0035 [0,0027; 0,0045] | 0,0032 [0,0023; 0,0044] | ≤ 0,0001 |
| GC split | 0,0072 [0,0062; 0,0083] | 0,0039 [0,0030; 0,0050] | 0,0025 [0,0018; 0,0034] | ≤ 0,0001 |
| GC join | 0,0070 [0,0059; 0,0080] | 0,0035 [0,0028; 0,0047] | 0,0025 [0,0018; 0,0034] | ≤ 0,0001 |

Drift doprava existuje i bez prvního a posledního slova řádku: zhruba polovina poziční informace o glyfech
(`rel_pos_mid`) zůstává. V ZL split u vnitřních tokenů od začátku ke konci řádku roste `a`
0,075 → 0,104, `d` 0,072 → 0,087, `t` 0,025 → 0,044 a `r` 0,045 → 0,053; klesají `e` 0,086 → 0,064
a `ee` 0,043 → 0,029. Dolů po odstavci se odlišuje hlavně řádek 0 (víc „other“, tj. mimo jiné p/f
gallows: 0,160 vs 0,119–0,130; méně `k`: 0,050 vs 0,068–0,070). Řádky 1 až 5+ jsou si podobné.
Generátor drift nemá.
**Kovariát:** korekce jemnými biny (`pcfine`: 10 binů `rel_pos` × `line_in_par` 0..5/6+, na stránce)
sníží plochu d = 2…10 reálu oproti `pcpage` jen o 5–20 % (např. ZL split kmen 0,081 → 0,071, afix
0,116 → 0,103) a dolní mez > 0 drží do d = 2–5 stejně jako u `pcpage`. Zbytková závislost tedy
driftem vysvětlená není. Stejně jako u `pcpage` ji ale vysvětluje složení řádku (míchání v řádku má
stejnou plochu: kmen 0,055–0,070, afix 0,072–0,095).

## Interpretace

Pozice v řádku a odstavci nese o tvaru tokenu víc informace než sekce, a to i o kmeni. Prefix a kmen
jsou nejvíc poziční, sufix nejvíc vázaný na písaře/jazyk, a vše platí shodně v ZL i GC a v obou
politikách mezer. Sekce nese o kmenech i afixech informaci nezávislou na pozici i na písaři,
několikanásobně větší, než dává lokálně kopírující generátor. Sekvenční závislost mezi tokeny se po
odečtení pozičních efektů omezuje na bezprostředního souseda (d = 1, hlavně hrana slovo–slovo). Vše
za d = 2 je stejně velké jako u nullů a vysvětluje ho složení řádků, ne pořadí. To odpovídá spíš
lokálnímu (generátorovému či „blábolovému“) než jazykovému profilu, jak jej formulovala spec. Pozor
ale, že reál se od kopírovacího nullu liší opačně, než spec předpokládala (ostrý vrchol, ne pomalý
pokles).

## Co výsledek neukazuje

* Neukazuje, že by text byl bez obsahu: sekční informace kmenů (nezávislá na písaři) je reálná a
  kopírovací null ji reprodukuje jen z ~20 %. Absence závislosti na d ≥ 3 na úrovni glyfů, afixů a
  kmenů s ≥ 20 výskyty nevylučuje závislosti mezi vzácnými slovy ani syntax, kterou naše jednotky
  nezachytí (BPE-64 a afixové rámce jsou hrubé).
* „Pozice > sekce“ u kmene může zčásti plynout z toho, že naučené afixy neodstraní poziční glyfy
  (p/f gallows v řádku 0, `-m`/`-g` na konci řádku) a ty pak zůstanou v kmeni. Tabulka tedy měří
  „kmen“ jako zbytek po greedy rozkladu, ne morfém. Poziční třída má 22 hodnot proti 7 sekcím.
  Korekce mícháním odstraní zkreslení konečným vzorkem, ne rozdíl v kapacitě proměnné.
* Písař a jazyk jsou zaměnitelní (I_lang < I_hand všude, ale ze stejných stránek), sekce s písařem také.
  Podmíněné členy je oddělují jen částečně.
* Chybí oba „silné“ nully spec (004 selfcitation, 005 šablona s tématem). Srovnání s generátorem
  je zatím jen proti autocopy 86, který žádný poziční mechanismus nemá, takže jeho nulová poziční
  informace je z konstrukce.
* Korekce `pc` je konzervativní (zamrazí ~41 % tokenů). Závěr „jen d ≤ 2“ je s ní bezpečný, jemné
  rozdíly ploch mezi zdroji ale nejsou.

## Další krok

1. Po sloučení 004/005 pustit `run.py decomp/mi --source …` pro `selfcitation_like` (má pravidla
   pro začátky řádků → test, zda poziční informace prefixu/kmene je reprodukovatelná) a pro šablonu
   005 (test sekční informace kmenů). Stačí přidat zdroje do `Context.corpus`.
2. Rozklad kmene s afixovými seznamy, které zahrnou poziční glyfy (p/f, `-m`, `-g`), a znovu
   I(kmen; pozice) vs I(kmen; sekce): rozhodne, zda je „poziční kmen“ artefakt parseru.
3. MI(d) mezi typy s nízkou frekvencí (kmeny < 20 výskytů) metodou sdílení kmenů na vzdálenost
   (návaznost na 006).
