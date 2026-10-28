# 011 – Rozklad vnitroslovní variace: je „paradigma" poziční jev?

Kód: `experiments/011_suffix_variation/run.py`, `summarize.py`; estimátory `voytools/variation.py`
(testy `tests/test_variation.py`). Deník a odchylky od spec: `experiments/011_suffix_variation/notes.md`.

## Co se počítalo

* **Sady:** ZL3b a GC2a, `P`, clean, split a join; dekompozice (a) naučené afixy, (b) + poziční glyfy (008).
  Kmeny s ≥ 20 tokeny (hlavní; ZL split 120 kmenů / 29 731 tokenů), kontrola ≥ 5 (`explained_min5`).
* **Kontroly 009** (přes `voytools.controls`, layout ZL split): Caesar, Celsus, Plinius, Culpeper (próza),
  Linnaeus (próza i položky). **Nully:** ts2020 calB (ZL split/join), šablona 005 a ac-tri 86 (ZL/GC split),
  30 seedů.
* **Rozklad:** vysvětlený podíl = (H(S|Z) − H(S|Z,X)) / H(S|Z) s korekcí podle spec – H(S|Z,X) po permutaci X
  uvnitř kmene (100× bod, 20× subsample/null) se odečítá. Interval = bod + 2 × kvantily odchylek ze 100 subsamplů
  80 % stránek (kontroly 30), nully kvantily přes seedy. Citlivost: křížově validovaný podíl (`explained_cv`,
  5 foldů, vyhlazený backoff).
* **Pozor na interpretaci korekce:** u X s mnoha hodnotami (stránka, „vše" ≈ 4 000 buněk kmen × X) je
  permutační korekce konzervativní (dolní mez), proto „vše" < „vše kromě stránky". Hlavní číslo pro „všechno
  dohromady" je proto **vše kromě stránky**; stránka zvlášť.

## (i) Kolik H(sufix | kmen) vysvětlí pozice, jazyk, písař, sekce (`decomposition.csv`, `decomposition.png`)

Dekompozice (a), split, bod [95% interval]; H(S|Z) v bitech.

| sada | H(S\|Z) | pozice (22 tř.) | místo v řádku | jazyk | písař | sekce | stránka | pozice+jazyk+sekce | vše kromě stránky |
|---|---|---|---|---|---|---|---|---|---|
| ZL | 1,817 | 0,026 [0,021; 0,030] | 0,025 | 0,058 [0,050; 0,067] | 0,068 [0,060; 0,076] | 0,053 [0,045; 0,063] | 0,082 [0,074; 0,090] | 0,080 | **0,079 [0,070; 0,088]** |
| GC | 1,804 | 0,022 [0,018; 0,027] | 0,020 | 0,049 | 0,059 | 0,049 | 0,068 | 0,065 | **0,065 [0,058; 0,072]** |

* **Prefix** (ZL / GC): pozice 0,046 / 0,043, písař 0,061 / 0,059, sekce 0,065 / 0,063, stránka 0,091 / 0,089,
  vše kromě stránky **0,097 / 0,095** (H(P|Z) = 2,22 / 2,16 bitu).
* **Join** stejně (sufix vše kromě stránky ZL 0,079, GC 0,068). **Dekompozice (b)** s pozičními glyfy: pozice
  vysvětlí o málo víc (ZL sufix 0,031, GC 0,034), celek stejný (0,080 / 0,071). **Kmeny ≥ 5**: rozdíl ≤ 0,005.
* Křížově validovaně je vše ještě menší: ZL sufix pozice 0,007, písař 0,049, vše kromě stránky 0,030.
* V bitech: pozice odebere ze 1,82 bitu 0,048 bitu, všechny kovariáty kromě stránky 0,146 bitu.
* Popisně (`top_stems.csv`): pozice posouvá jen proporce – např. kmen `k` má na začátku řádku `-aiin` 0,29,
  uvnitř 0,23, na konci `-y` 0,27 a `-am` 0,21; kmen `d` na začátku `-aiin` 0,32, uvnitř `-y` 0,57. Sada
  sufixů kmene je na všech místech stejně široká.

## (ii) Latina, ostatní kontroly a nully

Sufix, dekompozice (a), bod; „vše" = vše kromě stránky.

| zdroj | H(S\|Z) | pozice | písař | sekce | stránka | vše |
|---|---|---|---|---|---|---|
| Caesar | 0,849 | 0,003 | 0,007 | 0,015 | 0,033 | 0,004 |
| Celsus | 0,770 | −0,002 | 0,030 | 0,061 | 0,071 | 0,021 |
| Plinius | 0,944 | −0,003 | 0,023 | 0,044 | 0,043 | 0,015 |
| Culpeper | 0,723 | −0,001 | 0,003 | 0,010 | 0,033 | 0,002 |
| Linnaeus položky | 0,547 | 0,054 | 0,033 | 0,058 | 0,118 | 0,065 |
| ts2020 ZL | 1,664 | 0,015 | 0,019 | 0,026 | 0,113 | 0,023 |
| šablona 005 ZL | 2,893 | 0,000 | 0,001 | 0,002 | 0,003 | 0,001 |
| ac-tri 86 ZL | 1,658 | −0,000 | 0,024 | 0,030 | 0,131 | 0,010 |

* Latinská próza: pozice ≈ 0 (záporné hodnoty −0,002/−0,003 jsou v intervalu nuly), sekce (kniha) 1,5–6 %.
  Voynich má proti próze navíc poziční složku (2–3 %) a jazyk/písaře, ale **zbytek je i u něj > 90 %**.
* ts2020 vysvětlí pozicí 1,5 % (řádková pravidla), šablona nic – podle predikcí spec.
* Kopírovací generátory mají **silnější vazbu na stránku** než Voynich (ts2020 0,113, ac86 0,131 vs ZL 0,082;
  v bitech 0,19–0,22 vs 0,15), ale slabou na pozici/jazyk/písaře – variace u nich vzniká lokálním kopírováním.
* Jediná kontrola, kde pozice vysvětlí podobně jako u Voynichu, je Linnaeus v položkovém layoutu (0,054;
  řádek = položka katalogu).

## (iii) Má zbytek strukturu? (`residual_structure.csv`, `consistency.csv`)

Po podmínění kmen + poziční třída + jazyk + písař; bod [interval]; dekompozice (a), split.

| metrika | ZL | GC | ts2020 | ac86 ZL | šablona ZL | latina (rozsah) |
|---|---|---|---|---|---|---|
| I(S_t; S_{t−1}) bity | 0,007 [−0,001; 0,014] | 0,005 [−0,003; 0,012] | 0,014 | 0,006 | 0,000 | 0,009–0,013 |
| I(S_t; P_{t+1}) bity | **0,038 [0,028; 0,047]** | **0,034 [0,024; 0,042]** | 0,004 | 0,002 | 0,001 | 0,006–0,010 (Culpeper 0,036) |
| I(S; sekce \| ·) bity (podíl H) | 0,019 [0,012; 0,027] (1,6 %) | 0,019 (1,5 %) | 0,010 | 0,006 | 0,001 | 0,001–0,006 |
| modální sufix uvnitř stránky − očekávání | 0,023 [0,020; 0,028] | 0,023 | 0,054 | 0,062 | 0,002 | 0,009–0,017 |
| totéž uvnitř odstavce | 0,027 [0,023; 0,032] | 0,025 | 0,052 | 0,095 | 0,001 | 0,009–0,024 |

* **Shoda se sousedem:** ano, ale jen přes hranici slov dopředu – sufix tokenu nese informaci o prefixu
  následujícího tokenu (0,038 bitu) navíc k tomu, co dává pozice/jazyk/písař, ~10× víc než kopírovací nully
  a 4–6× víc než latinská próza. Jde o jev na švu slov (srov. okrajová MI glyfů ≈ 0,20 bitu, `VERIFICATION.md`),
  ne o shodu tvarů v gramatickém smyslu. Sufix ↔ sufix předchozího tokenu: nerozlišitelné od nuly (ZL, GC).
* **Konzistence ve stránce:** ano, malá (+2,3 p. b. nad očekáváním při stejné pozici/jazyku/písaři), ale
  **menší** než u kopírovacích nullů (+5–6 p. b.) a jen o něco větší než v latině – kopírování by dalo víc.
* **Sekce po podmínění:** malý zbytek 0,019 bitu (1,6 % H(S|Z)), 2× nad ts2020, nad latinou (kde sekce po
  podmínění téměř zmizí).
* Dekompozice (b) to nemění (I(S_t; P_{t+1}) 0,045 / 0,044; konzistence +0,022).

## Kontroly podezřelých hodnot (kritérium 3)

* Záporný vysvětlený podíl po korekci: jen u kontrol a nullů, kde je skutečná hodnota ≈ 0 (min. −0,004, vždy
  s intervalem přes nulu) – šum korekce, ne chyba. U Voynichu žádný.
* Stránka > 0,95: ne (surový podíl stránky 0,45, korigovaný 0,08; korekce odebere 0,37 bitu, tj. přeučení na
  207 stránkách je velké a uvedené).

## Odpověď

**(iv) Šířka paradigmatu není poziční jev:** pozice v řádku vysvětlí 2–3 % H(sufix | kmen) (prefix 4–5 %),
pozice, jazyk, písař a sekce dohromady 7–8 % (GC 6,5 %), stránka 8 % – přes 90 % variace tvaru při daném
kmeni zůstává (ZL i GC, split i join, obě dekompozice).

## Interpretace

Volba sufixu/prefixu při daném kmeni je z velké většiny nevysvětlená pozicí ani registrem, takže široké
a vyrovnané „paradigma" Voynichu (009) je vlastnost samotné volby tvaru, ne poziční konvence. Zbytek má jen
slabou strukturu: hlavně vazbu na začátek následujícího slova a malou konzistenci uvnitř stránky, která je
slabší než u kopírovacích generátorů; obsahová volba tvaru tím není prokázaná, kopírování jako zdroj variace
je ale nepravděpodobné.

## Co výsledek neukazuje

* Neukazuje, že zbylá variace je smysluplná – může jít o volnou (nepodmíněnou) volbu písaře, tj. právě to,
  co by dal lidský blábol (testuje 013).
* Kovariáty jsou hrubé (5 písařů, 7 sekcí, 22 pozičních tříd); jemnější kontext (např. konkrétní předchozí
  token) se zde netestoval kromě sousedních afixů.
* Permutační korekce je pro jemné X konzervativní; „vše" včetně stránky je dolní mez. Křížová validace dává
  ještě nižší čísla, takže závěr „malý podíl" na volbě estimátoru nezávisí.
* Afixy jsou naučené parserem; jiný parser by mohl přesunout část variace mezi kmen a afix (008: poziční
  glyfy ji mění jen o desetiny procentního bodu).

## Další krok

* 012: A/B jako registr – jazyk vysvětlí ~6 % H(sufix | kmen); otázka, zda je to systematická náhrada tvarů
  při stejných kmenech.
* 013: porovnat vysvětlený podíl a I(S_t; P_{t+1}) u lidského blábolu (pisatel = „písař").
* Zbytková vazba sufix → prefix dalšího slova: rozložit na konkrétní páry (je to jen fonotaktika švu, nebo
  shoda tvarů?).
