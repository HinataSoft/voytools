# 006 – Téma na dálku: sdílení kmenů mezi stránkami, popisky versus text

Kód `experiments/006_topic_at_distance/run.py`, rozhodnutí a odchylky od spec `notes.md`. Všechna čísla
pocházejí ze souborů v této složce. Hlavní nastavení: kosinus vektorů log(1+f) kmenů, vyřazeno
K = 30 nejčastějších kmenů, politika `split`; intervaly jsou 95% bootstrap po stránkách (200×), u nullů
2,5–97,5% kvantily replik; p = podíl replik nullu, které jsou stejně extrémní jako reál nebo extrémnější
((1+#)/(1+n); nejmenší možné p je 0,0099 pro 100 permutací a 0,0196 pro 50 auto-copy vzorků).

## Data

- 201 stránek s ≥ 40 čistými tokeny `P` (ZL split, GC split, GC join; ZL join 200), 20 100 dvojic,
  z toho 8 745 ze stejné sekce (`pairs.csv`).
- Zdroje (`runs.csv`, 2 204 řádků): reál, bootstrap po stránkách 200×, permutace stránek globálně / uvnitř
  sekce / uvnitř písaře, každá 100×, auto-copy config 86 z 003 (trigram), 50×. `voytools.selfcitation`
  (004) zatím neexistuje, takže druhý kopírovací null chybí.
- 30 vyřazených funkčních kmenů ZL split: `k o t aiin d a kee l e y r ke che s cth ain te ckh ee tee ch kch kche
  pche she ai am p al tch`; GC split: `h a am o e 9 8 c y hc hC s k 7 C an 1c K 1 hcc h1 H ap h1c kc cc az k1 ae 4`
  (`function_stems.csv`, 100 kmenů s počty pro každou kombinaci).

## Výsledky

### 1. Sdílení podle vzdálenosti a sekce (`curve.csv`, `curve.png`)

Podobnost klesá se vzdáleností monotónně, a to **ve stejné i v jiné sekci**. ZL split, stejná sekce:
0,188 [0,158; 0,226] při vzdálenosti 1 → 0,097 [0,089; 0,107] při 32+; jiná sekce: 0,232 [0,159; 0,285] →
0,123 [0,114; 0,133]. Dvojice z různých sekcí jsou v každém binu podobnější než dvojice ze stejné sekce.
Je to důsledek velikosti stránek: krátké herbářové stránky A (sekce H, 95 stránek) mají mezi sebou
nízký kosinus. Globální permutace dává plochou křivku (≈ 0,12), auto-copy klesá od vzdálenosti 1 (0,147)
k 0,054 ve stejné sekci.

### 2. Rozklad OLS (`ols.csv`)

Spec model `share ~ log(dist) + same_section + same_language + same_hand` dává pro ZL split záporný
koeficient sekce −0,032 [−0,041; −0,021] i písaře −0,035 [−0,057; −0,014]. Stejné znaménko má i auto-copy,
které o sekci nic neví: −0,035 a −0,029. Spec model je tedy zkreslený délkou stránek. Proto jsem
k němu přidal model `spec+size` (+ log n_i + log n_j a |log n_i − log n_j|; mimo spec, viz `notes.md`):

| cos, K=30, spec+size | log_dist | same_section | same_language | same_hand |
|---|---|---|---|---|
| reál ZL split | −0,0099 [−0,0139; −0,0057] | +0,0032 [−0,0098; +0,0194] | +0,0295 [+0,0189; +0,0440] | −0,0029 [−0,0216; +0,0110] |
| reál GC split | −0,0093 [−0,0120; −0,0067] | +0,0060 [−0,0043; +0,0191] | +0,0277 [+0,0157; +0,0419] | −0,0011 [−0,0161; +0,0118] |
| reál ZL join | −0,0091 | −0,0007 [−0,0120; +0,0159] | +0,0271 | −0,0042 |
| reál GC join | −0,0086 | +0,0049 [−0,0047; +0,0200] | +0,0229 | +0,0007 |
| auto-copy 86 ZL split (průměr) | −0,0026 (p 0,02) | +0,0020 | +0,0027 | −0,0010 |
| perm. v sekci ZL split | −0,0028 (p 0,01) | +0,0111 | −0,0002 | +0,0035 |
| perm. v písaři ZL split | −0,0012 (p 0,01) | −0,0019 | +0,0288 | +0,0010 |

- **Vzdálenostní složka** je v reálu 3–4× silnější než u auto-copy a než u permutace uvnitř sekce.
  Žádná replika žádného nullu nedosáhla reálné hodnoty (p ≤ 0,02 u všech nullů, ZL i GC, split i join).
- **Jazyková složka** (Currier A/B) je kladná. Permutace uvnitř písaře ji zachová (ZL split p 0,40),
  permutace uvnitř sekce ji zruší (p 0,01). Je to tedy rozdíl A/B v rámci sekcí (např. herbář A × B)
  a nedá se oddělit od písaře.
- **Sekční složka nad vzdálenost, velikost, jazyk a písaře je v plné sadě stránek nulová** (CI obsahuje 0
  ve všech čtyřech kombinacích).

### 3. Klíčový kontrast: stejná sekce, dist ≥ 16 minus jiná sekce, dist ≤ 3 (`contrast.csv`)

| cos, K=30 | surový | po odečtení velikosti | jen dvojice téhož písaře (surový) |
|---|---|---|---|
| ZL split | −0,096 [−0,142; −0,053] | −0,060 [−0,092; −0,020] | −0,092 [−0,139; −0,037] |
| ZL join | −0,086 [−0,125; −0,048] | −0,052 [−0,087; −0,014] | −0,083 [−0,127; −0,033] |
| GC split | −0,091 [−0,142; −0,049] | −0,045 [−0,073; −0,016] | −0,098 [−0,167; −0,048] |
| GC join | −0,088 [−0,135; −0,046] | −0,048 [−0,070; −0,020] | −0,094 [−0,150; −0,039] |

Kontrast je **záporný**: daleké stránky téže sekce sdílejí méně kmenů než blízké stránky různých sekcí.
Všechny nully leží výš (průměr po odečtení velikosti: globálně −0,002, v sekci −0,010, v písaři −0,010,
auto-copy −0,014; p ≤ 0,02 všude). Reál je tedy „lokálnější“ než všechny nully, včetně kopírovacího.

### 4. Citlivost na stránky `T`/`C`/`A` (`sensitivity_sections.csv`, mimo spec)

Stránky s kódem `T` (jen text) leží uvnitř jiných sekcí (f1r, f66r, f76r, f85r1 …), a proto tvoří většinu
dvojic „jiná sekce, ≤ 3“. Když se vyřadí `T`, `C` a `A` (zbývá 188 stránek, ZL join 187):

| cos, K=30, bez T/C/A | same_section (spec+size) | p vs perm. v písaři | p vs auto-copy | kontrast po odečtení velikosti |
|---|---|---|---|---|
| ZL split | +0,0216 [+0,0090; +0,0354] | 0,010 | 0,020 | −0,057 [−0,119; −0,001] |
| ZL join | +0,0147 [+0,0017; +0,0273] | 0,010 | 0,020 | −0,058 [−0,096; +0,000] |
| GC split | +0,0150 [+0,0032; +0,0277] | 0,010 | 0,020 | −0,040 [−0,083; +0,003] |
| GC join | +0,0189 [+0,0068; +0,0313] | 0,010 | 0,020 | −0,035 [−0,080; +0,011] |

Teprve pak se objeví malá kladná sekční složka, větší než u permutace v písaři i u auto-copy (u auto-copy
+0,003 až +0,004). Při vyřazení jen `T` je pořád v CI nuly. Vzdálenostní složka se nemění
(−0,007 až −0,010, p ≤ 0,02 vůči všem nullům).

### 5. Citlivost na vyřazení funkčních kmenů 0 / 30 / 100 (`ols.csv`, `contrast.csv`)

| ZL split, reál, spec+size | K=0 | K=30 | K=100 |
|---|---|---|---|
| cos log_dist | −0,0106 [−0,0145; −0,0065] | −0,0099 | −0,0030 [−0,0043; −0,0017] |
| cos same_section | +0,0132 [−0,0052; +0,0326] | +0,0032 | +0,0034 [−0,0011; +0,0098] |
| cos kontrast (surový) | −0,034 [−0,104; +0,031] | −0,096 | −0,022 [−0,047; +0,001] |
| cos kontrast po velikosti | +0,009 [−0,032; +0,052] | −0,060 | −0,013 [−0,033; +0,004] |
| jac log_dist | −0,0093 | −0,0056 | −0,0005 [−0,0015; +0,0003] |

Vzdálenostní složka je kosinem prokázána při všech K a Jaccardem při K=0 a 30. Při K=100 Jaccard
nestačí, protože kmenů s ≥ 2 výskyty na stránce zbývá málo. Záporný kontrast je výrazný jen při K=30.
Při K=0 (dominují funkční kmeny) a K=100 obsahují CI nulu. Sekční složka v plné sadě není při žádném K.
GC se chová stejně.

### 6. Popisky versus text (`labels.csv`, `labels.png`)

Popisky typů Lz, Lf, Lc, Ls, Ln, čisté tokeny, kmen podle afixů z P textu. Při K=30 zbývá v ZL split
462 tokenů ze 38 stránek. `same_page` = kmen se vyskytuje v P textu téže stránky; null = permutace
popisků mezi stránkami téže sekce (200×).

| K=30, same_page | ZL split | ZL join | GC split | GC join |
|---|---|---|---|---|
| všechny: reál / null | 0,048 / 0,048 (p 0,58) | 0,032 / 0,039 (p 0,90) | 0,066 / 0,049 (p 0,025) | 0,039 / 0,036 (p 0,38) |
| ok-/ot-: reál / null | 0,041 / 0,044 (p 0,73) | 0,035 / 0,043 (p 0,87) | 0,094 / 0,066 (p 0,030) | 0,041 / 0,037 (p 0,45) |
| ostatní: reál / null | 0,052 / 0,050 (p 0,49) | 0,031 / 0,036 (p 0,81) | 0,049 / 0,040 (p 0,18) | 0,038 / 0,035 (p 0,47) |
| Lf (rostliny, farm.): reál / null | 0,081 / 0,080 (p 0,58) | 0,073 / 0,079 (p 0,71) | 0,153 / 0,090 (p 0,010) | 0,082 / 0,061 (p 0,12) |

Zodiakové popisky (Lz, 220 v ZL split) nemají na svých stránkách P text, takže tam je podíl 0 v reálu
i v nullu. Proti textu P+C+R (`text == 'PCR'`) mají Lz 0,155 proti nullu 0,172 (p 0,84). Pro `same_section`
nevychází nic nad null (ZL split všechny p 0,68). Při K=0 a K=100 (`labels.csv`) je jediné p < 0,05 u GC
split K=0 (0,045) a ZL split K=100 (0,030; podíl 0,025 proti 0,013). Signál je nekonzistentní.

### 7. Tokeny ok-/ot- v běžném textu (`okot_profile.csv`, `initial_map.csv`)

Třídy podle počátku tokenu. Třídy v GC jsou odvozené zarovnáním lokusů a jsou čisté: 1 359 z 1 445 ZL ok-
odpovídá GC oh-/Ah-, 1 347 z 1 388 ZL ot- odpovídá GC ok-/Ak-.

| ZL split (GC split) | ok-/ot- | qo- | ch-/sh- | d- | ostatní | null (posun v řádku) |
|---|---|---|---|---|---|---|
| podíl prvních v řádku | 0,058 (0,059) | 0,098 (0,098) | 0,037 (0,035) | 0,194 (0,182) | 0,170 (0,158) | ≈ 0,11–0,12 |
| podíl posledních v řádku | 0,111 (0,106) | 0,064 (0,061) | 0,081 (0,079) | 0,205 (0,196) | 0,135 (0,130) | ≈ 0,11–0,12 |
| průměrná rel_pos | 0,554 (0,553) | 0,472 (0,469) | 0,491 (0,491) | 0,540 (0,537) | 0,490 (0,493) | ≈ 0,50 |
| délka (glyfy) | 5,06 (4,54) | 5,81 (5,20) | 4,07 (3,79) | 4,15 (3,30) | 4,10 (3,43) | – |

ok-/ot- se vyhýbají začátku řádku (0,058 [0,047; 0,070] proti 0,113 po posunu) a stojí spíš v druhé
polovině řádku. JSD jejich rozdělení pozic, předchůdců i následovníků od každé z tříd qo-, ch-/sh-, d-
a ostatních převyšuje hodnotu po posunu i dolní mez bootstrapu (např. JSD předchůdců proti d-
0,041 [0,030; 0,053] vs. 0,007 po posunu). Totéž ale platí pro každou dvojici tříd: d- má na krajích
řádku nadreprezentaci, ch-/sh- se začátku vyhýbá ještě víc (0,037). Nejblíž mají ok-/ot- k ch-/sh-
(JSD předchůdců 0,007 [0,005; 0,011] proti 0,004 po posunu).

### Podezřelé hodnoty (kritérium 4)

Žádná dvojice různých stránek nemá kosinus > 0,95. Maximum je 0,865 v ZL (f78r–f84v, K=0) a 0,812 v GC
(`pairs.csv`, `max_cos_k0` v `runs.csv`). Absolutní hodnota kontrastu nikde nepřesáhla 0,5 (maximum je
0,098).

## Odpovědi (kritérium 2)

| otázka | ZL | GC |
|---|---|---|
| (i) sekční složka nad vzdálenostní? | **ne** v plné sadě (+0,003, CI obsahuje 0); slabé **ano** po vyřazení T/C/A (+0,022 [0,009; 0,035]) | **ne** (+0,006, CI obsahuje 0); slabé **ano** bez T/C/A (+0,015 [0,003; 0,028]) |
| (ii) přežije kontrolu písaře? | v plné sadě není co kontrolovat; bez T/C/A **ano** (nad perm. v písaři, p 0,01). Jazyková složka A/B je ale písařská | stejně: bez T/C/A **ano** (p 0,01) |
| (iii) reprodukuje kopírovací null sekční složku? | **ne** (+0,002 až +0,004); nereprodukuje ani sílu vzdálenostního poklesu (−0,0026 vs −0,0099) ani záporný kontrast (−0,014 vs −0,060). Záporné koeficienty spec modelu (efekt velikosti) reprodukuje | **ne**, totéž |
| (iv) stopa popisků v textu téže stránky nad náhodu? | **ne** (0,048 vs 0,048, p 0,58; join p 0,90) | **nejednoznačně**: split ano (0,066 vs 0,049, p 0,025, hlavně Lf), join ne (p 0,38) |
| (v) chovají se ok-/ot- odlišně? | **ano** (vyhýbají se začátku řádku, pozdější pozice, delší) | **ano**, stejný profil |

## Interpretace

Sdílení kmenů mezi stránkami je dáno hlavně **blízkostí v rukopisu**: podobnost klesá se vzdáleností
3–4× strměji než u kopírovacího generátoru i u permutací. Sekce jako celek k tomu v plné sadě nic nepřidává.
Malá sekční složka se objeví jen po vyřazení nejednoznačně kódovaných stránek T/C/A a je pak nad
písařským i kopírovacím nullem. Predikci „téma“ ze spec (stejná sekce daleko > jiná sekce blízko) data
nepotvrzují: kontrast je záporný ve všech čtyřech kombinacích. Silný lokální pokles je slučitelný s pozvolným
posunem slovníku nebo pravopisu během psaní, s tématickou kontinuitou sousedních listů i s kopírováním
s delší pamětí, než má config 86. Popisky nemají v textu své stránky robustní stopu.

## Co výsledek **neukazuje**

- Že sekce *nenesou* obsah. Kód ilustrace je hrubá a místy sporná proxy tématu a téma v rámci sekce
  (konkrétní rostlina) se tu netestovalo. Netestovaly se ani dvojice herbář–farmacie se stejnou rostlinou.
- Že vzdálenostní pokles je obsahový. Rozlišit „posun slovníku během psaní“ (písař, čas, A→B) od
  „sousední stránky o podobném“ tato data neumí. Vzdálenost v souboru navíc není fyzické pořadí
  bifolií (přeházené listy, chybějící folia).
- Že ok-/ot- jsou jména nebo substantiva. Mají odlišný poziční profil, ale ten má každá počáteční třída.
  Posun v řádku je null jen pro pořadí, ne pro „jmennost“.
- Kopírovací null je jediný (config 86, provizorní, kritérium 1 v 003 nesplněno). Self-citation port (004)
  chybí, takže kopírování s delší pamětí nebo s cíleným výběrem zdroje zatím vyloučit nelze.
- Model spec bez velikosti stránky dává zavádějící znaménka (záporná sekce a písař i u auto-copy).
  Hlavní čísla jsou proto z modelu s velikostí, který spec nepředepisuje.

## Další krok

1. Zopakovat rozklad s **fyzickým pořadím** (bifolio, quire) a s permutacemi bifolií, aby se oddělila
   „blízkost v psaní“ od „blízkosti v knize“; přidat `strict` politiku.
2. Až bude 004, pustit `selfcitation_like` stejným kódem (`run.py pages --sources …` po doplnění zdroje)
   a ověřit, zda kopírování s pamětí přes stránky dá tak strmý pokles jako reál.
3. Pro téma jemněji než po sekcích: dvojice herbář–farmacie s toutéž rostlinou (seznam ze sekundární
   identifikace ilustrací) a podobnost očištěná o vzdálenost a velikost (rezidua modelu `spec+size`).
4. ok-/ot-: porovnat poziční profil s nullem, který zachová glyfovou třídu a pozici kmene (007).
