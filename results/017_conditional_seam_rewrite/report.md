# 017 – Dá lokální přepis pozici, šev a soudržnost slova zároveň?

Kód: `voytools/placerules.py` (přidáno: `ConditionalRules`, `learn_conditional_rules`, `apply_conditional_rules`;
R1–R3 z 014 beze změny; testy `tests/test_placerules.py`), `experiments/017_conditional_seam_rewrite/run.py`,
`summarize.py`. Rozhodnutí a interpretace spec: `experiments/017_conditional_seam_rewrite/notes.md`
(kopie v `results/017_conditional_seam_rewrite/notes.md`). Všechna čísla níže jsou z `results/017_conditional_seam_rewrite/`
(soubor v závorce). Jen ZL (ts2020 je Eva-only, jako v 014); split i join.

Sady (každá 30 seedů ts2020 config 21, tabulky učené z reálu s afixy b): `base` = ts2020; `R12` = +R1+R2;
`R123` = +R1+R2+R3 (014, přepočítáno touž cestou); `R12pA` = +R1+R2′+R3′a (prefix | poslední glyf předchozího, sufix);
`R12pB` = +R1+R2′+R3′b (… | první glyf kmene); `R12pC` = +R1+R2′+R3′c (… | sufix, hrubá pozice); `R12pA3` = R12pA + R3″
(první glyf kmene u tokenů bez prefixu | poslední glyf předchozího); `R12pA3L` = totéž s licencí typů (výsledný token
musí mít v reálu ≥ 2 výskyty). z-skóre: `z_seed` = (reál − průměr sady)/SD 30 seedů, `z_boot` = (reál − průměr sady)/SD
reálu ze 100 bootstrapů stránek, `z_strict` = to z obou s větší |z|; „dožene“ = |z| < 3 ve split **i** join.

## 1. Tabulky a cena přepisu (`table_support.csv`, `rewrite_cost.csv`, `rewrite_cost.png`)

Podpora tabulek reálu (buňka s ≥ 5 pozorováními; ZL split / join): R3 23/29 kontextů, 99,9 % tokenů; R3′a 171/304,
99,2 %; R3′b 204/413, 98,8 %; R3′c 368/768, 97,6 %; R3″ 16/25, 99,4 % (join: 165/288 atd., tytéž podíly). Koncové
glyfy pro R2′ (poměr konec/vnitřek > 5): `g`, `m`.

| sada | přepsáno tokenů (split / join) | z backoffu (podíl podmíněných) | zamítnuto licencí | typy (reál 6 796 / 7 729; ts2020 6 500 / 6 294) | nové typy proti ts2020 |
|---|---|---|---|---|---|
| R12 | 13,5 % / 14,3 % | – | – | 7 678 / 7 466 | 1 776 / 1 762 |
| R123 | 35,9 % / 35,3 % | – | – | 8 927 / 8 644 | 3 793 / 3 661 |
| R12pA | 45,2 % / 44,4 % | 1,7 % / 1,8 % | – | 9 029 / 8 771 | 4 469 / 4 314 |
| R12pB | 54,7 % / 53,6 % | 0,7 % / 1,1 % | – | 8 658 / 8 455 | 4 619 / 4 514 |
| R12pC | 47,8 % / 47,4 % | 3,7 % / 4,0 % | – | 9 160 / 8 917 | 4 741 / 4 609 |
| R12pA3 | 47,8 % / 46,6 % | 1,5 % / 1,7 % | – | 9 337 / 9 021 | 4 872 / 4 653 |
| R12pA3L | 26,3 % / 25,2 % | 1,0 % / 1,2 % | **19,1 % / 19,3 %** | **6 657 / 6 459** | 353 / 346 |

Podmíněný přepis mění 45–55 % tokenů (R3′a přepíše 11 212 prefixů vnitřních tokenů proti 8 185 u R3; R2′ 768;
R3″ 1 349), z backoffu jde ≤ 4 % (R3′c), takže tabulky nejsou „příliš jemné“ – ale každá sada bez licence přidá
2 400–2 900 typů (+37–44 %), tj. 4 300–4 900 typů, které ts2020 nemá. Licence zamítne pětinu přepisů (R1 1 588,
R2 624, R2′ 507, R3′a 5 382, R3″ 1 528 losů na seed), typy pak sedí na reálu (6 657 vs 6 796) – a s nimi zmizí
skoro vše ostatní (odd. 2).

## 2. Profil sad proti reálu (`augmented_profile.csv`, `metric_lists.csv`, `augmented.png`)

Voynich-specifické vlastnosti (SYNTHESIS odd. 1) a soudržnost slova, ZL split: reál / sada (z_seed; z_boot).
Join v závorce, kde se liší.

| metrika | reál | ts2020 | R123 (014) | R12pA (R3′a) | R12pB (R3′b) | R12pC (R3′c) | R12pA3 (+R3″) | R12pA3L (licence) |
|---|---|---|---|---|---|---|---|---|
| **vazba prefix–sufix** mi_prefix_suffix | 0,205 (join 0,183) | 0,182 (1,0; 3,6) | **0,094 (11,9; 16,7)** | **0,188 (3,5; 2,6)**; join 0,176 (1,1; 1,1) | 0,114 (8,0; 13,7) | 0,185 (3,7; 3,0); join 0,172 (1,8; 1,8) | 0,184 (4,4; 3,2); join 0,175 (1,2; 1,3) | 0,168 (2,2; 5,6) |
| šev: edge_mi_excess | 0,186 (join 0,153) | 0,006 (70; 24) | 0,126 (7,0; 8,0) | 0,119 (8,1; 9,0) | 0,155 (3,4; 4,2) | 0,120 (6,8; 8,9) | **0,161 (1,6; 3,3)**; join 0,141 (2,1; 1,5) | 0,057 (7,4; 17) |
| šev: bnd_edge_in (010) | 0,196 | 0,008 | 0,131 (8,6) | 0,123 (9,6) | 0,160 (4,2) | 0,124 (7,9) | 0,166 (2,0; 3,8); join (2,8; 2,1) | 0,060 (7,8) |
| sufix → prefix dalšího (011) | 0,038 | 0,004 (9,8; 1,9) | 0,024 (1,8; 0,8) | 0,021 (3,0; 1,0) | 0,025 (2,6; 0,7) | 0,021 (3,2; 0,9) | 0,023 (2,5; 0,9) | 0,010 (6,9; 1,6) |
| konec řádku 17×: end_enrich_max | 17,5 (join 16,9) | 8,0 (5,9; 4,1) | 6,5 (5,6; 4,8) | **19,9 (−2,6; −1,0)** | 19,9 | 19,9 | 19,7 (−2,5; −1,0) | 9,0 (4,9; 3,7) |
| `-m/-g` uvnitř řádku: end_mg_other | 0,0094 | 0,035 (−2,8; −20) | 0,035 | **0,0094 (0,1; 0,0)** | 0,0094 | 0,0094 | 0,0094 | 0,026 (−2,5; −14) |
| poziční vazba prefixu: prefix__I_pos | 0,134 (join 0,147) | 0,063 (12,8; 12,6) | 0,115 (4,9; 3,3) | 0,107 (5,5; 4,8) | 0,097 (7,0; 6,6) | 0,111 (6,0; 4,0) | 0,101 (6,7; 5,9) | 0,068 (9,1; 12) |
| (R12 pro srovnání) | | | R12: **0,164 (−1,1; −5,5)** | | | | | |
| prefix pos explained (011) | 0,046 | 0,026 (4,9; 3,2) | 0,038 (2,1; 1,1) | 0,039 (1,5; 1,0) | 0,052 (−1,6; −1,1) | 0,042 (1,4; 0,6) | 0,038 (1,8; 1,2) | 0,036 (1,7; 1,5) |
| poziční vazba sufixu: suffix__I_pos | 0,048 | 0,029 (3,8; 5,2) | 0,068 (−2,5; −5,5) | 0,074 (−3,1; −7,1) | 0,074 | 0,074 | 0,074 | 0,037 (2,5; 2,9) |
| začátek řádku: jsd_first_mid | 0,194 | 0,062 (21; 14) | 0,177 (3,1; 1,8) | **0,194 (0,0; −0,1)** | 0,168 (2,4; 2,8) | 0,201 (−0,7; −0,7) | 0,184 (1,0; 1,1) | 0,100 (5,2; 10) |
| konec řádku: jsd_last_mid | 0,107 | 0,048 (7,2; 9,7) | 0,150 (−1,7; −7,0) | 0,155 (−1,9; −7,9) | 0,156 | 0,199 (−3,4; −15) | 0,155 | 0,076 (2,1; 5,1) |
| záhlaví řádku: I_first_rest; head_rest | 0,057; 0,012 | 0,029 (4,9; 1,1); 0,007 (2,6; 1,7) | 0,023 (6,2; 1,3); 0,004 (7,4; 2,7) | 0,020 (8,9; 1,5); 0,003 (7,9; 2,9) | 0,016; 0,003 | 0,015; 0,003 | 0,020 (9,1; 1,5); 0,003 (8,0; 3,0) | 0,020; 0,004 |
| sekční vazba prefixu: prefix__I_sec | 0,085 | 0,012 (20; 9,3) | 0,006 (42; 10) | 0,005 (51; 10) | 0,004 | 0,004 | 0,004 (52; 10) | 0,007 (35; 9,9) |
| uzavřený řádek: bnd_*_ratio_line | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 | ≈ 0 |
| H(sufix \| kmen) | 1,82 (join 1,89) | 1,61 (1,9; 10) | 1,78 (0,4; 2,0) | **1,81 (0,0; 0,2)**; join 1,83 (1,0; 2,8) | 1,79 | 1,81 | 1,84 (−0,2; −1,0); join 1,85 (0,8; 2,0) | 1,71 (1,2; 5,6) |
| sufixů na kmen | 5,25 | 3,78 (9,5; 11) | 4,37 (5,7; 6,6) | 4,43 (4,7; 6,1) | 4,49 (4,5; 5,6) | 4,44 | 4,43 | 3,88 (8,4; 10) |
| cmi_excess | 0,088 | 0,114 (−1,3; −3,7) | 0,065 (2,3; 3,2) | 0,123 (−3,6; −5,0) | 0,053 (3,4; 4,9) | 0,122 (−3,4; −4,9) | 0,132 (−3,3; −6,2) | 0,138 (−2,7; −7,2) |
| H2 glyfů (extended) | 2,42 | 2,39 (0,5; 1,4) | 2,53 (−1,6; −4,8) | 2,56 (−2,3; −6,4) | 2,45 (−0,6; −1,5) | 2,57 | 2,60 (−2,9; −8,2) | 2,46 (−0,7; −1,9) |
| typy (calib) | 6 796 | 6 500 (0,3; 1,8) | 8 927 (−1,6; −13) | 9 029 (−1,7; −13) | 8 658 (−1,4; −11) | 9 160 (−1,7; −14) | 9 337 (−1,9; −15) | **6 657 (0,1; 0,8)** |
| hapax | 0,676 | 0,578 (9,0; 2,8) | 0,645 (2,0; 0,9) | 0,656 (1,1; 0,6) | 0,656 | 0,661 | 0,661 (0,8; 0,4) | 0,569 (7,7; 3,0) |

Počty metrik (65 × sada, kombinovaný verdikt split & join; `metric_lists.csv`): podle z_seed dožene R12pA 6 (+8 jen
v jedné politice), pokazí 4; R12pA3 dožene 8 (+8), pokazí 4; R123 dožene 4 (+8), pokazí 3; R12pA3L dožene 2 (+2),
pokazí 2. Podle z_strict: R12pA 8 (+4), R12pA3 9 (+6), R123 2 (+6). Bootstrap reálu je u typových a délkových metrik
přísnější (SD typů 168 proti SD seedů 1 126: +2 400 typů je z_boot −13, z_seed −1,7); u švu a pozice naopak přísnější
SD seedů.

## 3. Vlastnosti × sady (`augmented_lists.csv`, `minimal_unreachable.csv`)

„Dosaženo“ = všechny metriky vlastnosti |z| < 3 ve split i join. Podle **z_seed** (konvence 014):

| vlastnost | ts2020 | R12 | R123 | R12pA | R12pB | R12pC | R12pA3 | R12pA3L |
|---|---|---|---|---|---|---|---|---|
| poziční vazba prefixu | ne | jen split | přiblíží | přiblíží | přiblíží | přiblíží | přiblíží | ne |
| poziční vazba sufixu | ne | ne (přestřelí) | ne | ne | ne | ne | ne | jen split |
| začátek řádku | ne | **ano** | jen join | **ano** | jen split | **ano** | **ano** | přiblíží |
| konec řádku (jsd, -m, délka) | ne | jen join | jen join | jen join | přiblíží | přiblíží | jen join | jen split |
| konec řádku 17× | ne | ne | ne | **ano** | **ano** | **ano** | **ano** | ne |
| šev glyfový | ne | ne | přiblíží | přiblíží | jen join | přiblíží | **ano** | přiblíží |
| šev sufix → prefix | ne | ne | **ano** | přiblíží | **ano** | přiblíží | jen split | ne |
| uzavřený řádek | ano | ano | ano | ano | ano | ano | ano | ano |
| záhlaví řádku | ne | jen split | pokazí | pokazí | pokazí | pokazí | pokazí | pokazí |
| vazba prefix–sufix | ano | pokazí | pokazí | pokazí (z 3,5 split; join 1,1) | pokazí (8,0) | pokazí (3,7) | pokazí (4,4) | ano |
| sekční vazba prefixu | ne | ne | zhorší | zhorší | zhorší | zhorší | zhorší | zhorší |
| šířka paradigmatu | ne | ne | ne | přiblíží | přiblíží | ne | ne | ne |
| H(sufix \| kmen) vysvětleno (011) | ne | ne | ne | ne | ne | ne | ne | ne |
| typy | ano | ano | ano | ano | ano | ano | ano | ano |
| hapax | ne | ne | jen split | jen split | jen split | jen split | jen split | ne |
| cmi_excess | ano | ano | ano | pokazí | pokazí | pokazí | pokazí | pokazí |
| H2, délka | ano | ano | ano | ano | ano | ano | H2 pokazí | ano |

Podle **z_strict** navíc odpadnou: typy (všechny sady bez licence pokazí, z_boot −11 až −15), konec řádku (jsd_last_mid
přestřelí, z_boot −7), H2 a délka (z_boot −4 až −8), vazba prefix–sufix u R12pA3L (z_boot 5,6); šev glyfový u R12pA3
jen join (z_boot 3,3 split). Nejmenší kombinace, které žádná sada nedá najednou (`minimal_unreachable.csv`, z_seed):
samostatně nedosažitelné jsou **poziční vazba prefixu**, poziční vazba sufixu, konec řádku (jsd/délka), sekční vazba
prefixu, šířka paradigmatu, vysvětlený podíl H(sufix|kmen), hapax; z dosažitelných vlastností jsou nedosažitelné
dvojice **šev glyfový + vazba prefix–sufix**, šev glyfový + šev sufix→prefix, šev glyfový + záhlaví řádku, začátek
řádku + šev sufix→prefix, konec řádku 17× + záhlaví řádku, šev/začátek/17× + vazba prefix–sufix. Podle z_strict jsou
navíc samostatně nedosažitelné šev glyfový, vazba prefix–sufix a typy. Nejvíc vlastností najednou dá R12pA3: začátek
řádku + 17× + šev glyfový + uzavřený řádek (4; z_strict 3).

## 4. Šablona 005 (`template_check.csv`, ZL split, 30 seedů)

Šablona v konfiguraci 009 (`suffix_by=None`) **nemá** vazbu prefix–sufix z konstrukce: mi_prefix_suffix 0,007
(reál 0,205). R3′a ji přidá na 0,120 (z −9,1 seed / −12,8 boot; stále ne), šev dožene celý (edge_mi 0,183 vs 0,186,
z 0,16; R3″ přestřelí na 0,245), konec řádku 17× (17,6 vs 17,5), začátek řádku (jsd 0,178, z 1,8) a sufix → prefix
(0,077 vs 0,038, z −2,4); poziční vazba prefixu zůstane pod reálem (0,094 vs 0,134, z 9,8 – stejný deficit jako
u ts2020: 0,107), záhlaví řádku klesne (0,004 vs 0,057), sekční vazba prefixu klesne, typy +780 (8 122 vs 7 344;
reál 6 796), H(sufix|kmen) šablony 2,8 se nemění, sufixů na kmen 5,05 (z 1,5, reál 5,25). Přepsáno 40 % tokenů.
Obraz je tedy nezávislý na základním generátoru: podmíněný přepis dá šev, konec a začátek řádku, ale ne pozici
prefixu, záhlaví, sekci – a nenahradí paradigma generátoru.

## 5. Samokontrola (`selfcheck.csv`, `selfcheck_min_cell.csv`)

R12pA s tabulkami naučenými z ts2020 (seed + 1000) na ts2020 (seed), 30 seedů, z = (průměr − průměr ts2020)/SD ts2020:
šev (0,85 / 0,90), vazba prefix–sufix (0,00 / −0,04), sufix → prefix (−0,3), poziční vazba prefixu (0,85 / 0,71),
17× (0,74 / 1,38), H(sufix|kmen) (0,5 / 0,9), záhlaví (−1,0 / −0,7) jsou v mezích |z| < 1 nebo těsně nad; **mimo meze**
je 33 ze 130 metrik: hapax (4,3 / 6,5), jsd_first_mid (2,9 / 3,8), jsd_last_mid (2,0 / 4,8), suffix pos explained
(2,7 / 5,7), typy (1,0 / 2,1; +1 100), Zipf (join 2,5). Přepsáno je 23–26 % tokenů, z backoffu 0,7–1 %. Diagnostika
(`selfcheck_min_cell.csv`, 3 seedy): s prahem buňky 5 / 20 / 50 / 200 přepíše R3′a 17,9 / 17,3 / 16,5 / 15,1 %
tokenů a přidá 750 / 730 / 670 / 560 typů; R1 3,6 % a R2′ 1,7 % nezávisle na prahu. Přepis tedy neinjektuje nic ze
*směru* tabulek (šev, vazba, pozice jsou v mezích), ale **mechanismus sám zvyšuje hapax a typy a posouvá kraje řádku**,
protože dvě realizace ts2020 mají jiné podmíněné tabulky (generátor kopíruje vlastní text, tabulky jsou vzorkově
specifické) a minimální vazba pak přepíše 15 % tokenů i tam, kde by neměla. Část nárůstu typů a hapaxu v odd. 1–2
(řádově 700 z 2 500 typů) a posunu jsd_first/last je artefakt mechanismu, ne obsahu tabulek.

## 6. Odpovědi na kritéria

1. 8 sad × 30 seedů, split i join; z oběma způsoby; samokontrola **mimo meze** u 33/130 metrik (odd. 5) – u švu,
   vazby prefix–sufix, pozice a 17× v mezích.
2. (i) **Ano, téměř**: R3′a vrátí vazbu prefix–sufix z 0,094 (R3) na 0,188 / 0,176 (reál 0,205 / 0,183; z 3,5 / 1,1
   seed, 2,6 / 1,1 boot) a šev drží (0,119, tj. 64 % reálu – jako R3). Formálně „pokazí“ jen o práh ve split.
   (ii) R3″ dožene zbylou třetinu švu: 0,161 / 0,141 (z_seed 1,6 / 2,1; z_boot 3,3 / 1,5) – **ano podle seedů, ve split
   ne podle bootstrapu**; R2′ dožene 17×: 19,9 / 18,8 (z −2,6 / −1,8) a `-m` uvnitř řádku přesně (0,0094).
   (iii) Cena: 45–48 % tokenů přepsáno, 1,5–1,8 % z backoffu (R3′c 4 %), **+2 500 typů (+39 %)** a hapax 0,66 (reál
   0,68 – artefakt: z_boot typů −13 až −15); licence typů zamítne 19 % přepisů, drží typy (6 657) a vrátí se k ts2020:
   šev 0,057, pozice 0,068, 17× 9,0.
   (iv) **Ano.** Nejmenší nedosažitelné kombinace (z_seed): *samostatně* poziční vazba prefixu (všechny sady s R3 ji
   sníží na 0,10–0,115 vs 0,134: přepis prefixu podle švu vytlačí přepis podle pozice v témž slotu, R3′c s hrubou
   pozicí nepomůže), sekční vazba prefixu, šířka paradigmatu (sufixů na kmen 4,4 vs 5,25); *dvojice* šev glyfový +
   vazba prefix–sufix (R12pA3: 0,161 a z 4,4), šev + sufix→prefix, šev + záhlaví řádku. Na šabloně totéž: šev ano,
   pozice prefixu, záhlaví a sekce ne.
   **Věta: podpis (lokální přepis nestačí).** Podmíněný přepis dá najednou začátek řádku + konec řádku 17× + šev
   glyfový + uzavřený řádek, a za cenu +39 % typů; ale **poziční vazba prefixu + šev + vazba prefix–sufix + šířka
   paradigmatu zároveň** zůstává mimo dosah každé varianty, a jakmile přepis nesmí vytvářet nové typy, zmizí i šev
   a pozice.
3. Podezřelé: žádná sada nemá |z| < 3 u všech metrik; samokontrola mimo meze (odd. 5) – vysvětleno, ne chyba přepisu
   (tabulky ze stejného vzorku nemění nic: test `test_identity_tables_change_nothing`).
4. Testy (`tests/test_placerules.py`, 279 testů prochází), Python 3.10 / pandas 2.3 kompatibilní, bez scipy.

## Interpretace

Podmínění přepisu i na zbytku slova (R3′a) opravdu vrátí vazbu prefix–sufix, kterou R3 z 014 rozbil, a R3″ s R2′
dohoní šev a konec řádku – lokální glyfová pravidla s tabulkou podmíněnou sufixem tedy umí šev *a* soudržnost
prefix–sufix skoro najednou. Zaplatí to ale třemi věcmi, které se v žádné variantě nesejdou s dosaženým: prefix
nemůže zároveň nést švovou a poziční konvenci (poziční vazba prefixu klesne pod reál u každé sady s R3), přepis
musí vyrábět o 40 % víc typů, než má reál (a s licencí na existující slova mechanismus zkolabuje – reál dosahuje
švu i pozice *v rámci* svého slovníku, generátor s cizími kmeny ne), a šířku paradigmatu ani sekční vazbu prefixu
přepis nemění. Nejmenší nedosažitelná kombinace je tedy pozice prefixu + šev + prefix–sufix (+ šířka paradigmatu
a slovník), nezávisle na základním generátoru (ts2020 i šablona).

## Co výsledek neukazuje

* Neukazuje, že by *žádný* lokální mechanismus pozici, šev a soudržnost nedal: testovali jsme jednu rodinu
  (minimální vazba, tabulky z reálu, kmeny ts2020 beze změny). Mechanismus, který by volil *kmen* podle kontextu
  (ne jen afix), nebo generátor s vlastní podmíněnou gramatikou slova, je mimo rozsah.
* Samokontrola ukazuje, že minimální vazba sama o sobě přidává typy a hapax (15 % přepsaných tokenů i při tabulkách
  z téhož generátoru); nárůst typů v odd. 1 je proto horní odhad ceny, ne její čistá hodnota; „dohnaný“ hapax 0,66 je
  artefakt.
* Bootstrap stránek vychyluje typové metriky (duplicitní stránky): používáme jen jeho SD; u lokalitních metrik
  (loc_*) je bootstrap stránek nevhodný a jejich z_boot nehodnotíme.
* Jen ZL (ts2020 je Eva-only), šablona jen split; GC chybí, takže výsledek nemá kontrolu abecedou.
* Práh |z| < 3 je u metrik s malou SD seedů (mi_prefix_suffix SD 0,005) přísný: „pokazí“ u R3′a je rozdíl 0,017 bitu.

## Další krok

Dvě věci oddělují reál od každé varianty: (a) prefix nese pozici i šev zároveň – měřit v reálu
I(prefix; pozice | poslední glyf předchozího, sufix) a zjistit, zda je to jedna společná tabulka (prefix | pozice,
šev, sufix), kterou by další varianta R3‴ naučila, nebo zda ji reál dosahuje volbou *kmene*; (b) reál dosahuje
švu i pozice ve vlastním slovníku – zkusit přepis, který mění celý token na nejbližší existující tvar (editační
vzdálenost 1 v reálném slovníku) místo afixu, a měřit, kolik švu a pozice přežije bez nových typů. Dále 018–020.
