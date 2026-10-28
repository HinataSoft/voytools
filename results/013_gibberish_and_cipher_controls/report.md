# 013 – Chybějící referenční body: lidský blábol a šifra skutečného textu

Kód: `experiments/013_gibberish_and_cipher_controls/run.py`, `summarize.py`; načítání a vysázení
`voytools/controls.py` (`read_gibberish`, `lines_table`, `crop_template`, `stream_control`; testy
`tests/test_controls.py`). Rozhodnutí, formát dat a pravidla: `experiments/013_gibberish_and_cipher_controls/notes.md`.
Šifra Naibbe je zde **referenční bod** (jak vypadá smysl pod obfuskací na našich metrikách), ne hypotéza o rukopisu.

## Zdroje (`sources_used.csv`, `writers_inventory.csv`)

| zdroj | tokeny | řádky | stránky | sekce | typy | hapax | vysázení |
|---|---|---|---|---|---|---|---|
| blábol sdružený (38 pisatelů, seed 13) | 10 024 | 1 202 | 60 | 30 pisatelů | 6 859 | 0,851 | ZL split oříznutý na velikost blábolu, souvislý text |
| blábol ve vlastních řádcích | 10 032 | 1 617 | 88 | 38 | 6 861 | 0,851 | řádky a stránky pisatelů |
| Voynich ZL oříznutý | 10 024 | 1 202 | 60 | 7 | 2 826 | 0,678 | stejný layout jako blábol |
| Caesar / Celsus oříznuté | 10 024 | 1 202 | 60 | 7 | – | – | kontroly 009 ve stejném layoutu |
| Naibbe Caesar (ciphertext) | 34 268 z 80 877 | 4 128 | 207 | 7 bloků | 5 331 | 0,398 | ZL split, souvislý text |
| Naibbe Plinius (Greshko) | 34 268 z 34 764 | 4 128 | 207 | 7 bloků | 5 572 | 0,411 | totéž |
| Naibbe plaintext (respaced) | 34 268 | 4 128 | 207 | 7 bloků | 384 | 0,091 | totéž; tokeny = 1–2písmenné šifrovací jednotky |

Blábol: 38 pisatelů (README mluví o 42 dobrovolnících; v datech je 38 přepisů), 79–480 tokenů na pisatele
(průměr 264), 11–76 řádků, 1–4 stránky. Blábol má méně tokenů než Voynich (10 032 < 34 268), proto se srovnává
s Voynichem a prózou oříznutými na stejný počet (spec). Kontroly podezřelých hodnot: Naibbe plaintext 384 typů vs
ciphertext 5 572 (ne stejně – načtení v pořádku); hapax blábolu 0,851 (≥ 0,3, žádné opakování z tokenizace).

## Seznamy: které Voynich-specifické vlastnosti mají blábol a Naibbe (`feature_lists.csv`)

Pravidlo (spec ho nedefinuje, viz notes): f = (zdroj − R) / (Voynich − R), R = průměr smysluplné prózy na stejné
velikosti; ano f ≥ 0,5, částečně 0,2–0,5, ne < 0,2. Blábol se srovnává s Voynichem oříznutým (`voynich_crop`),
Naibbe s plným ZL split. Hodnoty: zdroj / Voynich / próza.

| vlastnost (SYNTHESIS odd. 1) | metrika | blábol | Naibbe (Plinius, Caesar) |
|---|---|---|---|
| poziční vazby afixů | I(prefix; pozice) | **ne**: 0,006 (vlastní řádky) / 0,130 / 0,001 | **ne**: 0,000 / 0,134 / 0,000 |
| sekční vazby afixů | I(sufix; sekce \| pozice) | **ano, silněji** (sekce = pisatel): 0,385 / 0,078 / 0,020 | **ne**: 0,000 |
| registr uvnitř písaře (analog) | podíl H(sufix\|kmen) vysvětlený písařem | **ano, silněji**: 0,353 (sdružený 0,295) / 0,061 / 0,016 | **ne**: 0,000 |
| šířka paradigmatu | sufixů na kmen; H(sufix\|kmen) | **ano**: 3,96 / 4,91 / 2,09; 1,89 / 1,77 / 0,84 | **ano**: 4,50 / 5,25 / 2,17; 1,99 / 1,82 / 0,80 |
| hapax | podíl typů | **ano, silněji**: 0,851 / 0,678 / 0,653 | **ne, opačně**: 0,411 / 0,676 / 0,610 |
| nízké pořadí tokenů | MI kmenů d = 1; token MI excess | kmeny **ano** (0,004 / 0,039 / 0,066), tokeny ne (0,021 / −0,004 / 0,023) | **ano** obojí: 0,029 / 0,051 / 0,168; 0,012 / 0,009 / 0,035 |
| okrajová MI | edge MI excess | **ne**: 0,013 / 0,199 / 0,052 | **ne**: 0,005 / 0,186 / 0,041 |
| konec řádku | obohacení koncového glyfu | **ne**: 2,1× / 16,4× / 1,2× (vlastní řádky) | nelze (řádky mechanické) |
| uzavřený řádek | poměr okrajové MI přes řádek / uvnitř | **ne**: 0,79 / 0,02 / 0,40 | nelze (řádky mechanické) |
| pokles sdílení se vzdáleností | `loc_log_dist` | nelze rozlišit (Voynich ≈ próza) | nelze rozlišit |
| paradigma není poziční (011) | podíl H(S\|Z) vysvětlený vším kromě stránky | ano, silněji (0,140 / 0,064 / 0,011) – ale pisatelem, ne pozicí | ne (0,000) |
| sufix → prefix dalšího slova (011) | I(S_t; P_{t+1} \| ·) | částečně: 0,011 / 0,022 / 0,006 | ne: 0,001 / 0,038 / 0,008 |
| BPE žlab | – | nepočítáno | nepočítáno |

**Odpověď (kritérium 2):**
* **Lidský blábol má:** šířku a vyrovnanost „paradigmat" (H(sufix | kmen) 1,89 bitu – jako Voynich, 2× próza),
  vysoký hapax (víc než Voynich), slabou posloupnostní vazbu kmenů a velmi silné registrové/sekční vazby afixů –
  ty ale nesou **pisatelé** (každý člověk má svůj „jazyk").
* **Šifra latiny (Naibbe) má:** šířku paradigmatu (vzniká z více šifrovacích tabulek pro jedno písmeno) a nízké
  pořadí tokenů (zhruba na úrovni Voynichu), nic víc.
* **Ani jeden nemá:** poziční vazby afixů (I(prefix; pozice) 0,13 bitu jen u Voynichu), okrajovou MI glyfů na
  švu slov (0,19 vs ≤ 0,013), konvenci konce řádku (16× obohacení vs 2×), uzavřený řádek (blábol ve vlastních
  řádcích 0,79 vs Voynich 0,02) a vazbu sufix → prefix dalšího slova v síle Voynichu. Naibbe navíc nemá ani
  sekční/registrové vazby (z konstrukce) a má opačný hapax (0,41).

## Vzdálenosti profilů (`distances.csv`, `profile_heatmap.png`)

RMS standardizovaných rozdílů 37 metrik 009 (SD přes všechny zdroje tabulky):

| pár | celkem | slovník | pozice | paradigma | lokalita | řádek |
|---|---|---|---|---|---|---|
| Voynich (oříznutý) ↔ blábol | **1,88** | 1,59 | 1,85 | 1,83 | 3,24 | 2,10 |
| Voynich (oříznutý) ↔ Caesar / Celsus (oříznuté) | 1,41 / 1,42 | 0,46 / 0,37 | 1,57 / 1,46 | 1,39 / 1,60 | 0,80 / 1,56 | 2,27 / 2,06 |
| Voynich ZL ↔ Naibbe Plinius / Caesar | **1,39 / 1,38** | 0,84 / 0,89 | 1,59 / 1,57 | 0,86 / 0,94 | 1,41 / 1,35 | 2,21 / 2,17 |
| Voynich ZL ↔ latinská próza (009) | 1,38–1,45 | 0,38–1,17 | 1,41–1,54 | 1,44–1,49 | 0,90–1,44 | 2,11–2,21 |
| Voynich ZL ↔ ts2020 | 0,89 | 0,32 | 0,95 | 0,60 | 1,24 | 1,26 |
| Naibbe Plinius ↔ Plinius (latina) | 0,84 | 1,64 | 0,23 | 1,31 | 0,77 | 0,34 |
| Naibbe Caesar ↔ Caesar | 0,65 | 0,96 | 0,09 | 1,19 | 0,56 | 0,12 |
| blábol ↔ Caesar / Celsus (oříznuté) | 1,60 / 1,45 | 1,21 / 1,29 | 1,97 / 1,77 | 1,57 / 1,47 | 2,48 / 1,84 | 1,05 / 1,16 |

* Blábol je od Voynichu **dál** než latinská próza stejné velikosti (1,88 vs 1,41) – hlavně lokalitou (pisatel =
  blok stránek) a pozicí.
* Naibbe je od Voynichu stejně daleko jako latinská próza (1,39 vs 1,38–1,45): přiblíží se v paradigmatu
  (0,86 vs ~1,45), ale ztratí slovník (hapax) a zůstane daleko v pozici a řádku. Nejbližší zůstává ts2020 (0,89).
* Šifra zachová z latiny pozici, řádek a posloupnost (0,09–0,34), mění slovník a paradigma (1,0–1,6).

## Po pisatelích a I(afix; pisatel) (`writers.csv`, `profile.csv`)

38 pisatelů vs 38 souvislých úseků ZL stejných velikostí (průměr 264 tokenů); průměr ± SD:

| | blábol | úseky Voynichu |
|---|---|---|
| hapax | 0,843 ± 0,105 | 0,795 ± 0,046 |
| H2 glyfů (uvnitř tokenu) | 2,86 ± 0,42 | 2,07 ± 0,22 |
| délka slova (glyfy) | 4,97 ± 0,87 | 4,00 ± 0,33 |
| sufixů na kmen (kmeny ≥ 5) | 2,21 ± 0,86 | 3,20 ± 0,49 |
| H(sufix \| kmen) | 0,72 ± 0,48 | 1,29 ± 0,24 |
| délka posledního slova řádku − ostatní | +0,41 ± 0,74 | +0,05 ± 0,38 |

* Pisatelé se mezi sebou liší **2–3× víc** než stejně velké úseky Voynichu (SD hapaxu, H2, délky, sufixů na kmen).
* **I(afix; pisatel) vs I(afix; písař):** pisatel vysvětlí 29,5 % H(sufix | kmen) (sdružený korpus; ve vlastním
  layoutu 35 %), písař u Voynichu na stejné velikosti 6,1 % (plný ZL 6,8 %). I(sufix; sekce) = I(sufix; pisatel)
  0,30 bitu vs Voynich I(sufix; sekce) 0,07–0,08. Výhrada: 38 pisatelů proti 5 písařům (víc kategorií; korekce
  permutací bias odstraní, ale ne rozdíl v jemnosti dělení).

## Interpretace

Lidský blábol reprodukuje vnitroslovní „bohatost" Voynichu (široké, vyrovnané sady tvarů, mnoho hapaxů, slabé
pořadí), ale ne jeho konvence vázané na místo – poziční vazby afixů, šev slov (okrajovou MI), konec řádku
a uzavřený řádek; a jeho „registry" jsou mnohem silnější a idiosynkratické (každý pisatel jinak), zatímco
Voynich má registr slabý a sdílený. Šifra latiny dává šířku paradigmatu mechanicky (víc zápisů téhož písmene),
ale nic z pozičních, řádkových ani registrových vlastností – ty tedy nejsou vedlejším produktem obfuskace typu
Naibbe.

## Co výsledek neukazuje

* Blábol je krátký (264 tokenů na pisatele, 10 000 celkem) a psaný v jiném prostředí (moderní pisatelé, latinka,
  zadání „piš nesmysl"); nevylučuje, že *jeden* trénovaný pisatel na 200 stránkách vyvine poziční a řádkové
  konvence. G&B sami píší, že vyšší struktura stránek a kvaternionů se na jejich vzorcích testovat nedá.
* Naibbe je jedna konkrétní šifra; jiná obfuskace může mít jiný profil. Vysázení do layoutu ZL je mechanické
  (řádky šifry nejsou řádky rukopisu), takže řádkové vlastnosti Naibbe nelze posoudit.
* Pravidlo ano/částečně/ne je konvence (f vůči rozdílu Voynich − próza); čísla a intervaly jsou ve
  `feature_lists.csv`.

## Další krok

* Pozice a šev slova jsou teď jediné vlastnosti, které nemá text, blábol, šifra ani generátor (ts2020 má polovinu
  pozice): rozložit I(prefix; pozice) a okrajovou MI na konkrétní glyfy/páry (SYNTHESIS: „rozklad okrajové MI na páry").
* Blábol jednoho pisatele na delším textu (pokud existuje korpus delší glosolálie) – test, zda poziční konvence
  vzniknou s délkou.
