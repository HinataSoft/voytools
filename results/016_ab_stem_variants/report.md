# 016 – Kmeny jen v B: glyfové varianty kmenů A, nebo nová slova?

Kód: `experiments/016_ab_stem_variants/run.py`; funkce `voytools/registers.py` (Levenshtein na jednotkách,
index symetrických smazání, inventář editací, substituce, sekční test; testy `tests/test_registers.py`).
Rozhodnutí a odchylky od spec: `experiments/016_ab_stem_variants/notes.md`.

## Co se počítalo

* ZL3b a GC2a, `P`, clean, split a join; výřezy **celkem** a **herbář**; čtyři reprezentace tokenu: kmen
  dekompozice (a), kmen (b) (+ poziční glyfy), **celý token** (`tok`) a celý token v **BPE jednotkách**
  (64 sloučení, `bpe`). Jednotky: Eva glyfy presetu `extended` (ZL), znaky v101 (GC), písmena (latina).
* Stejná velikost jako 012: N tokenů na stranu (celkem ZL 5 497 / 5 094, GC 5 692 / 5 389; herbář ≈ 1 600–1 750),
  celé stránky, **200 losů**. V každém losu: A × B (`AB`), dvě poloviny A (`AA`), dvě poloviny B (`BB`) a null tvaru
  (`SHUF`: B-only kmeny proti kmenům A s promíchanými glyfy při zachování délek). Oba směry: **B-only vs A**
  (`y_only`, hlavní) a A-only vs B (`x_only`).
* X-only kmen = ≥ 2 tokeny v X, 0 v Y; protějšek = nejbližší kmen Y (≥ 2 tokeny) podle Levenshteina, přesně do
  2 editací. Operace zapsány ve směru A → B, sousední vložení/smazání sloučena (`ins:ed@mid`).
* Nully: ts2020 config 21 (ZL split/join, 30 seedů) a šablona 005 (ZL/GC split, 30), A/B podle layoutu.
  Kontroly: latina **Plinius knihy 16–21 vs 23–30** (stejný jazyk, jiné téma) a **Caesar vs Celsus** (stejný jazyk,
  jiný registr), 200 losů, N = 7 392 / 17 134.
* Hodnoty: průměr [2,5 %; 97,5 %] přes losy; `A×B − poloviny` = rozdíl v témže losu proti průměru AA a BB.

## Kolik textu je „jen v B“ (`counterpart_excess.csv`, `share_only`)

Podíl tokenů strany B připadajících na B-only kmeny (celkem, split):

| | ZL kmen (a) | ZL token | GC kmen (a) | GC token | ts2020 token | šablona token | Plinius token | Caesar/Celsus token |
|---|---|---|---|---|---|---|---|---|
| A × B | 0,042 | **0,219** | 0,046 | **0,208** | 0,131 | 0,075 | 0,189 | 0,301 |
| poloviny | 0,018 | 0,054 | 0,023 | 0,057 | 0,124 | 0,072 | 0,098 | 0,066 |
| rozdíl [CI] | +0,024 [0,015; 0,037] | **+0,165 [0,125; 0,210]** | +0,024 [0,014; 0,033] | **+0,152 [0,118; 0,191]** | +0,007 [−0,013; 0,031] | +0,003 [−0,008; 0,020] | +0,091 [0,064; 0,119] | +0,235 [0,215; 0,251] |

Čtyřnásobek polovin u celých tokenů (ZL i GC, join stejně; herbář 0,27 / 0,24), v generátorech bez jazyka nula.
Latina: jiné téma 2×, jiný registr 4,5× polovin. A-only strana je menší (0,154 ZL, 0,148 GC tokenů A).

## (i) Mají B-only kmeny protějšek v A? (`counterpart_rates.csv`, `counterpart_excess.csv`, `counterpart_rates.png`)

Podíl B-only **tokenů**, jejichž kmen má v A protějšek ve vzdálenosti ≤ 1 (celkem; ≤ 2 je všude 0,98–1,00):

| sada | rep | A × B | poloviny (AA, BB) | A×B − poloviny [CI] | null tvaru SHUF |
|---|---|---|---|---|---|
| ZL split | kmen (a) | 0,927 [0,86; 0,98] | 0,939 | −0,011 [−0,076; 0,040] | 0,658 |
| ZL split | **token** | **0,836** [0,78; 0,89] | **0,921** | **−0,085 [−0,147; −0,026]** | 0,172 |
| ZL split | BPE | 0,935 | 0,955 | −0,021 [−0,063; 0,010] | 0,828 |
| ZL join | token | 0,807 | 0,904 | **−0,097 [−0,159; −0,035]** | 0,157 |
| GC split | kmen (a) | 0,889 | 0,957 | **−0,069 [−0,140; −0,005]** | 0,702 |
| GC split | **token** | **0,851** | **0,935** | **−0,084 [−0,140; −0,039]** | 0,194 |
| GC join | token | 0,826 | 0,927 | **−0,101 [−0,149; −0,056]** | 0,181 |
| ZL split herbář | token | 0,692 | 0,900 | **−0,208 [−0,314; −0,105]** | – |
| GC split herbář | token | 0,721 | 0,909 | **−0,189 [−0,284; −0,088]** | – |
| ts2020 ZL split | token | 0,731 | 0,737 | −0,006 [−0,070; 0,061] | 0,161 |
| šablona ZL split | token | 0,792 | 0,812 | −0,020 [−0,083; 0,039] | 0,327 |
| Plinius (téma) | token | 0,195 | 0,269 | **−0,074 [−0,115; −0,029]** | 0,039 |
| Caesar/Celsus (registr) | token | 0,202 | 0,307 | **−0,105 [−0,130; −0,077]** | 0,040 |

* **Podezřelá hodnota podle spec (krit. 3) nastala:** u kmenů (a), (b) i BPE mají poloviny téhož jazyka protějšek
  ≤ 1 v 0,93–0,96 (ZL i GC) – kmeny po odtržení afixů jsou 1–3 glyfy dlouhé (`k`, `a`, `at`, `o`), takže téměř
  každý kmen je od jiného jednu editaci; i promíchané glyfy dávají 0,66–0,83. **Na kmenech dekompozice se
  otázka rozhodnout nedá**; nosná je reprezentace **celý token** (poloviny 0,90–0,94, null tvaru 0,16–0,19).
* **Odpověď (i): protějšky *nejsou* nad vnitřní kontrolu – jsou pod ní.** B-only tokeny mají v A blízký tvar
  (≤ 1) v 81–85 % (ZL i GC, split i join), ale tokeny jen v jedné polovině téhož jazyka ho mají v druhé polovině
  v 90–94 %; rozdíl −0,08 až −0,10 s intervalem mimo nulu, v herbáři −0,19 až −0,21. Pokud by B byl A s jiným
  zápisem, B-only tvary by byly k A *blíž* než náhodné poloviny, ne dál. Generátory bez jazyka dávají nulu.
  V latině je obraz stejný: slova jen v jiné knize / jiném autorovi jsou od slovníku druhé strany dál než slova
  jen v jedné polovině (−0,07 / −0,11), a 80 % B-only tokenů nemá protějšek ≤ 1 vůbec.
* Opačný směr: A-only tokeny vs B mají protějšek ≤ 1 v 0,915 (ZL) / 0,920 (GC) ≈ poloviny (0,921 / 0,935,
  rozdíl −0,006 / −0,016, interval přes nulu). Asymetrie: tvary **jen v B** jsou od A dál, tvary jen v A od B ne
  (B obsahuje A + něco navíc; viz inventář).

## (ii) Inventář editací (`edit_inventory.csv`, `edit_coverage.csv`)

Operace dvojic B-only → A (zapsané A → B), celkem, inventář sdružený přes 20 losů; `ops_80` = počet operací,
které vysvětlí 80 % tokenů dvojic (celý skript v top-k); `Currier` = podíl tokenů na operacích složených
jen z `e`, `ee`, `d`, `q` nebo `ch`↔`sh`.

| sada | rep | n op. | ops 50 / 80 / 95 % | top-1 (podíl) | top-3 | Currier podíl |
|---|---|---|---|---|---|---|
| ZL split A×B | token | 283 | **4 / 30 / 94** | **`ins:d@mid` 0,46** | `sub:o>l@begin` 0,06, `ins:l@begin` 0,05 | 0,45 |
| ZL split poloviny | token | 258–304 | 30–32 / 92–96 / 175–193 | `ins:p@begin` 0,035–0,04 | `ins:e@mid`, `ins:d@begin` | 0,11–0,12 |
| ZL split A×B | kmen (a) | 134 | 6 / 24 / 64 | `ins:d@end` 0,27 | `ins:e@begin` 0,06, `ins:l@begin` 0,06 | 0,50 |
| ZL split poloviny | kmen (a) | 114–130 | 19–20 / 47–48 / 82–89 | `ins:o@end` / `ins:d@end` 0,04–0,07 | `ins:r@begin`, `ins:p@begin` | 0,15–0,21 |
| ZL split A×B | BPE | 585 | 20 / 77 / 258 | `sub:daiin>qokeedy@begin` 0,07 | `sub:daiin>chedy`, `ins:lk@begin` | 0,02 |
| ZL join A×B | token | 289 | 5 / 32 / 100 | `ins:d@mid` 0,45 | `ins:l@begin`, `sub:o>l@begin` | 0,43 |
| GC split A×B | token | 352 | **6 / 36 / 128** | **`ins:8@mid` 0,28** | `ins:7@mid` 0,16, `sub:o>e@begin` 0,05 | – |
| GC split poloviny | token | 381–395 | 33–37 / 113–116 / 210–252 | `sub:o>9@begin` 0,03 | `ins:c@mid`, `ins:g@begin` | – |
| GC split A×B | kmen (a) | 182 | 12 / 41 / 95 | `ins:8@end` 0,125 | `ins:c@begin` 0,07, `ins:7@end` 0,06 | – |
| Plinius A×B | token | 672 | 107 / 266 / 463 | `ins:r@mid` 0,04 | `ins:s@end`, `ins:d@begin` | – |
| Caesar/Celsus A×B | token | 849 | 123 / 304 / 538 | `ins:r@mid` 0,04 | `ins:c@begin`, `ins:b@mid` | – |

* Po losech (`runs.csv`): `ops_80` A×B token ZL 22 [16; 28], GC 29 [22; 38]; poloviny 47–58; ts2020 65, šablona
  69–84 (jejich protějšky jsou šum, inventář plochý).
* **Odpověď (ii): 80 % nese ~30 operací (po losech 22–29), ale jediná – vložení `d` uvnitř slova (ZL `chey → chedy`,
  `okey → okedy`; GC `ins:8@mid` = totéž ve v101, plus `ins:7@mid`) – nese 28–46 % protějšků, 5–10× víc než
  nejčastější operace mezi polovinami.** Ano, jsou to Currierovy znaky: `ed`/`-dy` dominuje, u kmenů (a)
  `ins:d@end` (= sufix `-dy` vzniklý z `-y`, který parser nechal v kmeni) a `ins:e@begin` / `ins:ee@begin`
  (`ed`, `eed`); Currierovy operace nesou 45–50 % tokenů dvojic (poloviny 11–21 %). Zbytek inventáře je dlouhý
  ohon (95 % až ~100 operací), tj. **jedna systematická operace + mnoho nesystematických**. Krit. 3 (jediná
  editace > 80 %) nenastalo. V BPE jednotkách `edy`/`dy` splývají do jedné jednotky, takže totéž vychází jako
  záměna `y>edy`, `daiin>chedy` (0,07) – inventář se rozpadá na drobné operace (ops 80 % = 77).

## (iii) Vrátí substituce sdílení na úroveň polovin? (`substitution_sharing.csv`, `substitution_sharing.png`)

Jaccard inventářů (≥ 2 tokeny) po přepsání B-only kmenů, jejichž skript leží v top-k operacích, na A protějšek
(celkem, split; poloviny týmž postupem):

| sada | rep | k = 0 | 1 | 3 | 5 | 10 | všechny ≤ 2 | poloviny k = 0 | poloviny všechny |
|---|---|---|---|---|---|---|---|---|---|
| ZL | token | 0,248 | 0,274 | 0,284 | 0,293 | 0,308 | 0,394 | **0,421–0,440** | 0,522–0,539 |
| GC | token | 0,236 | 0,253 | 0,269 | 0,278 | 0,292 | 0,381 | **0,415–0,424** | 0,512–0,524 |
| ZL | kmen (a) | 0,371 | 0,383 | 0,392 | 0,400 | 0,418 | 0,477 | **0,493–0,511** | 0,585–0,601 |
| GC | kmen (a) | 0,334 | 0,340 | 0,350 | 0,357 | 0,371 | 0,442 | **0,474–0,466** | 0,562–0,560 |
| ZL | BPE | 0,248 | 0,249 | 0,252 | 0,256 | 0,267 | 0,391 | 0,421–0,440 | 0,516–0,542 |
| ts2020 ZL | token | 0,292 | 0,297 | 0,305 | 0,312 | 0,325 | 0,455 | 0,291–0,313 | 0,452–0,467 |
| Plinius | token | 0,200 | 0,201 | 0,202 | 0,204 | 0,207 | 0,312 | 0,301–0,315 | 0,390–0,409 |
| Caesar/Celsus | token | 0,130 | 0,130 | 0,131 | 0,132 | 0,135 | 0,269 | 0,399–0,406 | 0,484–0,490 |

* **Odpověď (iii): žádný k ≤ 10 nevrátí sdílení na úroveň polovin.** Top-1 (`ins d`) přidá u tokenů 0,026 (ZL) /
  0,017 (GC) Jaccardu, top-10 0,06, zatímco mezera k polovinám je 0,17–0,19. Ani sloučení **všech** dvojic
  do 2 editací (≈ 300 operací, horní mez) na poloviny bez přepisu nedosáhne (0,39 vs 0,42–0,44; GC 0,38 vs 0,42)
  a poloviny stejným přepisem utečou na 0,52–0,54. U kmenů (a) sloučení všech dvojic dojde těsně pod poloviny
  (0,48 vs 0,49–0,51), ale to je artefakt krátkých kmenů (viz (i)). Latina: jiné téma / registr se chová stejně
  (top-10 nic, „všechny“ zdaleka ne).

## (iv) Zachovávají protějšky sekční profil? (`section_check.csv`)

Plný text, JSD profilu kmene v A vs profilu protějšku v B proti 500 permutacím přiřazení; referenční hladina =
**týž kmen** v A a v B. Profil `herbal_binary` (podíl tokenů v herbáři vs mimo něj):

| sada | rep | B-only ≤ 1: n, z, p | A-only ≤ 1: z, p | týž kmen: z, p |
|---|---|---|---|---|
| ZL split | token | 584, −1,15, 0,12 | +1,07, 0,84 | +0,54, 0,70 |
| ZL join | token | 578, +0,10, 0,54 | +0,74, 0,77 | +0,45, 0,64 |
| GC split | token | 798, −0,88, 0,19 | +0,13, 0,54 | +0,75, 0,77 |
| GC join | token | 738, −0,77, 0,23 | +0,44, 0,66 | +0,87, 0,81 |
| ZL split | kmen (a) | 177, −0,75, 0,25 | +0,73, 0,77 | −0,04, 0,48 |
| GC split | kmen (a) | 286, +0,59, 0,72 | −2,01, 0,02 | −0,14, 0,44 |

* Profil přes všechny sekce (`sections`) je ovládán tím, že A a B mají jiné sekce (A: herbář, farmacie; B: lázně,
  hvězdy, herbář, text), JSD ≈ 0,66–0,72 pro reálné i permutované dvojice, z ≈ 0 i u sdílených kmenů.
* **Odpověď (iv): nelze rozhodnout – test nemá sílu.** Ani týž kmen nemá v A a v B podobnější rozložení
  herbář / mimo herbář než náhodná dvojice kmenů (z −0,1 až +0,9 ve všech sadách), takže z protějšků (z −1,2 až
  +1,1, jedna izolovaná hodnota −2,0 u GC A-only (a) bez opakování v jiné sadě) nelze nic odvodit. To je v souladu
  s 006 (sdílení kmenů řídí blízkost, ne sekce): sekční rozložení kmene není mezi A a B stabilní vlastností.

## Nully a kontroly (shrnutí)

* ts2020 i šablona: `share_only`, protějšky ≤ 1, inventář i substituce A×B = poloviny (rozdíly v intervalech přes
  nulu); jejich „A/B“ je šum. Inventář je plochý (top-1 0,04–0,06, `ops_80` 65–84).
* Latinské kontroly: slovník jen v jiné knize (Plinius) nebo jiném autorovi (Caesar/Celsus) je od slovníku druhé
  strany **dál** než slovník jedné poloviny (−0,07 / −0,11 u ≤ 1), inventář plochý (top-1 ≤ 0,04, `ops_80`
  107–123), substituce nic nevrátí. Registr (Caesar/Celsus) zvětšuje B-only podíl víc (0,30) než téma (0,19).
* Join a split, ZL a GC: stejné znaménko i velikost u všech čtyř odpovědí.

## Interpretace

Kmeny (celé tvary) vyskytující se jen v B nejsou převážně glyfové varianty tvarů A: mají v A *méně* tvarů do jedné
editace než tvary jedné poloviny jazyka v druhé polovině, a přepis podle naučených operací – ani všech ~300 – sdílení
na úroveň polovin nevrátí. Existuje však jedna systematická operace, vložení `d` uvnitř slova (`-ey` → `-edy`, ve v101
`8`), která nese 28–46 % protějšků a 45–50 % tokenů dvojic připadá na Currierovy znaky; ta je ale sama o sobě málo
(Jaccard +0,02–0,03). Věta: **obojí, s převahou nového slovníku** – přibližně jeden z pěti B-tokenů (0,21–0,23)
je tvar jen v B (poloviny 0,05–0,06); z toho se asi 30–45 % dá vysvětlit jedinou záměnou zápisu (`ins d`, podíl
`ins d` × podíl s protějškem), zbytek jsou tvary bez systematického protějšku. Latinský text se stejným
jazykem a jiným tématem/registrem vykazuje tentýž vzorec (B-only slova dál než poloviny, plochý inventář), takže
nález je slučitelný i s „stejný systém, jiný obsah“, i s „dva pisatelské okruhy“; generátory bez jazyka ho nedávají.

## Co výsledek neukazuje

* Neříká, *co* B-only tvary jsou (obsah, jiný písař, jiný slovník téhož systému) – jen že nejsou jednoduchým
  přezápisem A. Jazyk a písař zůstávají zmatené (012).
* Kmeny dekompozice jsou pro Levenshtein ≤ 2 příliš krátké (poloviny > 0,9, krit. 3); nálezy (i) a (iii) stojí
  na celých tokenech a BPE, otázka „kmen vs afix“ se tím obchází, ne řeší.
* Substituce je licencovaná (přepisuje jen kmeny, jejichž skript sedí), ne plošná aplikace pravidla; je to horní mez
  pro hypotézu změny zápisu, pravidlo s kontextem (`-ey` → `-edy` jen po některých glyfech) by dalo méně, ne víc.
* Sekční test nemá sílu (ani sdílené kmeny), odpověď (iv) je „nevíme“, ne „ne“.
* Inventář operací závisí na jednotkách: v BPE jednotkách se `ins d` rozpadne na záměny celých koncovek; výroky
  o „jedné operaci“ platí pro glyfy presetu `extended` / znaky v101.

## Další krok

* Rozdělit B-only tvary na (1) `-edy`/`-dy` varianty (protějšek přes `ins d`) a (2) zbytek, a pro (2) spočítat
  profil 009/011 (hapax, pozice v řádku, sufixy): chovají se jako nová slova katalogu (jména/položky), nebo jako
  tvary s jinou morfologií? Totéž pro A-only (menší strana).
* Kontextové pravidlo `-ey → -edy`: po kterých glyfech B vkládá `d` a je to deterministické (změna zápisu) nebo
  pravděpodobnostní (jiný tvar téhož paradigmatu, 011)?
