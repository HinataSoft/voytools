# 013 – deník session

## Data

* **Blábol:** `git clone --depth 1 https://github.com/danielgaskell/voynich` do `../ext/gaskell` (cesta přes
  `VOYNICH_GIBBERISH`); `data/gibberish_transcriptions.zip` = 38 souborů `Gibberish - D{A,C}_NN.txt`
  (CRLF, UTF-8). Formát: řádky = řádky rukopisu pisatele, řádek `__________` = předěl stránky (G&B ho
  v `stats.py` používají jako „page marker"), některé řádky začínají odrážkou `*`. 38 pisatelů, 10 032 tokenů
  po čištění, 1 617 řádků, 88 stránek; na pisatele 79–480 tokenů (průměr 264), 11–76 řádků, 1–4 stránky
  (`writers_inventory.csv`). README uvádí 42 dobrovolníků – v zipu je 38 přepisů (ostatní zřejmě vyřazeny).
  Licence: modified MIT, citovat Gaskell & Bowern 2022 (CEUR WS, ICVM 2022).
* **Čištění** jako G&B `stats.py`: malá písmena, interpunkce (`string.punctuation`) a číslice smazány
  *uvnitř* tokenu (`co-nallol` → `conallol`), tokeny podle mezer. Glyf = písmeno (preset `None`).
  Implementace `voytools.controls.gibberish_words`, `parse_gibberish`, `read_gibberish` (zip i složka).
* **Naibbe:** z klonu `lrozanova/voynich-units` (`data/controls/naibbe/`, PROVENANCE.md): Caesar ciphertext
  (driver R&T, 80 877 tokenů, 9 985 řádků), Greshkův ciphertext Plinia XVI (34 764 tokenů) a k němu
  zarovnaný „respaced plaintext" (34 764 kousků po 1–2 písmenech – **to nejsou latinská slova**, ale
  šifrovací jednotky; latinský Plinius jako slova je kontrola `pliny_prose` z 009). Ciphertext je Eva →
  stejné presety jako Voynich (`extended` pro glyfové metriky).

## Vysázení

* **Blábol sdružený:** 38 pisatelů za sebou v náhodném pořadí (seed 13), jako souvislý text (prose) do layoutu
  ZL split **oříznutého na velikost blábolu** (`controls.crop_template`: stránky s konstantním krokem přes celý
  rukopis, 60 stránek, 10 024 tokenů); sekce = písař = většinový pisatel stránky (30 z 38 pisatelů je většinou
  na aspoň jedné stránce). Řádky jsou tu mechanické.
* **Blábol ve vlastním layoutu** (`gibberish_native`, `controls.lines_table`): vlastní řádky a stránky pisatelů
  – pro řádkové metriky 010 (konec řádku, uzavřený řádek). Afixy naučené na sdruženém korpusu.
* **Srovnání na stejné velikosti:** `voynich_crop` (ZL split ve stejném oříznutém layoutu) a `caesar_crop`,
  `celsus_crop` (kontroly 009 v oříznutém layoutu) – velikost slovníku, hapax a šířka paradigmatu na velikosti
  závisí, proto blábol srovnávám s Voynichem a prózou na ~10 000 tokenech, ne s plným ZL.
* **Naibbe:** text od začátku ořezaný na 34 268 tokenů, 7 souvislých bloků = sekce (knihy nejsou v souborech
  vyznačeny), souvislý text v plném layoutu ZL split.

## Metriky a srovnání

* Celý profil 009 (`M9.profile`, bod + 100 subsamplů), řádkové metriky 010 (`M10.line_profile`, bod + 100)
  a rozklad 011 (`M11.measure`, bod + 30; jen dekompozice (a)) – stejný kód jako v 009–011.
* Referenční zdroje převzaté z `results/009..011` (nepočítány znovu): ZL/GC split, ts2020, šablona, ac86,
  Caesar, Celsus, Plinius, Culpeper, Linnaeus, italská próza.
* **Vzdálenosti** metodou 009 (RMS standardizovaných rozdílů 37 metrik), SD přes všechny zdroje v tabulce
  (bez míchání uvnitř řádku); páry podle spec.
* **Seznamy vlastností:** pravidlo, které spec nedefinuje – podíl f = (zdroj − R) / (Voynich − R), R = průměr
  smysluplné prózy na stejné velikosti; ano f ≥ 0,5, částečně 0,2–0,5, ne < 0,2; f i na mezích intervalu
  zdroje. „A/B uvnitř písaře" nelze na blábolu ani Naibbe měřit přímo – analog: podíl H(sufix | kmen)
  vysvětlený pisatelem (011). BPE žlab není v profilech 009–011 → nepočítáno.
* **Per pisatel** (`writers.csv`): slovník, hapax, H1/H2 glyfů (uvnitř tokenu), délka slova, sufixů na kmen a
  H(sufix | kmen) pro kmeny s ≥ 5 tokeny, rozdíl délky posledního slova řádku; totéž pro 38 souvislých úseků ZL
  stejných velikostí (začátky rovnoměrně přes rukopis) – rozptyl mezi pisateli vs mezi úseky Voynichu.

## Průběh a uzavření

* Všech 24 běhů (8 zdrojů × 009/010/011) doběhlo bez chyb (`runs/p009|p010|p011/<zdroj>/`).
* Pravidlo seznamů dostalo dvě pojistky po prvním průchodu: Voynich ≈ próza (|V − R| < 2 SD Voynichu, např.
  `loc_log_dist`) → „nelze rozlišit"; řádkové metriky u zdrojů s mechanickými řádky se neposuzují.
* 011 je v této větvi hotové (není ještě v `main`); rozklad 011 jsem použil přímo (`M11.measure`), jak zadal
  uživatel.
* Stav: **hotovo**. Omezení dat (krátký blábol, mechanické řádky Naibbe) jsou popsána v reportu.
