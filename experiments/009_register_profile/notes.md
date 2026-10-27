# 009 – deník session

## Kontroly: co bylo v klonu `lrozanova/voynich-units` (commit z 27. 10. 2026, `--depth 1`)

| Spec | Nalezeno | Použito |
|---|---|---|
| latinský narativ (Caesar) | `data/voynich-units/latin/caesar.de_bello_gallico.part.1–3.tess` (knihy I–III, jen 17 765 slov – málo na 34 268 tokenů); `voynich_decipherment_repro_bundle/decipherment_attack_v5/latin_controls/caesar_la.txt` (Gutenberg #218, knihy I–IV) a `caesar_la2.txt` (knihy V–VIII), dohromady 51 370 tokenů | Gutenberg I–VIII (`caesar_prose`) |
| lékařská latina (Celsus) | `data/voynich-units/latin/celsus.de_medicina.part.1–5.tess` (knihy 1–5, 62 410 tokenů) | ano (`celsus_prose`) |
| botanická latina (Plinius) | `data/voynich-units/herbal/pliny_the_elder.naturalis_historia.part.4–6.*.tess` (knihy 16–30, 164 848 tokenů) | ano (`pliny_prose`) |
| anglický herbář (Culpeper) | `data/voynich-units/herbal/culpeper_en.txt` (Gutenberg #49513) | ano, dvakrát (`culpeper_prose`, `culpeper_items`) |
| italská próza | `data/voynich-units/italian/27608.txt` (Bollettino del Club Alpino Italiano 1895–96, 102 755 tokenů), `26166.txt` (Murri, 21 tis. slov – krátké) | `27608.txt` (`italian_prose`) |
| katalog (Linnaeus) | `voynich_decipherment_repro_bundle/decipherment_attack_v5/latin_controls/species1.txt` + `species2.txt` (*Species Plantarum*, třídy I–V, 56 505 tokenů těla) | ano, dvakrát (`linnaeus_prose`, `linnaeus_items`) |
| volitelně Boxer: `ellis_1854` | nalezeno (sparse checkout, commit 255bd63) – jen 10 565 slov, převážně zkratky receptů (`gr ij`, `misce`) | **ne** (krátké na layout 34 268 tokenů) |
| volitelně Boxer: `legacy/pliny` | `pliny_lat0.csv` (Perseus, celá NH) | ano, druhá kopie Plinia (`plinyP_prose`, stejné knihy 16–30) |

Další soubory v klonu (Naibbe, self-citation, čínské korpusy, `lm_*` jazykové modely) nejsou v zadání a nepoužity.

## Rozhodnutí

* **Tokenizace** (`voytools.controls.words`): malá písmena, slovo = maximální běh písmen (Unicode),
  číslice a interpunkce pryč; apostrof dělí (`dell'acqua` → `dell`, `acqua`); `æ`, `œ` zůstávají jako písmena.
  U Linnéa zůstávají zkratky citací (`roy lugdb`, `fl zeyl`, `bauh pin`) – jsou součástí katalogového registru.
* **Glyfy kontrol** = písmena (`glyph_preset=None`); u Voynichu ZL `extended` pro hranu a jednotky (jako 007),
  GC v101 znaky.
* **Sekce.** Spec chce přirozené oddíly „namapované na stránky podle pořadí“. Prosté oříznutí začátku dává
  2–4 oddíly (Plinius knihy 16–18, Celsus 1–3), zatímco ZL má 7 kódů `illustration`. Zvolil jsem
  **vyvážený výběr** (`controls.build_control`, `strategy='auto'`): řádky šablony se rozdělí na 7 bloků se
  stejným počtem tokenů, každý blok se naplní souvislým textem od začátku jednoho oddílu, oddíly se berou
  rovnoměrně rozložené po zdroji; příliš krátké sousední oddíly se slučují (`I MONOGYNIA+…`). Kde to
  nedá aspoň 4 sekce, použije se souvislý text od začátku zdroje (`sequential`, sekce = sloučené oddíly).
  Sekce stránky = většinový oddíl jejích řádků. Výsledek viz `results/009_register_profile/controls_used.csv`.
  - Linnaeus: oddíl = řád (třída × `MONOGYNIA` …), ne rod – rodů je ~340, tj. skoro jeden na stránku, což
    by sekční metriky nesrovnatelně nafouklo. Tělo od `_Classis I._` po první řádek hvězdiček.
  - Culpeper: herbář je abecední podle klíčového slova (`WILD TANSY` je pod T), oddíl = klouzavý (dolní)
    medián počátečních písmen nadpisů přes 9 hesel. Herbář od `AMARA DULCIS.` po `DIRECTIONS FOR MAKING`.
  - Italština: oddíly nejsou rozpoznatelné → 14 stejných bloků, blok 1 (titulní strana, kredity) vynechán.
  - Caesar: knihy (`COMMENTARIUS PRIMUS` … / `Liber V` …); vyvážení nevychází (knihy 3,6–11,5 tis.) →
    souvislý text knih I–VIII, 4 sekce.
* **Položkové vysázení** (`mode='items'`): každý řádek začíná novou položkou; kratší položka pokračuje
  dalšími, delší se na konci řádku utne (zbytek zahozen). Linnaeus: položka = odstavec (nadpis rodu,
  postranní jméno druhu `[_Indica._]`, diagnóza, synonymum, `Habitat`) – zdroj stačí jen na souvislé
  vysázení (4 sekce). Culpeper: odstavců (~2 000) je méně než řádků (4 128) → položka = věta/klauze
  (dělení za `. ; : ? !`), 6 sekcí. Culpeper „po heslech“ tedy nešel; použitá varianta je „po větách“.
* **Metadata kontrol**: stránky, odstavce, pozice v řádku, `language` a `hand` z reálu (spec); proto
  např. lokalitní model 006 obsahuje i `same_language`/`same_hand` – u kontrol jsou to jen kovariáty layoutu.
* **Afixy**: u každého korpusu `decompose_table` (tj. `learn_affixes`, 12 prefixů/sufixů) na celém textu;
  subsamply a nully téhož zdroje používají stejné seznamy (nully Voynichu reálné seznamy, jako 007).
* **Řádkové metriky (skupina 7)** – jednoduchá verze: třída = afixový rámec prefix|sufix s ≥ 20 výskyty
  (jednotka `affix` z 007); `first/last_top5_excess` = podíl 5 nejčastějších tříd na 1./posledním místě
  minus totéž uprostřed řádku; `jsd_first/last_mid`; `I_first_rest` = MI mezi třídou prvního tokenu a
  třídou každého dalšího tokenu řádku, korigovaná permutací prvních tříd mezi řádky.
* **Lokalita**: model `spec+size` z 006 (kosinus, K = 30 funkčních kmenů vyřazeno), koeficienty
  `log_dist`, `same_section` a kontrast po odečtení velikosti (`contrast_adj`); pořadí stránek ze souboru.
* **Shuffly**: bod 50 (MI tokenů, hrana, bias korekce, MI(d)), subsample/nully 50 pro MI tokenů a hranu,
  10 pro bias a MI(d), 20 pro `I_first_rest`.
* **Nully**: ts2020 calB jen Eva (ZL split + join); šablona 005 = plný model (rank 1 fáze B + α = 5,
  λ = 0,3 z `ablation_base`), ZL i GC split; ac-tri 86 = vzorky 8600–8629 jako 007 (generované znovu,
  cache v `scratch/`); míchání v řádcích ZL i GC split.
