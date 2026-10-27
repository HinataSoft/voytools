# 005 – Null „strukturovaný umělý jazyk" a robustnost experimentu 002

Kód: `experiments/005_structured_template_null/run.py`, generátor `voytools/templates.py` (`ParadigmTemplate`).
Rozhodnutí a odchylky od spec: `experiments/005_structured_template_null/notes.md`. Všechna čísla níže jsou ze
souborů v této složce: `best_configs.csv`, `grid_summary.csv`, `ablation.csv`, `ingredients.csv`,
`robustness.csv`, `robustness_verdicts.csv`, `ab_by_hand.csv`, `summary_info.json`, `runs.csv`
a po částech `grid/`, `stage2/`, `phaseB/`, `ablation/`, `robust/`, `ab/`, `real/`.

## Co bylo spočítáno

* **Data** jako v 002: ZL3b a GC2a, `P`, `clean`, politika `split` (hlavní), `join` (kontrola nejlepší konfigurace).
  Afixové seznamy se učí na reálném textu (`decompose_table`) a jsou stejné pro reál i všechny vzorky. Metriky 002
  jsou zkopírované beze změny. Na reálném textu dávají totéž co `results/002_paradigm_coherence/runs.csv`
  (max. rozdíl 4·10⁻⁶, jen zaokrouhlení CSV; `check_002.json`).
* **Šablona** (`ParadigmTemplate`): reálné inventáře prefixů a sufixů s reálnými frekvencemi. Kmeny jsou
  **umělé**, generuje je slotová/trigramová gramatika naučená na reálných kmenech, se Zipfovými vahami.
  Kmenů, které náhodou splývají s reálnými, je v ablacích ≈ 0,45–0,52 (`share_stems_real`); jsou to hlavně
  krátké řetězce. Každý kmen má pevnou sadu sufixů velikosti k, k se bere z reálných kmenů se stejnou
  frekvencí (`paradigm='fixed'`). Dále: téma α (Dirichlet po sekcích), lokalita stránky λ a volitelná vazba
  prefixu. Žádné kopírování tokenů.
* **Krok 1:** mřížka ze spec 216 konfigurací × 3 seedy na ZL i GC. Navíc **rozšíření** o 18 konfigurací
  (`zipf_s` 1,4/1,5/1,7, `n_stems` až 16 000), protože všechna nejlepší řešení ležela na hraně mřížky
  `zipf_s = 1,3` (odchylka od spec). Fáze B: 5 nejlepších × 50 seedů (u ZL 5 ze spec mřížky a 5 z rozšíření,
  celkem 10).
* **Krok 2:** ablace a „adice" (přidání jednoho ručního pravidla) plného modelu, 9 variant × 50 seedů × 2 soubory.
  Plný model = rank 1 fáze B plus α a λ, které cílové metriky neurčují. Ty jsou vybrané druhým krokem
  (sweep α × λ, 12 × 3 seedy) podle `jsd_stem_sections` a `stem_dispersion`. Výsledek u obou souborů: α = 5, λ = 0,3.
* **Krok 3:** reál = celý text + 200 podvzorků 80 % stránek bez vracení, slotový null 50×, `autocopy_536` 50×.
  Šest variant dekompozice: základ (14 sufixů, délka ≤ 4), 8 a 20 sufixů, délka ≤ 3, parser se sufixem napřed
  (`decompose(order='suffix_first')`), afixy = BPE jednotky (64 sloučení).
* **Krok 4:** JSD sufixových rozdělení mezi skupinami jazyk × písař, interval z podvzorků stránek, základní
  linie z permutace štítků stránek (200×). Totéž na šablonách (plný model, bez tématu, s pravidlem
  „sufixy podle jazyka").
* Celkem 3 190 běhů šablony (`summary_info.json`), 602 běhů robustnosti (každý vyhodnocený v 6 variantách dekompozice) a 804 běhů A/B.
* **Intervaly:** u šablon `z = (reál − průměr šablony) / SD šablony` přes 50 seedů. „Reprodukuje" znamená |z| < 3.
  U podvzorků se SD celého textu odhaduje jako SD podvzorku · √(0,8/0,2) = 2× (delete-d jackknife).
  Intervaly rozdílů jsou ± 1,96·√(SD_reál² + SD_null²).

## Krok 1 – kalibrace (kritérium 1)

| text | konfigurace (rank 1, fáze B) | loss | největší složka |
|---|---|---|---|
| ZL split | 8 000 kmenů, `zipf_s` 1,4, fixed, prefix none (id 2006) | **0,072** | n_stems_k 0,012, mean_suffixes 0,016 |
| GC split | 4 000 kmenů, `zipf_s` 1,3, fixed, prefix none (id 122) | **0,046** | hapax_share 0,022 |
| ZL join (kontrola, id 2006) | – | 0,123 | n_stems_k 0,041 |
| GC join (kontrola, id 194) | – | 0,060 | hapax_share 0,020 |

Kritérium 1 je **splněno**, ale u ZL až po rozšíření mřížky. Ve spec mřížce byla nejlepší ZL konfigurace
ve fázi B na loss 0,182 (id 204, `zipf_s` 1,3, prefix stem), tedy nad limitem 0,15. Chyběl hlavně podíl hapaxů
(složka 0,050 = relativní odchylka 22 %). Loss je prakticky plochá přes α a λ: cílové metriky téma ani
lokalitu neurčují (`grid_summary.csv`).

Napětí v kalibraci: šablona trefí `n_types`, `n_stems_k` i podíl kmenů s více sufixy, ale hapaxů má soustavně
méně než reál (ZL 0,610 vs 0,676, GC 0,621 vs 0,702; `ablation.csv`, varianta `best`, z = 10,6 / 12,4).
Reálný slovník má delší chvost jednorázových tvarů, než dává Zipf × pevné paradigma.

## Krok 2 – metrika × ingredience (kritérium 2)

Plný model (téma α = 5, lokalita λ = 0,3, pevné paradigma, bez vazby prefixu), z reálu proti 50 seedům.
Úplná tabulka je v `ingredients.csv` a `ablation.csv`, graf v `template_vs_real.png`.

| metrika 002 | reál ZL / GC | šablona ZL / GC | z ZL / GC | reprodukuje? | co ji rozbije / co ji spraví |
|---|---|---|---|---|---|
| frac_stems_multi_suffix | 0,912 / 0,839 | 0,906 / 0,844 | 0,5 / −0,3 | **ano** | rozbije `free_paradigm` (z −28 / −62) |
| frac_multi_suffix_5to9 | 0,832 / 0,777 | 0,835 / 0,747 | −0,1 / 1,0 | **ano** | `free_paradigm` (−31 / −39) |
| mean_suffixes_per_stem | 5,25 / 4,72 | 4,44 / 4,23 | 7,0 / 3,9 | ne (šablona o 10–15 % užší) | žádná adice |
| cmi_excess | 0,088 / 0,089 | 0,030 / 0,040 | 2,9 / 2,3 | **formálně ano**, ale průměr šablony je jen ⅓–½ reálu; SD šablony 0,02 je velké | `add_prefix_given_suffix` ho přestřelí (0,20 / 0,17, z −4,7 / −3,9) |
| mi_prefix_suffix | 0,205 / 0,176 | 0,011 / 0,010 | 18,7 / 22,6 | **ne** | spraví `add_prefix_given_suffix` (0,183 / 0,151, z 1,8 / 2,5), ale pak přestřelí cmi_excess |
| H_suffix_given_stem | 1,82 / 1,80 bitu | 2,78 / 2,69 | −9,1 / −8,6 | **ne**, žádná varianta | – |
| stem_dispersion | 0,459 / 0,481 | 0,500 / 0,491 | −2,5 / −0,6 | **ano** | rozbije `no_locality` (0,65; z −21 / −18), u ZL i `no_topic` (z −5,0) |
| jsd_stem_sections | 0,227 / 0,271 | 0,282 / 0,299 | −1,6 / −0,9 | **ano** | rozbije `no_topic` (0,12 / 0,15; z 19 / 26) |
| jsd_suffix_sections | 0,046 / 0,037 | 0,009 / 0,008 | 6,1 / 9,5 | ne | spraví `add_suffix_by_language` (0,033 / 0,024; z 2,1 / 3,0) |
| jsd_suffix_AB | 0,088 / 0,065 | 0,0006 / 0,0006 | 269 / 165 | **ne** | `add_suffix_by_language` dá 0,068 / 0,042 (z 4,4 / 7,1; asi ¾ posunu) |
| jaccard_stem_sections | 0,325 / 0,268 | 0,257 / 0,239 | 4,1 / 2,2 | jen GC | `no_locality` ho srovná (1,6 / −0,4), ale rozbije disperzi |
| productivity_slope | 0,309 / 0,328 | 0,372 / 0,426 | −4,0 / −6,2 | ne | jen `free_paradigm` (0,8 / 0,2), ta rozbije vše ostatní |
| n_stems_k, mean_len | – | – | |z| < 1,4 | ano (cílové) | `no_locality` posune n_stems_k (z 3,1 / 4,0) |
| hapax_share, n_types | – | – | 10,6 / 12,4; −4,3 / −4,7 | ne (cílové, viz krok 1) | – |

**Kolik ručních pravidel umělý systém potřebuje:**

1. Pevné paradigma kmene (podíl kmenů s více sufixy).
2. Téma po sekcích (`jsd_stem_sections`).
3. Lokalita stránky (disperze 0,46–0,48).

Tyto tři ingredience stačí na to, co 002 bralo jako „paradigma, ne šablona" a „rozptýlené, ne kopírované".

4. Vazba prefix–sufix (P(prefix | sufix)).
5. Frekvence sufixů podle jazyka A/B.

Pravidla 4 a 5 jsou nutná pro `mi_prefix_suffix`, `jsd_suffix_sections` a `jsd_suffix_AB`.

Ani s pěti pravidly (`add_all`) šablona nereprodukuje:
* `H_suffix_given_stem` – reálné kmeny mají jeden převažující sufix, šablona rozkládá tokeny po sadě rovnoměrněji;
* `mi_prefix_suffix` a `cmi_excess` současně – vazba přes P(prefix | sufix) dává správné nepodmíněné MI,
  ale 2× větší podmíněné;
* hapaxy a `n_types` současně (`add_all` hapax 0,48, z 23 / 35);
* sklon produktivity.

Šesté a sedmé pravidlo by tedy muselo být „nerovnoměrné váhy sufixů uvnitř kmene" a „vazba prefixu závislá
na kmeni i sufixu jinak než globálně". Každé z nich je zase jen parametr.

## Krok 3 – robustnost 002 (kritérium 3)

`robustness_verdicts.csv`: všech 7 dílčích testů platí v **každé** variantě (6 variant × ZL a GC).

* **(i) stabilní.** `cmi_excess` reálu je nad slotovým nullem a pod `autocopy_536` ve všech variantách.
  Například ZL: 8 sufixů 0,069 vs 0,032 / 0,163; 20 sufixů 0,102 vs 0,036 / 0,225; suffix-first 0,104 vs
  0,029 / 0,203; BPE 0,022 vs 0,008 / 0,076. `mi_prefix_suffix` reálu je nad oběma nully všude (ZL 0,128–0,260,
  BPE 0,035). Šíře paradigmatu (`mean_suffixes_per_stem`, `frac_multi_suffix_5to9`) leží všude mezi nully.
  Nejtěsněji to vychází u ZL `len3`: 5–9 kmenů 0,662 vs autocopy 0,539, interval rozdílu [0,022; 0,223].
* **(ii) stabilní, s výhradou z 002.** Rozdíl Jaccardů (afixy − kmeny) je kladný všude (0,61–0,72), ale ve všech
  variantách se neliší od slotového nullu (intervaly obsahují 0). Je to strukturální vlastnost uzavřeného
  seznamu afixů, ne nález. Informativní `jsd_stem_sections` leží všude mezi slotem a autocopy (ZL 0,19–0,34).
  `jsd_suffix_sections` je všude nad slotem.
* **(iii) stabilní.** `stem_dispersion` je 0,456–0,490 ve všech variantách, slot 0,685–0,697, autocopy 0,03–0,10.
* Velikosti efektů na variantě závisí. Podle očekávání rostou s počtem sufixů a BPE afixy je zmenšují asi 4×
  (jednotky BPE jsou delší a vzácnější). U GC je `len3` totožné se základem, protože naučené afixy GC mají
  nejvýš 2 znaky.

## Krok 4 – posun sufixů A/B versus písař (kritérium 4)

`ab_by_hand.csv`. JSD sufixů a v závorce přebytek nad permutací stránek [95% interval]:

| kontrast | ZL | GC |
|---|---|---|
| A vs B (vše) | 0,088 (0,087 [0,071; 0,102]) | 0,065 (0,063 [0,052; 0,075]) |
| A1 vs B2 (jiný jazyk i písař) | 0,103 (0,100 [0,086; 0,114]) | 0,086 (0,083 [0,072; 0,094]) |
| A1 vs B3 | 0,093 (0,088 [0,074; 0,102]) | 0,060 (0,056 [0,044; 0,067]) |
| **B2 vs B3 (stejný jazyk, jiný písař)** | 0,009 (0,007 [−0,001; 0,015]) | 0,013 (0,010 [0,003; 0,018]) |
| B2 vs B5 / B3 vs B5 (ruka 5: 838 tok.) | 0,020 / 0,037 (intervaly přebytku obsahují 0, resp. jsou na hraně) | 0,023 / 0,042 |
| **písař 3: A vs B** (733 / 11 102 tok.) | 0,068 (0,048 [0,034; 0,063]) | 0,063 (0,045 [0,027; 0,063]) |
| sekce uvnitř písaře A1 / B2 / B3 | 0,027 / 0,030 / 0,011 | 0,011 / 0,029 / 0,015 |

**Posun sufixů A/B kontrolu písaře přežívá.** Dva písaři jazyka B se liší asi 10× méně než A a B. Písař 3,
který píše oba jazyky, mezi svými A a B stránkami ukazuje přebytek 0,045–0,048, tedy asi polovinu celého
posunu A/B. Dolní mez intervalu je přitom daleko od nuly.

Šablona bez pravidla dává pro všechny kontrasty ≤ 0,008. S pravidlem „sufixy podle jazyka" dává A/B 0,068 / 0,042
a B2 vs B3 0,001, tedy vzorec reálu, ale jen proto, že je do ní vložený.

## Kontrola podezřelých hodnot

* `cmi_excess` > 1: žádný běh. `stem_dispersion` > 10: žádný běh.
* Jaccard = 1: `jaccard_suffix_AB` je 1 v 3 190 bězích šablony a `jaccard_suffix_sections` v 551. Je to
  strukturální vlastnost stejně jako v 002 (uzavřený inventář 15 sufixů; šablona ho dědí z reálu).
* Testy procházejí (`python -m pytest -q`, 199 prošlo, 1 přeskočen). Kód běží na Pythonu 3.13 s pandas 2.3.3.
  Nepoužívá scipy ani syntaxi 3.11+.

## Interpretace

Systém bez významu, který má paradigmata, tematickou volbu kmenů a stránkovou lokalitu **z konstrukce**,
reprodukuje hlavní pozitivní nálezy 002: podíl kmenů s více sufixy, rozptýlené (ne kopírované) výskyty kmenů
a tematické rozdělení kmenů po sekcích. Stačí na to tři ruční pravidla. Tyto nálezy jsou proto doklad
*systému*, ne *smyslu*, a směr 1 je v této podobě vyčerpaný.

Robustnost 002 vůči afixům a parseru je dobrá. Posun sufixů A/B je jazykový, ne písařský.

Co šablona ani s pěti pravidly nedává:
* nerovnoměrné rozdělení sufixů uvnitř kmene (H 1,8 vs 2,7–2,8 bitu);
* současně správné nepodmíněné i podmíněné MI prefix–sufix;
* podíl hapaxů.

To jsou konkrétní a testovatelné rozdíly. Umělý systém by je ale mohl získat dalšími parametry, takže nejsou
důkazem smyslu.

## Co výsledek neukazuje

* Neukazuje, že text je bez významu. Ukazuje jen, že metriky 002 smysl od dostatečně bohatého systému neodliší.
* Šablona je jedna konkrétní konstrukce (Zipf, pevné sady, globální váhy sufixů). Jiný generátor glosolálie
  (např. lidský, Gaskell & Bowern 2022) může vypadat jinak. Počet pravidel je horní odhad pro *tuto* rodinu.
* Kalibrace ZL vyžadovala `zipf_s` 1,4 mimo mřížku spec. Téma a lokalita byly naladěné na metriky, které pak
  šablona „reprodukuje" (`jsd_stem_sections`, `stem_dispersion`). U nich jde o důkaz *existence* parametrů,
  ne o nezávislou predikci.
* U `cmi_excess` je |z| < 3 jen díky velkému rozptylu šablony mezi seedy (různé inventáře kmenů). Průměr
  šablony je přitom ⅓–½ reálu.
* Písař 3 má v A jen 733 tokenů. Permutační linie s tím počítá, ale přebytek 0,045 je o jednom písaři a dvou
  sekcích.
* Robustnost se týká jen politiky `split`.

## Další krok

1. Do šablony přidat nerovnoměrné váhy sufixů uvnitř kmene (Dirichlet na kmen) a vazbu prefixu přes kmen × sufix.
   Pak ověřit, jestli zbývající rozdíly (H, MI vs CMI, hapaxy) zmizí. Pokud ano, směr 1 uzavřít.
2. Rozdíly, které šablona nedává bez explicitního pravidla (posun sufixů A/B uvnitř písaře 3, sufixové frekvence
   po sekcích), testovat v experimentu 007 proti pozici v řádku. Část posunu může být poziční.
3. Pro syntézu: metriky 002 přeřadit z „indicie smyslu" na „vlastnosti systému" a jako nevysvětlené vést jen
   H_suffix_given_stem, dvojici MI/CMI a A/B uvnitř písaře.
