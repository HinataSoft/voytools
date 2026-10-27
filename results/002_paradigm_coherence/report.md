# 002 – Paradigma, nebo šablona? Koherence kmenů a afixů

Kód: `experiments/002_paradigm_coherence/run.py`. Rozhodnutí a výklad spec: `experiments/002_paradigm_coherence/notes.md`.
Všechna čísla níže jsou ze souborů v této složce (`effects.csv`, `summary.csv`, `inventories.csv`,
`summary_info.json`, `runs.csv`; po částech `parts/<kombinace>/<zdroj>/`).

## Co bylo spočítáno

* **Data:** ZL3b (Eva) a GC2a (v101), `generic_type="P"`, `clean=True`; politiky `split`, `join`,
  `strict` (u `strict` vyřazeny tokeny s `,` uvnitř: ZL 2 104, GC 1 921). Tokenů: ZL 34 268 / 31 808 /
  29 704, GC 36 385 / 34 143 / 32 222.
* **Kombinace:** 9 (spec píše 8, ale 2 × 3 + 3 = 9): ZL × 3 politiky × {naučené, referenční} afixy,
  GC × 3 politiky × naučené. Všech 9 má kompletní tabulku (`summary_info.json` → `complete`).
* **Zdroje:** reál (bodový odhad) + bootstrap po stránkách 200×; slotová nezávislost 50×;
  autokopírování podle **obou** prvních řádků `results/001_autocopy_calibration/best_configs.csv`
  (`autocopy_573` a `autocopy_536`) 50× každý; míchání v řádcích 5×; stem-conditional shuffle =
  `cmi_shuffled` uvnitř každého běhu. Celkem 3 204 běhů.
* **Intervaly:** v `effects.csv` je rozdíl = reál (bodový odhad) − průměr nulového modelu; 95% interval
  = rozdíl ± 1,96·√(SD_boot² + SD_null²); `z_null` = rozdíl / SD_null. Níže „[a; b]“ = tento interval.

**Důležité pro robustnost:** naučené afixové seznamy ZL jsou ve všech třech politikách **množinově
totožné** s referenčními (`parts/ZL_*_learned/real/meta.json`). Kombinace `*_reference` proto dávají
přesně stejná čísla jako `*_learned` a závislost na volbě afixů tímto experimentem **testována
nebyla**. Robustnost tu ověřuje jen přechod na GC (jiná abeceda, jiné naučené seznamy).

## Čísla (ZL split / GC split; ostatní politiky v `effects.csv`, ve všech stejný směr)

| metrika | reál ZL | slot. nezáv. | autocopy_573 | autocopy_536 | reál GC | slot GC | ac573 GC | ac536 GC |
|---|---|---|---|---|---|---|---|---|
| cmi_excess (bit) | 0,088 | 0,035 | 0,226 | 0,196 | 0,089 | 0,041 | 0,239 | 0,194 |
| frac_stems_multi_suffix | 0,912 | 1,000 | 0,508 | 0,616 | 0,839 | 1,000 | 0,570 | 0,699 |
| frac_multi_suffix_5to9 | 0,832 | 0,999 | 0,349 | 0,470 | 0,777 | 0,999 | 0,410 | 0,552 |
| mean_suffixes_per_stem | 5,25 | 7,81 | 2,28 | 2,69 | 4,72 | 7,44 | 2,59 | 3,16 |
| H_suffix_given_stem (bit) | 1,82 | 3,13 | 1,36 | 1,53 | 1,80 | 3,06 | 1,54 | 1,70 |
| mi_prefix_suffix (bit) | 0,205 | 0,006 | 0,088 | 0,060 | 0,176 | 0,006 | 0,078 | 0,048 |
| n_stems_k | 308 | 311 | 947 | 851 | 417 | 403 | 933 | 828 |
| productivity_slope | 0,309 | 0,296 | 0,480 | 0,492 | 0,328 | 0,312 | 0,514 | 0,524 |
| jaccard_suffix_sections | 0,981 | 0,984 | 0,852 | 0,861 | 0,962 | 0,994 | 0,881 | 0,891 |
| jaccard_stem_sections | 0,325 | 0,355 | 0,115 | 0,126 | 0,268 | 0,323 | 0,117 | 0,129 |
| jsd_suffix_sections | 0,046 | 0,005 | 0,024 | 0,022 | 0,037 | 0,005 | 0,025 | 0,019 |
| jsd_stem_sections | 0,227 | 0,128 | 0,493 | 0,486 | 0,271 | 0,158 | 0,466 | 0,465 |
| jaccard_suffix_AB / jaccard_stem_AB | 1,00 / 0,405 | 1,00 / 0,471 | 0,99 / 0,173 | 0,99 / 0,204 | 1,00 / 0,369 | 1,00 / 0,437 | 1,00 / 0,167 | 1,00 / 0,200 |
| jsd_suffix_AB | 0,088 | 0,0003 | 0,0025 | 0,0022 | 0,065 | 0,0003 | 0,0025 | 0,0017 |
| stem_dispersion | 0,459 | 0,687 | 0,005 | 0,039 | 0,481 | 0,687 | 0,004 | 0,053 |

Stem-conditional shuffle: I(P;S|kmen) reálu 0,217 vs 0,129 po permutaci prefixů uvnitř kmene (ZL split),
rozdíl 0,088 [0,074; 0,101]; GC split 0,227 vs 0,138, rozdíl 0,089 [0,075; 0,103]; v ostatních
kombinacích 0,076–0,078 se stejně širokými intervaly.

Míchání v řádcích: všechny typové metriky jsou **přesně** stejné jako u reálu (max. |rozdíl| = 0 pro
n_stems_k, mean_suffixes_per_stem, frac_multi, H, MI, CMI, Jaccardy, productivity_slope;
`summary_info.json`). `cmi_excess` se liší jen o šum permutací (≤ 0,001). `stem_dispersion` se mění,
ale jen málo (max. |rozdíl| 0,007), protože míchání přesouvá tokeny jen uvnitř řádku.

## Odpovědi na tři otázky spec

**(i) Je `cmi_excess` reálu nad slotovou nezávislostí i nad autokopírováním?**
Nad slotovou nezávislostí **ano** – ZL ano, GC ano: rozdíl ZL split +0,053 [0,038; 0,068],
GC split +0,048 [0,033; 0,063], ve všech 9 kombinacích +0,040 až +0,053 s dolní mezí > 0,02.
Slotový null nedává přesně 0, ale 0,025–0,041: greedy re-parse slepeného slova vytváří trochu
závislosti sám. Nad autokopírováním **ne** – ZL ne, GC ne: oba generátory dávají `cmi_excess`
0,18–0,24, tedy 2–3× víc než reál (ZL split −0,138 [−0,159; −0,118] vs 573, −0,108 [−0,127; −0,089]
vs 536; GC obdobně). Podobně `mi_prefix_suffix`: reál 0,18–0,21 bitu, generátory 0,04–0,09, slot 0,006.
Kmeny reálu přitom nesou pestřejší paradigma než u generátoru: 5,3 vs 2,3–2,7 sufixu na kmen, podíl
kmenů s více sufixy 0,91 vs 0,51–0,62. I u kmenů s 5–9 tokeny je to 0,83 vs 0,35–0,47. Zároveň je
paradigma reálu užší než u nezávislé šablony (7,8 sufixu, podíl 1,00, H 3,1 vs 1,8 bitu).
Sklon produktivity reálu (0,31–0,34) se od slotového nullu neliší (ZL split +0,013 [−0,021; 0,048]).
U generátorů je výrazně vyšší (0,48–0,52), ale ty mají 2–3× více kmenů s ≥ 5 tokeny (850–950 vs 308).
Predikce „pomalejší růst než u šablony“ se tedy **nepotvrdila**. Reál se od šablony liší nižší
hladinou (průsečík 0,51 vs 1,06), ne sklonem.

**(ii) Jsou afixy mezi sekcemi stálejší než kmeny?**
Podle Jaccardu inventářů **ano** – ZL ano, GC ano: sufixy 0,96–0,98, kmeny 0,27–0,33. Totéž ale
platí pro **všechny** nulové modely (gap reálu 0,66–0,69, slot 0,63–0,67, autocopy 0,72–0,76;
intervaly rozdílů obsahují 0). Pro 15 sufixů (14 + prázdný) je Jaccard ≈ 1 strukturálně
(uzavřený seznam). Spec takovou hodnotu označuje jako podezřelou, viz níže. Informativnější jsou
frekvenční profily (JSD):
* Rozdělení kmenů se mezi sekcemi liší víc než u slotového nullu (0,227 vs 0,128, ZL split
  +0,099 [0,058; 0,140]; GC +0,113 [0,073; 0,153]), ale méně než u autokopírování (0,49).
  To odpovídá tematické závislosti kmenů bez lokálního driftu.
* Frekvence sufixů se mezi sekcemi mění také: JSD 0,046 vs 0,005 slot a 0,022–0,024 autocopy.
  Nejvíc se mění mezi jazyky A a B: `jsd_suffix_AB` 0,088 (ZL) / 0,065 (GC) vs ≤ 0,003 u
  všech nullů. Afixový **inventář** je stálý, afixové **frekvence** stálé nejsou. Posun mezi A a B
  je velký a žádný null ho nereprodukuje. Je ale spojený s písaři (Currier A/B ~ hand), takže ho
  nelze číst jako změnu tématu.

**(iii) Jsou varianty kmene rozprostřené (jazyk), nebo shluknuté (kopírování)?**
**Rozprostřené** vzhledem ke kopírování – ZL ano, GC ano: `stem_dispersion` reálu je 0,46–0,48
ve všech 9 kombinacích. Autocopy_573 dává 0,004–0,005, autocopy_536 0,03–0,06 (ZL split
+0,454 [0,415; 0,493] a +0,419 [0,373; 0,465]). Proti náhodnému rozmístění (slot 0,687 ≈ ln 2)
je reál mírně shluknutý (−0,23 [−0,27; −0,19]). To odpovídá kmenům vázaným na stránku nebo téma.
Neodpovídá to kopírování na vzdálenost desítek tokenů.

## Kontrola podezřelých hodnot (kritérium 3)

* `cmi_excess` > 1 bit: žádný běh. `stem_dispersion` > 10: žádný běh. Jaccard = 0: žádný.
* Jaccard = 1: `jaccard_suffix_AB` u reálu ve všech 9 kombinacích a u většiny bootstrapů i nullů.
  `jaccard_suffix_sections` = 1 ve 32 bootstrapech a 203 slotových vzorcích. Příčina je
  strukturální: inventář je 15 hodnot a při ≥ 3 výskytech se v A i B objeví všechny. Proto jsou
  přidány JSD. U stemového Jaccardu dominuje velikost sekce: dvojice s malou astronomickou sekcí A
  (271 tokenů, 27 kmenů ≥ 3) mají 0,06–0,12, ostatní 0,3–0,5 (`inventories.csv`).
* Bootstrap po stránkách je u typových metrik vychýlený (duplicitní stránky). Například
  `cmi_excess` bootstrapu ZL split má rozsah 2,5–97,5 % [0,133; 0,160] proti bodovému odhadu 0,088,
  `n_stems_k` [314; 394] proti 308 (`summary.csv`). Proto se pracuje s bodovým odhadem a jen se SD
  bootstrapu (jako v 001). I graf `paradigm.png` ukazuje reál jako bod ± 1,96 SD. Jako míra
  nejistoty je to hrubé.

## Interpretace

Kmeny se chovají víc jako lexikální jednotky s paradigmatem než jako výplň nezávislé šablony.
Kombinují se s omezenou sadou sufixů a v rámci kmene platí vazba prefix–sufix nad slotovou
nezávislostí. Varianty kmene jsou rozptýlené po celém textu, ne v lokálních shlucích, což kalibrované
autokopírování nereprodukuje. Platí to v obou abecedách a všech třech politikách mezer. Samotný
`cmi_excess` ale jazyk od generátoru neodlišuje (generátor ho má víc). Rozlišuje až kombinace
„hodně sufixů na kmen + rozptýlené výskyty kmene + silná nepodmíněná vazba prefix–sufix“.

## Co výsledek neukazuje

* Neukazuje, že afixy jsou morfémy ani že text nese význam. Pravidelný glosolalický nebo ručně
  vytvářený systém s omezenými kombinacemi a tematickou volbou kmenů by vzorec mohl reprodukovat.
  Takový model (Gaskell & Bowern 2022) tu nebyl.
* Autokopírování je jen hrubá kalibrace z 001 (loss ≈ 1, H2 a n_types netrefené). Má 2–3× víc kmenů
  s ≥ 5 výskyty, takže část rozdílů (cmi_excess, slope) může být efektem nekalibrované velikosti
  slovníku, ne mechanismu kopírování.
* Závislost na volbě afixů netestována (naučené = referenční seznamy u ZL). Slotový null závisí na
  greedy parseru (re-parse dává `cmi_excess` ≈ 0,03, ne 0).
* Rozdíly sufixových frekvencí A/B mohou být písařské, ne jazykové (A/B a písař jsou zaměněné).

## Další krok

1. Robustnost afixů: zopakovat s jinými seznamy (např. `n_suffixes` 8/20, `max_len` 3, BPE jednotky)
   a s parserem s preferencí sufixu.
2. Autocopy null s velikostí slovníku srovnanou s reálem (diag bod z 001, p_new 0,5) a
   `stem_dispersion` jako cílová metrika kalibrace.
3. Rozdělit `jsd_suffix_AB` podle písaře (`hand`) uvnitř jazyka a podle sekce uvnitř písaře.
4. Přidat null pro „systematickou glosolálii“ (gramatika s omezenými kombinacemi a tematickou volbou
   kmenů) a ověřit, jestli trojice metrik (i)–(iii) reál od něj odliší.
