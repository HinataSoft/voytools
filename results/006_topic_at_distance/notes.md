# 006 – deník session

## Rozhodnutí a výklad spec

- **Větev** `exp/006-topic-at-distance`, výsledky `results/006_topic_at_distance/`. Kód: `run.py`;
  obecné pomocné funkce v `voytools/stats.py` (`page_term_counts`, `page_similarity`, `distance_bins`,
  `ols_fit`) s testy v `tests/test_stats.py`. `voytools/nulls.py` beze změny.
- **Stránky:** P text `clean`, stránky s ≥ 40 tokeny (201 z 207 stránek s P textem v ZL i GC; 20 100 dvojic).
  Pořadí = pořadí prvního výskytu folia v souboru přes **všechny** lokusy (227 stránek), takže stránky bez
  P textu (zodiak) se do vzdálenosti počítají.
- **Kmeny:** `decompose_table` s afixy naučenými na reálném P textu dané kombinace (soubor × politika);
  všechny vzorky nullů dekomponovány stejnými seznamy.
- **Funkční kmeny:** K nejčastějších kmenů *daného vzorku* (u permutací = reál, u auto-copy vlastní seznam
  generovaného textu). Seznam 30 reálných: `function_stems.csv` (100 nejčastějších s počty).
- **Sekce** = kód ilustrace (`illustration`), písař = `hand` (Davis), jazyk = `language`.
- **Permutační nully** permutují *obsah* stránek mezi pozicemi; pozice si nechává metadata (pořadí, sekci,
  jazyk, písaře). Délka stránky se stěhuje s obsahem.
- **Bootstrap po stránkách:** 201 stránek s opakováním, dvojice všech vybraných pozic kromě dvojic téže
  původní stránky.
- **Auto-copy:** config 86 z 003 (`trigram`, `seed`, p_new 0,3, window 32, p_copy_exact 0,35, p_two_edits 0,5,
  recency 0, edit_type_probs 0,6/0,2/0,2), `autocopy_like` s rozvržením reálného P textu (stejné stránky,
  řádky, délky), seedy 10 000–10 049. `voytools.selfcitation` (004) neexistuje → druhý generátor vynechán.
- **Kontrolní politika:** `join` (vedle hlavní `split`), jak říká spec; `strict` se nepočítal.

## Odchylky a doplňky mimo spec

1. **Velikost stránky jako kovariát.** Kosinus log(1+f) i Jaccard rostou s počtem tokenů stránky a velikost
   stránek se mezi sekcemi/písaři silně liší (krátké herbářové stránky A × dlouhé B). Ve spec modelu proto
   vycházel koeficient `same_section` i `same_hand` *záporný* i u auto-copy, které o sekci nic neví.
   Přidán model `spec+size` (+ log n_i + log n_j a |log n_i − log n_j|) a kontrast reziduí po odečtení
   velikosti (`contrast_adj`). Spec model je v `ols.csv` beze změny (`model == 'spec'`).
2. **Kontrola písaře:** kromě kovariátu `same_hand` a nullu „permutace uvnitř písaře“ i model jen na
   dvojicích téhož písaře (`model == 'same_hand'`: log_dist + same_section) a kontrast jen v nich.
3. **Stránky `T`/`C`/`A`.** Kód `T` (jen text) mají stránky uvnitř jiných sekcí (f1r, f66r, f76r, f85r1 …),
   takže tvoří většinu dvojic „jiná sekce, vzdálenost ≤ 3“. Citlivost: `sensitivity_sections.csv`
   (vyřazeno nic / `T` / `T`, `C`, `A`), reál + bootstrap + všechny nully.
4. **Popisky i proti textu P+C+R** (`text == 'PCR'`), protože zodiakové stránky nemají P text a popisky `Lz`
   tam nemohou mít „hit“ (podíl 0 v reálu i nullu).
5. **v101 třídy počátečních glyfů** odvozeny zarovnáním lokusů ZL/GC se stejným počtem tokenů
   (`initial_map.csv`): Eva ok- = v101 oh-/Ah-, ot- = ok-/Ak-, qo- = 4o-, ch-/sh- = 1-/2-/3-/5-, d- = 8-/7-.
6. **ok-/ot- profil:** null = posun uvnitř řádku (`shuffle_within_lines`, 20×), intervaly = bootstrap po
   stránkách (200×). Glyfové třídy sousedů = stejné počáteční třídy + okraj řádku.

## Slepé uličky

- `pkill -f "run.py sections"` zabil vlastní shell (vzor se shodoval s příkazem). Běh zopakován.
