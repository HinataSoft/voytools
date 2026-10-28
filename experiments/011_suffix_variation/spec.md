# 011 – Rozklad vnitroslovní variace: je „paradigma“ poziční jev?

## Otázka
Voynich má 5,25 sufixu na kmen s vyrovnaným rozdělením (H(sufix | kmen) 1,8 bitu), zatímco v latině,
angličtině i italštině převažuje u kmene jeden tvar (≤ 0,9 bitu) a žádný generátor šířku netrefí (009).
Zároveň sufix nese informaci o pozici v řádku (0,05 bitu), o písaři/jazyku (0,09) a o sekci (0,07) (007).
Otázka: **kolik z variace sufixu při daném kmeni vysvětlí pozice v řádku, registr A/B, písař a sekce –
a co zbude?** Pokud zbude málo, je „paradigma“ poziční/registrový jev. Pokud zbude hodně, je zbytková
volba tvaru první vážný kandidát na obsah a testuje se dál, zda má strukturu (shoda se sousedy, vazba na
stránku, konzistence napříč rukopisem).

## Hypotézy a predikce
- **Poziční/registrový jev:** H(sufix | kmen, pozice, jazyk, písař) ≪ H(sufix | kmen); zbytek blízko
  korekci na konečný vzorek; ts2020 ukáže podobný podíl vysvětlený pozicí (má řádková pravidla).
- **Obsahová volba:** po podmínění zbude většina variace (jako u flexe v latině, kde pozice v lámaném řádku
  nevysvětlí nic) **a** zbytková volba má strukturu: sufix tokenu závisí na sufixu/prefixu souseda při daném
  kmeni (shoda), nebo je pro daný kmen konzistentní uvnitř stránky/odstavce víc než napříč rukopisem.
- **Kontroly:** latina – variace sufixu je gramatická, pozice v lámaném řádku ≈ 0 %, sekce (kniha) málo;
  ts2020 – pozice vysvětlí část, zbytek je náhodná náhrada; šablona – nic nevysvětlí nic (sufixy losovány).

## Data a filtry
ZL3b a GC2a, `P`, `clean`, split (hlavní) a join. Dekompozice: (a) naučené afixy (výchozí), (b) s pozičními
glyfy (008, varianta b). Kmeny s ≥ 20 tokeny (hlavní) a ≥ 5 (kontrola). Kontroly z 009 přes
`voytools/controls.py`: Celsus, Plinius, Caesar (latina), Culpeper (en), Linnaeus; nully ts2020 (30 seedů),
šablona 005 (30), ac-tri 86 (převzít/30).

## Metriky
1. **Rozklad podmíněné entropie** (bity, plug-in, korekce: hodnota po permutaci X uvnitř kmene, 100×):
   H(S | Z) a H(S | Z, X) pro S = sufix, Z = kmen, X ∈ {poziční třída (22 tříd z 007), rel_pos (5 binů),
   místo v řádku (4), line_in_par (3), jazyk A/B, písař, sekce, stránka} a kombinace {pozice + jazyk},
   {pozice + písař}, {pozice + jazyk + sekce}, {vše kromě stránky}, {vše}. Výstup: vysvětlený podíl
   = (H(S|Z) − H(S|Z,X)_korig.) / H(S|Z) s intervalem (100 subsamplů stránek). Totéž pro prefix.
   Implementovat jako `voytools/variation.py: conditional_entropy_decomposition(dfd, target, given, by, ...)`
   s testem (syntetická data se známou odpovědí).
2. **Zbytková struktura** (jen na reálu a nullech, po podmínění pozicí + jazykem + písařem):
   (i) shoda se sousedem: I(S_t ; S_{t−1} | Z_t, pozice) a I(S_t ; P_{t+1} | Z_t, pozice) s korekcí permutací
   sousedů uvnitř poziční třídy; (ii) konzistence: pro kmeny s ≥ 20 tokeny podíl tokenů nesoucích nejčastější
   sufix kmene **uvnitř stránky** vs napříč rukopisem (očekávání z permutace tokenů kmene mezi stránkami);
   (iii) vazba na sekci po pozici+jazyku+písaři: I(S; sekce | Z, pozice, jazyk, písař).
3. **Deskriptivně:** 10 kmenů s nejvíce tokeny – tabulka sufix × poziční třída (podíly), aby bylo vidět,
   jak variace vypadá (např. `-y` na začátku řádku, `-dy` uvnitř, `-m` na konci?).

## Nejistota
Reál: 100 subsamplů 80 % stránek (SD × 2). Nully: přes seedy. Permutační korekce 100×.

## Výstupy (`results/011_suffix_variation/`)
`decomposition.csv` (S/P × X × zdroj × sada: H, H_cond, korekce, vysvětlený podíl, interval),
`residual_structure.csv`, `consistency.csv`, `top_stems.csv`, `decomposition.png`, `runs.csv`/`meta.json`,
`report.md`, `notes.md`.

## Kritéria přijetí
1. Rozklad pro sufix i prefix, ZL i GC, split i join, obě dekompozice; kontroly a nully.
2. Report odpoví: (i) vysvětlený podíl H(sufix | kmen) pozicí, jazykem, písařem, sekcí a vším dohromady (ZL/GC);
   (ii) totéž u latiny a ts2020; (iii) má zbytek strukturu (shoda, konzistence v rámci stránky, sekce)? –
   ano/ne s čísly; (iv) jedna věta: je šířka paradigmatu poziční jev?
3. Podezřelé: vysvětlený podíl < 0 po korekci (chyba), > 0,95 u stránky (přeučení – stránka má moc hodnot;
   uvést korekci).
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
Entropie z tabulek; ~1–2 s na sadu × (13 X × 2 cíle × 2 dekompozice) × (1 + 100 subsamplů) × 4 sady ≈ 40 min;
nully a kontroly ≈ 30 min. Jedna session, jedna větev, PR.
