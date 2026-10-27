# Referenční hodnoty z literatury (pro reporty; nestahovat znovu)

## Rozanova & Temerev (2026), arXiv:2608.17096 – tabulka 7 („profil“)

Selekce: ZL3b, odstavcový text, nejisté mezery = oddělovač (jako naše `split`), Eva kompozity
sloučeny (cth, ckh, cph, cfh, ch, sh, iin, in, ee ≈ náš preset `extended`), 32 747 tokenů.

| Zdroj | h2 (bit) | D64 (bit) | pořadí tokenů (excess / H) | hapax (podíl typů) | edge MI (bit) | edit-1 poměr | diferenciál útoku | crossing |
|---|---|---|---|---|---|---|---|---|
| Voynich (ZL) | 2.69 | 1.05 | 0.79 % | 70 % | 0.197 | 1.12 | −0.15 až +0.08 | 2.5 % |
| Naibbe šifra latiny (Caesar) | 2.71 | 1.03 | 1.25 % | 41 % | 0.006 | 0.82 | +0.00 | 2.2 % |
| Naibbe šifra latiny (Plinius) | 2.71 | 1.04 | 1.06 % | 42 % | – | – | −0.01 | 1.9 % |
| Self-citation generátor (3 seedy) | 2.6–2.7 | 0.65–0.86 | 0.5–0.6 % | 59–60 % | 0.006–0.009 | 1.06–1.23 | −0.02 až +0.02 | – |
| Latinská próza | 3.44 | 1.24 | 4.00 % | 59 % | 0.033 | 0.43 | (index 0.99) | 15.5 % |
| Syntetické homofonní šifry | – | 0.8–1.1 | – | – | – | – | +0.13 až +0.52 | 3.8–12.4 % |

Definice: h2 = H(X_{t+1}|X_t) glyfů bez mezer; D64 = H(U) − H(U_{t+1}|U_t) po 64 BPE sloučeních
uvnitř tokenů; pořadí tokenů = shuffle-korigovaná MI sousedních tokenů (cap 2 000 typů) / H(T);
edge MI = shuffle-korigovaná MI posledního glyfu tokenu t a prvního glyfu tokenu t+1; edit-1 =
pozorovaný/zamíchaný podíl sousedních tokenů ve vzdálenosti ≤ 1 editace; crossing = podíl
skrytých hranic tokenů, které po smazání mezer protínají naučené BPE jednotky.

Další hodnoty z téhož článku: podmíněná entropie kontrol (Plinius 3.51, Celsus+Caesar 3.47,
Culpeper 3.52, italština 3.45; H(X) ≈ 4.0 všude); vzájemná informace sousedních tokenů jako
podíl H(T): botanická latina 2.02 %, ostatní prózy 3.59–9.72 %, Linnéův katalog 15.81 %; index
„vnitřnosti“ hranic: nejistá mezera 0.494, jistá 0.029, zlom řádku −0.036; AUC rozpoznání
nejisté mezery ze šířky mezery v obrázku 0.905; JSD prvního glyfu řádku vs. uprostřed:
0.528 (první řádky odstavců) vs 0.179 (ostatní).

## voytools (VERIFICATION.md, ZL3b, `split`, čisté P tokeny, seed 0)

h2 (extended, scope all) 2.689; adjacent-token MI excess 0.0755 bit = 0.90 % H(X);
edge MI excess 0.185 bit (observed 0.203); near-neighbour share 0.173 vs 0.146 (ratio 1.18);
hapax 0.676 (P), 0.686 (všechny lokusy); typy 6 796 (P).

## Další často citované hodnoty

- Zandbergen (voynich.nu): ZL ≈ 36 278 „dlouhých slov“ – odpovídá politice `join`.
- Lindemann & Bowern (2020): h2 Voynichu nižší než u všech srovnávaných jazyků; hodnota silně
  závisí na transliterační abecedě.
- Gaskell & Bowern (2022): lidský blábol reprodukuje Zipfa, délky slov i poziční preference;
  má nižší informační obsah a více opakování než jazyk.
- Timm & Schinner (2020): sousední tokeny ve vzdálenosti ≤ 1 editace častěji než náhodně;
  samotný rukopis: poměr ≈ 1.12–1.18 podle definice okna.
