# Verification of voytools 0.1.0

Numbers measured with `tools/make_verification.py` on the real transliteration files (Python 3.13.17, numpy 2.5.3, pandas 3.0.5).  Unless stated otherwise the statistics use paragraph text (generic type `P`), clean tokens (no unreadable `?`, no alternative readings, no high-ASCII characters), the default `split` space policy, seed 0 and 100 shuffles for the shuffle baselines.

## 1. Test suite

`pytest 9.1.1: 161 passed in 29 s with all real files present; 145 passed and 16 skipped in 5 s when the real files are not found`

## 2. Corpus size

| file | alphabet | IVTFF | pages | loci | tokens (all) | all, clean | P | P, clean | P0, clean | parser warnings |
|---|---|---|---|---|---|---|---|---|---|---|
| ZL3b-n.txt | Eva- | 2.0 | 227 | 5 385 | 39 026 | 37 978 | 35 075 | 34 268 | 33 205 | 0 |
| GC2a-n.txt | v101 | 2.0 | 226 | 5 367 | 40 701 | 40 283 | 36 703 | 36 385 | 35 273 | 0 |
| RF1b-er.txt | Eva- | 2.0 | 227 | 5 385 | 38 525 | 38 437 | 34 659 | 34 618 | 33 535 | 0 |
| IT2a-n.txt | EvaT | 2.0 | 225 | 5 215 | 37 919 | 37 759 | 34 491 | 34 411 | 33 383 | 0 |
| LSI_ivtff_0d.txt | Eva- | 1.5 | 225 | 5 216 | 37 919 | 37 759 | 34 491 | 34 411 | 33 383 | 0 |

The LSI file is interlinear (18 transcribers); the row shows transcriber H, the most complete one.

## 3. Vocabulary (paragraph text, clean tokens)

| file | tokens | types | hapax | hapax share | mean length (chars) | mean length (glyphs) |
|---|---|---|---|---|---|---|
| ZL3b-n.txt | 34 268 | 6 796 | 4 594 | 0.6760 | 4.9883 | 4.4597 |
| GC2a-n.txt | 36 385 | 8 351 | 5 859 | 0.7016 | 3.8736 | 3.8736 |
| RF1b-er.txt | 34 618 | 7 348 | 5 081 | 0.6915 | 5.0543 | 4.5140 |
| IT2a-n.txt | 34 411 | 7 069 | 4 791 | 0.6777 | 5.0659 | 4.5275 |
| LSI_ivtff_0d.txt | 34 411 | 7 069 | 4 791 | 0.6777 | 5.0659 | 4.5275 |

Hapax share = hapax legomena / types.  Glyph lengths use the `basic` Eva preset (identity for v101).

## 4. Glyph entropy of the space-less text (bits)

`H(X)` = unigram entropy `H1`, `H(X_{t+1} given X_t)` = conditional entropy `H2`, `H3` = the same given two preceding glyphs.  Scope `line`: glyph pairs do not cross line breaks; scope `all`: the whole text is one stream (this is what Rozanova & Temerev do).

### ZL3b-n.txt

| preset | scope | n_glyphs | n_symbols | H1 | H2 | H3 |
|---|---|---|---|---|---|---|
| basic | line | 152 826 | 32 | 3.8844 | 2.4452 | 2.2254 |
| basic | all | 152 826 | 32 | 3.8844 | 2.4925 | 2.2811 |
| extended | line | 139 009 | 35 | 3.9709 | 2.6435 | 2.4426 |
| extended | all | 139 009 | 35 | 3.9709 | 2.6894 | 2.4916 |

### GC2a-n.txt

| preset | scope | n_glyphs | n_symbols | H1 | H2 | H3 |
|---|---|---|---|---|---|---|
| none | line | 140 940 | 70 | 4.1533 | 2.8517 | 2.6156 |
| none | all | 140 940 | 70 | 4.1533 | 2.8973 | 2.6631 |

### Same measure on the other transcriptions (extended preset, scope `all`)

| file | H1 | H2 | H3 |
|---|---|---|---|
| ZL3b-n.txt | 3.9709 | 2.6894 | 2.4916 |
| RF1b-er.txt | 3.9738 | 2.7120 | 2.5167 |
| IT2a-n.txt | 3.9680 | 2.7071 | 2.5115 |
| LSI_ivtff_0d.txt | 3.9680 | 2.7071 | 2.5115 |

## 5. Mutual information

### Adjacent tokens (within lines; types capped at 2000 + `<other>`)

Excess = observed MI minus the mean of 100 within-line shuffles.  The *excess share* used in the comparison below is excess / H(X), where H(X) is the entropy of the (capped) token distribution, i.e. the fraction of the information of a token that is explained by its left neighbour; the column `excess / MI obs` is the share of the observed plug-in MI that survives the shuffle correction.

| file | token pairs | H(X) capped | MI obs | MI shuffled | excess (bits) | excess / H(X) | excess / MI obs | z |
|---|---|---|---|---|---|---|---|---|
| ZL3b-n.txt | 30 140 | 8.3896 | 3.2058 | 3.1304 | 0.0755 | 0.0090 | 0.0235 | 11.2299 |
| GC2a-n.txt | 32 257 | 8.1464 | 3.0326 | 2.9280 | 0.1047 | 0.0128 | 0.0345 | 17.0633 |
| RF1b-er.txt | 30 488 | 8.2972 | 3.1428 | 3.0682 | 0.0746 | 0.0090 | 0.0238 | 9.5825 |
| IT2a-n.txt | 30 294 | 8.3720 | 3.2012 | 3.1423 | 0.0589 | 0.0070 | 0.0184 | 8.1425 |
| LSI_ivtff_0d.txt | 30 294 | 8.3720 | 3.2012 | 3.1422 | 0.0590 | 0.0070 | 0.0184 | 8.0237 |

### Edge glyphs (last glyph of token t against first glyph of token t+1)

| file | pairs | MI obs | MI shuffled | excess (bits) | z |
|---|---|---|---|---|---|
| ZL3b-n.txt | 30 140 | 0.2027 | 0.0181 | 0.1846 | 209.4168 |
| GC2a-n.txt | 32 257 | 0.2661 | 0.0408 | 0.2253 | 175.1349 |
| RF1b-er.txt | 30 488 | 0.1969 | 0.0176 | 0.1793 | 177.7105 |
| IT2a-n.txt | 30 294 | 0.1896 | 0.0166 | 0.1730 | 160.3422 |
| LSI_ivtff_0d.txt | 30 294 | 0.1896 | 0.0166 | 0.1730 | 162.0021 |

### Neighbour similarity (window 5 tokens, edit distance <= 1, shuffle = within-page)

| file | share near | shuffled | ratio | z | share identical | identical shuffled |
|---|---|---|---|---|---|---|
| ZL3b-n.txt | 0.1726 | 0.1457 | 1.1843 | 17.3377 | 0.0443 | 0.0344 |
| GC2a-n.txt | 0.2023 | 0.1768 | 1.1445 | 16.6640 | 0.0352 | 0.0299 |

## 6. Positional profile of the first glyph (ZL3b)

Share of tokens in each position class that begin with the glyph (the eight glyphs with the highest share in any class; `parfirst_*` = lines that begin a paragraph).

| glyph | line_first | mid | line_last | parfirst_line_first | parfirst_mid | parfirst_line_last |
|---|---|---|---|---|---|---|
| p | 0.013 | 0.002 | 0.001 | 0.466 | 0.019 | 0.014 |
| o | 0.155 | 0.209 | 0.230 | 0.042 | 0.225 | 0.249 |
| t | 0.066 | 0.019 | 0.017 | 0.242 | 0.020 | 0.014 |
| ch | 0.039 | 0.201 | 0.146 | 0.001 | 0.123 | 0.090 |
| d | 0.187 | 0.076 | 0.170 | 0.017 | 0.063 | 0.156 |
| q | 0.149 | 0.163 | 0.083 | 0.037 | 0.181 | 0.096 |
| y | 0.174 | 0.036 | 0.039 | 0.017 | 0.036 | 0.040 |
| sh | 0.057 | 0.093 | 0.037 | 0.006 | 0.136 | 0.063 |
| tokens in class | 3 266 | 21 045 | 3 219 | 698 | 5 297 | 719 |

Jensen-Shannon divergence between classes of the first glyph (bits):

| class | line_first | mid | line_last | parfirst_line_first | parfirst_mid | parfirst_line_last |
|---|---|---|---|---|---|---|
| line_first | 0.0000 | 0.1860 | 0.1730 | 0.5320 | 0.1810 | 0.1410 |
| mid | 0.1860 | 0.0000 | 0.0450 | 0.6200 | 0.0260 | 0.0590 |
| line_last | 0.1730 | 0.0450 | 0.0000 | 0.6360 | 0.0800 | 0.0270 |
| parfirst_line_first | 0.5320 | 0.6200 | 0.6360 | 0.0000 | 0.5520 | 0.5760 |
| parfirst_mid | 0.1810 | 0.0260 | 0.0800 | 0.5520 | 0.0000 | 0.0510 |
| parfirst_line_last | 0.1410 | 0.0590 | 0.0270 | 0.5760 | 0.0510 | 0.0000 |

## 7. Prefixes and suffixes (top 15 by token count, up to 3 glyphs, at least 20 tokens)

### ZL3b-n.txt

Prefixes

| affix | n_glyphs | count | share | n_types | n_free | ext_share |
|---|---|---|---|---|---|---|
| o | 1 | 7 050 | 0.2057 | 1 299 | 597 | 0.2845 |
| ch | 1 | 5 540 | 0.1617 | 942 | 429 | 0.4653 |
| q | 1 | 5 232 | 0.1527 | 783 | 431 | 0.9778 |
| qo | 2 | 5 116 | 0.1493 | 689 | 407 | 0.5928 |
| d | 1 | 3 219 | 0.0939 | 506 | 319 | 0.6480 |
| qok | 3 | 3 033 | 0.0885 | 227 | 129 | 0.4517 |
| sh | 1 | 3 019 | 0.0881 | 495 | 234 | 0.5816 |
| che | 2 | 2 578 | 0.0752 | 353 | 240 | 0.2564 |
| da | 2 | 2 086 | 0.0609 | 186 | 75 | 0.5978 |
| ok | 2 | 2 006 | 0.0585 | 277 | 139 | 0.3903 |
| ot | 2 | 1 831 | 0.0534 | 295 | 170 | 0.3506 |
| a | 1 | 1 773 | 0.0517 | 243 | 113 | 0.4405 |
| she | 2 | 1 756 | 0.0512 | 214 | 145 | 0.2847 |
| y | 1 | 1 701 | 0.0496 | 494 | 317 | 0.3022 |
| cho | 2 | 1 529 | 0.0446 | 281 | 207 | 0.3179 |

Suffixes

| affix | n_glyphs | count | share | n_types | n_free | ext_share |
|---|---|---|---|---|---|---|
| y | 1 | 14 034 | 0.4095 | 2 648 | 499 | 0.4439 |
| dy | 2 | 6 229 | 0.1818 | 971 | 293 | 0.6410 |
| n | 1 | 5 673 | 0.1655 | 843 | 23 | 0.9804 |
| in | 2 | 5 562 | 0.1623 | 788 | 28 | 0.7073 |
| l | 1 | 5 366 | 0.1566 | 806 | 209 | 0.6034 |
| r | 1 | 5 059 | 0.1476 | 960 | 161 | 0.4649 |
| edy | 3 | 3 993 | 0.1165 | 379 | 59 | 0.2990 |
| iin | 3 | 3 934 | 0.1148 | 525 | 62 | 0.9174 |
| ey | 2 | 3 639 | 0.1062 | 582 | 74 | 0.4705 |
| ol | 2 | 3 238 | 0.0945 | 390 | 152 | 0.2329 |
| ar | 2 | 2 352 | 0.0686 | 432 | 198 | 0.2594 |
| or | 2 | 2 001 | 0.0584 | 341 | 146 | 0.2144 |
| al | 2 | 1 795 | 0.0524 | 328 | 162 | 0.2479 |
| eey | 3 | 1 712 | 0.0500 | 241 | 71 | 0.4340 |
| ain | 3 | 1 605 | 0.0468 | 253 | 130 | 0.3925 |

### GC2a-n.txt

Prefixes

| affix | n_glyphs | count | share | n_types | n_free | ext_share |
|---|---|---|---|---|---|---|
| o | 1 | 7 107 | 0.1953 | 1 523 | 681 | 0.2742 |
| 1 | 1 | 5 835 | 0.1604 | 1 080 | 474 | 0.4255 |
| 4 | 1 | 5 279 | 0.1451 | 949 | 496 | 0.9758 |
| 4o | 2 | 5 151 | 0.1416 | 845 | 461 | 0.5766 |
| 8 | 1 | 3 018 | 0.0829 | 513 | 327 | 0.6034 |
| 4oh | 3 | 2 970 | 0.0816 | 286 | 168 | 0.3162 |
| 1c | 2 | 2 483 | 0.0682 | 372 | 262 | 0.2573 |
| a | 1 | 2 323 | 0.0638 | 231 | 114 | 0.2970 |
| 2 | 1 | 2 317 | 0.0637 | 476 | 234 | 0.4998 |
| oh | 2 | 1 949 | 0.0536 | 329 | 176 | 0.3535 |
| ok | 2 | 1 869 | 0.0514 | 340 | 193 | 0.3392 |
| 8a | 2 | 1 821 | 0.0500 | 160 | 68 | 0.3921 |
| 9 | 1 | 1 696 | 0.0466 | 568 | 342 | 0.2919 |
| 1o | 2 | 1 560 | 0.0429 | 296 | 217 | 0.3090 |
| oe | 2 | 1 537 | 0.0422 | 291 | 230 | 0.2023 |

Suffixes

| affix | n_glyphs | count | share | n_types | n_free | ext_share |
|---|---|---|---|---|---|---|
| 9 | 1 | 13 831 | 0.3801 | 3 197 | 577 | 0.3225 |
| e | 1 | 5 852 | 0.1608 | 889 | 233 | 0.5819 |
| y | 1 | 4 864 | 0.1337 | 963 | 194 | 0.4918 |
| 89 | 2 | 4 460 | 0.1226 | 944 | 291 | 0.4955 |
| m | 1 | 3 796 | 0.1043 | 556 | 70 | 0.9341 |
| am | 2 | 3 546 | 0.0975 | 453 | 202 | 0.3029 |
| oe | 2 | 3 405 | 0.0936 | 425 | 160 | 0.1977 |
| ay | 2 | 2 392 | 0.0657 | 465 | 210 | 0.2115 |
| c9 | 2 | 2 352 | 0.0646 | 461 | 92 | 0.3474 |
| c89 | 3 | 2 210 | 0.0607 | 311 | 72 | 0.3661 |
| oy | 2 | 1 931 | 0.0531 | 352 | 146 | 0.2191 |
| ae | 2 | 1 894 | 0.0521 | 334 | 160 | 0.2339 |
| n | 1 | 1 687 | 0.0464 | 294 | 29 | 0.9769 |
| an | 2 | 1 648 | 0.0453 | 265 | 136 | 0.3823 |
| 79 | 2 | 1 572 | 0.0432 | 406 | 133 | 0.6081 |

Affix lists learned by `decompose_table` from ZL3b paragraph text (ranked by the number of attested stems, dominated affixes dropped):

- prefixes: `che cho ol qo ch d l o s sh t y`
- suffixes: `aiin ain al am ar chy dy ol or d l r s y`

Paradigm metrics with these lists (stems with >= 5 tokens):

| metric | ZL3b |
|---|---|
| n_tokens | 34 268 |
| n_stems | 2 345 |
| n_stems_k | 308 |
| token_coverage_k | 0.9199 |
| mean_suffixes_per_stem | 5.2500 |
| mean_prefixes_per_stem | 5.0714 |
| frac_stems_multi_suffix | 0.9123 |
| H_suffix_given_stem | 1.8170 |
| mi_prefix_suffix | 0.2052 |
| mi_prefix_suffix_shuffled | 0.0034 |
| cmi_prefix_suffix_given_stem | 0.2167 |
| cmi_shuffled | 0.1291 |
| cmi_excess | 0.0876 |

## 8. Comparison with published values

| quantity | published | voytools | verdict |
|---|---|---|---|
| glyph conditional entropy H(X_{t+1} given X_t), composites collapsed, ZL3b | about 2.69 (+-0.1 accepted) | 2.6894 (extended, scope all); 2.4925 (basic, scope all) | agrees |
| adjacent-token MI excess share = excess / H(X) (Rozanova & Temerev 2026) | < 1 % | 0.90 % (excess 0.0755 bits) | agrees |
| edge-glyph MI | about 0.2 bits | 0.2027 observed; 0.1846 excess over shuffles | agrees |

## 9. Command line interface

All commands run in-process through `voytools.cli.main` (the same code as `python -m voytools`); `mi` includes the neighbour-similarity test.

| file | command | exit code | seconds | output bytes |
|---|---|---|---|---|
| ZL3b-n.txt | summary | 0 | 1.88 | 4 327 |
| ZL3b-n.txt | entropy | 0 | 1.31 | 318 |
| ZL3b-n.txt | mi --shuffles 100 | 0 | 3.85 | 982 |
| ZL3b-n.txt | positions --permutations 100 | 0 | 1.53 | 4 920 |
| ZL3b-n.txt | affixes --top 15 | 0 | 0.96 | 2 596 |
| GC2a-n.txt | summary | 0 | 1.38 | 4 187 |
| GC2a-n.txt | entropy | 0 | 1.02 | 228 |
| GC2a-n.txt | mi --shuffles 100 | 0 | 4.15 | 982 |
| GC2a-n.txt | positions --permutations 100 | 0 | 1.37 | 4 878 |
| GC2a-n.txt | affixes --top 15 | 0 | 1.07 | 2 596 |

## 10. Timing

Parsing and building the tables of ZL3b takes 0.8 s; adjacent-token plus edge-glyph MI with 100 shuffles takes 1.7 s (GC2a: 1.7 s) and the neighbour-similarity test with 100 within-page shuffles 1.4 s (GC2a: 1.6 s), which together make up the time of the `mi` command in section 9.

## 11. Additional checks

### Token counts of ZL3b under different conventions

The task brief expected roughly 36-38 thousand tokens in total and about 32.7 thousand in paragraph text.  The counts depend on how the uncertain space `,` and unreadable tokens are treated:

| space policy | tokens | all loci | P | P0 |
|---|---|---|---|---|
| split | unfiltered | 39 026 | 35 075 | 33 976 |
| split | clean | 37 978 | 34 268 | 33 205 |
| join | unfiltered | 36 278 | 32 612 | 31 585 |
| join | clean | 35 234 | 31 808 | 30 817 |

The tokenizer of `voynich-units/analysis/reproduce_headlines.py` (regular expressions; commas and periods split; `<...>` tags and `{...}` ligatures deleted) finds 34,271 paragraph tokens in the same file, 34,175 of them without any of the characters `?*<>{}[]()|`; voytools finds 35,075 and 34,268.  The difference of 804 unfiltered tokens consists of 752 `<->` / `<~>` drawing breaks inside lines (voytools starts a new token there, as IVTFF 2.0 prescribes; the reference merges the two halves), 39 tokens in loci that the reference's pattern does not match (page `fRos`, `Pb` loci) and 13 further tokens (the remainder; ligature-only tokens, which the reference deletes, are the likely cause).  Numbers inside the brief's ranges are obtained with the `join` policy (unfiltered 36,278 tokens in all loci, 32,612 in paragraph text) and, with the default `split` policy, for clean tokens (37,978 in all loci, 33,205 in `P0` loci); the unfiltered `split` counts (39,026 and 35,075) lie above them.  The reference script itself splits on both `.` and `,`, i.e. it uses the same policy as the voytools default.

### Estimator check on a null text

After a within-line shuffle of the ZL3b paragraph text the vocabulary statistics are identical (True), the adjacent-token MI excess is -0.0034 bits (share -0.0004, z -0.46) and the edge-glyph MI excess is 0.0007 bits (z 0.64), as it must be.  The glyph entropy H(X_t+1|X_t) of the space-less text (extended preset, scope all) is 2.6894 for the real text and 2.7537 after the shuffle: it is not a bag-of-words statistic because glyph pairs at word boundaries depend on word order.

### Currier languages (paragraph text, extended preset, scope line)

| language | tokens | types | hapax share | H1 | H2 |
|---|---|---|---|---|---|
| A | 10 995 | 3 111 | 0.6863 | 3.9450 | 2.7555 |
| B | 22 742 | 4 654 | 0.6755 | 3.9554 | 2.4493 |

