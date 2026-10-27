# 009 – Registr: kde leží Voynich mezi narativem, katalogem a generátorem

## Otázka
Po 001–007 vypadá profil textu takto: hapaxy jako botanická latina, pořadí tokenů pod všemi prózami,
řádek jako jednotka, silné poziční konvence, žádná posloupnostní závislost za sousedem. To připomíná
**katalog, rejstřík nebo receptář** (položky, jména, formule), ne souvislou prózu. Otázka: **kterému
žánru skutečného textu je Voynich svým profilem nejblíž, a ve kterých vlastnostech se liší od všech
skutečných textů i od generátorů zároveň?** Experiment definuje „registrový profil“ a měří ho jednotně
na kontrolách, Voynichu a nullech.

## Kontroly (kladné)
Z klonu `git clone --depth 1 https://github.com/lrozanova/voynich-units` (adresář `data/voynich-units/`
a `voynich_decipherment_repro_bundle/`): latinský narativ (Caesar), lékařská latina (Celsus), botanická
latina (Plinius), anglický herbář (Culpeper), italská próza, a katalog (Linnaeus *Species Plantarum*),
jsou-li v klonu – uvést přesně, které soubory byly nalezeny a použity, a u chybějících to napsat.
Volitelně (jen pokud jde sparse checkout rychle, < 10 min): z `https://github.com/alexanderboxer/voynich-attack`
`corpora/latin/apothecary/ellis_1854` (lékárnický text) a `corpora/latin/legacy/pliny` jako druhá kopie.
Každá kontrola: malá písmena, interpunkce pryč, tokeny = slova; **oříznout na 34 268 tokenů**; „glyfy“ =
písmena (preset `None`); afixy naučit stejným algoritmem (`learn_affixes`), dekompozice stejná.

## Jak kontrolám dát strukturu Voynichu
Aby šly počítat poziční a stránkové metriky, každá kontrola se „vysází“ do layoutu ZL3b: posloupnost
délek řádků a stránek (počty tokenů na řádek a řádků na stránku) se převezme z reálné tabulky ZL split
(`generic_type="P"`), tj. kontrola = `as_table(lines, template)` s metadaty reálu kromě `illustration`.
**Sekce** kontrol: přirozené oddíly zdroje (knihy Plinia a Celsa, kapitoly Culpepera, rody u Linnéa),
namapované na stránky podle pořadí; u narativu kapitoly. Dvě varianty vysázení: (i) „próza“ – řádky
se lámou mechanicky; (ii) „položky“ – u Linnéa (a Culpepera, lze-li) řádek = jedna položka/heslo,
zkrácené/sloučené na délku řádku šablony. Uvést, která varianta se u které kontroly použila.

## Registrový profil (všechny metriky z `voytools`, stejný kód pro všechny zdroje)
1. Slovník: `hapax_share`, `n_types`, `top50_coverage`, Zipfův sklon (`vocabulary_summary`, `zipf_fit`).
2. Pořadí tokenů: `adjacent_token_mi` excess share (cap 2 000, 50 shufflů); hrana: `edge_glyph_mi`.
3. Pozice vs sekce (007): I(prefix / sufix / kmen / délka; pozice), I(·; sekce), I(·; sekce | pozice),
   korekce mícháním jako v 007 (kód z `experiments/007_…/run.py`).
4. Paradigma (002/004): `mean_suffixes_per_stem`, `frac_stems_multi_suffix`, `H_suffix_given_stem`,
   `mi_prefix_suffix`, `cmi_excess`, `stem_dispersion`, `jsd_stem_sections`.
5. Lokalita (006): koeficient `log_dist` a `same_section` modelu `spec+size`, kontrast (kód z 006).
6. Posloupnost (007): plocha MI(d) d = 2…10 s korekcí `pc` pro kmeny a afixy; hodnota při d = 1.
7. Řádek: podíl prvního/posledního tokenu řádku v nejčastějších třídách; I(první třída řádku; zbytek
   řádku) – jednoduchá verze toho, co 010 rozvede.

## Zdroje
Voynich ZL split a join, GC split; kontroly (5–7 textů × varianty vysázení); nully: ts2020 calB (30 seedů),
šablona 005 (30 seedů), ac-tri 86 (převzít z 003/007), míchání v řádcích. Nejistota: reál i kontroly
100 subsamplů 80 % stránek bez vracení (SD × 2 jako v 007); generátory přes seedy.

## Analýza
- `profile.csv`: metrika × zdroj (hodnota, interval). `profile_z.csv`: z-vzdálenost Voynichu od každé
  kontroly (SD z jejích subsamplů) a od každého nullu.
- **Vzdálenost profilů:** všechny metriky standardizovat (z vůči rozptylu přes zdroje), euklidovská
  vzdálenost Voynich–zdroj; pořadí zdrojů od nejbližšího; tytéž vzdálenosti s vyřazením skupiny metrik
  (slovník / pořadí / pozice / paradigma / lokalita) → které skupiny Voynich přibližují k prózám, které ke
  katalogu, které ke generátorům. Graf: `profile_radar.png` nebo heatmapa `profile_heatmap.png`.
- **Seznamy:** metriky, kde Voynich leží (a) uvnitř rozpětí skutečných textů, (b) mimo všechny skutečné
  texty i mimo nully (Voynich-specifické), (c) mimo skutečné texty, ale uvnitř nullů (generátorové).

## Výstupy (`results/009_register_profile/`)
`controls_used.csv` (soubor, zdroj, tokeny, sekce, vysázení), `profile.csv`, `profile_z.csv`,
`distances.csv`, `lists.csv`, grafy, `runs.csv`/`meta.json`, `report.md`, `notes.md`.

## Kritéria přijetí
1. Nejméně 5 skutečných kontrol včetně alespoň jedné katalogové/položkové (jinak uvést, co chybělo).
2. Profil kompletní pro všechny zdroje; seznamy (a)/(b)/(c); pořadí zdrojů podle vzdálenosti s uvedením,
   které skupiny metrik pořadí mění.
3. Report odpoví: je Voynich blíž položkovému textu než próze? Ve kterých metrikách je mimo všechny
   skutečné texty? Ve kterých mimo všechny nully? Podezřelé: kontrola s `hapax_share` > 0,9 (chyba
   tokenizace), kontrola s I(·; pozice) > Voynich u mechanicky lámané prózy (chyba vysázení).
4. Testy procházejí; Python 3.10 / pandas 2.3; bez scipy; žádné změny v `voytools` kromě případných
   pomocných funkcí pro načtení kontrol (`voytools/controls.py`, s testy).

## Rozpočet
Kontroly 7 × 2 varianty × profil (~30 s) + 100 subsamplů (~1 s každý na lehké metriky, MI(d) dražší) ≈
1 h; nully 60 vzorků ≈ 20 min. Celkem 1,5–2 h, jedna session, jedna větev, PR.
