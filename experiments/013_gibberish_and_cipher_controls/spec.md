# 013 – Chybějící referenční body: lidský blábol a šifra skutečného textu

## Otázka
Dva zdroje, na které se otázka smyslu nejvíc odvolává, v našich profilech dosud chybí: **lidský blábol**
(Gaskell & Bowern 2022: 42 dobrovolníků psalo nesmyslný text; data veřejně) a **šifra skutečného textu**
(Naibbe šifra latiny, Greshko 2025; vzorky v klonu Rozanové & Temereva). Otázka: **které z Voynich-specifických
vlastností (009/010/011) má lidský blábol a které má šifra latiny?** Blábol je null pro hypotézu „bez smyslu,
bez algoritmu“; šifra latiny je kalibrace, jak na našich metrikách vypadá *smysl pod obfuskací* – ne jako
hypotéza o rukopisu (vyloučeno v `CLAUDE.md`), ale jako referenční bod.

## Data
1. **Blábol:** `git clone --depth 1 https://github.com/danielgaskell/voynich` – složka `data/` (anonymizované
   přepisy; formát zjistit a popsat: řádky? stránky? délka na pisatele?). Licence: modified MIT, citovat
   článek. Vytvořit: (a) **sdružený korpus** (všichni pisatelé za sebou, pisatel = „sekce“ i „písař“, pořadí
   pisatelů náhodné s pevným seedem) vysázený do layoutu ZL (jako kontroly 009; pokud je tokenů méně než
   34 268, uvést a profil počítat na dostupném počtu a porovnávat s Voynichem oříznutým na stejný počet);
   (b) **po pisatelích** tam, kde metrika na krátký text stačí (slovník, H2, hapax, délky; řádkové metriky,
   má-li přepis řádky). Tokenizace: malá písmena, interpunkce pryč, glyf = písmeno (preset `None`); afixy
   naučené stejně jako u kontrol.
2. **Naibbe:** z klonu `lrozanova/voynich-units`: `data/controls/naibbe/naibbe_caesar_ciphertext_lines.txt`
   (Caesar) a `greshko_nathist_output_ciphertext.txt` (Plinius; plus `greshko_nathist_pre_encryption_respaced_plaintext.txt`
   pro srovnání „před šifrou“). Jsou v Eva → stejné presety glyfů jako u Voynichu. Vysázet do layoutu ZL;
   sekce = souvislé bloky (knihy nejsou v textu, uvést). Ořez na 34 268 tokenů.
3. **Srovnávací sady** (převzít z 009/010/011, nepočítat znovu, pokud jsou v `results/`): Voynich ZL/GC,
   ts2020, šablona, latinská próza, Linnaeus položky.

## Metriky
Celý registrový profil 009 (37 metrik, kód `experiments/009_register_profile/run.py`) + řádkové metriky
010 (hranice vs vnitřek, vertikální podobnost, konce řádků) + rozklad 011 (vysvětlený podíl H(sufix | kmen)
pozicí/„písařem“), pokud je 011 v `main` (jinak jen H(sufix | kmen) a sufixů na kmen). Pro blábol po
pisatelích navíc: rozptyl mezi pisateli u hapaxu, H2, sufixů na kmen.

## Analýza
- `profile.csv` jako v 009 se zdroji `gibberish_pooled`, `gibberish_<id>` (je-li), `naibbe_caesar`,
  `naibbe_pliny`, `naibbe_plaintext`.
- Vzdálenosti profilů (metoda 009) Voynich ↔ blábol, Voynich ↔ Naibbe, Naibbe ↔ latina (kolik šifra zachová),
  blábol ↔ próza; po skupinách metrik.
- **Seznamy:** pro každou Voynich-specifickou vlastnost z `results/SYNTHESIS.md` oddílu 1 (poziční a sekční
  vazby afixů, šířka a vyrovnanost paradigmat, hapax, nízké pořadí tokenů, okrajová MI, konec řádku, uzavřený
  řádek): má ji blábol? má ji Naibbe? (ano/ne/částečně s čísly a intervaly).
- Zvlášť: I(afix; pisatel) u blábolu vs I(afix; písař) u Voynichu – kolik „registru“ dělá prostě jiný člověk.

## Výstupy (`results/013_gibberish_and_cipher_controls/`)
`sources_used.csv` (soubory, tokeny, pisatelé, vysázení), `profile.csv`, `distances.csv`, `feature_lists.csv`,
`writers.csv`, `profile_heatmap.png`, `runs.csv`/`meta.json`, `report.md`, `notes.md`. Načítání do
`voytools/controls.py` (rozšířit, s testy).

## Kritéria přijetí
1. Oba zdroje načtené a popsané (počty tokenů, pisatelů, řádků); profil kompletní tam, kde délka stačí,
   jinak explicitně „nelze (krátký text)“.
2. `feature_lists.csv` pokrývá všechny vlastnosti oddílu 1 syntézy; report odpoví: které Voynich-specifické
   vlastnosti má lidský blábol, které šifra latiny, které ani jeden.
3. Podezřelé: Naibbe plaintext a ciphertext se stejným počtem typů (chyba načtení); blábol s `hapax_share`
   < 0,3 (opakování – zkontrolovat, zda není chyba tokenizace).
4. Testy, Python 3.10 / pandas 2.3, bez scipy.

## Rozpočet
Klonování + načtení 15 min; profily 3–5 zdrojů × (bod + 100 subsamplů) ≈ 40 min. Jedna session, jedna větev, PR.
