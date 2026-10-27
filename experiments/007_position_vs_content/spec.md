# 007 – Pozice versus obsah: co vysvětluje tvar slova a co zbývá po odečtení pozice

## Otázka
Tvar tokenu závisí na pozici v řádku a odstavci (Currier, Feaster, Stafford; 003: poziční JSD 0,19)
a slovník na sekci. V jazyce vysvětluje pozice hlavně *tvary* (afixy, funkční slova) a obsah hlavně
*kmeny*; v pozičním generátoru vysvětluje pozice všechno. Otázky: **(1) jak se dělí informace
o afixech a kmenech mezi pozici, sekci a písaře? (2) Zůstává po odečtení pozičních efektů závislost
mezi tokeny na vzdálenost 1–10, a jak rychle klesá?** Jazyk: pomalý pokles (syntax, pak téma);
generátor a blábol: lokální a rychle mizející.

## Hypotézy a predikce
- Kmen: informace ze sekce > z pozice; afixy: z pozice ≥ ze sekce. Kopírovací null: kmen ani afix
  nenesou sekční informaci nad rámec sousedství; šablonový null z 005 (s tématem): kmen sekční,
  afix žádnou poziční.
- Reziduální vzájemná informace (po podmínění poziční třídou) mezi jednotkami tokenů na vzdálenost
  d = 1…10: reál klesá pomaleji než kopírovací null (ten má jen okno) a než míchání v rámci poziční
  třídy (≈ 0). Vzhledem k 001/003 je známo, že na úrovni *identity tokenů* je pořadí slabé (0,9 %);
  proto měřit na úrovni (a) prvního/posledního glyfu, (b) BPE jednotek (64 sloučení), (c) afixových
  tříd, (d) kmenů s ≥ 20 výskyty + `<other>`.

## Data a filtry
ZL3b a GC2a, `P` text `clean`, politiky `split` a `join`. Dekompozice naučenými afixy (společné
seznamy). Poziční třídy jako v `positional_profile` (line_first/mid/line_last × parfirst) + `rel_pos`
binovaný na 5 úrovní + `line_in_par` (0, 1, 2+). Sekce = `illustration`, jazyk = `language`,
písař = `hand`.

## Metriky
1. **Rozklad informace** (plug-in entropie s korekcí mícháním; bity): pro cílovou proměnnou Y ∈
   {prefix, suffix, stem_class (top 200 + other), délka tokenu} spočítat H(Y), I(Y; pozice),
   I(Y; sekce), I(Y; písař), I(Y; sekce | pozice), I(Y; pozice | sekce), I(Y; sekce | písař);
   nejistota: 200 subsamplů po stránkách bez vracení (80 %) + permutační null (Y permutováno
   uvnitř stránky). Výsledek jako tabulka a graf „podíl H(Y) vysvětlený“.
2. **Křivka závislosti na vzdálenost:** pro jednotky U ∈ {poslední glyf→první glyf, BPE-64 jednotky,
   afixová třída, kmenová třída} spočítat I(U_t; U_{t+d}) pro d = 1…10 uvnitř odstavce (přes hranice
   řádků ano, přes hranice odstavců ne), odečíst průměr z (a) míchání tokenů uvnitř poziční třídy
   přes celý odstavec (zachová poziční efekty, zruší posloupnost) a (b) míchání uvnitř řádku. Totéž
   pro nully: `autocopy_like` config 86 (003) a `selfcitation_like` (pokud 004 sloučeno), šablonový
   null z 005 (pokud sloučeno) – každý 30 seedů. Shrnout polovinou poklesu (d, kde křivka klesne na
   50 % hodnoty při d = 1) a plochou pod křivkou d = 2…10.
3. **Rightward/downward drift** (Feaster): frekvence glyfových tříd jako funkce `rel_pos` a
   `line_in_par`; uvést jako kovariát, zda zbytková závislost po jeho odečtení přetrvá (varianta
   poziční třídy = jemné biny).

## Nulové modely
Míchání uvnitř poziční třídy (nové: `nulls.shuffle_within_position_class(df, seed)` – přidat do
`voytools/nulls.py` **jen tuto funkci**, s testem; pokud 004 změnilo `nulls.py`, nejdřív `git pull`
a rebase), míchání v řádcích, `autocopy_like` (003 config 86), `selfcitation_like` a šablona 005,
jsou-li v `main`. 30 seedů každý; reál 200 subsamplů.

## Výstupy (`results/007_position_vs_content/`)
`info_decomposition.csv`, `mi_distance.csv` (U × d × zdroj: MI, korigovaná MI, intervaly),
`drift.csv`, `runs.csv`/`meta.json`, `decomposition.png`, `mi_distance.png`, `report.md`, `notes.md`.

## Kritéria přijetí
1. Rozklad pro 4 cílové proměnné × ZL/GC × split/join, s intervaly; věta ke každé: pozice vs sekce
   vs písař, a zda to kopírovací null reprodukuje.
2. Křivky MI(d) pro 4 jednotky × zdroje; polovina poklesu a plocha s intervaly; explicitně říct,
   zda reál po odečtení pozičních efektů drží závislost dál než na d = 2 v ZL i GC.
3. Podezřelé: korigovaná MI < −0,01 bitu (chyba korekce), I(Y; pozice) > H(Y).
4. Testy procházejí; Python 3.10 / pandas 2.3; bez scipy.

## Rozpočet
Entropie a MI jsou tabulkové operace (numpy), ~1–3 s na sadu; 200 subsamplů × 8 sad + 4 nully × 30
seedů × 8 sad ≈ 1–1,5 h na 4 procesech. Jedna session, jedna větev, PR na konci.
