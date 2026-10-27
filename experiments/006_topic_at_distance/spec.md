# 006 – Téma na dálku: sdílení kmenů mezi stránkami a popisky versus text

## Otázka
Kopírování vysvětlí podobnost slov, která jsou blízko sebe; obsah by měl vytvářet podobnost mezi
místy, která spojuje téma, ne vzdálenost. Otázka: **sdílejí stránky stejné sekce víc kmenů, než
odpovídá jejich kodikologické vzdálenosti, a sdílejí vzdálené stránky téže sekce víc než sousední
stránky různých sekcí?** Druhá otázka: **chovají se popisky (zejména ok-/ot-) jako jména, která mají
stopu v běžném textu téže stránky?**

## Hypotézy a predikce
- **Téma (obsah):** sdílení kmenů mezi stránkami klesá se vzdáleností jen slabě a má silnou složku
  „stejná sekce“; dvojice (stejná sekce, daleko) sdílí víc než (jiná sekce, blízko). Po permutaci
  stránek uvnitř sekce se vzdálenostní složka ztratí, sekční zůstane.
- **Kopírování:** sdílení klesá monotónně se vzdáleností a sekce nehraje roli nad rámec vzdálenosti
  (sekce jsou souvislé bloky, proto je nutné obojí modelovat zároveň).
- **Písař:** část „sekčního“ efektu může být písařská (hand); testovat uvnitř písaře.
- **Popisky:** kmeny popisků se objevují v textu téže stránky častěji než v textu jiných stránek téže
  sekce (nad náhodu danou frekvencí); tokeny na ok-/ot- v běžném textu mají distribuci odlišnou od
  ostatních (pozice v řádku, sousedé), jak by odpovídalo substantivům/jménům.

## Data a filtry
ZL3b a GC2a; `P` text `clean` pro sdílení; `L` lokusy (popisky) pro druhou část, typy `Lz` (zodiak),
`Lf`, `Lc` (farmacie), `Ls` (hvězdy), `Ln` (nymfy). Politika `split` (hlavní), `join` (kontrola).
Dekompozice `decompose_table` s naučenými afixy; **kmen** = sloupec `stem`. Vyřadit 30
nejčastějších kmenů (funkční slova) – uvést seznam; citlivost: 0 / 30 / 100 vyřazených.

## Metriky
- Pro každou dvojici stránek (i, j) s ≥ 40 tokeny `P`: `share_ij` = kosinová podobnost vektorů
  frekvencí kmenů (log(1+f)) a Jaccard kmenů s ≥ 2 výskyty; `dist_ij` = |pořadí stránky i − j|
  v souboru (kodikologická vzdálenost); `same_section`, `same_language`, `same_hand`.
- Křivka sdílení podle vzdálenosti (biny 1, 2–3, 4–7, 8–15, 16–31, 32+), zvlášť pro stejná/jiná sekce.
- Rozklad (OLS, numpy): `share ~ log(dist) + same_section + same_language + same_hand` – koeficienty
  s intervaly z bootstrapu po stránkách (200) a z nullů.
- Klíčový kontrast: průměr `share` pro (stejná sekce, dist ≥ 16) minus (jiná sekce, dist ≤ 3).
- Popisky: pro každý popisek kmen → podíl popisků, jejichž kmen se vyskytuje v `P` textu téže stránky
  (`hit_same_page`), v textu jiné stránky téže sekce (`hit_same_section`), jinde; očekávání z frekvence
  kmene (permutace popisků mezi stránkami sekce, 200×). Zvlášť popisky začínající ok-/ot- a ostatní.
- Ok-/ot- v textu: rozdělení `rel_pos`, podíl na začátku/konci řádku, nejčastější předchůdci
  a následovníci (glyfové třídy), délka; porovnat s tokeny na `qo-`, `ch-/sh-`, `d-`.

## Nulové modely (stejný kód)
1. Reál s permutovanými stránkami **globálně** (ruší vzdálenost i sekci) – 100×.
2. Reál s permutovanými stránkami **uvnitř sekce** (ruší vzdálenost, zachová sekci) – 100×.
3. Reál s permutovanými stránkami uvnitř písaře – 100×.
4. Kopírovací generátor: `autocopy_like` s config 86 z 003 (`edit_model='trigram'`) – 50 seedů;
   pokud je k dispozici `voytools.selfcitation` (004 sloučeno), i `selfcitation_like` – 50 seedů.
5. Pro popisky: permutace popisků mezi stránkami téže sekce (200×).

## Výstupy (`results/006_topic_at_distance/`)
`pairs.csv` (dvojice stránek s metrikami – jen reál, obě abecedy), `curve.csv` (sdílení × vzdálenost ×
stejná/jiná sekce × zdroj), `ols.csv` (koeficienty s intervaly, zdroje), `contrast.csv`, `labels.csv`,
`okot_profile.csv`, `runs.csv`/`meta.json`, `curve.png`, `labels.png`, `report.md`, `notes.md`.

## Kritéria přijetí
1. Křivky a OLS pro reál a všech 4–5 nullů, ZL i GC, politiky split i join; kontrast „stejná sekce
   daleko vs. jiná sekce blízko“ s intervaly.
2. Report odpoví: (i) je sekční složka nad vzdálenostní? (ii) přežije kontrolu písaře? (iii) reprodukuje
   kopírovací null sekční složku? (iv) mají popisky stopu v textu téže stránky nad náhodu? (v) chovají
   se ok-/ot- tokeny odlišně? – každé ZL ano/ne, GC ano/ne.
3. Citlivost na vyřazení funkčních kmenů (0/30/100) uvedena.
4. Podezřelé: kosinus > 0,95 u dvojic různých stránek (duplicitní stránky?), kontrast > 0,5.

## Rozpočet
227 stránek → ~25 k dvojic; metriky vektorové, rychlé (< 1 s na kompletní sadu). Nully 450 běhů
≈ 20–40 min, generátory 100 vzorků ≈ 20 min. Celkem ≈ 1,5 h, jedna session, jedna větev, PR.
