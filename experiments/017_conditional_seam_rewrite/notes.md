# 017 – deník session

## Rozhodnutí a interpretace spec

* **Kód**: `voytools/placerules.py` doplněn jen přidáním (`ConditionalRules`, `learn_conditional_rules`,
  `apply_conditional_rules`, `CONDITIONAL_RULES`); R1–R3 (`PlaceRules`, `learn_place_rules`, `apply_place_rules`)
  beze změny. Sady `+R1+R2` a `+R1+R2+R3` jdou přes `apply_place_rules` (kód 014), ostatní přes
  `apply_conditional_rules`. Tabulky se učí z reálu s afixy (b) jako v 014 (jinak by `-m/-g`, `p-/f-` nebyly afixy).
* **Pořadí přepisu**: R1 (prefix krajních tokenů) → R2 (sufix posledního tokenu) → R2′ (poslední glyf tokenů uvnitř
  řádku) → přečtení kontextů z přepsaných částí (poslední glyf předchozího tokenu, vlastní sufix, první glyf kmene,
  hrubá pozice) → R3/R3′a/b/c (prefix vnitřních tokenů) → R3″ (první glyf kmene vnitřních tokenů bez prefixu).
  Kontext „poslední glyf předchozího tokenu“ se čte jednou před prefixovými pravidly (vektorově, jako 014); přepis
  prefixu ho změní jen u tokenů s prázdným kmenem i sufixem, R3″ u jednoglyfových kmenů – zanedbáno.
* **Minimální vazba** (014): token si afix nechá s pravděpodobností min(1, q(a)/p(a)), jinak losuje z přebytku
  (q − p)₊; p = rozdělení v generátoru **v dané buňce kontextu**, q = tabulka reálu pro tu buňku.
* **Backoff**: buňka tabulky reálu s < 5 pozorováními → méně podmíněná tabulka: R3′c → R3′a → R3 → nepodmíněné
  (prefixy vnitřních tokenů); R3′b → R3 → nepodmíněné; R3″ → nepodmíněné (první glyfy kmenů tokenů bez prefixu).
  Tokeny generátoru se seskupují podle (úroveň, klíč) a přepisují proti tabulce té úrovně. Podíl tokenů z backoffu
  je `share_backoff` v `runs.csv` / `rewrite_cost.csv`; podpora tabulek reálu v `table_support.csv`.
* **R2′**: množina koncových glyfů = glyfy s ≥ 20 výskyty na konci řádku a poměrem podíl(konec)/podíl(uvnitř) > 5
  (naučeno: `g`, `m`). U tokenů, které nejsou poslední v řádku, se glyf z množiny nechá s pravděpodobností
  min(1, q/p) (q = podíl mezi nekoncovými tokeny reálu, p = v generátoru), jinak se nahradí losem z rozdělení
  posledních glyfů nekoncových tokenů reálu bez koncové množiny. Glyf se mění v sufixu, má-li ho token, jinak v kmeni.
* **R3′c** hrubá pozice: první/vnitřní/poslední token × první řádek odstavce ano/ne (6 tříd; u vnitřních tokenů
  reálně jen vnitřní/poslední × 2).
* **Licence typů**: přepis přijat jen tehdy, když výsledný token má v reálu ≥ 2 výskyty; jinak nový los z přebytku
  (max 5×), pak původní část zůstává. Platí pro všechny kroky sady (`R12pA3L`). Zamítnuté losy = `n_rejected_*`.
* **z-skóre dvěma způsoby**: `z_seed` = (reál − průměr sady)/SD přes 30 seedů; `z_boot` = (reál − průměr sady)/SD
  reálu přes 100 **bootstrapů stránek** (tažení stránek s opakováním, kopie stránky dostane folio `f~k`, aby její
  řádky zůstaly vlastními jednotkami; pořadí stránek se přečísluje, aby lokalitní metriky viděly kopii ve
  vzdálenosti 1). `z_strict` = to z obou s větší |z|; verdikty používají `z_strict`. Bootstrap stránek *vychyluje*
  typové metriky (typy, hapax: duplicitní stránky snižují počet typů) – používáme jen jeho SD, ne průměr.
* **Vlastnost** = skupina metrik (`summarize.PROPERTIES`); „dosaženo“ = všechny metriky |z_strict| < 3 ve split
  **i** join. Verdikt proti základnímu generátoru = pravidlo 014 na největším |z| skupiny. Kalibrační cíl 004
  „loss“ není jedna metrika profilu – zastupují ho jeho složky (H2 glyfů, délka, typy, hapax, Zipf).
* **Minimální nedosažitelná kombinace**: podmnožina cílových vlastností (oddíl 1 syntézy bez A/B + soudržnost slova:
  prefix–sufix, sekční vazba prefixu, šířka paradigmatu, vysvětlený podíl H(sufix|kmen), typy, hapax), kterou žádná
  sada nedosáhne najednou a každá její vlastní podmnožina některou sadou dosažena je; hledáno do velikosti 4.
* **Samokontrola**: sada `+R1+R2′+R3′a` s tabulkami naučenými z ts2020 (seed + 1000) aplikovaná na ts2020 (seed),
  30 seedů; `z_mean` = (průměr samokontroly − průměr ts2020)/SD ts2020 (seedy párové: `paired_t`).
* **Šablona 005**: jen ZL_split – konfigurace šablony (`ablation_base`) existuje pro ZL_split a GC_split; join nelze.
  Základní šablona přepočítána touž cestou (sada `base`).
* **GC**: ts2020 je Eva-only (jako v 014) → jen ZL; uvedeno v reportu.

## Odchylky a poznámky

* Prostředí: Python 3.13 / pandas 3.0 / numpy 2.5 (cloud); kód drží 3.10 / pandas 2.3 kompatibilitu (žádný scipy).
* Plain ts2020 se počítá znovu touž cestou (`base`), ne přebírá z 009–011 jako v 014 – kvůli nákladovým sloupcům
  a shodné cestě; hodnoty se s 014 shodují v rámci SD.
