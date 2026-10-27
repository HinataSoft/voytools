# Repo na GitHubu a cloudové session – postup

## 1. Založení repa (jednorázově, na vašem počítači)

V `D:\Work\Voynich\analysis` (PowerShell nebo Git Bash):

```
git init -b main
git add .
git commit -m "voytools 0.1.0: IVTFF parser, token tables, stats, morphology, nulls, experiment runner"
```

Na GitHubu založte prázdné repo (např. `voynich-analysis`, privátní nebo veřejné – kód i data
jsou malé, ~4 MB) a připojte ho:

```
git remote add origin https://github.com/<user>/voynich-analysis.git
git push -u origin main
```

Poznámky:

* `.gitignore` už vynechává cache, `*.egg-info`, virtuální prostředí a velké mezivýsledky.
* `data/` obsahuje kopie transliterací z voynich.nu (viz `data/README.md`). U veřejného repa
  zvažte, zda je nechat, nebo je nahradit `data/download.sh`.
* `pyproject.toml` deklaruje licenci MIT; pokud chcete jinou, změňte ji před prvním pushem a
  přidejte soubor `LICENSE`.

## 2. Cloudové prostředí (claude.ai/code → Environments)

Doporučené nastavení prostředí, ve kterém poběží session nad tímto repem:

| Položka | Hodnota | Proč |
|---|---|---|
| Network access | **Trusted** (výchozí) | stačí pip + GitHub; data jsou v repu. `Custom` + `www.voynich.nu` jen pokud chcete `data/download.sh`. |
| Setup script | obsah `scripts/setup_env.sh` | nainstaluje balík s pytest a matplotlib (< 1 min) a prostředí se zacachuje |
| Env var `BASH_DEFAULT_TIMEOUT_MS` | `600000` | výchozí limit příkazu 10 min místo 2 |
| Env var `BASH_MAX_TIMEOUT_MS` | `3600000` | příkaz na pozadí smí běžet až hodinu |
| Env var `VOYTOOLS_DATA` | (nenastavovat) | testy i příklady najdou `data/` v repu samy |

Limity VM podle dokumentace: 4 vCPU, 16 GB RAM, 30 GB disk, Ubuntu 24.04, Python 3 s pip/uv/pytest.
Po několika minutách nečinnosti se VM pozastaví (soubory zůstanou); pokud je později zrušena,
běžící procesy zmizí – proto experimenty zapisují průběžně a commitují (viz `experiments/README.md`).

## 3. Spuštění experimentu

Z prohlížeče (claude.ai/code → New session → repo → prompt), z desktopové aplikace, nebo z terminálu
ve složce repa:

```
claude --cloud "Proveď experiment experiments/001_autocopy_calibration/spec.md. Pracuj na větvi exp/001-autocopy-calibration, výsledky do results/001_autocopy_calibration/, průběžně pushuj."
```

Paralelní běhy: každé session dejte disjunktní rozsah seedů, např.

```
claude --cloud "Experiment 002 podle spec.md, seedy seed_range(0, 2000, 4)[0], výsledky do results/002_paradigm_coherence/part0/"
claude --cloud "... seed_range(0, 2000, 4)[1] ... part1/"
```

a nakonec jedna session (nebo vy lokálně) spojí části:
`python -c "from voytools import merge_runs; merge_runs(['results/002_paradigm_coherence/part0', ...], 'results/002_paradigm_coherence')"`.

Běžící session lze řídit z terminálu: `claude -p "zpráva" --cloud <session-id>`.

## 4. Jak se výsledky dostanou zpět a jak do toho vstupuje Cowork

1. Session pushne větev `exp/...` (commit + případně PR).
2. Vy na Windows: `git fetch && git checkout exp/...` nebo merge do `main` a `git pull`.
3. Cowork (připojená složka `D:\Work\Voynich`) výsledky vidí okamžitě v `analysis\results\...`
   a může je interpretovat, porovnat s literaturou a napsat další `spec.md`.
4. Nový spec pushnete, další session ho vezme.

Cowork umí z veřejného repa i přímo `git pull` (pushovat ne), u privátního repa jde vše přes váš
checkout. Založit cloud session z Coworku nelze; zakládáte je vy.

## 5. Kontrola, že vše funguje (první session)

Prompt pro první session: *„Spusť `python -m pytest -q` a `python -m voytools summary data/ZL3b-n.txt`,
porovnej s VERIFICATION.md a nahlas rozdíly. Nic neměň.“* Očekávaný výsledek: 167 testů prošlo,
227 stran / 5 385 lokusů / 39 026 tokenů.
