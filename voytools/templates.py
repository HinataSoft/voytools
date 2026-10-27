"""Structured-template null: an artificial language with paradigms and topics by construction.

Experiment 002 found that Voynich stems take a restricted set of suffixes, that prefix and
suffix are strongly associated, that the occurrences of a stem are spread through the text and
that the stem distribution changes between sections.  :class:`ParadigmTemplate` builds a text
that has all of this *by construction* and carries no meaning (a "systematic glossolalia" or an
invented language with tables), so that one can ask which of those properties a meaningless
system reproduces and how many hand-made rules it needs.

Mechanism (one sample = one new artificial language plus one text in it)
-----------------------------------------------------------------------
* **Affixes** are the real prefix and suffix inventories with their real token frequencies
  (taken from a decomposed real table, see :func:`voytools.morphology.decompose_table`).
* **Stems** are *artificial*: ``n_stems`` distinct random words produced by a slot/trigram word
  grammar learned on the real stems (:meth:`voytools.wordgrammar.WordGrammar.sample`).  They are
  ranked in the order in which a grammar walk first produces them (a token-weighted walk first,
  so that the strings the grammar finds likely get the top ranks; a type-weighted walk fills the
  rest) and get Zipf weights ``rank ** -zipf_s``.  Short artificial stems inevitably coincide
  with real ones; the share is reported in ``attrs['template_info']``.
* **Paradigm**: with ``paradigm='fixed'`` every stem owns a suffix set of size ``k``; ``k`` is
  drawn from the real numbers of distinct suffixes of stems with about the same frequency
  (bins of ``log2`` count, using the stem's realised count in the sample), the members are drawn without
  replacement with the real suffix frequencies, and a token draws its suffix from the stem's
  set with the renormalised real frequencies.  ``paradigm='free'``: every token draws its
  suffix from the real marginal (slot independence).
* **Prefix**: ``prefix_dependence='none'`` draws from the real marginal; ``'stem'`` gives every
  stem a prefix set exactly like the suffix set; ``'suffix'`` draws from the real conditional
  distribution P(prefix | suffix); ``'stem_suffix'`` draws inside the stem's prefix set with
  weights P(prefix | suffix) (falling back to the set's marginal weights).
* **Topic**: ``topic_concentration`` α.  Each stem gets section weights
  θ ~ Dirichlet(α · K · base) (``K`` sections, ``base`` = token share of each section in the
  layout), and a token in section j draws its stem with probability ∝ w · θ_j.  ``α = inf``
  means no topic (θ = base); small α = strong topic.
* **Page locality** λ: every page gets a page vocabulary of ``page_vocab`` draws from the
  section distribution; with probability λ a token draws its stem from that page vocabulary
  (proportional to the draw counts), otherwise from the section distribution.  This makes stems
  clump on pages without copying tokens.
* **Suffix shift** ``suffix_by`` (default ``None``): the real suffix frequencies used for the
  sets and the token draws are taken separately for each value of this column of the layout
  (e.g. ``'language'``); this is the hand-made rule needed for a Currier A/B suffix shift.

The text has exactly the layout and metadata of the given template table
(:func:`voytools.nulls.as_table`); no token is ever copied from the real text or from earlier
in the generated text.  All randomness comes from the ``seed`` of :meth:`ParadigmTemplate.sample`.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Any

import numpy as np
import pandas as pd

from ._common import line_ids, resolve_preset
from .nulls import as_table
from .wordgrammar import WordGrammar, learn_word_grammar

__all__ = ["ParadigmTemplate", "PARADIGMS", "PREFIX_DEPENDENCE"]

PARADIGMS = ("fixed", "free")
PREFIX_DEPENDENCE = ("none", "stem", "suffix", "stem_suffix")
_MISSING = "?"


def _freq(values: list[str]) -> tuple[list[str], np.ndarray]:
    c = Counter(values)
    names = sorted(c)
    return names, np.asarray([c[n] for n in names], dtype=float)


def _k_by_bin(per_stem: pd.DataFrame, column: str) -> dict[int, np.ndarray]:
    """Observed numbers of distinct affixes of real stems, grouped by floor(log2(count))."""
    b = np.floor(np.log2(per_stem["count"].to_numpy(dtype=float))).astype(int)
    k = per_stem[column].to_numpy(dtype=int)
    return {int(v): k[b == v] for v in np.unique(b)}


def _validate(p: dict[str, Any]) -> None:
    if p["paradigm"] not in PARADIGMS:
        raise ValueError(f"paradigm must be one of {PARADIGMS}")
    if p["prefix_dependence"] not in PREFIX_DEPENDENCE:
        raise ValueError(f"prefix_dependence must be one of {PREFIX_DEPENDENCE}")
    if not (0.0 <= float(p["page_locality"]) <= 1.0):
        raise ValueError("page_locality must be in [0, 1]")
    if float(p["topic_concentration"]) <= 0:
        raise ValueError("topic_concentration must be > 0 (inf = no topic)")
    p["n_stems"] = int(p["n_stems"])


class ParadigmTemplate:
    """Generator of a meaningless structured text (see the module docstring).

    Build it with :meth:`learn` from a decomposed real table and call :meth:`sample`.
    Parameters (all keyword arguments of :meth:`learn`): ``n_stems`` (4000), ``zipf_s`` (1.1),
    ``paradigm`` ('fixed'), ``prefix_dependence`` ('none'), ``topic_concentration`` (inf),
    ``page_locality`` (0.0), ``page_vocab`` (20), ``suffix_by`` (None), ``section_column``
    ('illustration').
    """

    DEFAULTS: dict[str, Any] = {
        "n_stems": 4000, "zipf_s": 1.1, "paradigm": "fixed", "prefix_dependence": "none",
        "topic_concentration": math.inf, "page_locality": 0.0, "page_vocab": 20, "suffix_by": None,
        "section_column": "illustration",
    }

    def __init__(self, stats: dict[str, Any], grammars: tuple[WordGrammar, WordGrammar], params: dict[str, Any]):
        self.stats = stats
        self.grammars = grammars
        self.params = params
        self._pools: dict[int, list[str]] = {}

    # ------------------------------------------------------------------ learning

    @classmethod
    def learn(cls, dfd: pd.DataFrame, **params: Any) -> "ParadigmTemplate":
        """Learn affix frequencies, paradigm-size distributions and the stem grammar from ``dfd``.

        ``dfd`` is a decomposed real table (columns ``prefix``, ``stem``, ``suffix``, plus the
        ``suffix_by`` column if that option is used).  Unknown parameter names raise ``TypeError``.
        """
        unknown = set(params) - set(cls.DEFAULTS)
        if unknown:
            raise TypeError(f"unknown template parameters {sorted(unknown)}")
        p = {**cls.DEFAULTS, **params}
        _validate(p)
        preset = resolve_preset(dfd, "auto")
        pre, stem, suf = (dfd[c].astype(str).tolist() for c in ("prefix", "stem", "suffix"))
        st: dict[str, Any] = {"preset": preset}
        st["prefixes"], st["prefix_w"] = _freq(pre)
        st["suffixes"], st["suffix_w"] = _freq(suf)
        # P(prefix | suffix) as a (n_suffixes x n_prefixes) matrix
        pi = {a: i for i, a in enumerate(st["prefixes"])}
        si = {a: i for i, a in enumerate(st["suffixes"])}
        joint = np.zeros((len(si), len(pi)))
        np.add.at(joint, ([si[a] for a in suf], [pi[a] for a in pre]), 1.0)
        st["p_prefix_given_suffix"] = joint / joint.sum(axis=1, keepdims=True)
        # suffix frequencies per group (suffix_by)
        if p["suffix_by"] is not None:
            col = dfd[p["suffix_by"]].astype(object).where(dfd[p["suffix_by"]].notna(), _MISSING).astype(str)
            groups = {}
            for g, sub in pd.Series(suf).groupby(col.to_numpy()):
                c = Counter(sub.tolist())
                groups[str(g)] = np.asarray([c.get(a, 0) for a in st["suffixes"]], dtype=float)
            st["suffix_w_by"] = groups
        # paradigm sizes of real stems by frequency
        g = pd.DataFrame({"stem": stem, "p": pre, "s": suf}).groupby("stem")
        per = pd.DataFrame({"count": g.size(), "n_suffixes": g["s"].nunique(), "n_prefixes": g["p"].nunique()})
        st["k_suffix"] = _k_by_bin(per, "n_suffixes")
        st["k_prefix"] = _k_by_bin(per, "n_prefixes")
        st["real_stems"] = set(per.index.tolist())
        st["n_real_stems"] = len(per)
        g_tok = learn_word_grammar(stem, preset)
        g_type = learn_word_grammar(sorted(set(stem)), preset)
        return cls(st, (g_tok, g_type), p)

    def with_params(self, **params: Any) -> "ParadigmTemplate":
        """A copy with some parameters changed, sharing the learned statistics, the grammars and
        the stem-pool cache (cheap; used to run a parameter grid on one learned template)."""
        unknown = set(params) - set(self.DEFAULTS)
        if unknown:
            raise TypeError(f"unknown template parameters {sorted(unknown)}")
        p = {**self.params, **params}
        _validate(p)
        out = ParadigmTemplate(self.stats, self.grammars, p)
        out._pools = self._pools
        return out

    # ------------------------------------------------------------------ artificial stems

    def stem_pool(self, seed: int, n: int | None = None) -> list[str]:
        """The first ``n`` (default ``n_stems``) artificial stems for ``seed``, in rank order.

        Pools are cached per seed and extended on demand, so the stems of rank ``r`` are the same
        for every ``n_stems`` at a given seed.
        """
        n = self.params["n_stems"] if n is None else int(n)
        pool = self._pools.get(seed)
        if pool is not None and len(pool) >= n:
            return pool[:n]
        rng = np.random.default_rng([int(seed), 7_919])
        seen: dict[str, None] = {}
        for gram, patience in zip(self.grammars, (2_000, max(20_000, 20 * n))):
            misses = 0
            while len(seen) < n and misses < patience:
                w = gram.sample(rng)
                s = "".join(w) if w is not None else ""
                if s and s not in seen:
                    seen[s] = None
                    misses = 0
                else:
                    misses += 1
        pool = list(seen)
        if len(pool) < n:
            raise RuntimeError(f"the stem grammar produced only {len(pool)} distinct stems (< {n})")
        self._pools[seed] = pool
        return pool[:n]

    # ------------------------------------------------------------------ sampling

    @staticmethod
    def _sets(rng: np.random.Generator, counts: np.ndarray, kbins: dict[int, np.ndarray],
              weights: np.ndarray) -> dict[int, np.ndarray]:
        """Affix set (indices) of every stem with ``counts > 0``: size from the frequency bin of a
        real stem with the same count, members drawn without replacement by ``weights``."""
        keys = np.asarray(sorted(kbins))
        nz = int((weights > 0).sum())
        p = weights / weights.sum()
        out = {}
        for z in np.flatnonzero(counts > 0).tolist():
            b = int(np.floor(np.log2(counts[z])))
            key = keys[max(int(np.searchsorted(keys, b, side="right")) - 1, 0)]
            ks = kbins[int(key)]
            k = int(min(ks[int(rng.integers(len(ks)))], nz))
            out[z] = np.sort(rng.choice(len(weights), size=max(k, 1), replace=False, p=p))
        return out

    def sample(self, template_df: pd.DataFrame, seed: int = 0, pool_seed: int | None = None) -> pd.DataFrame:
        """A text in a new artificial language with the layout and metadata of ``template_df``.

        Returns a token table (see :func:`voytools.nulls.as_table`) with ``attrs['source'] =
        'template'`` and ``attrs['template_info']`` (number of stems, share of artificial stems
        that coincide with real stems, realised share of locality draws ...).  Deterministic for
        a given ``seed``.  ``pool_seed`` (default: ``seed``) selects the artificial stem inventory
        separately, so that several configurations can share the same stems (common random numbers).
        """
        p, st = self.params, self.stats
        rng = np.random.default_rng(int(seed))
        n = len(template_df)
        stems = self.stem_pool(int(seed if pool_seed is None else pool_seed))
        m = len(stems)
        w = np.arange(1, m + 1, dtype=float) ** (-float(p["zipf_s"]))
        w /= w.sum()

        # sections and topic: P(stem | section) as columns of an (m x K) matrix
        sec = template_df[p["section_column"]].astype(object).where(template_df[p["section_column"]].notna(), _MISSING)
        sec = sec.astype(str).to_numpy()
        sec_names, sec_codes = np.unique(sec, return_inverse=True)
        K = len(sec_names)
        base = np.bincount(sec_codes, minlength=K) / n
        alpha = float(p["topic_concentration"])
        if math.isinf(alpha):
            cond = np.repeat(w[:, None], K, axis=1)
        else:
            theta = rng.gamma(np.maximum(alpha * K * base, 1e-12), size=(m, K))
            theta /= np.maximum(theta.sum(axis=1, keepdims=True), 1e-300)
            cond = w[:, None] * theta
            cond /= cond.sum(axis=0, keepdims=True)

        # stems of all tokens
        z = np.empty(n, dtype=np.int64)
        for j in range(K):
            idx = np.flatnonzero(sec_codes == j)
            z[idx] = rng.choice(m, size=len(idx), p=cond[:, j])
        lam = float(p["page_locality"])
        n_local = 0
        if lam > 0:
            pages = template_df["folio"].astype(str).to_numpy()
            local = rng.random(n) < lam
            for pg in pd.unique(pages):
                idx = np.flatnonzero(pages == pg)
                j = int(np.bincount(sec_codes[idx]).argmax())
                vocab = rng.choice(m, size=int(p["page_vocab"]), p=cond[:, j])
                li = idx[local[idx]]
                z[li] = vocab[rng.integers(len(vocab), size=len(li))]
                n_local += len(li)

        counts = np.bincount(z, minlength=m)
        suffixes, prefixes = st["suffixes"], st["prefixes"]
        # suffix weights for every token (global or per suffix_by group)
        if p["suffix_by"] is not None:
            grp = template_df[p["suffix_by"]].astype(object).where(template_df[p["suffix_by"]].notna(), _MISSING)
            grp = grp.astype(str).to_numpy()
            wby = st["suffix_w_by"]
            grp_names = sorted(set(grp.tolist()))
            gw = np.stack([wby.get(g, st["suffix_w"]) for g in grp_names]) + 1e-12
            gcode = np.searchsorted(np.asarray(grp_names), grp)
        else:
            gw = st["suffix_w"][None, :] + 1e-12
            gcode = np.zeros(n, dtype=np.int64)
        gw = gw / gw.sum(axis=1, keepdims=True)

        s_idx = np.empty(n, dtype=np.int64)
        if p["paradigm"] == "free":
            for g in range(len(gw)):
                ii = np.flatnonzero(gcode == g)
                s_idx[ii] = rng.choice(len(suffixes), size=len(ii), p=gw[g])
        else:
            ssets = self._sets(rng, counts, st["k_suffix"], st["suffix_w"])
            order = np.lexsort((gcode, z))
            keys = z[order] * len(gw) + gcode[order]
            cuts = np.flatnonzero(np.diff(keys)) + 1
            for chunk in np.split(order, cuts):
                if not len(chunk):
                    continue
                zz, gg = int(z[chunk[0]]), int(gcode[chunk[0]])
                S = ssets[zz]
                q = gw[gg][S]
                s_idx[chunk] = S[rng.choice(len(S), size=len(chunk), p=q / q.sum())]

        pw = st["prefix_w"] / st["prefix_w"].sum()
        dep = p["prefix_dependence"]
        p_idx = np.empty(n, dtype=np.int64)
        if dep == "none":
            p_idx[:] = rng.choice(len(prefixes), size=n, p=pw)
        elif dep == "suffix":
            cps = st["p_prefix_given_suffix"]
            for s in np.unique(s_idx):
                ii = np.flatnonzero(s_idx == s)
                p_idx[ii] = rng.choice(len(prefixes), size=len(ii), p=cps[s])
        else:
            psets = self._sets(rng, counts, st["k_prefix"], st["prefix_w"])
            cps = st["p_prefix_given_suffix"]
            use_s = dep == "stem_suffix"
            order = np.lexsort((s_idx, z)) if use_s else np.argsort(z, kind="stable")
            keys = z[order] * len(suffixes) + (s_idx[order] if use_s else 0)
            cuts = np.flatnonzero(np.diff(keys)) + 1
            for chunk in np.split(order, cuts):
                if not len(chunk):
                    continue
                P = psets[int(z[chunk[0]])]
                q = cps[int(s_idx[chunk[0]])][P] if use_s else pw[P]
                if q.sum() <= 0:
                    q = pw[P]
                p_idx[chunk] = P[rng.choice(len(P), size=len(chunk), p=q / q.sum())]

        toks = [prefixes[a] + stems[b] + suffixes[c] for a, b, c in zip(p_idx.tolist(), z.tolist(), s_idx.tolist())]
        lengths = np.bincount(line_ids(template_df))
        lines, pos = [], 0
        for L in lengths.tolist():
            lines.append(toks[pos: pos + L])
            pos += L
        out = as_table(lines, template_df, glyph_preset=st["preset"])
        used = np.unique(z)
        out.attrs["source"] = "template"
        out.attrs["template_info"] = {
            "n_stems_used": int(len(used)),
            "share_stems_real": float(np.mean([stems[i] in st["real_stems"] for i in used])) if len(used) else float("nan"),
            "share_tokens_local": n_local / n if n else float("nan"),
            "mean_stem_len_chars": float(np.mean([len(stems[i]) for i in z])) if n else float("nan"),
        }
        return out
