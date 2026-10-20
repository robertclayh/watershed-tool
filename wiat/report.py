"""Figures and markdown reports.

Plots use a restrained, colourblind-validated palette. The treated basin and the
estimate are always the saturated series, and placebo/donor context is drawn in
neutral grey so the eye goes to the claim being tested.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, PLACEBO = "#0b0b0b", "#52514e", "#e4e3df", "#b9b8b3"
SERIES = (BLUE, ORANGE, AQUA)


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
        "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
        "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False,
        "lines.linewidth": 2.0, "figure.dpi": 110,
    })
    return plt


def _fmt(x: float, digits: int = 4) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    if np.isinf(x):
        return "−∞" if x < 0 else "∞"
    return f"{x:.{digits}g}"


# ----------------------------------------------------------------------------- plots
def plot_synthetic_control(res, ax=None, show_placebos: bool = True):
    """Treated vs synthetic trajectory (top) and gap with placebo gaps (bottom)."""
    plt = _plt()
    if ax is None:
        fig, ax = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True, height_ratios=[3, 2])
    else:
        fig = ax[0].figure
    times = res.times
    start = [t for t, p in zip(times, res.pre_mask) if not p][0]
    a0, a1 = ax
    a0.plot(times, res.y_treated.values, color=BLUE, label=f"treated ({res.treated})")
    a0.plot(times, res.y_synthetic.values, color=ORANGE, ls="--", label="synthetic control")
    a0.axvline(start, color=MUTED, lw=1, ls=":")
    a0.set_ylabel(getattr(res, "outcome", "outcome"))
    a0.legend(loc="best")
    a0.set_title(f"Synthetic control — ATT {_fmt(res.att)}  [{_fmt(res.ci[0])}, {_fmt(res.ci[1])}], "
                 f"p = {_fmt(res.p_value, 2)}")
    if show_placebos and res.placebo_gaps is not None:
        kept = res.placebo.index[res.placebo.kept] if res.placebo is not None else res.placebo_gaps.columns
        for c in kept:
            a1.plot(times, res.placebo_gaps[c].values, color=PLACEBO, lw=0.8, alpha=0.8)
        a1.plot([], [], color=PLACEBO, lw=0.8, label=f"placebo basins (n={len(kept)})")
    a1.plot(times, res.gap.values, color=BLUE, label="treated − synthetic")
    a1.axhline(0, color=MUTED, lw=0.8)
    a1.axvline(start, color=MUTED, lw=1, ls=":")
    a1.set_ylabel("gap")
    a1.legend(loc="best")
    fig.tight_layout()
    return fig


def plot_event_study(es, ax=None, title: str = "Event study (reference: t = −1)"):
    plt = _plt()
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 3.6))
    else:
        fig = ax.figure
    c = es.coefs
    x = c.index.to_numpy()
    if "pl_low" in c:
        ax.fill_between(x, c.pl_low, c.pl_high, color=PLACEBO, alpha=0.45, lw=0,
                        label=f"placebo range ({1 - es.alpha:.0%})")
    ax.errorbar(x, c.coef, yerr=[c.coef - c.ci_low, c.ci_high - c.coef], fmt="o", color=BLUE,
                ms=5, capsize=0, lw=1.5, label="estimate (cluster-robust CI)")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.axvline(-0.5, color=MUTED, lw=1, ls=":")
    ax.set_xlabel("periods relative to intervention")
    p = es.pretrend_placebo_p if np.isfinite(es.pretrend_placebo_p) else es.pretrend_p
    ax.set_title(f"{title} — pre-trend p = {_fmt(p, 2)}")
    ax.legend(loc="upper left")
    fig.tight_layout()
    return fig


def plot_donor_thinning(df: pd.DataFrame, truth: float | None = None, ax=None):
    """Estimates and intervals versus donor-pool size (output of ``donor_thinning``)."""
    plt = _plt()
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 3.8))
    else:
        fig = ax.figure
    g = df.groupby("k")
    series = []
    for col, color, off, name in (("did", BLUE, -0.12, "DiD (placebo CI)"), ("scm", ORANGE, 0.12, "SCM (placebo-inverted CI)")):
        series.append((g[f"{col}_estimate"].median(), g[f"{col}_ci_low"].median(), g[f"{col}_ci_high"].median(),
                       color, off, name))
    # One scale for both series, set by the typical finite bound so a single enormous
    # (but finite) interval does not flatten the rest. Anything beyond is clipped.
    vals = np.abs(np.concatenate([np.r_[e.values, lo.values, hi.values] for e, lo, hi, *_ in series]))
    vals = vals[np.isfinite(vals)]
    span = 1.3 * float(np.percentile(vals, 75)) if len(vals) else 1.0
    edge = span * 1.6
    for est, lo, hi, color, off, name in series:
        k = np.arange(len(est))
        bounded = np.isfinite(lo) & np.isfinite(hi)
        ax.vlines((k + off)[bounded], np.clip(lo[bounded], -edge, edge), np.clip(hi[bounded], -edge, edge),
                  color=color, lw=2)
        ax.vlines((k + off)[~bounded], -edge, edge, color=color, lw=1.2, ls=(0, (2, 2)))
        for kk in k[~bounded]:
            ax.annotate("∞", (kk + off, edge), ha="center", va="bottom", fontsize=8, color=color)
        ax.plot(k + off, est, "o", color=color, ms=6, label=name)
    ax.set_ylim(-edge * 1.12, edge * 1.12)
    ax.set_xticks(np.arange(len(est)), [str(v) for v in est.index])
    ax.set_xlabel("donor basins available (k nearest matches)")
    ax.axhline(0, color=MUTED, lw=0.8)
    if truth is not None:
        ax.axhline(truth, color=AQUA, lw=1.5, ls="--", label="true effect")
    ax.legend(loc="lower right")
    ax.set_title("Same study, shrinking donor pool (dashed = unbounded)")
    fig.tight_layout()
    return fig


SCARCITY_PANELS = (
    ("confident_correct", "P(confident, correct attribution)", (0, 1)),
    ("coverage", "Interval coverage (target 1 − α)", (0, 1)),
    ("ci_bounded", "Share of intervals that are finite", (0, 1.02)),
    ("rel_rmse", "RMSE / |true effect|", None),
    ("size", "False-positive rate under no effect", (0, 1)),
    ("median_rel_width", "Median CI width / |true effect|", None),
)


def plot_scarcity(summary: pd.DataFrame, estimator: str = "scm", alpha: float = 0.1,
                  panels=SCARCITY_PANELS, title: str | None = None, labels: dict | None = None):
    """Degradation curves: one panel per metric, one line per treated-basin extremity."""
    plt = _plt()
    s = summary[summary.estimator == estimator]
    exts = sorted(s.extremity.unique())
    ncol = 3
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(10.5, 3.1 * nrow))
    labels = labels or {0.0: "typical basin", 1.0: "atypical (1 SD)", 2.0: "extreme (2 SD)"}
    for ax, (metric, label, ylim) in zip(axes.ravel(), panels):
        for e, color in zip(exts, SERIES):
            d = s[s.extremity == e].sort_values("n_donors")
            ax.plot(d.n_donors, d[metric], "-o", color=color, ms=4, label=labels.get(e, f"extremity {e}"))
        if metric == "coverage":
            ax.axhline(1 - alpha, color=MUTED, ls=":", lw=1)
        if metric == "size":
            ax.axhline(alpha, color=MUTED, ls=":", lw=1)
        ax.axvline(1 / alpha - 1, color=PLACEBO, lw=6, alpha=0.35)
        ax.set_xscale("log")
        ticks = sorted(s.n_donors.unique())
        ax.set_xticks(ticks, [str(t) for t in ticks])
        ax.minorticks_off()
        if ylim:
            ax.set_ylim(*ylim)
        ax.set_title(label)
        ax.set_xlabel("donor basins")
    axes.ravel()[0].legend(loc="upper left")
    for ax in axes.ravel()[len(panels):]:
        ax.set_visible(False)
    fig.suptitle(title or f"Attribution quality vs donor scarcity — {estimator}", fontweight="bold")
    fig.tight_layout()
    return fig


def plot_estimator_comparison(summary: pd.DataFrame, metric: str, extremity: float,
                              estimators=("did_cluster", "did_placebo", "scm", "ascm"),
                              ylabel: str | None = None, ref: float | None = None, ax=None):
    plt = _plt()
    if ax is None:
        fig, ax = plt.subplots(figsize=(5.5, 3.4))
    else:
        fig = ax.figure
    names = {"did_cluster": "DiD, cluster SE (naive)", "did_placebo": "DiD, placebo inference",
             "scm": "Synthetic control", "ascm": "Augmented SCM"}
    colors = {"did_cluster": PLACEBO, "did_placebo": BLUE, "scm": ORANGE, "ascm": AQUA}
    s = summary[summary.extremity == extremity]
    for est in estimators:
        d = s[s.estimator == est].sort_values("n_donors")
        if len(d):
            ax.plot(d.n_donors, d[metric], marker="o", ms=4, color=colors[est], label=names[est],
                    ls="--" if est == "did_cluster" else "-")
    if ref is not None:
        ax.axhline(ref, color=MUTED, ls=":", lw=1)
    ax.set_xscale("log")
    ticks = sorted(s.n_donors.unique())
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlabel("donor basins")
    ax.set_ylabel(ylabel or metric)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------- markdown
def _table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(_fmt(v) if isinstance(v, (float, np.floating)) else str(v) for v in r) + " |")
    return "\n".join(out)


def report_markdown(rep) -> str:
    lines = [f"# Attribution report: `{rep.outcome}`", ""]
    lines += [f"**Verdict: {rep.verdict}**", ""]
    treated = ", ".join(f"`{b}` (from {s})" for b, s in rep.treat_start.items())
    lines += [f"- Treated basin(s): {treated}", f"- Donor basins used: {len(rep.donors)}",
              f"- Significance level: α = {rep.alpha}", ""]
    if rep.selection is not None:
        lines += ["## Donor pool", "", _table(rep.selection.summary().reset_index()), ""]
        if rep.selection.excluded:
            lines += [f"Excluded candidates: {len(rep.selection.excluded)} "
                      f"({', '.join(f'{k}: {v}' for k, v in list(rep.selection.excluded.items())[:6])})", ""]
    lines += ["## Failure diagnostics", ""]
    if rep.flags:
        order = {"fatal": 0, "major": 1, "minor": 2, "info": 3}
        for f in sorted(rep.flags, key=lambda f: order[f.severity]):
            lines.append(f"- **{f.severity.upper()} · {f.code}**: {f.message}")
    else:
        lines.append("- None raised.")
    lines += ["", "## Estimates", "", _table(rep.estimates()), ""]
    if rep.did is not None:
        lines += ["### DiD inference comparison", "",
                  "Different inference methods applied to the same point estimate. With few "
                  "treated basins, cluster-robust intervals are typically far too narrow.", "",
                  _table(rep.did.summary().reset_index()), ""]
    for b, s in rep.scm.items():
        w = s.weights[s.weights > 0.005].sort_values(ascending=False)
        lines += [f"### Synthetic control weights: `{b}`", "",
                  _table(w.rename("weight").reset_index().rename(columns={"index": "donor"})), "",
                  "Fit diagnostics: " + ", ".join(f"{k} = {_fmt(v, 3)}" for k, v in s.diagnostics.items()), ""]
    if rep.event_study is not None:
        lines += ["### Event study", "", _table(rep.event_study.coefs.reset_index()), ""]
    lines += ["## Reading this report", "",
              "- *Detection threshold*: the smallest effect this donor pool could distinguish from "
              "the spread of placebo effects in untreated basins. Effects below it are invisible to "
              "the design, not shown to be absent.",
              "- An unbounded (∞) interval means there are too few donors to rule out *any* "
              "effect size at the chosen α.", ""]
    return "\n".join(lines)


def write_report(rep, outdir: str | Path) -> Path:
    """Write ``report.md`` plus figures for a study into ``outdir``."""
    plt = _plt()
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    md = report_markdown(rep)
    figs = []
    for b, s in rep.scm.items():
        s.outcome = rep.outcome
        fig = plot_synthetic_control(s)
        name = f"scm_{b}.png"
        fig.savefig(out / name, bbox_inches="tight")
        plt.close(fig)
        figs.append(name)
    if rep.event_study is not None:
        fig = plot_event_study(rep.event_study)
        fig.savefig(out / "event_study.png", bbox_inches="tight")
        plt.close(fig)
        figs.append("event_study.png")
    if figs:
        md += "\n## Figures\n\n" + "\n".join(f"![{f}]({f})" for f in figs) + "\n"
    (out / "report.md").write_text(md, encoding="utf-8")
    return out / "report.md"
