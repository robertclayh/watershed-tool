"""End-to-end attribution study with explicit failure diagnostics.

:class:`AttributionStudy` chains donor matching, difference-in-differences, an
event-study pre-trend check and synthetic control. It then evaluates a fixed set of
**failure conditions** and turns them into one of four verdicts:

==============================  =====================================================
``ATTRIBUTED``                  Design is adequate, both estimators detect an effect
                                of the same sign.
``NO DETECTABLE EFFECT``        Design is adequate, nothing larger than the
                                detection threshold is present.
``INCONCLUSIVE``                Estimates exist, but identification is doubtful
                                (pre-trends, poor pre-fit, estimators disagree...).
``NOT ATTRIBUTABLE (DESIGN)``   The donor pool cannot support a comparative claim
                                at the chosen confidence level, whatever the data say.
==============================  =====================================================

The last category is the one practitioners most need to see. For many basins where
nature-based investment is most needed, the honest answer is that a comparative
design cannot attribute the change.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .did import DiDResult, EventStudyResult, did, event_study, staggered_att
from .matching import DonorSelection, balance_table, basin_signatures, select_donors
from .panel import TIME, UNIT, aggregate, balance
from .scarcity import donors_needed
from .synth import synthetic_control

FATAL, MAJOR, MINOR, INFO = "fatal", "major", "minor", "info"


@dataclass
class Flag:
    code: str
    severity: str
    message: str

    def __str__(self) -> str:
        return f"[{self.severity.upper()}] {self.code}: {self.message}"


@dataclass
class AttributionReport:
    outcome: str
    treat_start: dict
    alpha: float
    donors: list
    selection: DonorSelection | None
    balance: dict
    did: DiDResult | None
    event_study: EventStudyResult | None
    scm: dict
    staggered: object | None
    flags: list = field(default_factory=list)
    verdict: str = ""
    dropped_units: list = field(default_factory=list)

    def flags_by(self, severity: str) -> list:
        return [f for f in self.flags if f.severity == severity]

    def estimates(self) -> pd.DataFrame:
        rows = []
        if self.did is not None:
            rows.append({
                "estimator": f"DiD ({self.did.inference})", "basin": "all treated",
                "estimate": self.did.estimate, "ci_low": self.did.ci[0], "ci_high": self.did.ci[1],
                "p_value": self.did.p_value, "detection_threshold": self.did.detection_threshold,
            })
        for b, s in self.scm.items():
            rows.append({
                "estimator": f"SCM ({s.method})", "basin": b, "estimate": s.att,
                "ci_low": s.ci[0], "ci_high": s.ci[1], "p_value": s.p_value,
                "detection_threshold": np.nan,
            })
        if self.staggered is not None:
            st = self.staggered
            rows.append({
                "estimator": "Callaway-Sant'Anna", "basin": "all treated", "estimate": st.att,
                "ci_low": st.ci[0], "ci_high": st.ci[1], "p_value": np.nan,
                "detection_threshold": np.nan,
            })
        return pd.DataFrame(rows)

    def to_markdown(self) -> str:
        from .report import report_markdown

        return report_markdown(self)


def _to_period(value, freq: str, water_year_start_month: int):
    if freq is None:
        return value
    if isinstance(value, (int, np.integer)):
        return int(value)
    ts = pd.Timestamp(value)
    if freq == "water_year":
        ts = ts + pd.DateOffset(months=(12 - water_year_start_month + 1) % 12)
    # An intervention starting mid-year contaminates that year; the first fully
    # treated year is the next one. Starts on the first month count from that year.
    first_month = water_year_start_month if freq == "water_year" else 1
    on_boundary = pd.Timestamp(value).month == first_month and pd.Timestamp(value).day == 1
    return ts.year if on_boundary else ts.year + 1


class AttributionStudy:
    """Configure and run an attribution study.

    Parameters
    ----------
    panel
        Long panel, monthly (``date`` = month start) or already aggregated.
    outcome
        Outcome column, e.g. ``"et_mm"``, ``"water_yield_mm"``, ``"evaporative_index"``,
        ``"ndvi"``.
    treat_start
        ``{basin_id: intervention start}``. Dates are converted to the first fully
        treated period after aggregation.
    static
        Static basin attributes (index = basin id) used for matching.
    freq
        ``"year"``, ``"water_year"`` or ``None`` (use the panel as given).
    did_covariates
        Time-varying covariates for DiD, typically precipitation.
    match_on
        Columns of the matching feature table (signatures plus ``static``). Defaults
        to hydrologically meaningful ones, where available: mean precipitation, its
        seasonality, aridity, the outcome's baseline mean, elevation, slope and forest
        cover. Matching on many noisy features makes every basin look atypical.
    caliper, hard_caliper
        RMS and per-covariate limits on standardised differences between a donor and
        the treated basin (see :func:`wiat.matching.select_donors`).
    min_donors
        If matching yields fewer suitable donors, the nearest ``min_donors`` candidates
        are used anyway and the study is flagged. Attribution is still reported, never
        silently skipped.
    exclude
        ``{basin_id: reason}``, e.g. hydrologically connected basins.
    """

    def __init__(
        self,
        panel: pd.DataFrame,
        outcome: str,
        treat_start: Mapping,
        static: pd.DataFrame | None = None,
        unit: str = UNIT,
        time: str = TIME,
        freq: str | None = "year",
        water_year_start_month: int = 1,
        did_covariates: Sequence[str] = ("precip_mm",),
        match_on: Sequence[str] | None = None,
        caliper: float = 1.0,
        hard_caliper: float = 2.5,
        min_donors: int = 5,
        max_donors: int | None = None,
        exclude: Mapping | None = None,
        scm_method: str = "demeaned",
        alpha: float = 0.1,
        event_window: tuple[int, int] = (-5, 5),
        seed: int = 0,
    ):
        self.panel = panel
        self.outcome = outcome
        self.treat_start_raw = dict(treat_start)
        self.static = static
        self.unit, self.time = unit, time
        self.freq = freq
        self.wy = water_year_start_month
        self.did_covariates = [c for c in did_covariates if c in panel.columns and c != outcome]
        self.match_on = match_on
        self.caliper = caliper
        self.hard_caliper = hard_caliper
        self.min_donors = min_donors
        self.max_donors = max_donors
        self.exclude = dict(exclude or {})
        self.scm_method = scm_method
        self.alpha = alpha
        self.event_window = event_window
        self.seed = seed

    # ------------------------------------------------------------------ steps
    def _features(self) -> pd.DataFrame:
        first = min(self.treat_start_raw.values())
        p = self.panel
        is_dt = not np.issubdtype(p[self.time].dtype, np.number)
        baseline_end = pd.Timestamp(first) if is_dt else first
        variables = [c for c in ("precip_mm", "et_mm", "ndvi", "evi", self.outcome) if c in p]
        variables = list(dict.fromkeys(variables))
        feats = basin_signatures(
            p, baseline_end, variables=variables,
            trend_variables=[self.outcome] if self.outcome in p else (),
            unit=self.unit, time=self.time,
        )
        if self.static is not None:
            feats = feats.join(self.static.select_dtypes("number"), how="left")
        return feats

    def _match_covariates(self, feats: pd.DataFrame) -> list:
        default = ["precip_mm_mean", "precip_mm_seasonality", "aridity", f"{self.outcome}_mean",
                   "elevation_m", "slope_deg", "forest_frac"]
        return [c for c in dict.fromkeys(self.match_on or default) if c in feats.columns]

    def run(self) -> AttributionReport:
        flags: list[Flag] = []
        treated = list(self.treat_start_raw)

        # 1. Matching: who can credibly stand in for the treated basin(s)?
        feats = self._features()
        cov = self._match_covariates(feats)
        sel = select_donors(feats, treated, covariates=cov, caliper=self.caliper,
                            hard_caliper=self.hard_caliper, exclude=self.exclude, max_donors=self.max_donors)
        donors = list(sel.suitable)
        need = donors_needed(self.alpha)
        if len(donors) < self.min_donors:
            ranked = sel.distances[treated].min(axis=1).sort_values().index.tolist()
            fallback = [d for d in ranked if d not in donors][: self.min_donors - len(donors)]
            if fallback:
                flags.append(Flag(
                    "UNMATCHED_DONORS_USED", MAJOR,
                    f"only {len(donors)} candidate basin(s) pass the matching caliper "
                    f"(RMS {self.caliper} SD, max {self.hard_caliper} SD); "
                    f"added {len(fallback)} nearest unmatched basin(s) {fallback}. Comparisons rest "
                    "on basins that differ materially from the treated basin.",
                ))
            donors += fallback

        # 2. Panel at analysis frequency, balanced on what we need.
        if self.freq is not None:
            data = aggregate(self.panel, self.freq, unit=self.unit, time=self.time,
                             water_year_start_month=self.wy)
        else:
            data = self.panel
        starts = {u: _to_period(s, self.freq, self.wy) for u, s in self.treat_start_raw.items()}
        data = data[data[self.unit].isin(treated + donors)]
        data, dropped = balance(data, [self.outcome, *self.did_covariates], self.unit, self.time)
        if any(t in dropped for t in treated):
            raise ValueError(f"treated basin(s) {[t for t in treated if t in dropped]} have gaps in "
                             f"{self.outcome}; fill or shorten the study period")
        donors = [d for d in donors if d not in dropped]
        if dropped:
            flags.append(Flag("DONORS_DROPPED_FOR_GAPS", MINOR,
                              f"{len(dropped)} donor(s) dropped for incomplete records: {dropped[:8]}"))
        J = len(donors)

        # 3. Design adequacy: purely a function of the donor count.
        if J < need["significance_possible"]:
            flags.append(Flag(
                "DONOR_SCARCITY", FATAL,
                f"{J} usable donor basin(s). Placebo inference can never reach p <= {self.alpha} "
                f"with fewer than {need['significance_possible']} (min attainable p = {1/(J+1):.3f}).",
            ))
        elif J < need["stable_threshold"]:
            flags.append(Flag(
                "THIN_PLACEBO_DISTRIBUTION", MINOR,
                f"{J} usable donor basin(s). The detection threshold and interval width rest on the "
                f"most extreme one or two placebo basins (fewer than {need['stable_threshold']}), so "
                "one unusual control can make or break significance.",
            ))

        bal = {}
        for t in treated:
            bt = balance_table(feats, t, donors, covariates=cov)
            bal[t] = bt
            out_cov = bt.index[bt.outside_pool_range].tolist()
            if out_cov:
                flags.append(Flag(
                    "TREATED_OUTSIDE_DONOR_SUPPORT", MINOR,
                    f"{t} lies outside the donor pool's range on {out_cov}. Results extrapolate there.",
                ))

        # 4. Estimators.
        did_res = es = st = None
        if J >= 1:
            did_res = did(data, self.outcome, starts, covariates=self.did_covariates,
                          unit=self.unit, time=self.time, alpha=self.alpha, seed=self.seed)
            es = event_study(data, self.outcome, starts, window=self.event_window,
                             covariates=self.did_covariates, unit=self.unit, time=self.time,
                             alpha=self.alpha, seed=self.seed)
            if not es.pretrend_ok:
                p = es.pretrend_placebo_p if np.isfinite(es.pretrend_placebo_p) else es.pretrend_p
                flags.append(Flag("PRE_TREND", MAJOR,
                                  f"pre-intervention differences in trend (p = {p:.3f}). "
                                  "Parallel trends is doubtful, so DiD is not credible."))
            if len(set(starts.values())) > 1:
                if J >= 2:
                    st = staggered_att(data, self.outcome, starts, unit=self.unit, time=self.time,
                                       alpha=self.alpha, seed=self.seed)
                flags.append(Flag("STAGGERED_ADOPTION", MINOR,
                                  "multiple start dates. TWFE DiD can be biased; compare with the "
                                  "Callaway-Sant'Anna estimate."))
            if did_res.n_treated == 1:
                flags.append(Flag("SINGLE_TREATED_BASIN", INFO,
                                  "cluster-robust and wild-bootstrap inference are unreliable with "
                                  "one treated basin; placebo (Conley-Taber) inference is used."))

        scm = {}
        if J >= 2:
            for t in treated:
                s = synthetic_control(data, self.outcome, t, starts[t], donors=donors, unit=self.unit,
                                      time=self.time, method=self.scm_method, alpha=self.alpha)
                scm[t] = s
                d = s.diagnostics
                if d["pre_rmspe_normalized"] > 0.5:
                    flags.append(Flag("POOR_PRE_FIT", MAJOR,
                                      f"{t}: synthetic control misses the pre-period by "
                                      f"{d['pre_rmspe_normalized']:.2f} SD (pre-R² {d['pre_r2']:.2f})."))
                if d["frac_pre_outside_donor_range"] > 0.25:
                    flags.append(Flag("OUTSIDE_CONVEX_HULL", MAJOR,
                                      f"{t}: {d['frac_pre_outside_donor_range']:.0%} of pre-periods lie "
                                      "outside the donor range, which no weighting can reproduce."))
                n_disc = d.get("placebos_discarded", 0)
                if n_disc and s.min_attainable_p > self.alpha:
                    flags.append(Flag("PLACEBOS_DISCARDED", MAJOR,
                                      f"{t}: {n_disc} placebo basin(s) fit their own pre-period too poorly "
                                      f"to be informative and were discarded. Only {len(donors) - n_disc} "
                                      f"remain, so SCM cannot reach p <= {self.alpha}."))
                n_perfect = d.get("placebos_perfect_fit", 0)
                if n_perfect and not np.all(np.isfinite(s.ci)) and J >= need["significance_possible"]:
                    n_pre = int(np.sum(np.asarray(sorted(data[self.time].unique())) < starts[t]))
                    flags.append(Flag("SCM_OVERFIT", MINOR,
                                      f"{t}: {n_perfect} of {J} placebo basins reproduce their {n_pre}-period "
                                      "baseline exactly, so no treated gap can look extreme and the SCM "
                                      "interval is unbounded. Donors outnumber pre-periods: use fewer donors "
                                      "or rely on DiD."))
                if d["effective_donors"] < 2:
                    flags.append(Flag("CONCENTRATED_WEIGHTS", MINOR,
                                      f"{t}: synthetic control rests on ~{d['effective_donors']:.1f} "
                                      f"basin(s) (max weight {d['max_weight']:.2f}). Idiosyncratic "
                                      "shocks in that donor pass straight into the estimate."))

        # 5. Cross-estimator consistency.
        if did_res is not None and scm:
            scm_att = np.mean([s.att for s in scm.values()])
            disagree = np.sign(scm_att) != np.sign(did_res.estimate)
            lo, hi = did_res.ci
            if np.isfinite(lo) and np.isfinite(hi) and not (lo <= scm_att <= hi):
                disagree = True
            if disagree:
                flags.append(Flag("ESTIMATOR_DISAGREEMENT", MAJOR,
                                  f"DiD ({did_res.estimate:.4g}) and SCM ({scm_att:.4g}) disagree. The "
                                  "result depends on the identifying assumption chosen."))
            if np.isfinite(did_res.detection_threshold) and abs(did_res.estimate) < did_res.detection_threshold:
                flags.append(Flag("BELOW_DETECTION_THRESHOLD", INFO,
                                  f"|DiD estimate| {abs(did_res.estimate):.4g} is below what this donor "
                                  f"pool can distinguish from noise ({did_res.detection_threshold:.4g})."))

        n_pre = min(starts.values()) - min(data[self.time]) if self.freq else None
        if n_pre is not None and n_pre < 5:
            flags.append(Flag("SHORT_BASELINE", MINOR,
                              f"only {n_pre} pre-intervention periods; pre-trend and pre-fit checks are weak."))

        rep = AttributionReport(
            outcome=self.outcome, treat_start=starts, alpha=self.alpha, donors=donors,
            selection=sel, balance=bal, did=did_res, event_study=es, scm=scm, staggered=st,
            flags=flags, dropped_units=dropped,
        )
        rep.verdict = self._verdict(rep)
        return rep

    def _verdict(self, rep: AttributionReport) -> str:
        if rep.flags_by(FATAL) or rep.did is None:
            return "NOT ATTRIBUTABLE (DESIGN)"
        if rep.flags_by(MAJOR):
            return "INCONCLUSIVE"
        did_sig = rep.did.significant
        scm_sig = all(s.significant for s in rep.scm.values()) if rep.scm else False
        same_sign = all(np.sign(s.att) == np.sign(rep.did.estimate) for s in rep.scm.values())
        if did_sig and scm_sig and same_sign:
            return "ATTRIBUTED"
        if not did_sig and not any(s.significant for s in rep.scm.values()):
            return "NO DETECTABLE EFFECT"
        return "INCONCLUSIVE"
