"""Watershed Intervention Attribution Toolkit (wiat).

Attribute hydrologic change to nature-based interventions with difference-in-
differences and synthetic control on remote-sensing basin panels. The toolkit also
reports where and why that attribution fails.
"""

from .attribution import AttributionReport, AttributionStudy, Flag
from .did import did, event_study, staggered_att
from .matching import balance_table, basin_signatures, hydrologically_connected, select_donors
from .panel import add_hydrologic_outcomes, aggregate, balance, deseasonalize
from .scarcity import donor_thinning, donors_needed, run_scarcity_experiment, summarize_scarcity
from .simulate import simulate_watersheds
from .synth import synthetic_control

__version__ = "0.1.0"

__all__ = [
    "AttributionReport",
    "AttributionStudy",
    "Flag",
    "add_hydrologic_outcomes",
    "aggregate",
    "balance",
    "balance_table",
    "basin_signatures",
    "deseasonalize",
    "did",
    "donor_thinning",
    "donors_needed",
    "event_study",
    "hydrologically_connected",
    "run_scarcity_experiment",
    "select_donors",
    "simulate_watersheds",
    "staggered_att",
    "summarize_scarcity",
    "synthetic_control",
]
