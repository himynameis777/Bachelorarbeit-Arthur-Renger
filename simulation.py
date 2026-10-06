"""Calibrate the control system or run the final habitat-loss experiment.

This script runs only 0% habitat loss, disables predator reintroduction, and
records whether both populations persist. It is intentionally separate from
the final experiment so unstable parameter combinations can be rejected first.

Examples (macOS):
    python simulation.py --quick
    python simulation.py --automated-calibration
    python simulation.py --validate-candidates candidate_configurations_v2.csv
    python simulation.py --random-runs 2000 --repetitions 1 --max-ticks 2000
    python simulation.py --habitat-experiment --habitat-loss-values 0,25,50,75,100 \
        --repetitions 2 --max-ticks 3000 --record-every 10

The eight sampled ranges follow the methodology. The sampled screen is
reproducible for a given --master-seed and avoids materializing the complete
parameter space. Completed runs are checkpointed to CSV, including when the
user interrupts the program. The explicit candidate-validation mode selects
configurations from persistence and repeated bounded oscillations; it reports
the older endpoint-size and CV criteria as diagnostics but does not use them as
its primary gate.

NETLOGO_HOME can be supplied as an environment variable if NetLogo is installed
somewhere other than /Applications/NetLogo 7.0.3.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import math
import os
import platform
import random
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# NOTE: `pandas` and `pynetlogo` are intentionally imported inside `main()`
# to avoid triggering heavy/native imports at module import time (which can
# fail on systems missing developer tools). This file provides lightweight
# environment checks before attempting those imports.


MASTER_SEED = 24680
SCRIPT_DIR = Path(__file__).resolve().parent
# Use the renamed NetLogo model in the `code/` folder
DEFAULT_MODEL = SCRIPT_DIR / "abm.nlogox"
DEFAULT_NETLOGO_HOME = Path("/Applications/NetLogo 7.0.3")
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "baseline_outputs"
DEFAULT_AUTOMATED_OUTPUT_DIR = SCRIPT_DIR / "automated_calibration_outputs"
DEFAULT_CANDIDATE_OUTPUT_DIR = SCRIPT_DIR / "calibration_v2_validation"
DEFAULT_HABITAT_OUTPUT_DIR = SCRIPT_DIR / "habitat_loss_outputs"

# Configuration selected by the screening/calibration/validation workflow.
FINAL_CONFIGURATION_ID = 275
FINAL_CONFIGURATION = {
    "n_preds": 10,
    "prey-reproduction-chance": 4,
    "predator-reproduction-chance": 7,
    "max-prey-hungry": 24,
    "max-pred-hungry": 75,
    "predator-vision": 6,
    "prey-vision": 6,
    "crowding-hunger-cost": 4,
    "max-local-prey-density": 3,
}


def integer_values(value: str) -> list[int]:
    """Parse comma-separated integers or an inclusive start:stop:step range."""
    try:
        if ":" in value:
            parts = [int(item.strip()) for item in value.split(":")]
            if len(parts) != 3:
                raise ValueError
            start, stop, step = parts
            if step == 0 or (stop - start) * step < 0:
                raise ValueError
            endpoint = stop + (1 if step > 0 else -1)
            values = list(range(start, endpoint, step))
        else:
            values = [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Expected comma-separated integers or start:stop:step, "
            f"received: {value!r}"
        ) from exc
    if not values:
        raise argparse.ArgumentTypeError("At least one value is required.")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calibrate predator-prey parameters under 0% habitat loss or run "
            "the final habitat-loss experiment."
        )
    )
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument(
        "--netlogo-home",
        type=Path,
        default=Path(os.environ.get("NETLOGO_HOME", DEFAULT_NETLOGO_HOME)),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--max-ticks", type=int, default=2000)
    parser.add_argument("--record-every", type=int, default=20)
    parser.add_argument("--master-seed", type=int, default=MASTER_SEED)
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=1,
        help="Save all completed runs after every N runs (default: 1).",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help=(
            "Run a diverse 20-configuration pilot with one repetition and "
            "500 ticks."
        ),
    )
    parser.add_argument(
        "--quick-configurations",
        type=int,
        default=20,
        help="Number of diverse configurations used by --quick (default: 20).",
    )
    parser.add_argument(
        "--random-runs",
        type=int,
        default=None,
        help=(
            "Sample N diverse parameter configurations reproducibly. "
            "Use --repetitions 1 for a small first-stage screen."
        ),
    )
    parser.add_argument(
        "--automated-calibration",
        action="store_true",
        help=(
            "Run a reproducible two-stage calibration: a diverse screen, "
            "automatic candidate selection, and multi-seed validation."
        ),
    )
    parser.add_argument(
        "--validate-candidates",
        type=Path,
        default=None,
        help=(
            "Run validation only for the exact parameter configurations in "
            "the supplied CSV file."
        ),
    )
    parser.add_argument(
        "--habitat-experiment",
        action="store_true",
        help=(
            "Bypass baseline calibration and run configuration 275 across "
            "the requested habitat-loss levels."
        ),
    )
    parser.add_argument(
        "--habitat-loss-values",
        type=integer_values,
        default=list(range(0, 101, 5)),
        help=(
            "Habitat-loss percentages for --habitat-experiment, as comma-"
            "separated values or an inclusive start:stop:step range "
            "(default: 0:100:5)."
        ),
    )
    parser.add_argument(
        "--screen-configurations",
        type=int,
        default=300,
        help="Configurations in the automated first-stage screen (default: 300).",
    )
    parser.add_argument(
        "--screen-ticks",
        type=int,
        default=1000,
        help="Maximum ticks in each screening run (default: 1000).",
    )
    parser.add_argument(
        "--validation-ticks",
        type=int,
        default=2000,
        help="Maximum ticks in each validation run (default: 2000).",
    )
    parser.add_argument(
        "--validation-repetitions",
        type=int,
        default=10,
        help="Deterministic validation seeds per candidate (default: 10).",
    )
    parser.add_argument(
        "--validation-stability-window",
        type=int,
        default=1000,
        help=(
            "Final ticks used to assess bounded oscillations during validation "
            "(default: 1000)."
        ),
    )
    parser.add_argument(
        "--validation-record-every",
        type=int,
        default=10,
        help="Ticks between validation records (default: 10).",
    )
    parser.add_argument(
        "--max-validation-configurations",
        type=int,
        default=30,
        help="Maximum screen candidates sent to validation (default: 30).",
    )
    parser.add_argument(
        "--screen-min-final-preys",
        type=int,
        default=50,
        help="Minimum prey at the end of a screening run (default: 50).",
    )
    parser.add_argument(
        "--screen-min-final-predators",
        type=int,
        default=2,
        help="Minimum predators at the end of a screening run (default: 2).",
    )
    parser.add_argument(
        "--screen-max-cv",
        type=float,
        default=2.0,
        help="Maximum population CV allowed in screen candidates (default: 2.0).",
    )
    parser.add_argument(
        "--screen-max-relative-trend",
        type=float,
        default=1.0,
        help=(
            "Maximum absolute relative trend for either population in the "
            "screening stability window (default: 1.0)."
        ),
    )
    parser.add_argument(
        "--min-validation-persistence-rate",
        type=float,
        default=0.8,
        help="Minimum persistence rate for final selection (default: 0.8).",
    )
    parser.add_argument(
        "--min-validation-stability-rate",
        type=float,
        default=0.6,
        help=(
            "Minimum legacy stable-dynamics rate for the two-stage automated "
            "mode (default: 0.6; diagnostic only in candidate validation)."
        ),
    )
    parser.add_argument(
        "--min-validation-acceptance-rate",
        type=float,
        default=0.0,
        help=(
            "Optional minimum strict acceptance rate for final selection. "
            "The default 0.0 does not impose the endpoint-size filter."
        ),
    )
    parser.add_argument(
        "--min-validation-oscillation-rate",
        type=float,
        default=0.6,
        help=(
            "Minimum fraction of validation runs with repeated bounded "
            "oscillations in both populations (default: 0.6)."
        ),
    )
    parser.add_argument(
        "--min-validation-cycles",
        type=int,
        default=2,
        help=(
            "Minimum meaningful cycles per population in the final stability "
            "window (default: 2)."
        ),
    )
    parser.add_argument(
        "--min-cycle-relative-amplitude",
        type=float,
        default=0.1,
        help=(
            "Minimum peak-to-trough change relative to the population mean "
            "for a reversal to count (default: 0.1)."
        ),
    )
    parser.add_argument("--n-preds", type=int, default=10)
    parser.add_argument(
        "--prey-reproduction-values",
        type=integer_values,
        default=list(range(1, 6)),
        help="Prey reproduction percentages (default: 1:5:1).",
    )
    parser.add_argument(
        "--predator-reproduction-values",
        type=integer_values,
        default=list(range(5, 26)),
        help="Predator reproduction percentages (default: 5:25:1).",
    )
    parser.add_argument(
        "--max-prey-hungry-values",
        type=integer_values,
        default=list(range(20, 81)),
        help="Maximum prey hunger values (default: 20:80:1).",
    )
    parser.add_argument(
        "--max-pred-hungry-values",
        type=integer_values,
        default=list(range(20, 81)),
        help="Maximum predator hunger values (default: 20:80:1).",
    )
    parser.add_argument(
        "--predator-vision-values",
        type=integer_values,
        default=list(range(3, 11)),
        help="Predator vision values (default: 3:10:1).",
    )
    parser.add_argument(
        "--prey-vision-values",
        type=integer_values,
        default=list(range(3, 11)),
        help="Prey vision values (default: 3:10:1).",
    )
    parser.add_argument(
        "--crowding-hunger-cost-values",
        type=integer_values,
        default=list(range(1, 5)),
        help="Crowding hunger costs (default: 1:4:1).",
    )
    parser.add_argument(
        "--max-local-prey-density-values",
        type=integer_values,
        default=list(range(3, 9)),
        help="Maximum local prey densities (default: 3:8:1).",
    )
    parser.add_argument(
        "--stability-window",
        type=int,
        default=500,
        help="Number of final ticks used to assess stability (default: 500).",
    )
    parser.add_argument("--min-final-preys", type=int, default=50)
    parser.add_argument("--max-final-preys", type=int, default=500)
    parser.add_argument("--min-final-predator-ratio", type=float, default=0.5)
    parser.add_argument("--max-final-predator-ratio", type=float, default=2.0)
    parser.add_argument("--max-prey-cv", type=float, default=0.5)
    parser.add_argument("--max-predator-cv", type=float, default=0.75)
    parser.add_argument(
        "--max-relative-window-trend",
        type=float,
        default=0.5,
        help=(
            "Maximum absolute fitted population change across the stability "
            "window, divided by its mean (default: 0.5)."
        ),
    )
    return parser.parse_args()


def set_variable(netlogo: pyNetLogo.NetLogoLink, name: str, value) -> None:
    if isinstance(value, bool):
        rendered = "true" if value else "false"
    else:
        rendered = value
    netlogo.command(f"set {name} {rendered}")


def report_number(netlogo: pyNetLogo.NetLogoLink, reporter: str, fallback=math.nan):
    try:
        return float(netlogo.report(reporter))
    except Exception:
        return fallback


def parameter_dimensions(args: argparse.Namespace) -> list[tuple[str, list[int]]]:
    """Return the eight parameter dimensions specified in the methodology."""
    return [
        ("prey-reproduction-chance", args.prey_reproduction_values),
        ("predator-reproduction-chance", args.predator_reproduction_values),
        ("max-prey-hungry", args.max_prey_hungry_values),
        ("max-pred-hungry", args.max_pred_hungry_values),
        ("predator-vision", args.predator_vision_values),
        ("prey-vision", args.prey_vision_values),
        ("crowding-hunger-cost", args.crowding_hunger_cost_values),
        ("max-local-prey-density", args.max_local_prey_density_values),
    ]


def parameter_space_size(args: argparse.Namespace) -> int:
    size = 1
    for _, values in parameter_dimensions(args):
        size *= len(values)
    return size


def full_parameter_grid(args: argparse.Namespace) -> list[dict]:
    """Build a complete grid only when the configured space is manageable."""
    size = parameter_space_size(args)
    if size > 100_000:
        raise ValueError(
            f"The complete parameter grid contains {size:,} configurations. "
            "Use --quick or --random-runs N to sample it without exhausting memory."
        )

    dimensions = parameter_dimensions(args)
    names = [name for name, _ in dimensions]
    grid = itertools.product(*(values for _, values in dimensions))
    return [
        {"n_preds": args.n_preds, **dict(zip(names, combination))}
        for combination in grid
    ]


def diverse_parameter_sample(
    args: argparse.Namespace, sample_size: int, seed: int
) -> list[dict]:
    """Select a reproducible discrete Latin-hypercube-style sample."""
    if sample_size < 1:
        raise ValueError("Sample size must be at least 1.")
    space_size = parameter_space_size(args)
    if sample_size > space_size:
        raise ValueError(
            f"Requested {sample_size:,} configurations from a space of "
            f"only {space_size:,}."
        )

    dimensions = parameter_dimensions(args)
    generator = random.Random(seed)
    sampled_columns: dict[str, list[int]] = {}

    for name, values in dimensions:
        ordered_values = sorted(set(values))
        if sample_size == 1:
            indices = [(len(ordered_values) - 1) // 2]
        else:
            indices = [
                round(index * (len(ordered_values) - 1) / (sample_size - 1))
                for index in range(sample_size)
            ]
        column = [ordered_values[index] for index in indices]
        generator.shuffle(column)
        sampled_columns[name] = column

    configurations: list[dict] = []
    seen: set[tuple[int, ...]] = set()
    for row in range(sample_size):
        config = {
            "n_preds": args.n_preds,
            **{name: sampled_columns[name][row] for name, _ in dimensions},
        }
        key = tuple(config[name] for name, _ in dimensions)
        if key not in seen:
            configurations.append(config)
            seen.add(key)

    while len(configurations) < sample_size:
        config = {
            "n_preds": args.n_preds,
            **{name: generator.choice(values) for name, values in dimensions},
        }
        key = tuple(config[name] for name, _ in dimensions)
        if key not in seen:
            configurations.append(config)
            seen.add(key)

    return configurations


def fixed_parameters() -> dict:
    return {
        "n_preys": 100,
        "herd-radius": 4,
        "herd-chance": 30,
        "min-habitat-loss-radius": 4,
        "max-habitat-loss-radius": 10,
        "fixed-habitat-loss-percent": 0,
        "predator-reintroduction?": False,
        "reintroduced-predators": 2,
    }


def slope_per_tick(ticks: list[int], values: list[float]) -> float:
    """Return the least-squares population slope per simulation tick."""
    if len(values) < 2 or len(ticks) != len(values):
        return math.nan
    x_mean = sum(ticks) / len(ticks)
    y_mean = sum(values) / len(values)
    denominator = sum((tick - x_mean) ** 2 for tick in ticks)
    if denominator == 0:
        return 0.0
    return sum(
        (tick - x_mean) * (value - y_mean)
        for tick, value in zip(ticks, values)
    ) / denominator


def slope_per_record(values: list[float]) -> float:
    """Legacy index-based slope retained for downstream file compatibility."""
    return slope_per_tick(list(range(len(values))), values)


def coefficient_of_variation(values: list[float]) -> float:
    if len(values) < 2:
        return math.nan
    # local import to avoid heavy module import at top-level
    import pandas as pd

    series = pd.Series(values, dtype="float64")
    mean = series.mean()
    return float(series.std(ddof=1) / mean) if mean != 0 else math.nan


def relative_window_trend(
    ticks: list[int], values: list[float], window: int
) -> float:
    """Absolute fitted change across a window relative to mean population."""
    if len(values) < 2:
        return math.nan
    mean = sum(values) / len(values)
    if mean == 0:
        return math.nan
    return abs(slope_per_tick(ticks, values)) * window / mean


def count_meaningful_cycles(
    values: list[float], minimum_relative_amplitude: float
) -> int:
    """Count repeated, non-trivial cycles after three-point smoothing.

    A complete cycle requires three alternating meaningful extrema, such as
    peak-trough-peak. Reversals smaller than the configured fraction of the
    series mean are ignored so small stochastic fluctuations are not counted
    as population cycles.
    """
    if len(values) < 5:
        return 0
    mean = sum(values) / len(values)
    if mean <= 0:
        return 0
    threshold = mean * minimum_relative_amplitude
    smoothed = [
        (values[index - 1] + values[index] + values[index + 1]) / 3
        for index in range(1, len(values) - 1)
    ]
    candidates: list[tuple[str, float]] = []
    for index in range(1, len(smoothed) - 1):
        previous_value = smoothed[index - 1]
        value = smoothed[index]
        next_value = smoothed[index + 1]
        if value > previous_value and value >= next_value:
            candidates.append(("peak", value))
        elif value < previous_value and value <= next_value:
            candidates.append(("trough", value))

    meaningful: list[tuple[str, float]] = []
    for kind, value in candidates:
        if not meaningful:
            meaningful.append((kind, value))
            continue
        previous_kind, previous_value = meaningful[-1]
        if kind == previous_kind:
            more_extreme = (
                kind == "peak" and value > previous_value
            ) or (
                kind == "trough" and value < previous_value
            )
            if more_extreme:
                meaningful[-1] = (kind, value)
            continue
        if abs(value - previous_value) >= threshold:
            meaningful.append((kind, value))

    return max(0, (len(meaningful) - 1) // 2)


def oscillation_metrics(
    rows: list[dict],
    max_ticks: int,
    stability_window: int,
    minimum_relative_amplitude: float,
    minimum_cycles: int,
) -> dict:
    stability_start = max(0, max_ticks - stability_window)
    late_rows = [row for row in rows if row["tick"] >= stability_start]
    prey_cycles = count_meaningful_cycles(
        [float(row["preys"]) for row in late_rows], minimum_relative_amplitude
    )
    predator_cycles = count_meaningful_cycles(
        [float(row["predators"]) for row in late_rows],
        minimum_relative_amplitude,
    )
    return {
        "prey_meaningful_cycles_stability_window": prey_cycles,
        "predator_meaningful_cycles_stability_window": predator_cycles,
        "oscillatory_dynamics": (
            prey_cycles >= minimum_cycles and predator_cycles >= minimum_cycles
        ),
    }


def run_one(
    netlogo: pyNetLogo.NetLogoLink,
    config_id: int,
    params: dict,
    repetition: int,
    seed: int,
    max_ticks: int,
    record_every: int,
    args: argparse.Namespace,
) -> tuple[dict, list[dict]]:
    all_params = {**fixed_parameters(), **params}
    for name, value in all_params.items():
        set_variable(netlogo, name, value)

    netlogo.command(f"random-seed {seed}")
    netlogo.command("setup")

    rows: list[dict] = []

    def collect() -> None:
        rows.append(
            {
                "config_id": config_id,
                "repetition": repetition,
                "seed": seed,
                "tick": int(netlogo.report("ticks")),
                "preys": int(netlogo.report("count preys")),
                "predators": int(netlogo.report("count predators")),
                "mean_prey_hunger": report_number(
                    netlogo, "ifelse-value any? preys [mean [prey-hungry] of preys] [0]"
                ),
                "mean_predator_hunger": report_number(
                    netlogo,
                    "ifelse-value any? predators [mean [pred-hungry] of predators] [0]",
                ),
                "prey_deaths_predation": report_number(
                    netlogo, "prey-deaths-predation"
                ),
                "prey_deaths_starvation": report_number(
                    netlogo, "prey-deaths-starvation"
                ),
                "prey_deaths_crowding": report_number(
                    netlogo, "prey-deaths-crowding"
                ),
                "actual_habitat_loss_percent": report_number(
                    netlogo, "habitat-loss-percent", fallback=0.0
                ),
                "habitat_loss_applied": int(
                    bool(netlogo.report("habitat-loss-applied?"))
                ),
                "trapped_preys": report_number(
                    netlogo,
                    "count preys with [[no-food?] of patch-here]",
                    fallback=0.0,
                ),
                **{name.replace("-", "_").replace("?", ""): value
                   for name, value in all_params.items()},
            }
        )

    collect()
    extinct_population = "none"

    while int(netlogo.report("ticks")) < max_ticks:
        netlogo.command("go")
        current_tick = int(netlogo.report("ticks"))

        prey_alive = bool(netlogo.report("any? preys"))
        predator_alive = bool(netlogo.report("any? predators"))

        if current_tick % record_every == 0:
            collect()

        # Required because NetLogo's stop exits one `go` call, not this Python loop.
        if not prey_alive or not predator_alive:
            extinct_population = "preys" if not prey_alive else "predators"
            if not rows or rows[-1]["tick"] != current_tick:
                collect()
            break

    final_tick = int(netlogo.report("ticks"))
    final_preys = int(netlogo.report("count preys"))
    final_predators = int(netlogo.report("count predators"))
    if not rows or rows[-1]["tick"] != final_tick:
        collect()

    stability_start = max(0, max_ticks - args.stability_window)
    late_rows = [r for r in rows if r["tick"] >= stability_start]
    late_ticks = [r["tick"] for r in late_rows]
    late_preys = [r["preys"] for r in late_rows]
    late_predators = [r["predators"] for r in late_rows]

    survived = (
        final_tick >= max_ticks and final_preys > 0 and final_predators > 0
    )
    plausible_final_sizes = (
        args.min_final_preys <= final_preys <= args.max_final_preys
        and args.min_final_predator_ratio * params["n_preds"]
        <= final_predators
        <= args.max_final_predator_ratio * params["n_preds"]
    )

    prey_cv = coefficient_of_variation(late_preys)
    predator_cv = coefficient_of_variation(late_predators)
    prey_slope = slope_per_tick(late_ticks, late_preys)
    predator_slope = slope_per_tick(late_ticks, late_predators)
    prey_relative_trend = relative_window_trend(
        late_ticks, late_preys, args.stability_window
    )
    predator_relative_trend = relative_window_trend(
        late_ticks, late_predators, args.stability_window
    )

    expected_late_records = args.stability_window // record_every + 1
    enough_stability_data = len(late_rows) >= max(3, expected_late_records - 1)
    stable_dynamics = (
        survived
        and enough_stability_data
        and math.isfinite(prey_cv)
        and math.isfinite(predator_cv)
        and math.isfinite(prey_relative_trend)
        and math.isfinite(predator_relative_trend)
        and prey_cv <= args.max_prey_cv
        and predator_cv <= args.max_predator_cv
        and prey_relative_trend <= args.max_relative_window_trend
        and predator_relative_trend <= args.max_relative_window_trend
    )

    summary = {
        "config_id": config_id,
        "repetition": repetition,
        "seed": seed,
        **{name.replace("-", "_"): value for name, value in params.items()},
        "final_tick": final_tick,
        "final_preys": final_preys,
        "final_predators": final_predators,
        "survived_to_max_ticks": survived,
        "extinct_population": extinct_population,
        "stability_window_ticks": args.stability_window,
        "stability_records": len(late_rows),
        "prey_cv_stability_window": prey_cv,
        "predator_cv_stability_window": predator_cv,
        "prey_slope_per_tick_stability_window": prey_slope,
        "predator_slope_per_tick_stability_window": predator_slope,
        "prey_relative_trend_stability_window": prey_relative_trend,
        "predator_relative_trend_stability_window": predator_relative_trend,
        # Retained so existing analysis code can still read the old columns.
        "prey_cv_last_500": prey_cv,
        "predator_cv_last_500": predator_cv,
        "prey_slope_per_record_last_500": slope_per_record(late_preys),
        "predator_slope_per_record_last_500": slope_per_record(late_predators),
        "stable_dynamics": stable_dynamics,
        "plausible_final_sizes": plausible_final_sizes,
        "accepted_run": survived and plausible_final_sizes and stable_dynamics,
        "predator_reintroduction_count": report_number(
            netlogo, "predator-reintroduction-count", fallback=-1
        ),
        "actual_habitat_loss_percent": report_number(
            netlogo, "habitat-loss-percent", fallback=0.0
        ),
        "prey_extinct": int(final_preys == 0),
        "predator_extinct": int(final_predators == 0),
        "extinction_tick": report_number(
            netlogo,
            "ifelse-value extinction-tick = nobody [-1] [extinction-tick]",
            fallback=-1,
        ),
        "prey_deaths_predation": report_number(
            netlogo, "prey-deaths-predation", fallback=0.0
        ),
        "prey_deaths_starvation": report_number(
            netlogo, "prey-deaths-starvation", fallback=0.0
        ),
        "prey_deaths_crowding": report_number(
            netlogo, "prey-deaths-crowding", fallback=0.0
        ),
        "trapped_preys_final": report_number(
            netlogo,
            "count preys with [[no-food?] of patch-here]",
            fallback=0.0,
        ),
    }
    return summary, rows


def aggregate_configurations(run_summary: pd.DataFrame) -> pd.DataFrame:
    parameter_columns = [
        "n_preds",
        "prey_reproduction_chance",
        "predator_reproduction_chance",
        "max_prey_hungry",
        "max_pred_hungry",
        "predator_vision",
        "prey_vision",
        "crowding_hunger_cost",
        "max_local_prey_density",
    ]
    # Import pandas locally to avoid top-level import side-effects
    import pandas as pd

    return (
        run_summary.groupby(["config_id", *parameter_columns], as_index=False)
        .agg(
            n_runs=("repetition", "size"),
            persistence_rate=("survived_to_max_ticks", "mean"),
            stability_rate=("stable_dynamics", "mean"),
            acceptance_rate=("accepted_run", "mean"),
            median_final_preys=("final_preys", "median"),
            median_final_predators=("final_predators", "median"),
            median_prey_cv=("prey_cv_stability_window", "median"),
            median_predator_cv=("predator_cv_stability_window", "median"),
            median_prey_relative_trend=(
                "prey_relative_trend_stability_window", "median"
            ),
            median_predator_relative_trend=(
                "predator_relative_trend_stability_window", "median"
            ),
        )
        .sort_values(
            [
                "acceptance_rate",
                "persistence_rate",
                "stability_rate",
                "median_predator_cv",
                "median_predator_relative_trend",
            ],
            ascending=[False, False, False, True, True],
        )
    )


def save_outputs(
    pd,
    run_rows: list[dict],
    time_rows: list[dict],
    output_dir: Path,
) -> None:
    """Atomically checkpoint completed runs so interruptions do not lose work."""
    if not run_rows:
        return

    run_summary = pd.DataFrame(run_rows)
    time_series = pd.DataFrame(time_rows)
    config_ranking = aggregate_configurations(run_summary)

    outputs = {
        "baseline_run_summary.csv": run_summary,
        "baseline_timeseries.csv": time_series,
        "baseline_configuration_ranking.csv": config_ranking,
    }
    for filename, frame in outputs.items():
        destination = output_dir / filename
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        frame.to_csv(temporary, index=False)
        temporary.replace(destination)


def format_duration(seconds: float) -> str:
    seconds = max(0, round(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours:d}h {minutes:02d}m"
    if minutes:
        return f"{minutes:d}m {seconds:02d}s"
    return f"{seconds:d}s"


PARAMETER_OUTPUT_TO_NETLOGO = {
    "n_preds": "n_preds",
    "prey_reproduction_chance": "prey-reproduction-chance",
    "predator_reproduction_chance": "predator-reproduction-chance",
    "max_prey_hungry": "max-prey-hungry",
    "max_pred_hungry": "max-pred-hungry",
    "predator_vision": "predator-vision",
    "prey_vision": "prey-vision",
    "crowding_hunger_cost": "crowding-hunger-cost",
    "max_local_prey_density": "max-local-prey-density",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json_write(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def deterministic_run_seed(
    master_seed: int,
    stage: str,
    params: dict,
    repetition: int,
) -> int:
    """Derive a run seed independent of execution order and resume state."""
    ordered = [
        str(params["n_preds"]),
        *(str(params[name]) for name, _ in parameter_dimensions_from_params(params)),
    ]
    payload = "|".join([str(master_seed), stage, *ordered, str(repetition)])
    digest = hashlib.sha256(payload.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % 2_147_483_646 + 1


def parameter_dimensions_from_params(params: dict) -> list[tuple[str, int]]:
    """Return parameter values in a fixed order for seed derivation."""
    return [
        ("prey-reproduction-chance", params["prey-reproduction-chance"]),
        ("predator-reproduction-chance", params["predator-reproduction-chance"]),
        ("max-prey-hungry", params["max-prey-hungry"]),
        ("max-pred-hungry", params["max-pred-hungry"]),
        ("predator-vision", params["predator-vision"]),
        ("prey-vision", params["prey-vision"]),
        ("crowding-hunger-cost", params["crowding-hunger-cost"]),
        ("max-local-prey-density", params["max-local-prey-density"]),
    ]


def params_from_output_row(row) -> dict:
    return {
        netlogo_name: int(row[output_name])
        for output_name, netlogo_name in PARAMETER_OUTPUT_TO_NETLOGO.items()
    }


def load_candidate_configurations(pd, path: Path, args: argparse.Namespace):
    """Load and validate an explicit, reproducible candidate parameter grid."""
    if not path.is_file():
        raise FileNotFoundError(f"Candidate configuration CSV not found: {path}")
    source = pd.read_csv(path)
    required_columns = ["config_id", *PARAMETER_OUTPUT_TO_NETLOGO.keys()]
    missing = [column for column in required_columns if column not in source.columns]
    if missing:
        raise ValueError(
            "Candidate CSV is missing required columns: " + ", ".join(missing)
        )
    if source.empty:
        raise ValueError("Candidate CSV contains no configurations.")
    if source[required_columns].isna().any().any():
        raise ValueError("Candidate CSV contains missing required values.")

    normalized = source[required_columns].copy()
    for column in required_columns:
        numeric = pd.to_numeric(normalized[column], errors="raise")
        if not ((numeric % 1) == 0).all():
            raise ValueError(f"Candidate column {column!r} must contain integers.")
        normalized[column] = numeric.astype("int64")

    if normalized["config_id"].duplicated().any():
        duplicates = sorted(
            normalized.loc[normalized["config_id"].duplicated(), "config_id"]
            .unique()
            .tolist()
        )
        raise ValueError(f"Duplicate candidate configuration IDs: {duplicates}")
    if (normalized["config_id"] < 1).any():
        raise ValueError("Candidate configuration IDs must be positive integers.")
    if (normalized["n_preds"] < 1).any():
        raise ValueError("Candidate n_preds values must be at least 1.")

    allowed_values = {
        name.replace("-", "_"): set(values)
        for name, values in parameter_dimensions(args)
    }
    for column, allowed in allowed_values.items():
        invalid = sorted(set(normalized[column]) - allowed)
        if invalid:
            raise ValueError(
                f"Candidate column {column!r} contains values outside the "
                f"configured methodological range: {invalid}"
            )

    normalized = normalized.sort_values("config_id").reset_index(drop=True)
    configurations = [
        (int(row["config_id"]), params_from_output_row(row))
        for _, row in normalized.iterrows()
    ]
    return configurations, normalized


def configuration_table(pd, configurations: list[tuple[int, dict]]):
    return pd.DataFrame(
        [
            {
                "config_id": config_id,
                **{
                    name.replace("-", "_"): value
                    for name, value in params.items()
                },
            }
            for config_id, params in configurations
        ]
    )


def validate_or_write_parameter_grid(
    pd,
    configurations: list[tuple[int, dict]],
    destination: Path,
) -> None:
    expected = configuration_table(pd, configurations).sort_values("config_id")
    if destination.exists():
        existing = pd.read_csv(destination).sort_values("config_id")
        if list(existing.columns) != list(expected.columns) or not existing.reset_index(
            drop=True
        ).equals(expected.reset_index(drop=True)):
            raise ValueError(
                f"Existing parameter grid does not match this command: {destination}. "
                "Use a new --output-dir to preserve reproducibility."
            )
        return
    expected.to_csv(destination, index=False)


def load_stage_checkpoint(pd, stage_dir: Path) -> tuple[list[dict], list[dict]]:
    summary_path = stage_dir / "baseline_run_summary.csv"
    timeseries_path = stage_dir / "baseline_timeseries.csv"
    if not summary_path.exists() and not timeseries_path.exists():
        return [], []
    if not summary_path.exists() or not timeseries_path.exists():
        raise ValueError(
            f"Incomplete checkpoint in {stage_dir}. Both the run summary and "
            "time series are required."
        )
    summary = pd.read_csv(summary_path)
    timeseries = pd.read_csv(timeseries_path)
    if summary.duplicated(["config_id", "repetition"]).any():
        raise ValueError(f"Duplicate completed runs found in {summary_path}.")
    return summary.to_dict("records"), timeseries.to_dict("records")


def run_calibration_stage(
    pd,
    netlogo,
    stage_name: str,
    configurations: list[tuple[int, dict]],
    repetitions: int,
    max_ticks: int,
    args: argparse.Namespace,
    stage_dir: Path,
    stability_window: int | None = None,
    record_every: int | None = None,
):
    """Run or resume one calibration stage and return its three data frames."""
    stage_dir.mkdir(parents=True, exist_ok=True)
    validate_or_write_parameter_grid(
        pd, configurations, stage_dir / "baseline_parameter_grid.csv"
    )

    stage_args = copy.copy(args)
    requested_window = (
        args.stability_window if stability_window is None else stability_window
    )
    stage_args.stability_window = min(requested_window, max_ticks)
    stage_record_every = args.record_every if record_every is None else record_every
    pd.DataFrame(
        [
            {
                "stage": stage_name,
                "master_seed": args.master_seed,
                "configurations": len(configurations),
                "repetitions": repetitions,
                "max_ticks": max_ticks,
                "record_every": stage_record_every,
                "stability_window": stage_args.stability_window,
                "max_prey_cv": args.max_prey_cv,
                "max_predator_cv": args.max_predator_cv,
                "max_relative_window_trend": args.max_relative_window_trend,
            }
        ]
    ).to_csv(stage_dir / "baseline_run_settings.csv", index=False)

    run_rows, time_rows = load_stage_checkpoint(pd, stage_dir)
    completed_keys = {
        (int(row["config_id"]), int(row["repetition"])) for row in run_rows
    }
    expected_ids = {config_id for config_id, _ in configurations}
    unexpected_ids = {key[0] for key in completed_keys} - expected_ids
    if unexpected_ids:
        raise ValueError(
            f"Checkpoint for {stage_name} contains unexpected configuration IDs: "
            f"{sorted(unexpected_ids)}"
        )

    total = len(configurations) * repetitions
    completed = len(completed_keys)
    newly_completed = 0
    started_at = time.monotonic()
    print(
        f"\n{stage_name.upper()}: {len(configurations)} configurations x "
        f"{repetitions} repetitions = {total} runs; max ticks={max_ticks}"
    )
    if completed:
        print(f"Resuming from {completed}/{total} completed runs.")

    try:
        for config_id, params in configurations:
            for repetition in range(repetitions):
                key = (config_id, repetition)
                if key in completed_keys:
                    continue

                seed = deterministic_run_seed(
                    args.master_seed, stage_name, params, repetition
                )
                summary, series = run_one(
                    netlogo,
                    config_id,
                    params,
                    repetition,
                    seed,
                    max_ticks,
                    stage_record_every,
                    stage_args,
                )
                summary.update(
                    oscillation_metrics(
                        series,
                        max_ticks,
                        stage_args.stability_window,
                        args.min_cycle_relative_amplitude,
                        args.min_validation_cycles,
                    )
                )
                bounded_trend = (
                    bool(summary["survived_to_max_ticks"])
                    and finite_and_at_most(
                        summary["prey_relative_trend_stability_window"],
                        args.max_relative_window_trend,
                    )
                    and finite_and_at_most(
                        summary["predator_relative_trend_stability_window"],
                        args.max_relative_window_trend,
                    )
                )
                no_reintroduction = (
                    float(summary["predator_reintroduction_count"]) == 0.0
                )
                summary["bounded_trend_dynamics"] = bounded_trend
                summary["no_predator_reintroduction"] = no_reintroduction
                summary["bounded_oscillation"] = (
                    bounded_trend
                    and bool(summary["oscillatory_dynamics"])
                    and no_reintroduction
                )
                summary["calibration_stage"] = stage_name
                for record in series:
                    record["calibration_stage"] = stage_name
                run_rows.append(summary)
                time_rows.extend(series)
                completed_keys.add(key)
                completed += 1
                newly_completed += 1

                if newly_completed % args.checkpoint_every == 0:
                    save_outputs(pd, run_rows, time_rows, stage_dir)

                elapsed = time.monotonic() - started_at
                estimated_remaining = (
                    elapsed / newly_completed * (total - completed)
                    if newly_completed
                    else 0
                )
                print(
                    f"[{stage_name} {completed}/{total}] config={config_id} "
                    f"rep={repetition} tick={summary['final_tick']} "
                    f"prey={summary['final_preys']} "
                    f"predators={summary['final_predators']} "
                    f"stable={summary['stable_dynamics']} "
                    f"ETA={format_duration(estimated_remaining)}"
                )
    finally:
        save_outputs(pd, run_rows, time_rows, stage_dir)

    summary = pd.DataFrame(run_rows).sort_values(["config_id", "repetition"])
    timeseries = pd.DataFrame(time_rows).sort_values(
        ["config_id", "repetition", "tick"]
    )
    ranking = aggregate_configurations(summary)
    return summary, timeseries, ranking


def finite_and_at_most(value, maximum: float) -> bool:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(numeric) and numeric <= maximum


def make_screen_decisions(pd, summary, args: argparse.Namespace):
    decisions = summary.copy()
    eligible_values: list[bool] = []
    reasons: list[str] = []
    max_trends: list[float] = []
    max_cvs: list[float] = []

    for _, row in decisions.iterrows():
        failures: list[str] = []
        if not bool(row["survived_to_max_ticks"]):
            failures.append("population extinction before screen endpoint")
        if int(row["final_preys"]) < args.screen_min_final_preys:
            failures.append("final prey below screen minimum")
        if int(row["final_predators"]) < args.screen_min_final_predators:
            failures.append("final predators below screen minimum")

        prey_cv_ok = finite_and_at_most(
            row["prey_cv_stability_window"], args.screen_max_cv
        )
        predator_cv_ok = finite_and_at_most(
            row["predator_cv_stability_window"], args.screen_max_cv
        )
        prey_trend_ok = finite_and_at_most(
            row["prey_relative_trend_stability_window"],
            args.screen_max_relative_trend,
        )
        predator_trend_ok = finite_and_at_most(
            row["predator_relative_trend_stability_window"],
            args.screen_max_relative_trend,
        )
        if not prey_cv_ok:
            failures.append("prey CV above screen maximum or undefined")
        if not predator_cv_ok:
            failures.append("predator CV above screen maximum or undefined")
        if not prey_trend_ok:
            failures.append("prey relative trend above screen maximum or undefined")
        if not predator_trend_ok:
            failures.append(
                "predator relative trend above screen maximum or undefined"
            )

        trend_values = [
            float(row["prey_relative_trend_stability_window"]),
            float(row["predator_relative_trend_stability_window"]),
        ]
        cv_values = [
            float(row["prey_cv_stability_window"]),
            float(row["predator_cv_stability_window"]),
        ]
        max_trends.append(
            max(trend_values) if all(math.isfinite(v) for v in trend_values) else math.inf
        )
        max_cvs.append(
            max(cv_values) if all(math.isfinite(v) for v in cv_values) else math.inf
        )
        eligible_values.append(not failures)
        reasons.append("eligible" if not failures else "; ".join(failures))

    decisions["screen_max_observed_relative_trend"] = max_trends
    decisions["screen_max_observed_cv"] = max_cvs
    decisions["screen_eligible"] = eligible_values
    decisions["screen_decision_reason"] = reasons
    decisions["selected_for_validation"] = False

    eligible_indices = (
        decisions.loc[decisions["screen_eligible"]]
        .sort_values(
            [
                "screen_max_observed_relative_trend",
                "screen_max_observed_cv",
                "config_id",
            ]
        )
        .head(args.max_validation_configurations)
        .index
    )
    decisions.loc[eligible_indices, "selected_for_validation"] = True
    return decisions.sort_values("config_id")


def make_validation_decisions(
    pd,
    ranking,
    run_summary,
    args: argparse.Namespace,
    require_legacy_stability: bool = True,
):
    oscillation_summary = (
        run_summary.groupby("config_id", as_index=False)
        .agg(
            oscillation_rate=("oscillatory_dynamics", "mean"),
            bounded_trend_rate=("bounded_trend_dynamics", "mean"),
            bounded_oscillation_rate=("bounded_oscillation", "mean"),
            no_reintroduction_rate=("no_predator_reintroduction", "mean"),
            median_prey_meaningful_cycles=(
                "prey_meaningful_cycles_stability_window",
                "median",
            ),
            median_predator_meaningful_cycles=(
                "predator_meaningful_cycles_stability_window",
                "median",
            ),
        )
    )
    decisions = ranking.merge(oscillation_summary, on="config_id", how="left")
    selected: list[bool] = []
    reasons: list[str] = []
    for _, row in decisions.iterrows():
        failures: list[str] = []
        if float(row["persistence_rate"]) < args.min_validation_persistence_rate:
            failures.append("persistence rate below threshold")
        if (
            require_legacy_stability
            and float(row["stability_rate"])
            < args.min_validation_stability_rate
        ):
            failures.append("stability rate below threshold")
        if float(row["acceptance_rate"]) < args.min_validation_acceptance_rate:
            failures.append("strict acceptance rate below threshold")
        if (
            float(row["bounded_oscillation_rate"])
            < args.min_validation_oscillation_rate
        ):
            failures.append("bounded-oscillation rate below threshold")
        selected.append(not failures)
        reasons.append("selected" if not failures else "; ".join(failures))
    decisions["selected_configuration"] = selected
    decisions["validation_decision_reason"] = reasons
    return decisions


def command_version(command: list[str]) -> str:
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        output = (result.stdout + "\n" + result.stderr).strip()
        return output
    except Exception as exc:
        return f"unavailable: {exc}"


def reproducibility_protocol(args: argparse.Namespace) -> dict:
    """Return every input that must remain fixed when resuming a run."""
    protocol = {
        "mode": (
            "habitat_loss_experiment"
            if args.habitat_experiment
            else "candidate_validation"
            if args.validate_candidates is not None
            else "automated_calibration"
            if args.automated_calibration
            else "single_stage"
        ),
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "model_sha256": sha256_file(args.model.resolve()),
        "master_seed": args.master_seed,
        "fixed_parameters": fixed_parameters(),
        "sampled_parameter_values": {
            name: list(values) for name, values in parameter_dimensions(args)
        },
        "initial_predators": args.n_preds,
        "record_every": args.record_every,
        "stability_window": args.stability_window,
        "run_acceptance_thresholds": {
            "min_final_preys": args.min_final_preys,
            "max_final_preys": args.max_final_preys,
            "min_final_predator_ratio": args.min_final_predator_ratio,
            "max_final_predator_ratio": args.max_final_predator_ratio,
            "max_prey_cv": args.max_prey_cv,
            "max_predator_cv": args.max_predator_cv,
            "max_relative_window_trend": args.max_relative_window_trend,
        },
        "screen": {
            "configurations": args.screen_configurations,
            "repetitions": 1,
            "max_ticks": args.screen_ticks,
            "minimum_final_preys": args.screen_min_final_preys,
            "minimum_final_predators": args.screen_min_final_predators,
            "maximum_cv": args.screen_max_cv,
            "maximum_relative_trend": args.screen_max_relative_trend,
            "maximum_validation_candidates": args.max_validation_configurations,
        },
        "validation": {
            "repetitions": args.validation_repetitions,
            "max_ticks": args.validation_ticks,
            "record_every": args.validation_record_every,
            "stability_window": args.validation_stability_window,
            "minimum_persistence_rate": args.min_validation_persistence_rate,
            "minimum_stability_rate": args.min_validation_stability_rate,
            "minimum_oscillation_rate": args.min_validation_oscillation_rate,
            "minimum_cycles_per_population": args.min_validation_cycles,
            "minimum_cycle_relative_amplitude": (
                args.min_cycle_relative_amplitude
            ),
            "minimum_strict_acceptance_rate": (
                args.min_validation_acceptance_rate
            ),
        },
    }
    if args.habitat_experiment:
        protocol["habitat_loss_experiment"] = {
            "configuration_id": FINAL_CONFIGURATION_ID,
            "configuration": FINAL_CONFIGURATION,
            "habitat_loss_values": list(args.habitat_loss_values),
            "repetitions": args.repetitions,
            "max_ticks": args.max_ticks,
            "record_every": args.record_every,
            "treatment_tick": 2000,
            "predator_reintroduction": False,
            "paired_seeds_across_treatments": True,
        }
    if args.validate_candidates is not None:
        candidate_path = args.validate_candidates.resolve()
        protocol["candidate_file"] = {
            "path": str(candidate_path),
            "sha256": sha256_file(candidate_path),
        }
    return protocol


def manifest_payload(args, pd, jpype, pyNetLogo) -> dict:
    java_executable = Path(os.environ.get("JAVA_HOME", "")) / "bin" / "java"
    java_command = str(java_executable) if java_executable.is_file() else "java"
    protocol = reproducibility_protocol(args)
    protocol_sha256 = hashlib.sha256(
        json.dumps(protocol, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return {
        "schema_version": 1,
        "status": "running",
        "started_utc": utc_now(),
        "command": [sys.executable, *sys.argv],
        "arguments": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
        },
        "files": {
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "model": str(args.model.resolve()),
            "model_sha256": sha256_file(args.model.resolve()),
        },
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "pandas_version": getattr(pd, "__version__", "unknown"),
            "jpype_version": getattr(jpype, "__version__", "unknown"),
            "pynetlogo_version": getattr(pyNetLogo, "__version__", "unknown"),
            "java_version": command_version([java_command, "-version"]),
            "netlogo_home": str(args.netlogo_home.resolve()),
        },
        "reproducibility_protocol": protocol,
        "reproducibility_protocol_sha256": protocol_sha256,
    }


def save_habitat_outputs(
    pd,
    run_rows: list[dict],
    time_rows: list[dict],
    output_dir: Path,
) -> None:
    """Atomically save final-experiment checkpoints and treatment summaries."""
    if not run_rows:
        return

    run_summary = (
        pd.DataFrame(run_rows)
        .sort_values(["fixed_habitat_loss_percent", "repetition"])
        .reset_index(drop=True)
    )
    time_series = (
        pd.DataFrame(time_rows)
        .sort_values(["fixed_habitat_loss_percent", "repetition", "tick"])
        .reset_index(drop=True)
    )
    treatment_summary = (
        run_summary.groupby("fixed_habitat_loss_percent", as_index=False)
        .agg(
            completed_runs=("repetition", "size"),
            prey_extinction_rate=("prey_extinct", "mean"),
            predator_extinction_rate=("predator_extinct", "mean"),
            mean_actual_habitat_loss=("actual_habitat_loss_percent", "mean"),
            median_final_preys=("final_preys", "median"),
            median_final_predators=("final_predators", "median"),
            median_final_tick=("final_tick", "median"),
        )
        .sort_values("fixed_habitat_loss_percent")
    )

    outputs = {
        "habitat_run_summary.csv": run_summary,
        "habitat_timeseries.csv": time_series,
        "habitat_treatment_summary.csv": treatment_summary,
    }
    for filename, frame in outputs.items():
        destination = output_dir / filename
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        frame.to_csv(temporary, index=False)
        temporary.replace(destination)


def run_habitat_experiment(args, pd, jpype, pyNetLogo, jvm_path: Path) -> None:
    """Run the selected baseline configuration across habitat-loss treatments."""
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "habitat_experiment_manifest.json"
    manifest = manifest_payload(args, pd, jpype, pyNetLogo)

    summary_path = output_dir / "habitat_run_summary.csv"
    series_path = output_dir / "habitat_timeseries.csv"
    if manifest_path.exists():
        existing_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing_manifest.get("reproducibility_protocol_sha256") != manifest.get(
            "reproducibility_protocol_sha256"
        ):
            raise ValueError(
                "The existing output directory uses a different model, seed, "
                "configuration, treatment grid, or run protocol. Use a new "
                "--output-dir instead of mixing experiments."
            )
        manifest["started_utc"] = existing_manifest.get(
            "started_utc", manifest["started_utc"]
        )
        resume_events = list(existing_manifest.get("resume_events_utc", []))
        resume_events.append(utc_now())
        manifest["resume_events_utc"] = resume_events
    elif summary_path.exists() or series_path.exists():
        raise ValueError(
            "Habitat-loss checkpoints exist without a matching manifest. "
            "Use a new --output-dir."
        )
    atomic_json_write(manifest_path, manifest)

    configuration_row = {
        "config_id": FINAL_CONFIGURATION_ID,
        **{name.replace("-", "_"): value for name, value in FINAL_CONFIGURATION.items()},
    }
    configuration_path = output_dir / "habitat_configuration.csv"
    expected_configuration = pd.DataFrame([configuration_row])
    if configuration_path.exists():
        existing_configuration = pd.read_csv(configuration_path)
        if not existing_configuration.equals(expected_configuration):
            raise ValueError(
                "Saved habitat configuration does not match configuration 275. "
                "Use a new --output-dir."
            )
    else:
        expected_configuration.to_csv(configuration_path, index=False)

    paired_seeds = {
        repetition: deterministic_run_seed(
            args.master_seed,
            "habitat_loss_config_275",
            FINAL_CONFIGURATION,
            repetition,
        )
        for repetition in range(args.repetitions)
    }
    seed_table = pd.DataFrame(
        [
            {"repetition": repetition, "seed": seed}
            for repetition, seed in paired_seeds.items()
        ]
    )
    seed_path = output_dir / "habitat_seed_table.csv"
    if seed_path.exists():
        if not pd.read_csv(seed_path).equals(seed_table):
            raise ValueError(
                "Saved seed table does not match this command. Use a new "
                "--output-dir."
            )
    else:
        seed_table.to_csv(seed_path, index=False)

    if summary_path.exists() != series_path.exists():
        raise ValueError(
            "Only one habitat checkpoint file is present. Use a new output "
            "directory or restore the matching checkpoint."
        )
    if summary_path.exists():
        run_frame = pd.read_csv(summary_path)
        time_frame = pd.read_csv(series_path)
        run_rows = run_frame.to_dict("records")
        time_rows = time_frame.to_dict("records")
    else:
        run_rows = []
        time_rows = []

    completed_keys = {
        (int(row["fixed_habitat_loss_percent"]), int(row["repetition"]))
        for row in run_rows
    }
    expected_keys = {
        (loss, repetition)
        for loss in args.habitat_loss_values
        for repetition in range(args.repetitions)
    }
    unexpected_keys = completed_keys - expected_keys
    if unexpected_keys:
        raise ValueError(
            "Existing habitat checkpoints contain runs outside this design: "
            f"{sorted(unexpected_keys)[:10]}"
        )

    total = len(expected_keys)
    completed = len(completed_keys)
    newly_completed = 0
    started_at = time.monotonic()
    netlogo = None
    interrupted = False
    run_error: Exception | None = None

    print(
        f"Habitat-loss experiment: configuration {FINAL_CONFIGURATION_ID}; "
        f"{len(args.habitat_loss_values)} treatments x {args.repetitions} "
        f"paired repetitions = {total} runs"
    )
    print("NetLogo treatment tick: 2000; experiment end: tick 3000")
    if completed:
        print(f"Resuming from {completed}/{total} completed runs.")

    try:
        netlogo = pyNetLogo.NetLogoLink(
            gui=False,
            netlogo_home=str(args.netlogo_home),
            jvm_path=str(jvm_path),
        )
        netlogo.load_model(str(args.model))

        for habitat_loss in args.habitat_loss_values:
            for repetition in range(args.repetitions):
                key = (habitat_loss, repetition)
                if key in completed_keys:
                    continue

                parameters = {
                    **FINAL_CONFIGURATION,
                    "fixed-habitat-loss-percent": habitat_loss,
                }
                summary, series = run_one(
                    netlogo,
                    FINAL_CONFIGURATION_ID,
                    parameters,
                    repetition,
                    paired_seeds[repetition],
                    args.max_ticks,
                    args.record_every,
                    args,
                )
                summary["experiment_stage"] = "final_habitat_loss"
                summary["treatment_tick"] = 2000
                for record in series:
                    record["experiment_stage"] = "final_habitat_loss"
                    record["treatment_tick"] = 2000

                run_rows.append(summary)
                time_rows.extend(series)
                completed_keys.add(key)
                completed += 1
                newly_completed += 1

                if newly_completed % args.checkpoint_every == 0:
                    save_habitat_outputs(pd, run_rows, time_rows, output_dir)

                elapsed = time.monotonic() - started_at
                estimated_remaining = (
                    elapsed / newly_completed * (total - completed)
                    if newly_completed
                    else 0
                )
                print(
                    f"[{completed}/{total}] loss={habitat_loss}% "
                    f"rep={repetition} seed={paired_seeds[repetition]} "
                    f"tick={summary['final_tick']} prey={summary['final_preys']} "
                    f"predators={summary['final_predators']} "
                    f"prey_extinct={bool(summary['prey_extinct'])} "
                    f"ETA={format_duration(estimated_remaining)}"
                )
    except KeyboardInterrupt:
        interrupted = True
        print("\nInterrupted by user; saving completed habitat-loss runs.")
    except Exception as exc:
        run_error = exc
    finally:
        save_habitat_outputs(pd, run_rows, time_rows, output_dir)
        if netlogo is not None:
            try:
                netlogo.kill_workspace()
            except Exception as exc:
                print(f"WARNING: NetLogo workspace did not close cleanly: {exc}")
                if run_error is None and not interrupted:
                    run_error = exc

    manifest["completed_runs"] = len(completed_keys)
    manifest["expected_runs"] = total
    if run_error is not None:
        manifest["status"] = "failed"
        manifest["failed_utc"] = utc_now()
        manifest["error"] = f"{type(run_error).__name__}: {run_error}"
    elif interrupted:
        manifest["status"] = "interrupted"
        manifest["interrupted_utc"] = utc_now()
    else:
        manifest["status"] = "completed"
        manifest["completed_utc"] = utc_now()
    atomic_json_write(manifest_path, manifest)

    if run_error is not None:
        raise run_error
    if interrupted:
        print(f"Saved partial habitat-loss outputs to {output_dir}")
        sys.exit(130)
    print(f"\nSaved completed habitat-loss experiment to {output_dir}")


def run_automated_calibration(args, pd, jpype, pyNetLogo, jvm_path: Path) -> None:
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "calibration_manifest.json"
    manifest = manifest_payload(args, pd, jpype, pyNetLogo)
    checkpoint_paths = [
        output_dir / "screen" / "baseline_run_summary.csv",
        output_dir / "validation" / "baseline_run_summary.csv",
    ]
    if manifest_path.exists():
        existing_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        previous_signature = existing_manifest.get(
            "reproducibility_protocol_sha256"
        )
        current_signature = manifest["reproducibility_protocol_sha256"]
        if previous_signature != current_signature:
            raise ValueError(
                "The existing output directory was created with a different "
                "model, script, parameter range, seed, or selection protocol. "
                "Use a new --output-dir instead of mixing calibration designs."
            )
        manifest["started_utc"] = existing_manifest.get(
            "started_utc", manifest["started_utc"]
        )
        resume_events = list(existing_manifest.get("resume_events_utc", []))
        resume_events.append(utc_now())
        manifest["resume_events_utc"] = resume_events
    elif any(path.exists() for path in checkpoint_paths):
        raise ValueError(
            "Calibration checkpoints exist without a reproducibility manifest. "
            "Use a new --output-dir so incompatible results are not combined."
        )
    atomic_json_write(manifest_path, manifest)

    sampled = diverse_parameter_sample(
        args, args.screen_configurations, args.master_seed
    )
    screen_configurations = list(enumerate(sampled, start=1))
    netlogo = None
    try:
        netlogo = pyNetLogo.NetLogoLink(
            gui=False,
            netlogo_home=str(args.netlogo_home),
            jvm_path=str(jvm_path),
        )
        netlogo.load_model(str(args.model))

        screen_summary, _, _ = run_calibration_stage(
            pd,
            netlogo,
            "screen",
            screen_configurations,
            repetitions=1,
            max_ticks=args.screen_ticks,
            args=args,
            stage_dir=output_dir / "screen",
        )
        screen_decisions = make_screen_decisions(pd, screen_summary, args)
        screen_decisions.to_csv(
            output_dir / "screen_configuration_decisions.csv", index=False
        )

        validation_rows = screen_decisions.loc[
            screen_decisions["selected_for_validation"]
        ]
        validation_configurations = [
            (int(row["config_id"]), params_from_output_row(row))
            for _, row in validation_rows.iterrows()
        ]
        manifest["screen_completed_utc"] = utc_now()
        manifest["screen_results"] = {
            "completed_runs": int(len(screen_summary)),
            "eligible_configurations": int(
                screen_decisions["screen_eligible"].sum()
            ),
            "sent_to_validation": len(validation_configurations),
        }
        atomic_json_write(manifest_path, manifest)

        if not validation_configurations:
            empty_columns = list(aggregate_configurations(screen_summary).columns) + [
                "selected_configuration",
                "validation_decision_reason",
            ]
            pd.DataFrame(columns=empty_columns).to_csv(
                output_dir / "validation_configuration_decisions.csv", index=False
            )
            pd.DataFrame(columns=empty_columns).to_csv(
                output_dir / "selected_configurations.csv", index=False
            )
            manifest["status"] = "completed_no_screen_candidates"
            manifest["completed_utc"] = utc_now()
            atomic_json_write(manifest_path, manifest)
            print(
                "\nCalibration completed, but no configurations passed the "
                "predefined screen. Review screen_configuration_decisions.csv."
            )
            return

        validation_summary, _, validation_ranking = run_calibration_stage(
            pd,
            netlogo,
            "validation",
            validation_configurations,
            repetitions=args.validation_repetitions,
            max_ticks=args.validation_ticks,
            args=args,
            stage_dir=output_dir / "validation",
            stability_window=args.validation_stability_window,
            record_every=args.validation_record_every,
        )
        validation_decisions = make_validation_decisions(
            pd, validation_ranking, validation_summary, args
        )
        validation_decisions.to_csv(
            output_dir / "validation_configuration_decisions.csv", index=False
        )
        selected = validation_decisions.loc[
            validation_decisions["selected_configuration"]
        ].copy()
        selected.to_csv(output_dir / "selected_configurations.csv", index=False)

        manifest["validation_completed_utc"] = utc_now()
        manifest["validation_results"] = {
            "validated_configurations": len(validation_configurations),
            "completed_runs": int(
                len(validation_configurations) * args.validation_repetitions
            ),
            "selected_configurations": int(len(selected)),
        }
        manifest["status"] = "completed"
        manifest["completed_utc"] = utc_now()
        atomic_json_write(manifest_path, manifest)

        print("\nAutomated calibration completed.")
        print(f"Selected configurations: {len(selected)}")
        print(f"Results: {output_dir}")
    except KeyboardInterrupt:
        manifest["status"] = "interrupted"
        manifest["interrupted_utc"] = utc_now()
        atomic_json_write(manifest_path, manifest)
        print("\nInterrupted. Completed runs were checkpointed and can be resumed.")
        raise
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["failed_utc"] = utc_now()
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        atomic_json_write(manifest_path, manifest)
        raise
    finally:
        if netlogo is not None:
            try:
                netlogo.kill_workspace()
            except Exception as exc:
                print(f"WARNING: NetLogo workspace did not close cleanly: {exc}")


def run_candidate_validation(args, pd, jpype, pyNetLogo, jvm_path: Path) -> None:
    """Validate an explicit candidate grid without repeating parameter search."""
    configurations, normalized_candidates = load_candidate_configurations(
        pd, args.validate_candidates, args
    )
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "calibration_manifest.json"
    manifest = manifest_payload(args, pd, jpype, pyNetLogo)
    checkpoint_path = output_dir / "validation" / "baseline_run_summary.csv"

    if manifest_path.exists():
        existing_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing_manifest.get("reproducibility_protocol_sha256") != manifest.get(
            "reproducibility_protocol_sha256"
        ):
            raise ValueError(
                "The existing output directory uses a different script, model, "
                "candidate grid, seed, or validation protocol. Use a new "
                "--output-dir instead of mixing designs."
            )
        manifest["started_utc"] = existing_manifest.get(
            "started_utc", manifest["started_utc"]
        )
        resume_events = list(existing_manifest.get("resume_events_utc", []))
        resume_events.append(utc_now())
        manifest["resume_events_utc"] = resume_events
    elif checkpoint_path.exists():
        raise ValueError(
            "Validation checkpoints exist without a reproducibility manifest. "
            "Use a new --output-dir so incompatible results are not combined."
        )
    atomic_json_write(manifest_path, manifest)

    candidate_copy = output_dir / "candidate_configurations.csv"
    if candidate_copy.exists():
        existing_candidates = pd.read_csv(candidate_copy)
        if not existing_candidates.equals(normalized_candidates):
            raise ValueError(
                "The saved candidate grid differs from the supplied candidate "
                "file. Use a new --output-dir."
            )
    else:
        normalized_candidates.to_csv(candidate_copy, index=False)

    netlogo = None
    try:
        netlogo = pyNetLogo.NetLogoLink(
            gui=False,
            netlogo_home=str(args.netlogo_home),
            jvm_path=str(jvm_path),
        )
        netlogo.load_model(str(args.model))
        validation_summary, _, validation_ranking = run_calibration_stage(
            pd,
            netlogo,
            "validation_v2",
            configurations,
            repetitions=args.validation_repetitions,
            max_ticks=args.validation_ticks,
            args=args,
            stage_dir=output_dir / "validation",
            stability_window=args.validation_stability_window,
            record_every=args.validation_record_every,
        )
        validation_decisions = make_validation_decisions(
            pd,
            validation_ranking,
            validation_summary,
            args,
            require_legacy_stability=False,
        )
        validation_decisions.to_csv(
            output_dir / "validation_configuration_decisions.csv", index=False
        )
        selected = validation_decisions.loc[
            validation_decisions["selected_configuration"]
        ].copy()
        selected.to_csv(output_dir / "selected_configurations.csv", index=False)

        manifest["validation_completed_utc"] = utc_now()
        manifest["validation_results"] = {
            "candidate_configurations": len(configurations),
            "completed_runs": int(
                len(configurations) * args.validation_repetitions
            ),
            "selected_configurations": int(len(selected)),
        }
        manifest["status"] = "completed"
        manifest["completed_utc"] = utc_now()
        atomic_json_write(manifest_path, manifest)

        print("\nCandidate validation completed.")
        print(f"Selected configurations: {len(selected)}")
        print(f"Results: {output_dir}")
    except KeyboardInterrupt:
        manifest["status"] = "interrupted"
        manifest["interrupted_utc"] = utc_now()
        atomic_json_write(manifest_path, manifest)
        print("\nInterrupted. Completed runs were checkpointed and can be resumed.")
        raise
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["failed_utc"] = utc_now()
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        atomic_json_write(manifest_path, manifest)
        raise
    finally:
        if netlogo is not None:
            try:
                netlogo.kill_workspace()
            except Exception as exc:
                print(f"WARNING: NetLogo workspace did not close cleanly: {exc}")


def main() -> None:
    args = parse_args()

    selected_modes = sum(
        [
            bool(args.quick),
            args.random_runs is not None,
            bool(args.automated_calibration),
            args.validate_candidates is not None,
            bool(args.habitat_experiment),
        ]
    )
    if selected_modes > 1:
        raise ValueError(
            "Use only one of --quick, --random-runs, "
            "--automated-calibration, --validate-candidates, or "
            "--habitat-experiment."
        )
    if args.repetitions < 1:
        raise ValueError("--repetitions must be at least 1.")
    if args.n_preds < 1:
        raise ValueError("--n-preds must be at least 1.")
    if args.max_ticks < 1:
        raise ValueError("--max-ticks must be at least 1.")
    if args.record_every < 1:
        raise ValueError("--record-every must be at least 1.")
    if args.checkpoint_every < 1:
        raise ValueError("--checkpoint-every must be at least 1.")
    if args.quick_configurations < 1:
        raise ValueError("--quick-configurations must be at least 1.")
    if args.random_runs is not None and args.random_runs < 1:
        raise ValueError("--random-runs must be at least 1.")
    if args.screen_configurations < 1:
        raise ValueError("--screen-configurations must be at least 1.")
    if args.screen_ticks < 1 or args.validation_ticks < 1:
        raise ValueError("Screen and validation ticks must be at least 1.")
    if args.validation_repetitions < 1:
        raise ValueError("--validation-repetitions must be at least 1.")
    if args.validation_stability_window < 1:
        raise ValueError("--validation-stability-window must be at least 1.")
    if args.validation_record_every < 1:
        raise ValueError("--validation-record-every must be at least 1.")
    if args.max_validation_configurations < 1:
        raise ValueError("--max-validation-configurations must be at least 1.")
    if args.screen_min_final_preys < 1 or args.screen_min_final_predators < 1:
        raise ValueError("Screen population minima must be at least 1.")
    if args.screen_max_cv < 0 or args.screen_max_relative_trend < 0:
        raise ValueError("Screen CV and trend thresholds cannot be negative.")
    for option_name in [
        "min_validation_persistence_rate",
        "min_validation_stability_rate",
        "min_validation_acceptance_rate",
        "min_validation_oscillation_rate",
    ]:
        value = getattr(args, option_name)
        if not 0 <= value <= 1:
            raise ValueError(f"--{option_name.replace('_', '-')} must be in [0, 1].")
    if args.min_validation_cycles < 1:
        raise ValueError("--min-validation-cycles must be at least 1.")
    if not 0 < args.min_cycle_relative_amplitude <= 1:
        raise ValueError("--min-cycle-relative-amplitude must be in (0, 1].")
    if any(value < 0 or value > 100 for value in args.habitat_loss_values):
        raise ValueError("--habitat-loss-values must be between 0 and 100.")
    if len(set(args.habitat_loss_values)) != len(args.habitat_loss_values):
        raise ValueError("--habitat-loss-values cannot contain duplicates.")
    args.habitat_loss_values = sorted(args.habitat_loss_values)
    if args.habitat_experiment and args.max_ticks != 3000:
        raise ValueError(
            "The final protocol requires --max-ticks 3000 because habitat "
            "loss is applied by NetLogo at tick 2000."
        )

    configured_space_size = parameter_space_size(args)
    if args.habitat_experiment:
        parameter_grid = []
        if args.output_dir == DEFAULT_OUTPUT_DIR:
            args.output_dir = DEFAULT_HABITAT_OUTPUT_DIR
        run_mode = "final habitat-loss experiment"
    elif args.validate_candidates is not None:
        parameter_grid = []
        if args.output_dir == DEFAULT_OUTPUT_DIR:
            args.output_dir = DEFAULT_CANDIDATE_OUTPUT_DIR
        run_mode = "explicit candidate validation"
    elif args.automated_calibration:
        parameter_grid = diverse_parameter_sample(
            args, args.screen_configurations, args.master_seed
        )
        if args.output_dir == DEFAULT_OUTPUT_DIR:
            args.output_dir = DEFAULT_AUTOMATED_OUTPUT_DIR
        run_mode = "automated two-stage calibration"
    elif args.quick:
        args.repetitions = 1
        args.max_ticks = 500
        args.stability_window = min(args.stability_window, args.max_ticks)
        parameter_grid = diverse_parameter_sample(
            args,
            args.quick_configurations,
            args.master_seed,
        )
        if args.output_dir == DEFAULT_OUTPUT_DIR:
            args.output_dir = SCRIPT_DIR / "pilot_outputs"
        run_mode = "diverse 500-tick pilot"
    elif args.random_runs is not None:
        if (
            args.random_runs >= configured_space_size
            and configured_space_size <= 100_000
        ):
            parameter_grid = full_parameter_grid(args)
        else:
            parameter_grid = diverse_parameter_sample(
                args, args.random_runs, args.master_seed
            )
        run_mode = "diverse parameter screen"
    else:
        parameter_grid = full_parameter_grid(args)
        run_mode = "full parameter grid"

    if args.automated_calibration and args.stability_window > args.screen_ticks:
        raise ValueError(
            "--stability-window cannot exceed --screen-ticks."
        )
    if args.validate_candidates is None and not args.automated_calibration:
        if not 1 <= args.stability_window <= args.max_ticks:
            raise ValueError(
                "--stability-window must be between 1 and --max-ticks."
            )
    if (
        args.validate_candidates is not None or args.automated_calibration
    ) and args.validation_stability_window > args.validation_ticks:
        raise ValueError(
            "--validation-stability-window cannot exceed --validation-ticks."
        )

    # Quick environment checks before importing heavy modules
    print(f"Default model path: {args.model}")
    print(f"Model exists: {args.model.exists()}")
    print(f"NETLOGO_HOME: {args.netlogo_home}")
    print(f"NETLOGO_HOME exists: {args.netlogo_home.exists()}")
    print(f"Mode: {run_mode}")
    print(f"Configured parameter space: {configured_space_size:,} combinations")
    if args.habitat_experiment:
        print(f"Selected configuration: {FINAL_CONFIGURATION_ID}")
        print(f"Habitat-loss levels: {args.habitat_loss_values}")
        print(
            f"Experiment: {len(args.habitat_loss_values)} treatments x "
            f"{args.repetitions} paired repetitions = "
            f"{len(args.habitat_loss_values) * args.repetitions} runs; "
            f"habitat loss at tick 2000; end at tick {args.max_ticks}"
        )
    elif args.validate_candidates is not None:
        print(f"Candidate file: {args.validate_candidates}")
        print(
            f"Validation: {args.validation_repetitions} repetitions x "
            f"{args.validation_ticks} ticks; assessment window="
            f"{args.validation_stability_window}; record every="
            f"{args.validation_record_every}"
        )
    elif args.automated_calibration:
        print(
            f"Screen: {args.screen_configurations} configurations x 1 repetition "
            f"x {args.screen_ticks} ticks"
        )
        print(
            f"Validation cap: {args.max_validation_configurations} configurations "
            f"x {args.validation_repetitions} repetitions x "
            f"{args.validation_ticks} ticks"
        )
    else:
        print(
            f"Configurations: {len(parameter_grid)}; repetitions: {args.repetitions}; "
            f"total runs: {len(parameter_grid) * args.repetitions}; "
            f"maximum ticks: {args.max_ticks}"
        )

    if not args.model.is_file():
        raise FileNotFoundError(f"NetLogo model not found: {args.model}")
    if not args.netlogo_home.exists():
        raise FileNotFoundError(
            f"NetLogo installation not found: {args.netlogo_home}. "
            "Pass --netlogo-home or set NETLOGO_HOME."
        )

    # Attempt to import heavy dependencies now and provide friendly errors
    try:
        import pandas as pd
    except Exception as e:
        print("ERROR: failed to import pandas.")
        print("Reason:", e)
        print("Try: pip3 install pandas")
        sys.exit(1)

    # Validate Java without starting it manually. pyNetLogo must start the JVM
    # itself so that it can supply the required NetLogo classpath.
    try:
        import jpype
        jvm_path = Path(os.environ["JAVA_HOME"]) / "lib" / "server" / "libjvm.dylib"
        print(f"JAVA_HOME: {os.environ.get('JAVA_HOME', '(not set)')}")
        print(f"JPype JVM path: {jvm_path}")
        print(f"JVM library exists: {jvm_path.is_file()}")
        if not jvm_path.is_file():
            raise FileNotFoundError(f"JPype selected a missing JVM library: {jvm_path}")
    except Exception as exc:
        print("ERROR: JPype could not locate a usable JVM.")
        print("Reason:", exc)
        print("Select Java 17 or 21 in the terminal before running this script.")
        print('Example: export JAVA_HOME=$(/usr/libexec/java_home -v 17)')
        sys.exit(1)

    try:
        import pynetlogo as pyNetLogo
    except Exception as e:
        print("ERROR: failed to import pynetlogo.")
        print("Reason:", e)
        print("Install it in the active environment with: python -m pip install pynetlogo")
        sys.exit(1)

    if args.habitat_experiment:
        run_habitat_experiment(args, pd, jpype, pyNetLogo, jvm_path)
        return

    if args.validate_candidates is not None:
        run_candidate_validation(args, pd, jpype, pyNetLogo, jvm_path)
        return

    if args.automated_calibration:
        run_automated_calibration(args, pd, jpype, pyNetLogo, jvm_path)
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    parameter_grid_rows = [
        {
            "config_id": config_id,
            **{name.replace("-", "_"): value for name, value in params.items()},
        }
        for config_id, params in enumerate(parameter_grid, start=1)
    ]
    pd.DataFrame(parameter_grid_rows).to_csv(
        args.output_dir / "baseline_parameter_grid.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "mode": run_mode,
                "master_seed": args.master_seed,
                "configured_parameter_space": configured_space_size,
                "configurations": len(parameter_grid),
                "repetitions": args.repetitions,
                "fixed_initial_predators": args.n_preds,
                "max_ticks": args.max_ticks,
                "record_every": args.record_every,
                "stability_window": args.stability_window,
                "min_final_preys": args.min_final_preys,
                "max_final_preys": args.max_final_preys,
                "min_final_predator_ratio": args.min_final_predator_ratio,
                "max_final_predator_ratio": args.max_final_predator_ratio,
                "max_prey_cv": args.max_prey_cv,
                "max_predator_cv": args.max_predator_cv,
                "max_relative_window_trend": args.max_relative_window_trend,
            }
        ]
    ).to_csv(args.output_dir / "baseline_run_settings.csv", index=False)

    random_generator = random.Random(args.master_seed)

    try:
        netlogo = pyNetLogo.NetLogoLink(
            gui=False,
            netlogo_home=str(args.netlogo_home),
            jvm_path=str(jvm_path),
        )
    except Exception as exc:
        print("ERROR: Java was found, but pyNetLogo could not start NetLogo.")
        print("Reason:", exc)
        print(
            "Confirm that Python and Java use the same architecture "
            "(both arm64 or both x86_64)."
        )
        raise
    run_rows: list[dict] = []
    time_rows: list[dict] = []
    interrupted = False
    run_error: Exception | None = None
    started_at = time.monotonic()

    try:
        netlogo.load_model(str(args.model))
        total = len(parameter_grid) * args.repetitions
        completed = 0

        for config_id, params in enumerate(parameter_grid, start=1):
            for repetition in range(args.repetitions):
                completed += 1
                seed = random_generator.randint(1, 2_147_483_647)
                summary, series = run_one(
                    netlogo,
                    config_id,
                    params,
                    repetition,
                    seed,
                    args.max_ticks,
                    args.record_every,
                    args,
                )
                run_rows.append(summary)
                time_rows.extend(series)

                if completed % args.checkpoint_every == 0:
                    save_outputs(pd, run_rows, time_rows, args.output_dir)

                elapsed = time.monotonic() - started_at
                estimated_remaining = (
                    elapsed / completed * (total - completed) if completed else 0
                )
                print(
                    f"[{completed}/{total}] config={config_id} rep={repetition} "
                    f"tick={summary['final_tick']} prey={summary['final_preys']} "
                    f"predators={summary['final_predators']} "
                    f"stable={summary['stable_dynamics']} "
                    f"accepted={summary['accepted_run']} "
                    f"ETA={format_duration(estimated_remaining)}"
                )
    except KeyboardInterrupt:
        interrupted = True
        print("\nInterrupted by user; saving all completed runs before exit.")
    except Exception as exc:
        run_error = exc
    finally:
        try:
            netlogo.kill_workspace()
        except Exception as exc:
            print(f"WARNING: NetLogo workspace did not close cleanly: {exc}")
            if run_error is None and not interrupted:
                run_error = exc
        finally:
            save_outputs(pd, run_rows, time_rows, args.output_dir)

    if run_error is not None:
        raise run_error
    if interrupted:
        print(f"Saved partial calibration outputs to {args.output_dir}")
        sys.exit(130)

    config_ranking = aggregate_configurations(pd.DataFrame(run_rows))

    print("\nTop baseline configurations:")
    print(config_ranking.head(10).to_string(index=False))
    print(f"\nSaved calibration outputs to {args.output_dir}")


if __name__ == "__main__":
    main()
