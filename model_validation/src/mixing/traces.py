"""How fast a chain's joint density trace forgets where it was.

The joint density is the no-drift instrument (`src/field/joint_density.h`), which makes it the one
scalar worth taking an autocorrelation time of. Why that is the measurement, and what it said, is
in `mixing_cost/findings.md`.
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass

import numpy as np
import pandas as pd


def autocorrelation(trace, max_lag: int | None = None) -> np.ndarray:
    """The normalised autocorrelation of a trace, lag 0 first.

    Computed through the Fourier transform, which is the same estimator as the
    direct sum with the biased (divide by `n`) normalisation: the tail lags are
    shrunk towards zero, which is what keeps the sum below well behaved.
    """
    x = np.asarray(trace, dtype=float)
    if x.ndim != 1:
        raise ValueError(f"A trace is one dimensional, but this one has shape {x.shape}.")
    if x.size < 2:
        raise ValueError(f"An autocorrelation needs at least two samples, but got {x.size}.")

    centred = x - x.mean()
    if max_lag is None:
        max_lag = x.size - 1
    max_lag = min(max_lag, x.size - 1)

    # Padded to at least 2n so that the circular correlation the transform
    # computes is the linear one.
    size = 1 << int(2 * x.size - 1).bit_length()
    spectrum = np.fft.rfft(centred, size)
    autocovariance = np.fft.irfft(spectrum * np.conjugate(spectrum), size)[: max_lag + 1]

    if autocovariance[0] <= 0.0:
        raise ValueError("The trace does not vary, so it has no autocorrelation.")
    return autocovariance / autocovariance[0]


def integrated_autocorrelation_time(trace) -> float:
    """Geyer's initial monotone sequence estimator of the autocorrelation time.

    `tau = 1 + 2 * sum over positive lags of rho`, which for an independent trace
    is 1 and for a sticky one is the number of iterations one independent draw
    costs. The sum is not taken to the end: the tail lags are noise, and adding
    them makes the estimate diverge. Geyer's rule pairs consecutive lags --
    `Gamma_m = rho_2m + rho_2m+1`, which is positive and decreasing for a
    reversible chain -- and stops at the first pair that is not.
    """
    rho = autocorrelation(trace)
    n_pairs = (rho.size - 1) // 2
    if n_pairs == 0:
        raise ValueError(
            f"An autocorrelation time needs at least three samples, but got {rho.size}."
        )

    pairs = rho[: 2 * n_pairs].reshape(n_pairs, 2).sum(axis=1)
    nonpositive = np.flatnonzero(pairs <= 0.0)
    cut = int(nonpositive[0]) if nonpositive.size else pairs.size
    # The monotone half of the rule: a pair larger than one before it is noise,
    # so the kept sequence is made non-increasing before it is summed.
    kept = np.minimum.accumulate(pairs[:cut]) if cut else np.empty(0)

    # rho_0 = 1 sits inside the first pair, so subtracting one leaves 2 * sum
    # over the positive lags. A trace that is anti-correlated at lag one gives a
    # value below 1, which is a real answer and not clamped away.
    return -1.0 + 2.0 * float(kept.sum())


def effective_sample_size(trace) -> float:
    """How many independent draws the trace is worth."""
    x = np.asarray(trace, dtype=float)
    return x.size / integrated_autocorrelation_time(x)


def half_split_drift(trace) -> float:
    """How far the second half of a trace sits from the first, in trace standard deviations.

    An autocorrelation time describes the distribution a chain has reached. A
    chain still climbing out of its initialisation has not reached one, and the
    time it reports then measures the climb. So the two are reported side by
    side: a drift of a standard deviation or more says to read the time next to
    it as a lower bound rather than as an answer.
    """
    x = np.asarray(trace, dtype=float)
    if x.size < 2:
        raise ValueError(f"A drift needs at least two samples, but got {x.size}.")
    sd = float(np.std(x))
    if sd <= 0.0:
        raise ValueError("The trace does not vary, so there is nothing to scale a drift by.")
    first, second = np.array_split(x, 2)
    return float(second.mean() - first.mean()) / sd


@dataclass(frozen=True)
class TraceSummary:
    """One chain's joint density trace, as the report prints it."""

    n_samples: int
    mean: float
    sd: float
    acf_lag_one: float
    autocorrelation_time: float
    effective_sample_size: float
    drift: float


def summarise_trace(trace) -> TraceSummary:
    x = np.asarray(trace, dtype=float)
    return TraceSummary(
        n_samples=int(x.size),
        mean=float(np.mean(x)),
        sd=float(np.std(x, ddof=1)),
        acf_lag_one=float(autocorrelation(x, max_lag=1)[1]),
        autocorrelation_time=integrated_autocorrelation_time(x),
        effective_sample_size=effective_sample_size(x),
        drift=half_split_drift(x),
    )


def _read_trace_frame(path: str | pathlib.Path, burn_in_rows: int) -> pd.DataFrame:
    """`*_joint_density.txt`, with the burn-in rows dropped.

    The trace is written from the first iteration on. Nothing in the writer knows
    the chain has not started yet -- unlike the parameter traces, which stattools
    withholds until `--writeBurnin` says otherwise -- so the rows the burn-in
    wrote are dropped here, by count.
    """
    frame = pd.read_csv(path, sep="\t")
    if burn_in_rows >= len(frame):
        raise ValueError(
            f"{path} holds {len(frame)} rows, which is not more than the "
            f"{burn_in_rows} burn-in rows to drop."
        )
    return frame.iloc[burn_in_rows:]


def read_joint_density_factors(
    path: str | pathlib.Path, burn_in_rows: int = 0
) -> dict[str, np.ndarray]:
    """Every column of the trace, in file order.

    The columns are the point of the file (`src/field/joint_density.h`): a single
    total says a chain drifts, the columns say which factor is dragging it. That
    matters here, because the factors mix at very different speeds -- the two
    node-state factors carry the phylogenetic parameters, and the data factor is
    a function of the field and the data-source parameters alone.
    """
    frame = _read_trace_frame(path, burn_in_rows)
    return {name: frame[name].to_numpy(dtype=float) for name in frame.columns}


