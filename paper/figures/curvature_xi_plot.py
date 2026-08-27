"""
Numerical plot of the logarithmic curvature of Xi on the critical line.

We use the zero-product expression

    Xi(t) = Xi(0) prod_n (1 - t^2/gamma_n^2),

where gamma_n are the positive ordinates of the non-trivial zeros on the
critical line.

For this exploratory plot we approximate

    N(t) = -2 d^2/dt^2 log |Xi(t)|

by summing the contributions of the first N positive zeros:

    N(t) ≈ sum_n 4 (gamma_n^2 + t^2) / (gamma_n^2 - t^2)^2.

Zeros with gamma_n > 50 still contribute on [0, 50]; these are the "tails"
included in the computation. This figure is only numerical motivation and is
not used as a formal proof.
"""

from pathlib import Path

import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt


# ----- Parameters ---------------------------------------------------------

T_MIN = 0.0
T_MAX = 50.0
NUM_T_POINTS = 8000
NUM_ZEROS = 500
SCRIPT_DIR = Path(__file__).resolve().parent

OUTPUT = SCRIPT_DIR / "curvature_xi_0_50.png"

# ----- Computation --------------------------------------------------------

def riemann_zero_ordinates(n_zeros: int) -> np.ndarray:
    """Return the first n_zeros positive ordinates gamma_n."""
    gammas = []
    for n in range(1, n_zeros + 1):
        z = mp.zetazero(n)
        gammas.append(float(mp.im(z)))
    return np.array(gammas)


def curvature_from_zeros(t: np.ndarray, gammas: np.ndarray) -> np.ndarray:
    """
    Compute the truncated logarithmic curvature

        N(t) = sum_gamma 4 (gamma^2 + t^2) / (gamma^2 - t^2)^2.

    The expression has poles at the zeros gamma. For plotting purposes we do
    not remove them here; we later clip the vertical axis.
    """
    values = np.zeros_like(t)

    for gamma in gammas:
        g2 = gamma * gamma
        values += 4.0 * (g2 + t * t) / (g2 - t * t) ** 2

    return values


def main() -> None:
    mp.mp.dps = 50

    print(f"Computing the first {NUM_ZEROS} zero ordinates...")
    gammas = riemann_zero_ordinates(NUM_ZEROS)

    print("Computing curvature values...")
    t = np.linspace(T_MIN, T_MAX, NUM_T_POINTS)
    N = curvature_from_zeros(t, gammas)

    # Clip only for visualization: the true function has poles at zeros.
    y_max = np.percentile(N[np.isfinite(N)], 98.5)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 5))
    plt.plot(t, N, linewidth=0.8)
    plt.ylim(0, y_max)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$N(t)$")
    plt.title(r"Logarithmic curvature of $\Xi(t)$ on $0 \leq t \leq 50$")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUTPUT, dpi=300)

    print(f"Figure saved to: {OUTPUT.resolve()}")


if __name__ == "__main__":
    main()