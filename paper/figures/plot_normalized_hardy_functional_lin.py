import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT = SCRIPT_DIR / "normalized_hardy_functional_t500.png"

mp.mp.dps = 30

def Xi(t):
    t = mp.mpf(t)
    s = mp.mpf("0.5") + 1j*t
    v = (
        mp.mpf("0.5") * s * (s - 1)
        * mp.power(mp.pi, -s/2)
        * mp.gamma(s/2) * mp.zeta(s)
    )
    return mp.re(v)

N0 = -2 * mp.diff(Xi, mp.mpf("0"), 2) / Xi(0)

def lambda0(t):
    t = mp.mpf(t)
    return (
        2*(t*t - mp.mpf("0.25")) / (t*t + mp.mpf("0.25"))**2
        + mp.mpf("0.25")
        * mp.re(mp.polygamma(1, mp.mpf("0.25") + 0.5j*t))
    )

def B_N0(t):
    t = mp.mpf(t)
    z = mp.siegelz(t, derivative=0)
    zp = mp.siegelz(t, derivative=1)
    zpp = mp.siegelz(t, derivative=2)
    return zp*zp - z*zpp + (lambda0(t) - N0/2)*z*z

t = np.linspace(0.02, 500.0, 2000)

B = []

for i, x in enumerate(t):
    B.append(float(B_N0(x)))

    if i % 100 == 0:
        print(f"{i}/{len(t)}  t={x:.2f}")

B = np.array(B)

fig = plt.figure(figsize=(12, 6.5))
ax = fig.add_subplot(111)
ax.plot(t, B)
ax.set_xlim(0.02, 500.0)
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$\mathcal{B}_{\mathcal{N}(0)}[Z](t)$")
ax.set_title(r"Normalized Hardy functional $\mathcal{B}_{\mathcal{N}(0)}[Z](t)$")
ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
fig.savefig(OUTPUT, dpi=240)
plt.close(fig)

print("N(0) =", mp.nstr(N0, 30))
