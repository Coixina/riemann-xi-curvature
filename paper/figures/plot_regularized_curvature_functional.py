import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT = SCRIPT_DIR / "regularized_curvature_functional.png"

mp.mp.dps = 50

def Xi(t):
    t = mp.mpf(t)
    s = mp.mpf("0.5") + 1j*t
    v = (mp.mpf("0.5") * s * (s-1)
         * mp.power(mp.pi, -s/2)
         * mp.gamma(s/2) * mp.zeta(s))
    return mp.re(v)

N0 = -2 * mp.diff(Xi, mp.mpf("0"), 2) / Xi(0)

def A(t):
    t = mp.mpf(t)
    return (mp.mpf("0.5") * (t*t + mp.mpf("0.25"))
            * mp.power(mp.pi, -mp.mpf("0.25"))
            * abs(mp.gamma(mp.mpf("0.25") + 0.5j*t)))

def lambda0(t):
    t = mp.mpf(t)
    return (2*(t*t-mp.mpf("0.25"))/(t*t+mp.mpf("0.25"))**2
            + mp.mpf("0.25") *
              mp.re(mp.polygamma(1, mp.mpf("0.25") + 0.5j*t)))

def H_N0(t):
    t = mp.mpf(t)
    z = mp.siegelz(t, derivative=0)
    zp = mp.siegelz(t, derivative=1)
    zpp = mp.siegelz(t, derivative=2)
    B = zp*zp - z*zpp + (lambda0(t)-N0/2)*z*z
    return A(t)**2 * B

t = np.linspace(0.02, 50.0, 900)
H = np.array([float(H_N0(x)) for x in t])

fig = plt.figure(figsize=(12, 6.5))
ax = fig.add_subplot(111)
ax.plot(t, H)
ax.set_yscale("log")
ax.set_xlim(0.02, 50.0)
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$H_{\mathcal{N}(0)}(t)$")
ax.set_title(r"Regularized curvature functional $H_{\mathcal{N}(0)}(t)$")
ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
fig.savefig(OUTPUT, dpi=240)
plt.close(fig)

print("N(0) =", mp.nstr(N0, 30))
