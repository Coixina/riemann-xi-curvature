# Rigorous Certification of the Logarithmic Curvature of the Riemann Ξ Function

This repository contains the software, rigorous finite certificate, and
large-t constant verifier accompanying the manuscript

> **The Global Minimum and Large-Height Growth of the Logarithmic Curvature of the Riemann Xi Function**

The computational part of the proof uses Arb/FLINT ball arithmetic. The
finite-range certificate is based on direct rigorous evaluation of the
Z-function; the large-t verifier rigorously checks the explicit
numerical inequalities used in the analytic argument.

The repository provides:

- rigorous finite-interval certification using Arb/FLINT ball arithmetic;
- a supplied finite certificate with an independent verifier;
- rigorous verification of the explicit constants used in the large-t
  argument;
- an automated test suite covering the certification and verification
  machinery;
- the LaTeX sources of the associated manuscript.

The mathematical proof itself is given in the manuscript.

---

## 1. Proof architecture

The proof is divided into three ranges.

1. **Initial analytic range $0<t\le 2\pi.$**
   This range is handled analytically by a Hadamard-product argument and
   requires no computer-assisted certification.

2. **Finite certified range  $2\pi\le t\le 882\pi.$**
   On this compact interval, rigorous Arb/FLINT interval arithmetic
   certifies $\mathcal N(t)>\frac{1}{10}.$

3. **Large-t analytic range $t\ge 882\pi.$**
   An explicit Riemann-Siegel argument proves the stronger bound
   $\mathcal N(t)>\frac{1}{2}.$
   The executable `verify_large_t_constants` rigorously checks the explicit
   numerical inequalities entering this part of the proof and the subsequent
   height-dependent estimates.

The finite certifier in this repository addresses the second range. The
large-t verifier supports the explicit-constant computations in the third
range.

---

## 2. Quantity certified on the finite interval

The finite computation certifies a lower bound $\kappa$ for the logarithmic
curvature of the Riemann $\Xi$-function.

On the finite interval, the problem is reduced to rigorously proving
positivity of the Hardy-Z bilinear expression

```math
B_\kappa[Z](t)
=
Z'(t)^2-Z(t)Z''(t)+\lambda_\kappa(t)Z(t)^2.
```

where

$$
\lambda_\kappa(t)
    =
    \lambda_0(t)-\frac{\kappa}{2}.
$$

The supplied certificate uses

$$
\boxed{\kappa=\frac{1}{10}}.
$$

All evaluations of $Z$, $Z'$, and $Z''$ in the finite certifier are
performed directly with rigorous Arb/FLINT ball arithmetic. No
floating-point sampling is used in the certification.

---

## 3. Finite certification strategy

For a positive integer parameter $N$, the finite domain is

$$
[ 2\pi ,2\pi N^2 ].
$$

For the supplied certificate,

$$
N=21,
    \qquad
    2\pi N^2=882\pi.
$$

The interval is subdivided adaptively into exact symbolic dyadic intervals.
Arb balls are used only for rigorous analytic evaluation.

For each candidate interval, the certifier attempts to prove

$$
Z'^2-ZZ''+\lambda_\kappa Z^2>0
$$

throughout the complete interval. If the enclosure is inconclusive, the
interval is bisected exactly and its two children are tested recursively.

A successful run terminates only after the complete original interval has
been covered by rigorously certified terminal leaves.

---

## 4. Certificate and trust model

A successful certification run writes the terminal dyadic partition to a
plain-text certificate.

The expensive adaptive search used to discover this partition is **not**
part of the trusted certificate.

The independent verifier:

1. reads the stored certificate;
2. checks with exact integer arithmetic that the terminal dyadic leaves form
   a complete partition of the original interval;
3. rejects gaps, overlaps, duplicate leaves, inconsistent counts, malformed
   certificate data, and insufficient verification precision;
4. independently recomputes the rigorous Hardy-Z inequality on every
   stored leaf.

In particular, `verify_certificate` does **not** use the adaptive
`interval_engine` that generated the certificate.

If a stored leaf becomes inconclusive when checked at a different precision,
the verifier may refine that leaf locally. This refinement is used only for
the analytic verification and does not alter the exact partition stored in
the certificate.

This separation between **certificate generation** and **certificate
verification** is intentional. Verification does not require reproducing the
adaptive search that discovered the partition.

---

## 5. Repository structure

The principal files and directories are:

```text
.
├── CMakeLists.txt
├── LICENSE
├── README.md
│
├── include/                  # public headers for the finite certifier
├── src/                      # Hardy-only production implementation
│
├── certifier/
│   ├── rh_certifier.c        # adaptive finite-certificate generator
│   └── verify_certificate.c  # independent finite-certificate verifier
│
├── certificates/
│   └── finite_N21_kappa_1_10.txt
│
├── verify/
│   ├── verify_large_t_constants.c
│   ├── verify_large_t_constants_output.txt
│   └── verify_finite_N21_kappa_1_10.log
│
├── tests/                    # CTest unit/integration tests and fixtures
│
└── paper/
    └── figures/              # figures and code/data used to generate them
```

---

## 6. Software requirements

The code requires:

- a C compiler with C11 support;
- CMake >= 3.15;
- FLINT 3.x with Arb/Acb support.

The reference build and computations were performed on GNU/Linux using GCC
and FLINT 3.1.3. The finite certificate reported below was generated on a
system with an Intel Core i7-1165G7 processor (4 cores / 8 threads) and
16 GB of DDR4 memory.

On Debian/Ubuntu systems, CMake searches for `arb.h` in the usual FLINT
include locations and for the FLINT shared library automatically.

---

## 7. Building

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

The principal executables are then

```text
build/rh_certifier
build/verify_certificate
build/verify_large_t_constants
```

The CMake configuration prints the detected FLINT include directory and
library.

---

## 8. Automated test suite

After building, run

```bash
ctest --test-dir build --output-on-failure
```

The suite tests, among other things:

- exact interval geometry;
- the $\lambda_\kappa$ model;
- rigorous Hardy jets;
- the bilinear certifier;
- the adaptive interval engine;
- verifier local refinement;
- certificate-fixture generation;
- successful verification at several precisions;
- rejection of gaps, overlaps, duplicate leaves, and inconsistent counts;
- rejection of insufficient verification precision;
- the implemented large-t constant checks.

The current CMake configuration defines 18 tests.

---

## 9. Quick verification of the supplied results

It is **not necessary to regenerate the finite certificate** in order to
verify the supplied computation.

### 9.1 Large-t constants

Run

```bash
./build/verify_large_t_constants 256
```

A successful run ends with

```text
ALL IMPLEMENTED CURRENT LARGE-t RIGOROUS CHECKS PASSED
```

A reference output is included as

```text
verify/verify_large_t_constants_output.txt
```

### 9.2 Finite certificate

Verify the supplied finite certificate with

```bash
./build/verify_certificate \
    certificates/finite_N21_kappa_1_10.txt
```

The certificate records the precision at which it was generated. If no
verification precision is supplied, the verifier uses the stored generator
precision.

The same exact partition may also be checked at a higher precision, for
example

```bash
./build/verify_certificate \
    certificates/finite_N21_kappa_1_10.txt \
    384
```

or

```bash
./build/verify_certificate \
    certificates/finite_N21_kappa_1_10.txt \
    512
```

At higher precision the verifier is allowed to perform additional local
dyadic refinement when necessary.

Verification of the complete certificate is computationally intensive and
may take several hours, depending on the machine.

---

## 10. Supplied finite certificate

The certificate distributed with the repository is

```text
certificates/finite_N21_kappa_1_10.txt
```

Its header is

```text
RH_FINITE_CERTIFICATE 1
large_t_N 21
kappa_num 1
kappa_den 10
generator_precision 256
leaf_count 00000000000001914979
BEGIN_LEAVES
```

Thus the certificate corresponds to

```text
N                   = 21
interval            = [2*pi, 882*pi]
kappa               = 1/10
generator precision = 256 bits
terminal leaves     = 1,914,979
```

The certificate file is approximately 20 MB.

### SHA-256

The SHA-256 digest of the supplied certificate is

```text
3272f10d33533257e12dab5ff498969efc1525756825717e4f6fb4d6fee1a045
```

It can be checked with

```bash
sha256sum certificates/finite_N21_kappa_1_10.txt
```

The digest identifies the exact certificate file associated with the
reported computation. It is **not** a mathematical verification of the
certificate; mathematical verification is performed by
`verify_certificate`.

---

## 11. Regenerating the finite certificate

Independent verification does not require regeneration.

For users who wish to repeat the adaptive search itself, the supplied
certificate can be regenerated with

```bash
./build/rh_certifier \
    21 \
    1/10 \
    certificates/finite_N21_kappa_1_10.txt \
    256 \
    40 \
    10000
```

The command-line interface is

```text
rh_certifier N kappa certificate [precision] [max_depth] [progress_every]
```

where `kappa` may be entered as an exact rational or in decimal notation.

The default optional parameters are

```text
precision       = 256 bits
maximum depth   = 40
progress report = every 10000 evaluations
```

During a long computation, progress output reports both the number of
analytic evaluations and the fraction of the original interval already
rigorously covered.

Full regeneration is computationally intensive. The supplied
$N=21$, $\kappa=1/10$, 256-bit certificate was generated in
approximately 2 h 50 min on the reference system described above. This
timing is only indicative; runtime depends strongly on the processor, system
load, and other running processes.

The supplied certificate is provided specifically so that independent
verification does not depend on reproducing the adaptive search path.

**Important:** regenerating directly to
`certificates/finite_N21_kappa_1_10.txt` will replace the supplied certificate.
To preserve the release copy, use another output filename, for example

```bash
./build/rh_certifier \
    21 \
    1/10 \
    certificates/finite_N21_kappa_1_10_regenerated.txt \
    256 \
    40 \
    10000
```

and compare the resulting certificate as appropriate.

---

## 12. Large-t constant verifier

The executable `verify_large_t_constants` uses rigorous Arb/FLINT ball
arithmetic to verify explicit numerical inequalities entering the
large-t proof.

At the reference precision,

```bash
./build/verify_large_t_constants 256
```

it checks, among other quantities:

- the Hurwitz-zeta bounds used in the upper-endpoint tail;
- the current upper non-stationary remainder at $N=21$;
- the lower non-stationary remainder coefficient assembly;
- the Stirling/Cauchy scalar bound and associated positivity checks;
- the global Hardy-real $C_1^Z,C_2^Z$ jet bounds;
- conversion of those jet bounds into derivative-remainder constants;
- the transition corrections at $t=882\pi$;
- $L_{1/10}(20)<0$ and $L_{1/10}(21)>1.5318$;
- the exact-to-truncated perturbation bound at $t=882\pi$;
- the first-cell margin used in the uniform large-height theorem;
- the $\kappa=1/2$ first-cell master-majorant check;
- the explicit height-dependent cell checks for $N=22,\ldots,33$;
- the auxiliary inequalities used from $N=34$ onward.

Higher precision may be requested independently, for example

```bash
./build/verify_large_t_constants 384
```

or

```bash
./build/verify_large_t_constants 512
```

This program does not establish the large-t theorem by numerical
sampling. It rigorously verifies the explicit numerical inequalities and
constant comparisons used by the analytic proof.

---

## 13. Reproducibility levels

There are three distinct levels of computational reproduction.

### A. Automated implementation tests

```bash
ctest --test-dir build --output-on-failure
```

These test individual components, integration paths, certificate rejection
cases, and the large-t verifier.

### B. Independent verification of the supplied certificate

```bash
./build/verify_certificate \
    certificates/finite_N21_kappa_1_10.txt
```

This recomputes the rigorous analytic certification on the stored exact
partition without repeating the adaptive search.

### C. Full regeneration

```bash
./build/rh_certifier \
    21 1/10 \
    certificates/finite_N21_kappa_1_10_regenerated.txt \
    256 40 10000
```

This repeats the adaptive search from the original interval.

The distinction is important: the certificate is a compact witness of a
successful search, while the independent verifier rechecks the mathematical
content of that witness without trusting the search path.

---

## 14. Implementation boundary

The production finite-certification path is deliberately small and
Hardy-only:

```text
exact dyadic geometry
        ↓
direct rigorous Hardy-Z interval jet
        ↓
rigorous lambda_kappa interval
        ↓
bilinear positivity test
        ↓
adaptive exact subdivision
        ↓
plain-text certificate
```

No Riemann-Siegel approximation or external remainder budget is used by the
finite certifier.

The Riemann-Siegel machinery belongs to the separate analytic large-t
argument in the manuscript; `verify_large_t_constants` checks its explicit
numerical constants.

---

## 15. Figure sources

The `paper/figures/` directory contains the figures used in the manuscript,
together with the code and data needed to generate them.

The manuscript PDF and LaTeX sources are intentionally not distributed in
this repository. The repository is intended to provide the computational
artifacts supporting the paper: source code, rigorous certificates,
independent verifiers, tests, and reproducible figure material.

---

## 16. License

The source code in this repository is released under the MIT License. See
`LICENSE` for details.

This project uses FLINT, which is distributed under the GNU Lesser General
Public License (LGPL), version 3 or later.
