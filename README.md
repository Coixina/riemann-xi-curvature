# Riemann Xi Curvature — rigorous certification code

Computational artifacts accompanying the manuscript

> **The Global Minimum and Large-Height Growth of the Logarithmic Curvature of the Riemann Xi Function**

The mathematical proof is given in the manuscript. This repository contains the code and rigorous computational artifacts used for:

- the finite interval $2\pi \le t \le 882\pi$;
- verification of the explicit constants in the large-height argument;
- reproducible generation of the numerical figures.

The computations use Arb/FLINT ball arithmetic.

## Requirements

- C compiler with C11 support
- CMake >= 3.15
- FLINT 3.x with Arb/Acb support

Reference computations were performed with FLINT 3.1.3 on GNU/Linux.

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

## Run the tests

```bash
ctest --test-dir build --output-on-failure
```

## Verify the large-height constants

```bash
./build/verify_large_t_constants 256
```

A successful run ends with

```text
ALL IMPLEMENTED CURRENT LARGE-t RIGOROUS CHECKS PASSED
```

A reference output is provided in

[`verify/verify_large_t_constants_output.txt`](verify/verify_large_t_constants_output.txt).


## Verify the supplied finite certificate

```bash
./build/verify_certificate \
    certificates/finite_N21_kappa_1_10.txt
```

The certificate covers

$$
2\pi \le t \le 882\pi,
\qquad
\kappa=\frac{1}{10}.
$$

with 256-bit generator precision and 1,914,979 terminal dyadic leaves.

Its SHA-256 digest is

```text
3272f10d33533257e12dab5ff498969efc1525756825717e4f6fb4d6fee1a045
```

Check it with

```bash
sha256sum certificates/finite_N21_kappa_1_10.txt
```

The verifier checks the exact dyadic partition with integer arithmetic and independently recomputes the rigorous analytic inequality on every stored leaf. Regenerating the adaptive search is not required for verification.

## Regenerate the finite certificate

Optional:

```bash
./build/rh_certifier \
    21 \
    1/10 \
    certificates/finite_N21_kappa_1_10_regenerated.txt \
    256 \
    40 \
    10000
```

Full regeneration is substantially more expensive than verification.

## Main files

- [`certifier/rh_certifier.c`](certifier/rh_certifier.c)
- [`certifier/verify_certificate.c`](certifier/verify_certificate.c)
- [`certificates/finite_N21_kappa_1_10.txt`](certificates/finite_N21_kappa_1_10.txt)
- [`verify/verify_large_t_constants.c`](verify/verify_large_t_constants.c)
- [`verify/verify_large_t_constants_output.txt`](verify/verify_large_t_constants_output.txt)
- [`verify/verify_finite_N21_kappa_1_10.log`](verify/verify_finite_N21_kappa_1_10.log)
- [`paper/figures/`](paper/figures/)

The manuscript PDF and LaTeX sources are not distributed in this repository.

## License

The source code is released under the MIT License. See `LICENSE`.

FLINT is distributed under the GNU Lesser General Public License (LGPL), version 3 or later.
