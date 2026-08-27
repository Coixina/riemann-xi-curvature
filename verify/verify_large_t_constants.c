/*
 * verify_large_t_constants.c
 *
 * Rigorous/reproducible checks for the scalar constants used in the
 * large-t theorem for
 *
 *   H_kappa(t) = Xi'(t)^2 - Xi(t) Xi''(t) - (kappa/2) Xi(t)^2.
 *
 * This program does NOT prove the structural Poisson/intertwining/Mittag-Leffler
 * identities. It verifies the scalar inequalities that enter after those analytic
 * identities have been proved.
 *
 * Target: FLINT 3.x (Arb integrated in FLINT).
 *
 * Checks currently included:
 *   1. Hurwitz-zeta sums used in the upper-endpoint tail:
 *        zeta(2,4/5) < 2.300
 *        zeta(3,4/5) < 2.216
 *   2. Upper non-stationary remainder at the first cell N=21:
 *        near part < 4.854
 *        far part  < 0.0091
 *        total     < 4.863
 *   3. Lower non-stationary remainder coefficient assembly:
 *        critical + Euler + non-stationary
 *        <= (4.23, 18.20, 42.61)
 *   4. Stirling/Cauchy scalar constant at t=200 and exact positivity
 *      of the shifted integer polynomials proving
 *        1/t^2 <= lambda_0^(3) <= 2/t^2  (t>=200).
 *   5. Direct certification of the global Hardy-real C1^Z,C2^Z jet
 *      bounds on |p|<=1.
 *   6. Conversion of those jet bounds into the derivative-remainder
 *      constants 0.954 and 2.000.
 *   7. Fixed-kappa transition corrections at t=882*pi:
 *        |R_tr| < 0.00708, |R_Q| < 0.00243.
 *   8. Current discrete majorant:
 *        L_{1/10}(20) < 0, L_{1/10}(21) > 1.5318.
 *   9. Exponential intertwining defect and exact-to-truncated
 *      perturbation at t=882*pi:
 *        |E_exp| negligible, perturbation < 0.00885.
 *  10. Uniform large-height theorem first-cell margin:
 *        B_{1/10}[Z] > 1.5134, hence > 1.50.
 *  11. Height-dependent master majorant at kappa=1/2:
 *        B_{1/2}[Z] > 0.0336 at N=21.
 *  12. Explicit quadratic bound checks for cells N=22,...,33,
 *      together with the N>=34 scalar inequalities used in Cor. 3.31.
 *
 * No numerical grid in p is used for the C1^Z,C2^Z bounds.
 */

#include <stdio.h>
#include <stdlib.h>

#include "arb.h"
#include "acb.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#define DEFAULT_PREC 256

static int failures = 0;

static void
set_q(arb_t x, slong num, slong den, slong prec)
{
    arb_set_si(x, num);
    arb_div_si(x, x, den, prec);
}

static void
pow_q(arb_t y, const arb_t x, slong num, slong den, slong prec)
{
    arb_t e;
    arb_init(e);
    set_q(e, num, den, prec);
    arb_pow(y, x, e, prec);
    arb_clear(e);
}

static void
print_ball(const char *name, const arb_t x)
{
    printf("  %-34s = ", name);
    arb_printn(x, 30, 0);
    printf("\n");
}

static void
check_lt(const char *name, const arb_t lhs, const arb_t rhs)
{
    printf("[%-52s] ", name);
    if (arb_lt(lhs, rhs))
    {
        printf("PASS\n");
    }
    else
    {
        printf("FAIL\n");
        print_ball("lhs", lhs);
        print_ball("rhs", rhs);
        failures++;
    }
}

static void
check_gt(const char *name, const arb_t lhs, const arb_t rhs)
{
    printf("[%-52s] ", name);
    if (arb_gt(lhs, rhs))
    {
        printf("PASS\n");
    }
    else
    {
        printf("FAIL\n");
        print_ball("lhs", lhs);
        print_ball("rhs", rhs);
        failures++;
    }
}


static void
verify_hurwitz_sums(slong prec)
{
    arb_t s, a, z, bound;

    arb_init(s);
    arb_init(a);
    arb_init(z);
    arb_init(bound);

    printf("\n=== Hurwitz sums ===\n");

    set_q(a, 4, 5, prec);

    arb_set_ui(s, 2);
    arb_hurwitz_zeta(z, s, a, prec);
    set_q(bound, 2300, 1000, prec);
    print_ball("zeta(2,4/5)", z);
    check_lt("zeta(2,4/5) < 2.300", z, bound);

    arb_set_ui(s, 3);
    arb_hurwitz_zeta(z, s, a, prec);
    set_q(bound, 2216, 1000, prec);
    print_ball("zeta(3,4/5)", z);
    check_lt("zeta(3,4/5) < 2.216", z, bound);

    arb_clear(s);
    arb_clear(a);
    arb_clear(z);
    arb_clear(bound);
}



static int
shifted_poly_has_strictly_positive_coeffs(const slong *coeff, slong len,slong shift)
{
    fmpz_poly_t f, g;
    fmpz_t c, a;
    slong k;
    int ok = 1;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_init(c);
    fmpz_init(a);

    for (k = 0; k < len; k++)
        fmpz_poly_set_coeff_si(f, k, coeff[k]);

    fmpz_set_si(c, shift);
    fmpz_poly_taylor_shift(g, f, c);

    for (k = 0; k < fmpz_poly_length(g); k++)
    {
        fmpz_poly_get_coeff_fmpz(a, g, k);
        if (fmpz_sgn(a) <= 0)
        {
            ok = 0;
            break;
        }
    }

    fmpz_clear(c);
    fmpz_clear(a);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);

    return ok;
}

static void
verify_stirling_lambda0(slong prec)
{
    /*
       P_-(x): numerator proving lambda0^(3) - 1/x > 0.
       P_+(x): numerator proving 2/x - lambda0^(3) > 0.
       Coefficients are in increasing degree order.
    */
    static const slong Pm[8] = {
        -105, 6159, -873208, 5647824,
        -5644800, -2983680, -645120, 1290240
    };
    static const slong Pp[8] = {
        210, -3219, 908488, -5412624,
        6585600, 5241600, 3655680, 430080
    };

    arb_t t, q, tmp, bound;

    arb_init(t);
    arb_init(q);
    arb_init(tmp);
    arb_init(bound);

    printf("\n=== Stirling / lambda_0 ===\n");

    printf("[%-52s] %s\n",
           "P_-(x+40000) has positive coefficients",
           shifted_poly_has_strictly_positive_coeffs(Pm, 8, 40000)
           ? "PASS" : "FAIL");
    if (!shifted_poly_has_strictly_positive_coeffs(Pm, 8, 40000))
        failures++;

    printf("[%-52s] %s\n",
           "P_+(x+40000) has positive coefficients",
           shifted_poly_has_strictly_positive_coeffs(Pp, 8, 40000)
           ? "PASS" : "FAIL");
    if (!shifted_poly_has_strictly_positive_coeffs(Pp, 8, 40000))
        failures++;

    /*
       Verify at t=200:
       t^2 * (32/105) * (8/(4t-1))^7 < 1.23e-10.
       The analytic proof separately establishes that the left-hand side
       decreases for t>=200.
    */
    arb_set_ui(t, 200);

    arb_mul_ui(tmp, t, 4, prec);
    arb_sub_ui(tmp, tmp, 1, prec);      /* 4t - 1 */

    arb_set_ui(q, 8);
    arb_div(q, q, tmp, prec);           /* 8/(4t-1) */
    arb_pow_ui(q, q, 7, prec);

    arb_mul_ui(q, q, 32, prec);
    arb_div_ui(q, q, 105, prec);

    arb_mul(q, q, t, prec);
    arb_mul(q, q, t, prec);             /* multiply by t^2 */

    set_q(bound, 123, 1000000000000L, prec); /* 1.23e-10 */
    print_ball("scaled Cauchy remainder at 200", q);
    check_lt("scaled Cauchy remainder < 1.23e-10", q, bound);

    arb_clear(t);
    arb_clear(q);
    arb_clear(tmp);
    arb_clear(bound);
}


/* ============================================================
   Direct C1/C2 global jet certification
   ============================================================ */

/*
   Gabcke global derivative bound for |p| <= 1.

   For r = 2n:
       |F^(2n)(p)| <= (2n)!/(2^n n!) * pi^n.

   For r = 2n+1:
       |F^(2n+1)(p)| <= 2^(n+1) n! pi^n.
*/
static void
build_F_derivative_bound(arb_t bound, slong r, slong prec)
{
    arb_t pi, pi_power;

    arb_init(pi);
    arb_init(pi_power);
    arb_const_pi(pi, prec);

    if ((r % 2) == 0)
    {
        ulong n = (ulong) (r / 2);
        fmpz_t num, den;

        fmpz_init(num);
        fmpz_init(den);

        fmpz_fac_ui(num, 2 * n);
        fmpz_fac_ui(den, n);

        arb_pow_ui(pi_power, pi, n, prec);
        arb_set_fmpz(bound, num);
        arb_mul(bound, bound, pi_power, prec);
        arb_div_fmpz(bound, bound, den, prec);
        arb_mul_2exp_si(bound, bound, -(slong) n);

        fmpz_clear(den);
        fmpz_clear(num);
    }
    else
    {
        ulong n = (ulong) ((r - 1) / 2);
        fmpz_t fac;

        fmpz_init(fac);
        fmpz_fac_ui(fac, n);

        arb_pow_ui(pi_power, pi, n, prec);
        arb_set_fmpz(bound, fac);
        arb_mul(bound, bound, pi_power, prec);
        arb_mul_2exp_si(bound, bound, (slong) (n + 1));

        fmpz_clear(fac);
    }

    arb_clear(pi_power);
    arb_clear(pi);
}


/*
   Direct Hardy-real C1^Z/C2^Z global jet certification.

   In the notation of the theorem,
       C1^Z(p) = F^(3)(p)/(12*pi^2),
       C2^Z(p) = F^(6)(p)/(288*pi^4) + F^(2)(p)/(16*pi^2).

   These are the coefficients after passage to Hardy's real function.
   In particular, the pre-Hardy term -i F/(96*pi) does not occur in C2^Z.

   The bounds below use Gabcke's global derivative bounds for F on
   |p| <= 1 and outward-rounded Arb arithmetic.
*/
static void
verify_C1C2_global_jets(slong prec)
{
    arb_t pi, pi2, pi4;
    arb_t F2, F3, F4, F5, F6, F7, F8;
    arb_t C10, C11, C12, C20, C21, C22;
    arb_t den, term, bound;

    arb_init(pi); arb_init(pi2); arb_init(pi4);
    arb_init(F2); arb_init(F3); arb_init(F4); arb_init(F5);
    arb_init(F6); arb_init(F7); arb_init(F8);
    arb_init(C10); arb_init(C11); arb_init(C12);
    arb_init(C20); arb_init(C21); arb_init(C22);
    arb_init(den); arb_init(term); arb_init(bound);

    printf("\n=== Hardy-real C1^Z,C2^Z global jet bounds on |p| <= 1 ===\n");

    arb_const_pi(pi, prec);
    arb_mul(pi2, pi, pi, prec);
    arb_mul(pi4, pi2, pi2, prec);

    build_F_derivative_bound(F2, 2, prec);
    build_F_derivative_bound(F3, 3, prec);
    build_F_derivative_bound(F4, 4, prec);
    build_F_derivative_bound(F5, 5, prec);
    build_F_derivative_bound(F6, 6, prec);
    build_F_derivative_bound(F7, 7, prec);
    build_F_derivative_bound(F8, 8, prec);

    /* C1^Z and its first two p-derivatives. */
    arb_mul_ui(den, pi2, 12, prec);
    arb_div(C10, F3, den, prec);
    arb_div(C11, F4, den, prec);
    arb_div(C12, F5, den, prec);

    print_ball("sup |C1^Z|", C10);
    set_q(bound, 106103296, 1000000000L, prec);
    check_lt("||C1^Z|| < 0.106103296", C10, bound);

    print_ball("sup |(C1^Z)'|", C11);
    set_q(bound, 250000001, 1000000000L, prec);
    check_lt("||(C1^Z)'|| < 0.250000001", C11, bound);

    print_ball("sup |(C1^Z)''|", C12);
    set_q(bound, 1333333334, 1000000000L, prec);
    check_lt("||(C1^Z)''|| < 1.333333334", C12, bound);

    /*
       C2^Z and its first two p-derivatives:
         F^(6)/(288*pi^4) + F^(2)/(16*pi^2),
         F^(7)/(288*pi^4) + F^(3)/(16*pi^2),
         F^(8)/(288*pi^4) + F^(4)/(16*pi^2).
    */
    arb_mul_ui(den, pi4, 288, prec);
    arb_div(C20, F6, den, prec);
    arb_div(C21, F7, den, prec);
    arb_div(C22, F8, den, prec);

    arb_mul_ui(den, pi2, 16, prec);

    arb_div(term, F2, den, prec);
    arb_add(C20, C20, term, prec);

    arb_div(term, F3, den, prec);
    arb_add(C21, C21, term, prec);

    arb_div(term, F4, den, prec);
    arb_add(C22, C22, term, prec);

    print_ball("sup |C2^Z|", C20);
    set_q(bound, 36474, 1000000L, prec);
    check_lt("||C2^Z|| < 0.036474", C20, bound);

    print_ball("sup |(C2^Z)'|", C21);
    set_q(bound, 185681, 1000000L, prec);
    check_lt("||(C2^Z)'|| < 0.185681", C21, bound);

    print_ball("sup |(C2^Z)''|", C22);
    set_q(bound, 552084, 1000000L, prec);
    check_lt("||(C2^Z)''|| < 0.552084", C22, bound);

    arb_clear(bound); arb_clear(term); arb_clear(den);
    arb_clear(C22); arb_clear(C21); arb_clear(C20);
    arb_clear(C12); arb_clear(C11); arb_clear(C10);
    arb_clear(F8); arb_clear(F7); arb_clear(F6); arb_clear(F5);
    arb_clear(F4); arb_clear(F3); arb_clear(F2);
    arb_clear(pi4); arb_clear(pi2); arb_clear(pi);
}


static void
verify_derivative_remainder_conversion(slong prec)
{
    /*
       The value remainder is taken directly from Gabcke at order K=2:
           |Z-P2| < 0.011 t^(-7/4).

       Here the derivative reconstruction is recomputed from the
       rounded Hardy-real C1^Z,C2^Z jet bounds used in the theorem.
    */
    arb_t pi, twopi, T0, C1, C1p, C1pp, C2, C2p, C2pp;
    arb_t d1, d2, dd1, dd2, tmp, tmp2, sum, bound;

    arb_init(pi); arb_init(twopi); arb_init(T0);
    arb_init(C1); arb_init(C1p); arb_init(C1pp);
    arb_init(C2); arb_init(C2p); arb_init(C2pp);
    arb_init(d1); arb_init(d2); arb_init(dd1); arb_init(dd2);
    arb_init(tmp); arb_init(tmp2); arb_init(sum); arb_init(bound);

    printf("\n=== C1^Z,C2^Z derivative-remainder conversion ===\n");

    arb_const_pi(pi, prec);
    arb_mul_ui(twopi, pi, 2, prec);
    arb_set_ui(T0, 200);

    set_q(C1, 106103296, 1000000000L, prec);
    set_q(C1p, 250000001, 1000000000L, prec);
    set_q(C1pp, 1333333334, 1000000000L, prec);

    set_q(C2, 36474, 1000000L, prec);
    set_q(C2p, 185681, 1000000L, prec);
    set_q(C2pp, 552084, 1000000L, prec);

    /* First derivative C1^Z, r=3/2. */
    pow_q(tmp, twopi, 3, 4, prec);
    arb_mul(tmp, tmp, C1, prec);
    arb_mul_ui(tmp, tmp, 3, prec);
    arb_div_ui(tmp, tmp, 4, prec);
    arb_sqrt(tmp2, T0, prec);
    arb_div(tmp, tmp, tmp2, prec);

    pow_q(tmp2, twopi, 1, 4, prec);
    arb_mul(tmp2, tmp2, C1p, prec);
    arb_add(d1, tmp, tmp2, prec);

    /* First derivative C2^Z, r=5/2. */
    pow_q(tmp, twopi, 5, 4, prec);
    arb_mul(tmp, tmp, C2, prec);
    arb_mul_ui(tmp, tmp, 5, prec);
    arb_div_ui(tmp, tmp, 4, prec);
    arb_div(tmp, tmp, T0, prec);

    pow_q(tmp2, twopi, 3, 4, prec);
    arb_mul(tmp2, tmp2, C2p, prec);
    arb_sqrt(sum, T0, prec);
    arb_div(tmp2, tmp2, sum, prec);
    arb_add(d2, tmp, tmp2, prec);

    set_q(bound, 418140, 1000000L, prec);
    print_ball("C1^Z first-derivative coefficient", d1);
    check_lt("C1^Z first-derivative coeff < 0.418140", d1, bound);

    set_q(bound, 54374, 1000000L, prec);
    print_ball("C2^Z first-derivative coefficient", d2);
    check_lt("C2^Z first-derivative coeff < 0.054374", d2, bound);

    /* Second derivative C1^Z, r=3/2. */
    pow_q(tmp, twopi, 3, 4, prec);
    arb_mul(tmp, tmp, C1, prec);
    arb_mul_ui(tmp, tmp, 21, prec);
    arb_div_ui(tmp, tmp, 16, prec);
    arb_div(tmp, tmp, T0, prec);
    arb_set(dd1, tmp);

    pow_q(tmp, twopi, 1, 4, prec);
    arb_mul(tmp, tmp, C1p, prec);
    arb_mul_ui(tmp, tmp, 5, prec);
    arb_div_ui(tmp, tmp, 2, prec);
    arb_sqrt(tmp2, T0, prec);
    arb_div(tmp, tmp, tmp2, prec);
    arb_add(dd1, dd1, tmp, prec);

    pow_q(tmp, twopi, -1, 4, prec);
    arb_mul(tmp, tmp, C1pp, prec);
    arb_add(dd1, dd1, tmp, prec);

    /* Second derivative C2^Z, r=5/2. */
    pow_q(tmp, twopi, 5, 4, prec);
    arb_mul(tmp, tmp, C2, prec);
    arb_mul_ui(tmp, tmp, 45, prec);
    arb_div_ui(tmp, tmp, 16, prec);
    arb_sqrt(tmp2, T0, prec);
    arb_mul(tmp2, tmp2, T0, prec);
    arb_div(tmp, tmp, tmp2, prec);
    arb_set(dd2, tmp);

    pow_q(tmp, twopi, 3, 4, prec);
    arb_mul(tmp, tmp, C2p, prec);
    arb_mul_ui(tmp, tmp, 7, prec);
    arb_div_ui(tmp, tmp, 2, prec);
    arb_div(tmp, tmp, T0, prec);
    arb_add(dd2, dd2, tmp, prec);

    pow_q(tmp, twopi, 1, 4, prec);
    arb_mul(tmp, tmp, C2pp, prec);
    arb_sqrt(tmp2, T0, prec);
    arb_div(tmp, tmp, tmp2, prec);
    arb_add(dd2, dd2, tmp, prec);

    set_q(bound, 914892, 1000000L, prec);
    print_ball("C1^Z second-derivative coefficient", dd1);
    check_lt("C1^Z second-derivative coeff < 0.914892", dd1, bound);

    set_q(bound, 75063, 1000000L, prec);
    print_ball("C2^Z second-derivative coefficient", dd2);
    check_lt("C2^Z second-derivative coeff < 0.075063", dd2, bound);

    set_q(sum, 481, 1000, prec);
    arb_add(sum, sum, d1, prec);
    arb_add(sum, sum, d2, prec);
    set_q(bound, 954, 1000, prec);
    print_ball("final E' coefficient", sum);
    check_lt("0.481 + C1^Z' + C2^Z' < 0.954", sum, bound);

    set_q(sum, 101, 100, prec);
    arb_add(sum, sum, dd1, prec);
    arb_add(sum, sum, dd2, prec);
    set_q(bound, 2, 1, prec);
    print_ball("final E'' coefficient", sum);
    check_lt("1.01 + C1^Z'' + C2^Z'' < 2.000", sum, bound);

    arb_clear(pi); arb_clear(twopi); arb_clear(T0);
    arb_clear(C1); arb_clear(C1p); arb_clear(C1pp);
    arb_clear(C2); arb_clear(C2p); arb_clear(C2pp);
    arb_clear(d1); arb_clear(d2); arb_clear(dd1); arb_clear(dd2);
    arb_clear(tmp); arb_clear(tmp2); arb_clear(sum); arb_clear(bound);
}

static void
harmonic_number(arb_t H, ulong N, slong prec)
{
    arb_t term;
    ulong n;

    arb_init(term);
    arb_zero(H);

    for (n = 1; n <= N; n++)
    {
        arb_set_ui(term, 1);
        arb_div_ui(term, term, n, prec);
        arb_add(H, H, term, prec);
    }

    arb_clear(term);
}

static void
compute_A_ell(arb_t A, arb_t ell, const arb_t t, slong prec)
{
    arb_t pi, T;

    arb_init(pi);
    arb_init(T);

    arb_const_pi(pi, prec);
    arb_mul_ui(T, pi, 2, prec);
    arb_div(T, t, T, prec);              /* T = t/(2*pi) */

    pow_q(A, T, 1, 4, prec);             /* A = T^(1/4) */
    arb_log(ell, T, prec);
    arb_div_ui(ell, ell, 2, prec);       /* ell = (1/2) log T */

    arb_clear(pi);
    arb_clear(T);
}

static void
compute_Eexp(arb_t out, const arb_t t, slong prec)
{
    arb_t pi, T, q, x;

    arb_init(pi);
    arb_init(T);
    arb_init(q);
    arb_init(x);

    arb_const_pi(pi, prec);
    arb_mul_ui(T, pi, 2, prec);
    arb_div(T, t, T, prec);              /* T = t/(2*pi) */

    pow_q(q, T, 3, 4, prec);
    arb_mul(q, q, pi, prec);
    arb_mul_ui(q, q, 8, prec);

    arb_log(x, T, prec);
    arb_add(x, x, pi, prec);
    arb_mul(q, q, x, prec);

    arb_mul(x, pi, t, prec);
    arb_neg(x, x);
    arb_exp(x, x, prec);
    arb_mul(out, q, x, prec);

    arb_clear(pi);
    arb_clear(T);
    arb_clear(q);
    arb_clear(x);
}

static void
compute_E0_DE_Ttr(arb_t E0, arb_t DE, arb_t Ttr,
                  const arb_t t, slong prec)
{
    arb_t A, ell, ep1, shape, tp, term, tmp;

    arb_init(A); arb_init(ell); arb_init(ep1); arb_init(shape);
    arb_init(tp); arb_init(term); arb_init(tmp);

    compute_A_ell(A, ell, t, prec);

    arb_add_ui(ep1, ell, 1, prec);
    arb_mul(shape, ep1, ep1, prec);
    arb_add_ui(shape, shape, 1, prec);

    /* E0 = 9.6 A(ell+1)t^-5/4
          + 10.05 A t^-7/4
          + 0.055 A((ell+1)^2+1)t^-7/4
          + 0.9216 t^-5/2
          + 0.02211 t^-7/2. */
    arb_zero(E0);

    pow_q(tp, t, -5, 4, prec);
    arb_mul(term, A, ep1, prec);
    arb_mul(term, term, tp, prec);
    set_q(tmp, 96, 10, prec);
    arb_mul(term, term, tmp, prec);
    arb_add(E0, E0, term, prec);

    pow_q(tp, t, -7, 4, prec);
    arb_mul(term, A, tp, prec);
    set_q(tmp, 1005, 100, prec);
    arb_mul(term, term, tmp, prec);
    arb_add(E0, E0, term, prec);

    arb_mul(term, A, shape, prec);
    arb_mul(term, term, tp, prec);
    set_q(tmp, 55, 1000, prec);
    arb_mul(term, term, tmp, prec);
    arb_add(E0, E0, term, prec);

    pow_q(tp, t, -5, 2, prec);
    set_q(tmp, 9216, 10000, prec);
    arb_mul(term, tp, tmp, prec);
    arb_add(E0, E0, term, prec);

    pow_q(tp, t, -7, 2, prec);
    set_q(tmp, 2211, 100000, prec);
    arb_mul(term, tp, tmp, prec);
    arb_add(E0, E0, term, prec);

    /* DE = 0.11(2*pi)^(-1/4)t^-3/2 + 0.000121 t^-7/2. */
    arb_zero(DE);
    arb_const_pi(tmp, prec);
    arb_mul_ui(tmp, tmp, 2, prec);
    pow_q(tmp, tmp, -1, 4, prec);
    pow_q(tp, t, -3, 2, prec);
    arb_mul(term, tmp, tp, prec);
    set_q(tmp, 11, 100, prec);
    arb_mul(term, term, tmp, prec);
    arb_add(DE, DE, term, prec);

    pow_q(tp, t, -7, 2, prec);
    set_q(tmp, 121, 1000000, prec);
    arb_mul(term, tp, tmp, prec);
    arb_add(DE, DE, term, prec);

    /* Ttr = 7.632 A(ell+1)t^-5/4
            + 8.044 A t^-7/4
            + 0.044 A((ell+1)^2+1)t^-7/4. */
    arb_zero(Ttr);

    pow_q(tp, t, -5, 4, prec);
    arb_mul(term, A, ep1, prec);
    arb_mul(term, term, tp, prec);
    set_q(tmp, 7632, 1000, prec);
    arb_mul(term, term, tmp, prec);
    arb_add(Ttr, Ttr, term, prec);

    pow_q(tp, t, -7, 4, prec);
    arb_mul(term, A, tp, prec);
    set_q(tmp, 8044, 1000, prec);
    arb_mul(term, term, tmp, prec);
    arb_add(Ttr, Ttr, term, prec);

    arb_mul(term, A, shape, prec);
    arb_mul(term, term, tp, prec);
    set_q(tmp, 44, 1000, prec);
    arb_mul(term, term, tmp, prec);
    arb_add(Ttr, Ttr, term, prec);

    arb_clear(A); arb_clear(ell); arb_clear(ep1); arb_clear(shape);
    arb_clear(tp); arb_clear(term); arb_clear(tmp);
}

static void
compute_GN(arb_t G, ulong n, slong prec)
{
    arb_t N, B, tmp, P, sqrtN, two_over_N;

    arb_init(N); arb_init(B); arb_init(tmp); arb_init(P);
    arb_init(sqrtN); arb_init(two_over_N);

    arb_set_ui(N, n);
    arb_add_ui(tmp, N, 2, prec);
    arb_set_ui(two_over_N, 2);
    arb_div(two_over_N, two_over_N, N, prec);
    arb_add(tmp, tmp, two_over_N, prec);
    arb_log(B, tmp, prec);

    arb_mul(P, B, B, prec);
    set_q(tmp, 423, 100, prec);
    arb_mul(P, P, tmp, prec);

    set_q(tmp, 1820, 100, prec);
    arb_mul(tmp, tmp, B, prec);
    arb_add(P, P, tmp, prec);

    set_q(tmp, 4261, 100, prec);
    arb_add(P, P, tmp, prec);

    arb_sqrt(sqrtN, N, prec);
    arb_div(G, P, sqrtN, prec);

    arb_clear(N); arb_clear(B); arb_clear(tmp); arb_clear(P);
    arb_clear(sqrtN); arb_clear(two_over_N);
}

static void
compute_L_fixed(arb_t L, ulong n, slong k_num, slong k_den, slong prec)
{
    arb_t N, logN, H, G, tmp;

    arb_init(N); arb_init(logN); arb_init(H); arb_init(G); arb_init(tmp);

    arb_set_ui(N, n);
    arb_log(logN, N, prec);
    harmonic_number(H, n, prec);
    compute_GN(G, n, prec);

    arb_pow_ui(L, logN, 3, prec);
    arb_mul_ui(L, L, 4, prec);
    arb_div_ui(L, L, 3, prec);

    set_q(tmp, k_num, k_den, prec);
    arb_mul(tmp, tmp, H, prec);
    arb_sub(L, L, tmp, prec);

    set_q(tmp, 4863, 1000, prec);
    arb_sub(L, L, tmp, prec);
    arb_sub(L, L, G, prec);

    set_q(tmp, 1, 1000, prec);
    arb_sub(L, L, tmp, prec);

    arb_clear(N); arb_clear(logN); arb_clear(H); arb_clear(G); arb_clear(tmp);
}

static void
verify_upper_residual_current(slong prec)
{
    arb_t N, near, far, total, x, y, z, pi, pi2, alpha, tmp, bound;

    arb_init(N); arb_init(near); arb_init(far); arb_init(total);
    arb_init(x); arb_init(y); arb_init(z); arb_init(pi); arb_init(pi2);
    arb_init(alpha); arb_init(tmp); arb_init(bound);

    printf("\n=== Upper non-stationary remainder (current N=21 constants) ===\n");

    arb_set_ui(N, 21);
    arb_const_pi(pi, prec);
    arb_mul(pi2, pi, pi, prec);
    set_q(alpha, 1098, 1000, prec);      /* alpha_N < 1.098 */

    /* near = 34.84 * [0.147*2.300/N + 0.055*2.216]. */
    set_q(x, 147, 1000, prec);
    set_q(y, 2300, 1000, prec);
    arb_mul(x, x, y, prec);
    arb_div(x, x, N, prec);

    set_q(y, 55, 1000, prec);
    set_q(z, 2216, 1000, prec);
    arb_mul(y, y, z, prec);
    arb_add(x, x, y, prec);
    set_q(y, 3484, 100, prec);
    arb_mul(near, y, x, prec);

    set_q(bound, 4854, 1000, prec);
    print_ball("near upper endpoint at N=21", near);
    check_lt("near upper endpoint < 4.854", near, bound);

    /* Far part: equations (97)--(100) of the current paper. */
    arb_zero(far);

    /* (19/(16*pi^2))*(9.82/N)*(2/N + 4/N^2) */
    set_q(x, 19, 16, prec);
    arb_div(x, x, pi2, prec);
    set_q(y, 982, 100, prec);
    arb_div(y, y, N, prec);
    arb_mul(x, x, y, prec);
    arb_set_ui(y, 2);
    arb_div(y, y, N, prec);
    arb_set_ui(z, 4);
    arb_mul(tmp, N, N, prec);
    arb_div(z, z, tmp, prec);
    arb_add(y, y, z, prec);
    arb_mul(x, x, y, prec);
    arb_add(far, far, x, prec);

    /* (4.5*alpha/(4*pi^2))*2.62*(2/N^2 + 8/N^3) */
    set_q(x, 45, 10, prec);
    arb_mul(x, x, alpha, prec);
    arb_div_ui(x, x, 4, prec);
    arb_div(x, x, pi2, prec);
    set_q(y, 262, 100, prec);
    arb_mul(x, x, y, prec);
    arb_mul(tmp, N, N, prec);
    arb_set_ui(y, 2);
    arb_div(y, y, tmp, prec);
    arb_mul(tmp, tmp, N, prec);
    arb_set_ui(z, 8);
    arb_div(z, z, tmp, prec);
    arb_add(y, y, z, prec);
    arb_mul(x, x, y, prec);
    arb_add(far, far, x, prec);

    /* (alpha/(2*pi^2))*2.62*(2/N^2 + 8/N^3) */
    arb_set(x, alpha);
    arb_div_ui(x, x, 2, prec);
    arb_div(x, x, pi2, prec);
    set_q(z, 262, 100, prec);
    arb_mul(x, x, z, prec);
    arb_mul(tmp, N, N, prec);
    arb_set_ui(y, 2);
    arb_div(y, y, tmp, prec);
    arb_mul(tmp, tmp, N, prec);
    arb_set_ui(z, 8);
    arb_div(z, z, tmp, prec);
    arb_add(y, y, z, prec);
    arb_mul(x, x, y, prec);
    arb_add(far, far, x, prec);

    /* (3*alpha^2/(4*pi^2))*0.89*N*(8/(3N^3)+16/N^4) */
    arb_mul(x, alpha, alpha, prec);
    arb_mul_ui(x, x, 3, prec);
    arb_div_ui(x, x, 4, prec);
    arb_div(x, x, pi2, prec);
    set_q(y, 89, 100, prec);
    arb_mul(x, x, y, prec);
    arb_mul(x, x, N, prec);

    arb_mul(tmp, N, N, prec);
    arb_mul(tmp, tmp, N, prec);          /* N^3 */
    arb_set_ui(y, 8);
    arb_div_ui(y, y, 3, prec);
    arb_div(y, y, tmp, prec);
    arb_mul(tmp, tmp, N, prec);          /* N^4 */
    arb_set_ui(z, 16);
    arb_div(z, z, tmp, prec);
    arb_add(y, y, z, prec);
    arb_mul(x, x, y, prec);
    arb_add(far, far, x, prec);

    set_q(bound, 91, 10000, prec);
    print_ball("far upper endpoint at N=21", far);
    check_lt("far upper endpoint < 0.0091", far, bound);

    arb_add(total, near, far, prec);
    set_q(bound, 4863, 1000, prec);
    print_ball("total upper remainder at N=21", total);
    check_lt("total upper remainder < 4.863", total, bound);

    arb_clear(N); arb_clear(near); arb_clear(far); arb_clear(total);
    arb_clear(x); arb_clear(y); arb_clear(z); arb_clear(pi); arb_clear(pi2);
    arb_clear(alpha); arb_clear(tmp); arb_clear(bound);
}

static void
verify_lower_remainder_rounding_current(slong prec)
{
    arb_t x, y, z, sum, bound;

    arb_init(x); arb_init(y); arb_init(z); arb_init(sum); arb_init(bound);

    printf("\n=== Lower remainder coefficient assembly ===\n");

    /* B_H^2: exact rational identity
    2.108 + 1.324 + 0.798 = 4.230 = 4.23. */
    if (2108 + 1324 + 798 == 4230)
        printf("[%-52s] PASS\n", "2.108 + 1.324 + 0.798 = 4.230");
    else
    {
        printf("[%-52s] FAIL\n", "2.108 + 1.324 + 0.798 = 4.230");
        failures++;
    }

    /* B_H: 9.706 + 5.297 + 3.193 = 18.196 < 18.20. */
    set_q(x, 9706, 1000, prec);
    set_q(y, 5297, 1000, prec);
    set_q(z, 3193, 1000, prec);
    arb_add(sum, x, y, prec);
    arb_add(sum, sum, z, prec);
    set_q(bound, 1820, 100, prec);
    print_ball("B_H coefficient sum", sum);
    check_lt("18.196 < 18.20", sum, bound);

    /* constant: 23.292 + 12.049 + 7.262 = 42.603 < 42.61. */
    set_q(x, 23292, 1000, prec);
    set_q(y, 12049, 1000, prec);
    set_q(z, 7262, 1000, prec);
    arb_add(sum, x, y, prec);
    arb_add(sum, sum, z, prec);
    set_q(bound, 4261, 100, prec);
    print_ball("constant coefficient sum", sum);
    check_lt("42.603 < 42.61", sum, bound);

    arb_clear(x); arb_clear(y); arb_clear(z); arb_clear(sum); arb_clear(bound);
}

static void
verify_transition_corrections_882pi(slong prec)
{
    arb_t pi, t, A, ell, ep1, shape, e0, e1, e2;
    arb_t S0, S1, S2, Rtr, RQ, tmp, bound;
    arb_t q0, q1, q2;

    arb_init(pi); arb_init(t); arb_init(A); arb_init(ell); arb_init(ep1);
    arb_init(shape); arb_init(e0); arb_init(e1); arb_init(e2);
    arb_init(S0); arb_init(S1); arb_init(S2); arb_init(Rtr); arb_init(RQ);
    arb_init(tmp); arb_init(bound); arb_init(q0); arb_init(q1); arb_init(q2);

    printf("\n=== Transition corrections at t=882*pi ===\n");

    arb_const_pi(pi, prec);
    arb_mul_ui(t, pi, 882, prec);
    compute_A_ell(A, ell, t, prec);
    arb_add_ui(ep1, ell, 1, prec);
    arb_mul(shape, ep1, ep1, prec);
    arb_add_ui(shape, shape, 1, prec);

    pow_q(tmp, t, -7, 4, prec);
    set_q(bound, 11, 1000, prec);
    arb_mul(e0, bound, tmp, prec);
    pow_q(tmp, t, -5, 4, prec);
    set_q(bound, 954, 1000, prec);
    arb_mul(e1, bound, tmp, prec);
    pow_q(tmp, t, -7, 4, prec);
    arb_mul_ui(e2, tmp, 2, prec);

    arb_mul_ui(S0, A, 4, prec);
    arb_mul(S1, S0, ep1, prec);
    arb_mul(S2, S0, shape, prec);

    /* Rtr <= 2 S1 e1 + S0 e2 + S2 e0 + 0.1 S0 e0. */
    arb_zero(Rtr);
    arb_mul(tmp, S1, e1, prec); arb_mul_ui(tmp, tmp, 2, prec); arb_add(Rtr, Rtr, tmp, prec);
    arb_mul(tmp, S0, e2, prec); arb_add(Rtr, Rtr, tmp, prec);
    arb_mul(tmp, S2, e0, prec); arb_add(Rtr, Rtr, tmp, prec);
    arb_mul(tmp, S0, e0, prec); set_q(bound, 1, 10, prec); arb_mul(tmp, tmp, bound, prec); arb_add(Rtr, Rtr, tmp, prec);

    set_q(bound, 708, 100000, prec);
    print_ball("R_tr majorant at 882*pi", Rtr);
    check_lt("R_tr < 0.00708", Rtr, bound);

    set_q(q0, 219339, 1000000, prec);
    set_q(q1, 3349, 1000000, prec);
    set_q(q2, 4201, 100000000, prec);

    /* RQ <= q1^2 + q0*q2 + 0.05*q0^2. */
    arb_mul(RQ, q1, q1, prec);
    arb_mul(tmp, q0, q2, prec); arb_add(RQ, RQ, tmp, prec);
    arb_mul(tmp, q0, q0, prec); set_q(bound, 5, 100, prec); arb_mul(tmp, tmp, bound, prec); arb_add(RQ, RQ, tmp, prec);

    set_q(bound, 243, 100000, prec);
    print_ball("R_Q majorant at 882*pi", RQ);
    check_lt("R_Q < 0.00243", RQ, bound);

    arb_clear(pi); arb_clear(t); arb_clear(A); arb_clear(ell); arb_clear(ep1);
    arb_clear(shape); arb_clear(e0); arb_clear(e1); arb_clear(e2);
    arb_clear(S0); arb_clear(S1); arb_clear(S2); arb_clear(Rtr); arb_clear(RQ);
    arb_clear(tmp); arb_clear(bound); arb_clear(q0); arb_clear(q1); arb_clear(q2);
}

static void
verify_current_first_cell(slong prec)
{
    arb_t L20, L21, pi, t, Eexp, perturb, A, ell, ep1, shape;
    arb_t e0, e1, e2, P0, P1, P2, tmp, bound, theorem_margin;
    arb_t E0, DE, Ttr, GN, H, Cstar, half_margin, tq;

    arb_init(L20); arb_init(L21); arb_init(pi); arb_init(t); arb_init(Eexp);
    arb_init(perturb); arb_init(A); arb_init(ell); arb_init(ep1); arb_init(shape);
    arb_init(e0); arb_init(e1); arb_init(e2); arb_init(P0); arb_init(P1); arb_init(P2);
    arb_init(tmp); arb_init(bound); arb_init(theorem_margin);
    arb_init(E0); arb_init(DE); arb_init(Ttr); arb_init(GN); arb_init(H);
    arb_init(Cstar); arb_init(half_margin); arb_init(tq);

    printf("\n=== Current first-cell checks: N=20/21 and t=882*pi ===\n");

    compute_L_fixed(L20, 20, 1, 10, prec);
    compute_L_fixed(L21, 21, 1, 10, prec);

    arb_zero(bound);
    print_ball("L_{1/10}(20)", L20);
    check_lt("L_{1/10}(20) < 0", L20, bound);

    set_q(bound, 15318, 10000, prec);
    print_ball("L_{1/10}(21)", L21);
    check_gt("L_{1/10}(21) > 1.5318", L21, bound);

    arb_const_pi(pi, prec);
    arb_mul_ui(t, pi, 882, prec);
    compute_Eexp(Eexp, t, prec);
    print_ball("E_exp(882*pi)", Eexp);
    set_q(bound, 1, 1000000000000L, prec); /* vastly weaker than actual */
    check_lt("E_exp(882*pi) < 1e-12", Eexp, bound);

    /* Fixed-kappa absolute perturbation from Lemma 3.22, using |lambda_1/10|<0.05. */
    compute_A_ell(A, ell, t, prec);
    arb_add_ui(ep1, ell, 1, prec);
    arb_mul(shape, ep1, ep1, prec);
    arb_add_ui(shape, shape, 1, prec);

    arb_mul_ui(P0, A, 5, prec);
    arb_mul(P1, P0, ep1, prec);
    arb_mul(P2, P0, shape, prec);

    pow_q(tmp, t, -7, 4, prec); set_q(bound, 11, 1000, prec); arb_mul(e0, bound, tmp, prec);
    pow_q(tmp, t, -5, 4, prec); set_q(bound, 954, 1000, prec); arb_mul(e1, bound, tmp, prec);
    pow_q(tmp, t, -7, 4, prec); arb_mul_ui(e2, tmp, 2, prec);

    arb_zero(perturb);
    arb_mul(tmp, P1, e1, prec); arb_mul_ui(tmp, tmp, 2, prec); arb_add(perturb, perturb, tmp, prec);
    arb_mul(tmp, P0, e2, prec); arb_add(perturb, perturb, tmp, prec);
    arb_mul(tmp, P2, e0, prec); arb_add(perturb, perturb, tmp, prec);
    arb_mul(tmp, e1, e1, prec); arb_add(perturb, perturb, tmp, prec);
    arb_mul(tmp, e0, e2, prec); arb_add(perturb, perturb, tmp, prec);
    arb_mul(tmp, P0, e0, prec); arb_mul_ui(tmp, tmp, 2, prec);
    arb_mul_2exp_si(tmp, tmp, -1);       /* multiply by 1/2 */
    arb_div_ui(tmp, tmp, 10, prec);      /* total factor 0.05 */
    arb_add(perturb, perturb, tmp, prec);
    arb_mul(tmp, e0, e0, prec);
    arb_mul_2exp_si(tmp, tmp, -1);
    arb_div_ui(tmp, tmp, 10, prec);
    arb_add(perturb, perturb, tmp, prec);

    set_q(bound, 885, 100000, prec);
    print_ball("absolute perturbation at 882*pi", perturb);
    check_lt("absolute perturbation < 0.00885", perturb, bound);

    /* Theorem 3.23 first-cell margin: L21 - .00708 - .00243 - .00885 - Eexp. */
    arb_set(theorem_margin, L21);
    set_q(tmp, 708, 100000, prec); arb_sub(theorem_margin, theorem_margin, tmp, prec);
    set_q(tmp, 243, 100000, prec); arb_sub(theorem_margin, theorem_margin, tmp, prec);
    set_q(tmp, 885, 100000, prec); arb_sub(theorem_margin, theorem_margin, tmp, prec);
    arb_sub(theorem_margin, theorem_margin, Eexp, prec);

    set_q(bound, 15134, 10000, prec);
    print_ball("Theorem 3.23 first-cell margin", theorem_margin);
    check_gt("first-cell margin > 1.5134", theorem_margin, bound);

    /* Corollary 3.29, kappa=1/2: formula (153) at N=21. */
    compute_E0_DE_Ttr(E0, DE, Ttr, t, prec);
    compute_GN(GN, 21, prec);
    harmonic_number(H, 21, prec);

    arb_set(Cstar, GN);
    set_q(tmp, 4863, 1000, prec); arb_add(Cstar, Cstar, tmp, prec);
    arb_add(Cstar, Cstar, Ttr, prec);
    set_q(tq, 2408, 100000, prec); arb_add(Cstar, Cstar, tq, prec);
    arb_set(tmp, DE); arb_mul_2exp_si(tmp, tmp, -1); arb_add(Cstar, Cstar, tmp, prec);

    arb_set_ui(tmp, 21); arb_log(tmp, tmp, prec); arb_pow_ui(half_margin, tmp, 3, prec);
    arb_mul_ui(half_margin, half_margin, 4, prec); arb_div_ui(half_margin, half_margin, 3, prec);
    arb_set(tmp, H); arb_mul_2exp_si(tmp, tmp, -1); arb_sub(half_margin, half_margin, tmp, prec);
    set_q(tmp, 1, 1000, prec); arb_sub(half_margin, half_margin, tmp, prec);
    arb_sub(half_margin, half_margin, Eexp, prec);
    arb_sub(half_margin, half_margin, E0, prec);
    arb_sub(half_margin, half_margin, Cstar, prec);

    set_q(bound, 336, 10000, prec);
    print_ball("kappa=1/2 first-cell margin", half_margin);
    check_gt("B_{1/2}[Z] first-cell majorant > 0.0336", half_margin, bound);

    arb_clear(L20); arb_clear(L21); arb_clear(pi); arb_clear(t); arb_clear(Eexp);
    arb_clear(perturb); arb_clear(A); arb_clear(ell); arb_clear(ep1); arb_clear(shape);
    arb_clear(e0); arb_clear(e1); arb_clear(e2); arb_clear(P0); arb_clear(P1); arb_clear(P2);
    arb_clear(tmp); arb_clear(bound); arb_clear(theorem_margin);
    arb_clear(E0); arb_clear(DE); arb_clear(Ttr); arb_clear(GN); arb_clear(H);
    arb_clear(Cstar); arb_clear(half_margin); arb_clear(tq);
}

static void
verify_quadratic_bound_cells(slong prec)
{
    static const slong milli_bounds[12] = {
        1011, 2910, 4770, 6592, 8378, 10131,
        11853, 13543, 15205, 16839, 18446, 20028
    };
    ulong N;
    arb_t pi, t, logN, logNp1, H, G, E0, DE, Ttr, Eexp;
    arb_t kplus, Q, Cstar, B, tmp, tq, bound, one;
    arb_t F, Fp, x;

    arb_init(pi); arb_init(t); arb_init(logN); arb_init(logNp1);
    arb_init(H); arb_init(G); arb_init(E0); arb_init(DE); arb_init(Ttr);
    arb_init(Eexp); arb_init(kplus); arb_init(Q); arb_init(Cstar); arb_init(B);
    arb_init(tmp); arb_init(tq); arb_init(bound); arb_init(one);
    arb_init(F); arb_init(Fp); arb_init(x);

    printf("\n=== Corollary 3.31 explicit quadratic-bound checks ===\n");
    arb_const_pi(pi, prec);
    arb_one(one);
    set_q(tq, 2408, 100000, prec);

    for (N = 22; N <= 33; N++)
    {
        /* left endpoint t = 2*pi*N^2 */
        arb_set_ui(t, N);
        arb_mul(t, t, t, prec);
        arb_mul(t, t, pi, prec);
        arb_mul_ui(t, t, 2, prec);

        arb_set_ui(logN, N);
        arb_log(logN, logN, prec);
        arb_set_ui(logNp1, N + 1);
        arb_log(logNp1, logNp1, prec);

        harmonic_number(H, N, prec);
        compute_GN(G, N, prec);
        compute_E0_DE_Ttr(E0, DE, Ttr, t, prec);
        compute_Eexp(Eexp, t, prec);

        /* kappa_N^+ = (2/25) log^2(N+1). */
        arb_mul(kplus, logNp1, logNp1, prec);
        arb_mul_ui(kplus, kplus, 2, prec);
        arb_div_ui(kplus, kplus, 25, prec);

        /* Q <= max(1,kappa_N^+). */
        if (arb_lt(kplus, one))
            arb_set(Q, one);
        else
            arb_set(Q, kplus);

        arb_set(Cstar, G);
        set_q(tmp, 4863, 1000, prec); arb_add(Cstar, Cstar, tmp, prec);
        arb_add(Cstar, Cstar, Ttr, prec);
        arb_add(Cstar, Cstar, tq, prec);
        arb_set(tmp, DE); arb_mul_2exp_si(tmp, tmp, -1); arb_add(Cstar, Cstar, tmp, prec);

        /* B = 4/3 log^3 N - kplus H - .001 - Eexp - E0 - Q*Cstar */
        arb_pow_ui(B, logN, 3, prec);
        arb_mul_ui(B, B, 4, prec); arb_div_ui(B, B, 3, prec);
        arb_mul(tmp, kplus, H, prec); arb_sub(B, B, tmp, prec);
        set_q(tmp, 1, 1000, prec); arb_sub(B, B, tmp, prec);
        arb_sub(B, B, Eexp, prec);
        arb_sub(B, B, E0, prec);
        arb_mul(tmp, Q, Cstar, prec); arb_sub(B, B, tmp, prec);

        set_q(bound, milli_bounds[N - 22], 1000, prec);
        {
            char name[80];
            snprintf(name, sizeof(name), "cell N=%lu lower bound > %.3f",
                     (unsigned long) N, milli_bounds[N - 22] / 1000.0);
            check_gt(name, B, bound);
        }
    }

    /* N=34 auxiliary numerical inequalities from the published proof. */
    compute_GN(G, 34, prec);
    set_q(bound, 27823, 1000, prec);
    print_ball("G_34", G);
    check_lt("G_34 < 27.823", G, bound);

    arb_set_ui(t, 34);
    arb_mul(t, t, t, prec);
    arb_mul(t, t, pi, prec);
    arb_mul_ui(t, t, 2, prec);
    compute_E0_DE_Ttr(E0, DE, Ttr, t, prec);
    compute_Eexp(Eexp, t, prec);

    /* remaining C* terms = Ttr + TQ + DE/2 < 0.028 */
    arb_set(tmp, Ttr);
    arb_add(tmp, tmp, tq, prec);
    arb_set(bound, DE); arb_mul_2exp_si(bound, bound, -1);
    arb_add(tmp, tmp, bound, prec);
    set_q(bound, 28, 1000, prec);
    print_ball("Ttr + TQ + DE/2 at N=34", tmp);
    check_lt("remaining C* terms at N=34 < 0.028", tmp, bound);

    arb_add(tmp, E0, Eexp, prec);
    set_q(bound, 9, 1000, prec);
    print_ball("E0 + Eexp at N=34", tmp);
    check_lt("E0 + Eexp at N=34 < 0.009", tmp, bound);

    /*
       F(x) = 4/3 x^3 - 0.01
              - (2/25)(x+1/34)^2(x+843/25).
       Check F(log 34)>20 and F'(log 34)>27.
    */
    arb_set_ui(x, 34);
    arb_log(x, x, prec);

    arb_pow_ui(F, x, 3, prec);
    arb_mul_ui(F, F, 4, prec); arb_div_ui(F, F, 3, prec);
    set_q(tmp, 1, 100, prec); arb_sub(F, F, tmp, prec);

    arb_set(tmp, x);
    set_q(bound, 1, 34, prec); arb_add(tmp, tmp, bound, prec);
    arb_mul(tmp, tmp, tmp, prec);
    arb_set(bound, x);
    set_q(Q, 843, 25, prec); arb_add(bound, bound, Q, prec);
    arb_mul(tmp, tmp, bound, prec);
    arb_mul_ui(tmp, tmp, 2, prec); arb_div_ui(tmp, tmp, 25, prec);
    arb_sub(F, F, tmp, prec);

    arb_set_ui(bound, 20);
    print_ball("F(log 34)", F);
    check_gt("F(log 34) > 20", F, bound);

    /* F' = 94/25 x^2 - 57424/10625 x - 57349/361250. */
    arb_mul(Fp, x, x, prec);
    set_q(tmp, 94, 25, prec); arb_mul(Fp, Fp, tmp, prec);
    set_q(tmp, 57424, 10625, prec); arb_mul(tmp, tmp, x, prec); arb_sub(Fp, Fp, tmp, prec);
    set_q(tmp, 57349, 361250, prec); arb_sub(Fp, Fp, tmp, prec);

    arb_set_ui(bound, 27);
    print_ball("F'(log 34)", Fp);
    check_gt("F'(log 34) > 27", Fp, bound);

    arb_clear(pi); arb_clear(t); arb_clear(logN); arb_clear(logNp1);
    arb_clear(H); arb_clear(G); arb_clear(E0); arb_clear(DE); arb_clear(Ttr);
    arb_clear(Eexp); arb_clear(kplus); arb_clear(Q); arb_clear(Cstar); arb_clear(B);
    arb_clear(tmp); arb_clear(tq); arb_clear(bound); arb_clear(one);
    arb_clear(F); arb_clear(Fp); arb_clear(x);
}

int
main(int argc, char **argv)
{
    slong prec = DEFAULT_PREC;

    if (argc >= 2)
    {
        long p = strtol(argv[1], NULL, 10);
        if (p < 64)
        {
            fprintf(stderr, "precision must be at least 64 bits\n");
            return 2;
        }
        prec = (slong) p;
    }

    printf("============================================================\n");
    printf("LARGE-t CONSTANT VERIFIER -- CURRENT ARTICLE CONSTANTS\n");
    printf("Arb/FLINT rigorous ball arithmetic\n");
    printf("precision = %ld bits\n", (long) prec);
    printf("paper interface: N=21, t=882*pi\n");
    printf("============================================================\n");

    verify_hurwitz_sums(prec);
    verify_upper_residual_current(prec);
    verify_lower_remainder_rounding_current(prec);
    verify_stirling_lambda0(prec);
    verify_C1C2_global_jets(prec);
    verify_derivative_remainder_conversion(prec);
    verify_transition_corrections_882pi(prec);
    verify_current_first_cell(prec);
    verify_quadratic_bound_cells(prec);

    printf("\n============================================================\n");
    if (failures == 0)
    {
        printf("ALL IMPLEMENTED CURRENT LARGE-t RIGOROUS CHECKS PASSED\n");
        printf("============================================================\n");
        return 0;
    }
    else
    {
        printf("%d CHECK(S) FAILED\n", failures);
        printf("============================================================\n");
        return 1;
    }
}
