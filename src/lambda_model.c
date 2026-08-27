#include "lambda_model.h"

#include "acb.h"


void
lambda_0(
    arb_t result,
    const arb_t t,
    slong prec
)
{
    arb_t t2;
    arb_t tmp;
    arb_t numerator;
    arb_t denominator;
    arb_t imag_part;

    acb_t z;
    acb_t order;
    acb_t psi1;

    arb_init(t2);
    arb_init(tmp);
    arb_init(numerator);
    arb_init(denominator);
    arb_init(imag_part);

    acb_init(z);
    acb_init(order);
    acb_init(psi1);

    /*
     * Rational contribution:
     *
     *   2 (t^2 - 1/4) / (t^2 + 1/4)^2.
     */

    arb_mul(t2, t, t, prec);

    /*
     * numerator = 2 (t^2 - 1/4).
     *
     * Build the rational constant exactly before Arb rounding.
     */
    arb_mul_ui(numerator, t2, 4, prec);
    arb_sub_ui(numerator, numerator, 1, prec);
    arb_div_ui(numerator, numerator, 2, prec);

    /*
     * denominator = (t^2 + 1/4)^2.
     */
    arb_mul_ui(tmp, t2, 4, prec);
    arb_add_ui(tmp, tmp, 1, prec);
    arb_div_ui(tmp, tmp, 4, prec);
    arb_mul(denominator, tmp, tmp, prec);

    arb_div(result, numerator, denominator, prec);

    /*
     * z = 1/4 + i t/2.
     */
    arb_one(acb_realref(z));
    arb_div_ui(
        acb_realref(z),
        acb_realref(z),
        4,
        prec
    );

    arb_set(imag_part, t);
    arb_div_ui(imag_part, imag_part, 2, prec);
    arb_set(acb_imagref(z), imag_part);

    /*
     * order = 1, so acb_polygamma evaluates
     *
     *   psi^(1)(z),
     *
     * the trigamma function.
     */
    acb_one(order);

    acb_polygamma(
        psi1,
        order,
        z,
        prec
    );

    /*
     * lambda_0(t)
     *   = rational contribution
     *     + (1/4) Re psi^(1)(z).
     */
    arb_set(tmp, acb_realref(psi1));
    arb_div_ui(tmp, tmp, 4, prec);
    arb_add(result, result, tmp, prec);

    acb_clear(psi1);
    acb_clear(order);
    acb_clear(z);

    arb_clear(imag_part);
    arb_clear(denominator);
    arb_clear(numerator);
    arb_clear(tmp);
    arb_clear(t2);
}


void
lambda_kappa(
    arb_t result,
    const arb_t t,
    const arb_t kappa,
    slong prec
)
{
    arb_t half_kappa;

    arb_init(half_kappa);

    lambda_0(
        result,
        t,
        prec
    );

    arb_div_ui(
        half_kappa,
        kappa,
        2,
        prec
    );

    arb_sub(
        result,
        result,
        half_kappa,
        prec
    );

    arb_clear(half_kappa);
}