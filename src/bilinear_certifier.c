#include "bilinear_certifier.h"

#include "hardy_jet.h"
#include "lambda_model.h"


/*
 * Build
 *
 *     B_kappa[Z]
 *       =
 *     Z'^2
 *       - Z Z''
 *       + lambda_kappa Z^2
 *
 * using rigorous Arb arithmetic on the complete input
 * interval.
 */
static void
bilinear_build_value(
    arb_t result,
    const hardy_jet_t *jet,
    const arb_t lambda,
    slong prec
)
{
    arb_t tmp;
    arb_t z_squared;

    arb_init(tmp);
    arb_init(z_squared);


    /*
     * Z'^2.
     */
    arb_mul(
        result,
        jet->dt,
        jet->dt,
        prec
    );


    /*
     * - Z Z''.
     */
    arb_mul(
        tmp,
        jet->value,
        jet->d2t,
        prec
    );

    arb_sub(
        result,
        result,
        tmp,
        prec
    );


    /*
     * + lambda Z^2.
     */
    arb_mul(
        z_squared,
        jet->value,
        jet->value,
        prec
    );

    arb_mul(
        tmp,
        lambda,
        z_squared,
        prec
    );

    arb_add(
        result,
        result,
        tmp,
        prec
    );


    arb_clear(z_squared);
    arb_clear(tmp);
}


/*
 * Extract a rigorous directed lower bound
 *
 *     certified_lower
 *       <=
 *     inf(value).
 */
static int
bilinear_build_certified_lower(
    arf_t certified_lower,
    const arb_t value,
    slong prec
)
{
    if (!arb_is_finite(value))
        return 0;

    arb_get_lbound_arf(
        certified_lower,
        value,
        prec
    );

    return
        arf_is_finite(
            certified_lower
        );
}


void
bilinear_certificate_init(
    bilinear_certificate_t *certificate
)
{
    arb_init(
        certificate->value
    );

    arb_init(
        certificate->lambda
    );

    arf_init(
        certificate->certified_lower
    );

    certificate->status =
        BILINEAR_STATUS_INVALID;
}


void
bilinear_certificate_clear(
    bilinear_certificate_t *certificate
)
{
    arf_clear(
        certificate->certified_lower
    );

    arb_clear(
        certificate->lambda
    );

    arb_clear(
        certificate->value
    );
}


int
bilinear_certify(
    bilinear_certificate_t *certificate,
    const arb_t t,
    const arb_t kappa,
    slong prec
)
{
    hardy_jet_t jet;

    int success = 0;


    certificate->status =
        BILINEAR_STATUS_INVALID;

    arb_zero(
        certificate->value
    );

    arb_zero(
        certificate->lambda
    );

    arf_zero(
        certificate->certified_lower
    );


    if (prec <= 0)
        return 0;

    if (!arb_is_finite(t))
        return 0;

    if (!arb_is_finite(kappa))
        return 0;


    hardy_jet_init(
        &jet
    );


    /*
     * Direct rigorous Hardy-Z interval evaluation.
     */
    if (!hardy_jet_build(
            &jet,
            t,
            prec))
    {
        goto cleanup;
    }


    /*
     * lambda_kappa
     *   =
     * lambda_0 - kappa/2.
     */
    lambda_kappa(
        certificate->lambda,
        t,
        kappa,
        prec
    );

    if (!arb_is_finite(
            certificate->lambda))
    {
        goto cleanup;
    }


    /*
     * Build the exact regularized Hardy functional.
     */
    bilinear_build_value(
        certificate->value,
        &jet,
        certificate->lambda,
        prec
    );

    if (!arb_is_finite(
            certificate->value))
    {
        goto cleanup;
    }


    /*
     * Extract a rigorous lower bound.
     */
    if (!bilinear_build_certified_lower(
            certificate->certified_lower,
            certificate->value,
            prec))
    {
        goto cleanup;
    }


    if (arf_sgn(
            certificate->certified_lower) > 0)
    {
        certificate->status =
            BILINEAR_STATUS_CERTIFIED;
    }
    else
    {
        certificate->status =
            BILINEAR_STATUS_INCONCLUSIVE;
    }


    success = 1;


cleanup:

    hardy_jet_clear(
        &jet
    );

    return success;
}