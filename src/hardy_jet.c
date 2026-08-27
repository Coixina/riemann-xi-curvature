#include "hardy_jet.h"

#include "acb.h"
#include "acb_dirichlet.h"


#define HARDY_JET_LENGTH 3


void
hardy_jet_init(
    hardy_jet_t *jet
)
{
    arb_init(
        jet->value
    );

    arb_init(
        jet->dt
    );

    arb_init(
        jet->d2t
    );
}


void
hardy_jet_clear(
    hardy_jet_t *jet
)
{
    arb_clear(
        jet->d2t
    );

    arb_clear(
        jet->dt
    );

    arb_clear(
        jet->value
    );
}


int
hardy_jet_build(
    hardy_jet_t *jet,
    const arb_t t,
    slong prec
)
{
    slong j;

    acb_t t_acb;
    acb_ptr coeff;

    int success = 0;


    if (prec <= 0)
        return 0;

    if (!arb_is_finite(t))
        return 0;


    acb_init(
        t_acb
    );

    coeff =
        _acb_vec_init(
            HARDY_JET_LENGTH
        );


    /*
     * Embed the real Arb input into Acb.
     */
    arb_set(
        acb_realref(t_acb),
        t
    );

    arb_zero(
        acb_imagref(t_acb)
    );


    /*
     * FLINT returns the Taylor coefficients
     *
     *     Z(t+x)
     *       =
     *     c_0 + c_1 x + c_2 x^2 + ...
     *
     * through order 2.
     */
    acb_dirichlet_hardy_z(
        coeff,
        t_acb,
        NULL,
        NULL,
        HARDY_JET_LENGTH,
        prec
    );


    /*
     * The exact Hardy-Z function and all of its derivatives
     * are real on the real axis.
     *
     * Before projecting an Acb coefficient to its real part,
     * require the computed imaginary enclosure to contain zero.
     */
    for (j = 0;
         j < HARDY_JET_LENGTH;
         ++j)
    {
        if (!acb_is_finite(
                coeff + j))
        {
            goto cleanup;
        }

        if (!arb_contains_zero(
                acb_imagref(coeff + j)))
        {
            goto cleanup;
        }
    }


    /*
     * Z(t) = c_0.
     */
    arb_set(
        jet->value,
        acb_realref(coeff + 0)
    );


    /*
     * Z'(t) = c_1.
     */
    arb_set(
        jet->dt,
        acb_realref(coeff + 1)
    );


    /*
     * Z''(t) = 2 c_2.
     */
    arb_mul_ui(
        jet->d2t,
        acb_realref(coeff + 2),
        2,
        prec
    );


    if (!arb_is_finite(jet->value) ||
        !arb_is_finite(jet->dt) ||
        !arb_is_finite(jet->d2t))
    {
        goto cleanup;
    }


    success = 1;


cleanup:

    _acb_vec_clear(
        coeff,
        HARDY_JET_LENGTH
    );

    acb_clear(
        t_acb
    );

    return success;
}