#include <stdio.h>

#include "arb.h"
#include "arf.h"

#include "bilinear_certifier.h"


static int failures = 0;


static void
check(
    int condition,
    const char *message
)
{
    printf(
        "[%-72s] %s\n",
        message,
        condition ? "PASS" : "FAIL"
    );

    if (!condition)
        failures++;
}


int
main(void)
{
    const slong prec = 256;

    bilinear_certificate_t cert0;
    bilinear_certificate_t certk;

    arb_t t;
    arb_t kappa;

    arf_t lower;
    arf_t upper;


    printf(
        "============================================================\n"
        "BILINEAR CERTIFIER TESTS\n"
        "Hardy-only rigorous interval certification\n"
        "precision = %ld bits\n"
        "============================================================\n",
        prec
    );


    bilinear_certificate_init(
        &cert0
    );

    bilinear_certificate_init(
        &certk
    );

    arb_init(t);
    arb_init(kappa);

    arf_init(lower);
    arf_init(upper);


    /*
     * ========================================================
     * Point certification
     * ========================================================
     */

    printf(
        "\n"
        "=== Point certification ===\n"
    );

    arb_set_si(
        t,
        50
    );

    arb_zero(
        kappa
    );

    check(
        bilinear_certify(
            &cert0,
            t,
            kappa,
            prec
        ),
        "bilinear computation succeeds at t = 50"
    );

    check(
        cert0.status ==
            BILINEAR_STATUS_CERTIFIED,
        "t = 50 is rigorously certified for kappa = 0"
    );

    check(
        arf_sgn(
            cert0.certified_lower
        ) > 0,
        "certified lower bound is positive"
    );

    check(
        arb_is_finite(
            cert0.value
        ),
        "bilinear functional enclosure is finite"
    );

    check(
        arb_is_finite(
            cert0.lambda
        ),
        "lambda enclosure is finite"
    );


    /*
     * ========================================================
     * Genuine interval certification
     * ========================================================
     */

    printf(
        "\n"
        "=== Genuine interval certification ===\n"
    );

    arf_set_si(
        lower,
        50
    );

    arf_set_d(
        upper,
        50.01
    );

    arb_set_interval_arf(
        t,
        lower,
        upper,
        prec
    );

    arb_zero(
        kappa
    );

    check(
        bilinear_certify(
            &cert0,
            t,
            kappa,
            prec
        ),
        "bilinear computation succeeds on [50,50.01]"
    );

    check(
        cert0.status ==
            BILINEAR_STATUS_CERTIFIED ||
        cert0.status ==
            BILINEAR_STATUS_INCONCLUSIVE,
        "successful interval computation has a valid logical status"
    );


    /*
     * ========================================================
     * Kappa dependence
     * ========================================================
     */

    printf(
        "\n"
        "=== Kappa dependence ===\n"
    );

    arb_set_si(
        t,
        50
    );

    arb_zero(
        kappa
    );

    check(
        bilinear_certify(
            &cert0,
            t,
            kappa,
            prec
        ),
        "kappa = 0 computation succeeds"
    );


    arb_set_ui(
        kappa,
        1
    );

    arb_div_ui(
        kappa,
        kappa,
        8,
        prec
    );

    check(
        bilinear_certify(
            &certk,
            t,
            kappa,
            prec
        ),
        "kappa = 1/8 computation succeeds"
    );

    /*
     * Since
     *
     *   lambda_kappa
     *     =
     *   lambda_0 - kappa/2,
     *
     * positive kappa lowers the bilinear functional by
     *
     *   (kappa/2) Z^2.
     *
     * Hence the point value for kappa = 1/8 cannot exceed
     * the kappa = 0 value.
     */
    check(
        arf_cmp(
            certk.certified_lower,
            cert0.certified_lower
        ) <= 0,
        "positive kappa does not increase the certified lower bound"
    );


    /*
     * ========================================================
     * Invalid input
     * ========================================================
     */

    printf(
        "\n"
        "=== Invalid input ===\n"
    );

    arb_indeterminate(
        t
    );

    arb_zero(
        kappa
    );

    check(
        !bilinear_certify(
            &cert0,
            t,
            kappa,
            prec
        ),
        "indeterminate t is rejected"
    );


    arb_set_si(
        t,
        50
    );

    arb_indeterminate(
        kappa
    );

    check(
        !bilinear_certify(
            &cert0,
            t,
            kappa,
            prec
        ),
        "indeterminate kappa is rejected"
    );


    arb_zero(
        kappa
    );

    check(
        !bilinear_certify(
            &cert0,
            t,
            kappa,
            0
        ),
        "nonpositive precision is rejected"
    );


    /*
     * ========================================================
     * Summary
     * ========================================================
     */

    printf(
        "\n"
        "============================================================\n"
    );

    if (failures == 0)
    {
        printf(
            "ALL BILINEAR CERTIFIER TESTS PASSED\n"
        );
    }
    else
    {
        printf(
            "%d BILINEAR CERTIFIER TEST(S) FAILED\n",
            failures
        );
    }

    printf(
        "============================================================\n"
    );


    arf_clear(upper);
    arf_clear(lower);

    arb_clear(kappa);
    arb_clear(t);

    bilinear_certificate_clear(
        &certk
    );

    bilinear_certificate_clear(
        &cert0
    );

    return failures == 0 ? 0 : 1;
}