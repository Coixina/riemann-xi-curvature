#include <stdio.h>

#include "arb.h"

#include "hardy_jet.h"


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


static int
arb_overlaps_negated(
    const arb_t x,
    const arb_t y
)
{
    arb_t minus_y;

    int result;

    arb_init(
        minus_y
    );

    arb_neg(
        minus_y,
        y
    );

    result =
        arb_overlaps(
            x,
            minus_y
        );

    arb_clear(
        minus_y
    );

    return result;
}


int
main(void)
{
    const slong prec_low  = 128;
    const slong prec_high = 256;

    hardy_jet_t jet_zero;
    hardy_jet_t jet_pos;
    hardy_jet_t jet_neg;
    hardy_jet_t jet_low;
    hardy_jet_t jet_high;
    hardy_jet_t jet_interval;

    arb_t t;
    arb_t minus_t;
    arb_t interval;

    arf_t lower;
    arf_t upper;


    printf(
        "============================================================\n"
        "HARDY JET TESTS\n"
        "Direct rigorous FLINT/Arb Hardy-Z evaluation\n"
        "============================================================\n"
    );


    hardy_jet_init(&jet_zero);
    hardy_jet_init(&jet_pos);
    hardy_jet_init(&jet_neg);
    hardy_jet_init(&jet_low);
    hardy_jet_init(&jet_high);
    hardy_jet_init(&jet_interval);

    arb_init(t);
    arb_init(minus_t);
    arb_init(interval);

    arf_init(lower);
    arf_init(upper);


    /*
     * ========================================================
     * t = 0
     * ========================================================
     */

    printf(
        "\n"
        "=== Evaluation at t = 0 ===\n"
    );

    arb_zero(
        t
    );

    check(
        hardy_jet_build(
            &jet_zero,
            t,
            prec_high
        ),
        "Hardy jet can be evaluated at t = 0"
    );

    check(
        arb_is_finite(jet_zero.value) &&
        arb_is_finite(jet_zero.dt) &&
        arb_is_finite(jet_zero.d2t),
        "all three components at t = 0 are finite"
    );

    /*
     * Hardy Z is even:
     *
     *     Z(-t) = Z(t),
     *
     * hence
     *
     *     Z'(0) = 0.
     */
    check(
        arb_contains_zero(
            jet_zero.dt
        ),
        "Z'(0) enclosure contains zero"
    );


    /*
     * ========================================================
     * Parity
     * ========================================================
     */

    printf(
        "\n"
        "=== Hardy-Z parity ===\n"
    );

    arb_set_si(
        t,
        10
    );

    arb_neg(
        minus_t,
        t
    );

    check(
        hardy_jet_build(
            &jet_pos,
            t,
            prec_high
        ),
        "Hardy jet can be evaluated at t = 10"
    );

    check(
        hardy_jet_build(
            &jet_neg,
            minus_t,
            prec_high
        ),
        "Hardy jet can be evaluated at t = -10"
    );

    /*
     * Z is even.
     */
    check(
        arb_overlaps(
            jet_pos.value,
            jet_neg.value
        ),
        "Z(-10) overlaps Z(10)"
    );

    /*
     * Z' is odd.
     */
    check(
        arb_overlaps_negated(
            jet_pos.dt,
            jet_neg.dt
        ),
        "Z'(-10) overlaps -Z'(10)"
    );

    /*
     * Z'' is even.
     */
    check(
        arb_overlaps(
            jet_pos.d2t,
            jet_neg.d2t
        ),
        "Z''(-10) overlaps Z''(10)"
    );


    /*
     * ========================================================
     * Precision consistency
     * ========================================================
     */

    printf(
        "\n"
        "=== Precision consistency ===\n"
    );

    arb_set_si(
        t,
        100
    );

    check(
        hardy_jet_build(
            &jet_low,
            t,
            prec_low
        ),
        "128-bit Hardy jet at t = 100 succeeds"
    );

    check(
        hardy_jet_build(
            &jet_high,
            t,
            prec_high
        ),
        "256-bit Hardy jet at t = 100 succeeds"
    );

    check(
        arb_overlaps(
            jet_low.value,
            jet_high.value
        ),
        "128-bit and 256-bit Z enclosures overlap"
    );

    check(
        arb_overlaps(
            jet_low.dt,
            jet_high.dt
        ),
        "128-bit and 256-bit Z' enclosures overlap"
    );

    check(
        arb_overlaps(
            jet_low.d2t,
            jet_high.d2t
        ),
        "128-bit and 256-bit Z'' enclosures overlap"
    );

    /*
     * ========================================================
     * Interval-width diagnostic
     * ========================================================
     */

    printf(
        "\n"
        "=== Interval-width diagnostic ===\n"
    );

    {
        slong p;

        arb_t eps;
        arb_t upper_ball;

        arb_init(eps);
        arb_init(upper_ball);

        for (p = 0; p <= 20; ++p)
        {
            /*
             * interval = [100, 100 + 2^(-p)].
             */
            arb_set_si(
                t,
                100
            );

            arb_one(
                eps
            );

            arb_mul_2exp_si(
                eps,
                eps,
                -p
            );

            arb_add(
                upper_ball,
                t,
                eps,
                prec_high
            );

            arb_union(
                interval,
                t,
                upper_ball,
                prec_high
            );

            printf(
                "width = 2^(-%-2ld) : %s\n",
                p,
                hardy_jet_build(
                    &jet_interval,
                    interval,
                    prec_high
                )
                    ? "SUCCESS"
                    : "FAIL"
            );
        }

        arb_clear(upper_ball);
        arb_clear(eps);
    }


        /*
     * ========================================================
     * Genuine interval input
     * ========================================================
     */

    printf(
        "\n"
        "=== Genuine interval input ===\n"
    );


    /*
     * A relatively wide interval may be too pessimistic for
     * the direct interval evaluation to remain usable.
     *
     * This is not a loss of rigor.  The adaptive certifier will
     * simply subdivide such an interval.
     */
    arf_set_si(
        lower,
        100
    );

    arf_set_si(
        upper,
        101
    );

    arb_set_interval_arf(
        interval,
        lower,
        upper,
        prec_high
    );

    check(
        !hardy_jet_build(
            &jet_interval,
            interval,
            prec_high
        ),
        "wide interval [100,101] is rejected when the enclosure is unusable"
    );


    /*
     * Now use the narrower interval
     *
     *     [100, 100.5].
     *
     * The width diagnostic above showed that this interval is
     * successfully handled by the direct Hardy evaluator.
     */
    arf_set_si(
        lower,
        100
    );

    arf_set_d(
        upper,
        100.5
    );

    arb_set_interval_arf(
        interval,
        lower,
        upper,
        prec_high
    );

    check(
        hardy_jet_build(
            &jet_interval,
            interval,
            prec_high
        ),
        "Hardy jet succeeds on the interval [100,100.5]"
    );

    check(
        arb_is_finite(jet_interval.value) &&
        arb_is_finite(jet_interval.dt) &&
        arb_is_finite(jet_interval.d2t),
        "successful interval jet components are finite"
    );


    /*
     * Point evaluations at both endpoints must be contained
     * in the corresponding interval enclosures.
     */
    arb_set_si(
        t,
        100
    );

    check(
        hardy_jet_build(
            &jet_low,
            t,
            prec_high
        ),
        "endpoint jet at t = 100 succeeds"
    );

    check(
        arb_contains(
            jet_interval.value,
            jet_low.value
        ),
        "interval Z enclosure contains Z(100)"
    );

    check(
        arb_contains(
            jet_interval.dt,
            jet_low.dt
        ),
        "interval Z' enclosure contains Z'(100)"
    );

    check(
        arb_contains(
            jet_interval.d2t,
            jet_low.d2t
        ),
        "interval Z'' enclosure contains Z''(100)"
    );


    arb_set_d(
        t,
        100.5
    );

    check(
        hardy_jet_build(
            &jet_high,
            t,
            prec_high
        ),
        "endpoint jet at t = 100.5 succeeds"
    );

    check(
        arb_contains(
            jet_interval.value,
            jet_high.value
        ),
        "interval Z enclosure contains Z(100.5)"
    );

    check(
        arb_contains(
            jet_interval.dt,
            jet_high.dt
        ),
        "interval Z' enclosure contains Z'(100.5)"
    );

    check(
        arb_contains(
            jet_interval.d2t,
            jet_high.d2t
        ),
        "interval Z'' enclosure contains Z''(100.5)"
    );


    /*
     * Point evaluations at both endpoints must be contained
     * in the corresponding interval enclosures.
     */
    arb_set_si(
        t,
        100
    );

    check(
        hardy_jet_build(
            &jet_low,
            t,
            prec_high
        ),
        "endpoint jet at t = 100 succeeds"
    );

    check(
        arb_contains(
            jet_interval.value,
            jet_low.value
        ),
        "interval Z enclosure contains Z(100)"
    );

    check(
        arb_contains(
            jet_interval.dt,
            jet_low.dt
        ),
        "interval Z' enclosure contains Z'(100)"
    );

    check(
        arb_contains(
            jet_interval.d2t,
            jet_low.d2t
        ),
        "interval Z'' enclosure contains Z''(100)"
    );


    arb_set_si(
        t,
        101
    );

    check(
        hardy_jet_build(
            &jet_high,
            t,
            prec_high
        ),
        "endpoint jet at t = 101 succeeds"
    );

    check(
        arb_contains(
            jet_interval.value,
            jet_high.value
        ),
        "interval Z enclosure contains Z(101)"
    );

    check(
        arb_contains(
            jet_interval.dt,
            jet_high.dt
        ),
        "interval Z' enclosure contains Z'(101)"
    );

    check(
        arb_contains(
            jet_interval.d2t,
            jet_high.d2t
        ),
        "interval Z'' enclosure contains Z''(101)"
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

    check(
        !hardy_jet_build(
            &jet_low,
            t,
            prec_high
        ),
        "indeterminate input is rejected"
    );

    arb_zero(
        t
    );

    check(
        !hardy_jet_build(
            &jet_low,
            t,
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
            "ALL HARDY JET TESTS PASSED\n"
        );
    }
    else
    {
        printf(
            "%d HARDY JET TEST(S) FAILED\n",
            failures
        );
    }

    printf(
        "============================================================\n"
    );


    arf_clear(upper);
    arf_clear(lower);

    arb_clear(interval);
    arb_clear(minus_t);
    arb_clear(t);

    hardy_jet_clear(&jet_interval);
    hardy_jet_clear(&jet_high);
    hardy_jet_clear(&jet_low);
    hardy_jet_clear(&jet_neg);
    hardy_jet_clear(&jet_pos);
    hardy_jet_clear(&jet_zero);

    return failures == 0 ? 0 : 1;
}