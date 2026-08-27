#include <stdio.h>

#include "arb.h"
#include "arf.h"
#include "fmpz.h"

#include "interval_geometry.h"


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
arb_contains_si_value(
    const arb_t x,
    slong value
)
{
    arb_t expected;

    int result;

    arb_init(expected);

    arb_set_si(
        expected,
        value
    );

    result =
        arb_contains(
            x,
            expected
        );

    arb_clear(expected);

    return result;
}


int
main(void)
{
    const slong prec = 256;

    interval_domain_t domain;

    interval_leaf_t root;
    interval_leaf_t leaf;
    interval_leaf_t left;
    interval_leaf_t right;
    interval_leaf_t invalid_leaf;

    arb_t t_left;
    arb_t t_right;
    arb_t interval;

    arb_t left_ball;
    arb_t right_ball;

    arf_t a;
    arf_t b;

    arf_t lo;
    arf_t hi;

    fmpz_t index;

    printf(
        "============================================================\n"
        "INTERVAL GEOMETRY TESTS\n"
        "Exact dyadic geometry on an arbitrary frozen interval\n"
        "precision = %ld bits\n"
        "============================================================\n",
        prec
    );

    interval_domain_init(
        &domain
    );

    interval_leaf_init(
        &root
    );

    interval_leaf_init(
        &leaf
    );

    interval_leaf_init(
        &left
    );

    interval_leaf_init(
        &right
    );

    interval_leaf_init(
        &invalid_leaf
    );

    arb_init(t_left);
    arb_init(t_right);
    arb_init(interval);

    arb_init(left_ball);
    arb_init(right_ball);

    arf_init(a);
    arf_init(b);

    arf_init(lo);
    arf_init(hi);

    fmpz_init(index);


    /*
     * ========================================================
     * Domain validation
     * ========================================================
     */

    printf(
        "\n"
        "=== Frozen domain validation ===\n"
    );

    check(
        !interval_domain_is_valid(
            &domain
        ),
        "freshly initialized zero-width domain is invalid"
    );

    arf_set_si(
        a,
        2
    );

    arf_set_si(
        b,
        10
    );

    check(
        interval_domain_set_arf(
            &domain,
            a,
            b
        ),
        "exact ARF domain [2,10] is accepted"
    );

    check(
        interval_domain_is_valid(
            &domain
        ),
        "stored domain is valid"
    );

    check(
        arf_cmp(
            domain.left,
            a
        ) == 0,
        "stored left endpoint is exactly 2"
    );

    check(
        arf_cmp(
            domain.right,
            b
        ) == 0,
        "stored right endpoint is exactly 10"
    );

    check(
        !interval_domain_set_arf(
            &domain,
            b,
            a
        ),
        "reversed domain is rejected"
    );

    check(
        !interval_domain_set_arf(
            &domain,
            a,
            a
        ),
        "zero-width domain is rejected"
    );


    /*
     * ========================================================
     * Arb -> frozen outward domain
     * ========================================================
     */

    printf(
        "\n"
        "=== Outward freezing from Arb endpoint enclosures ===\n"
    );

    /*
     * left_ball = [2,3].
     */
    arf_set_si(
        lo,
        2
    );

    arf_set_si(
        hi,
        3
    );

    arb_set_interval_arf(
        left_ball,
        lo,
        hi,
        prec
    );

    /*
     * right_ball = [9,10].
     */
    arf_set_si(
        lo,
        9
    );

    arf_set_si(
        hi,
        10
    );

    arb_set_interval_arf(
        right_ball,
        lo,
        hi,
        prec
    );

    check(
        interval_domain_set_arb(
            &domain,
            left_ball,
            right_ball,
            prec
        ),
        "Arb endpoint enclosures can be frozen"
    );

    arf_set_si(
        a,
        2
    );

    arf_set_si(
        b,
        10
    );

    check(
        arf_cmp(
            domain.left,
            a
        ) <= 0,
        "frozen left endpoint is an outward lower bound"
    );

    check(
        arf_cmp(
            domain.right,
            b
        ) >= 0,
        "frozen right endpoint is an outward upper bound"
    );


    /*
     * Restore the exact simple domain [2,10].
     */
    check(
        interval_domain_set_arf(
            &domain,
            a,
            b
        ),
        "exact domain [2,10] restored"
    );


    /*
     * ========================================================
     * Root leaf
     * ========================================================
     */

    printf(
        "\n"
        "=== Root leaf ===\n"
    );

    interval_leaf_set_root(
        &root
    );

    check(
        interval_leaf_is_valid(
            &root
        ),
        "root leaf (0,0) is valid"
    );

    check(
        root.depth == 0,
        "root depth is zero"
    );

    check(
        fmpz_is_zero(
            root.index
        ),
        "root index is zero"
    );

    check(
        interval_leaf_endpoints(
            t_left,
            t_right,
            &domain,
            &root,
            prec
        ),
        "root endpoints can be evaluated rigorously"
    );

    check(
        arb_contains_si_value(
            t_left,
            2
        ),
        "root left endpoint contains 2"
    );

    check(
        arb_contains_si_value(
            t_right,
            10
        ),
        "root right endpoint contains 10"
    );


    /*
     * ========================================================
     * Arbitrary leaf
     * ========================================================
     */

    printf(
        "\n"
        "=== Arbitrary dyadic leaf ===\n"
    );

    /*
     * On [2,10], leaf (d,k)=(2,1) is
     *
     *     [4,6].
     */
    fmpz_set_ui(
        index,
        1
    );

    check(
        interval_leaf_set(
            &leaf,
            2,
            index
        ),
        "leaf (2,1) is accepted"
    );

    check(
        interval_leaf_is_valid(
            &leaf
        ),
        "leaf (2,1) is valid"
    );

    check(
        interval_leaf_endpoints(
            t_left,
            t_right,
            &domain,
            &leaf,
            prec
        ),
        "leaf (2,1) endpoints can be evaluated"
    );

    check(
        arb_contains_si_value(
            t_left,
            4
        ),
        "leaf (2,1) left endpoint contains 4"
    );

    check(
        arb_contains_si_value(
            t_right,
            6
        ),
        "leaf (2,1) right endpoint contains 6"
    );

    check(
        interval_leaf_arb(
            interval,
            &domain,
            &leaf,
            prec
        ),
        "leaf (2,1) can be converted to one Arb enclosure"
    );

    check(
        arb_contains_si_value(
            interval,
            4
        ),
        "leaf Arb enclosure contains its left endpoint"
    );

    check(
        arb_contains_si_value(
            interval,
            5
        ),
        "leaf Arb enclosure contains its interior midpoint"
    );

    check(
        arb_contains_si_value(
            interval,
            6
        ),
        "leaf Arb enclosure contains its right endpoint"
    );


    /*
     * ========================================================
     * Exact symbolic bisection
     * ========================================================
     */

    printf(
        "\n"
        "=== Exact symbolic bisection ===\n"
    );

    check(
        interval_leaf_bisect(
            &left,
            &right,
            &leaf
        ),
        "leaf (2,1) bisects successfully"
    );

    check(
        left.depth == 3,
        "left child has depth 3"
    );

    check(
        right.depth == 3,
        "right child has depth 3"
    );

    check(
        fmpz_cmp_ui(
            left.index,
            2
        ) == 0,
        "left child index is 2"
    );

    check(
        fmpz_cmp_ui(
            right.index,
            3
        ) == 0,
        "right child index is 3"
    );

    /*
     * left child = [4,5].
     */
    check(
        interval_leaf_endpoints(
            t_left,
            t_right,
            &domain,
            &left,
            prec
        ),
        "left child endpoints can be evaluated"
    );

    check(
        arb_contains_si_value(
            t_left,
            4
        ),
        "left child begins at 4"
    );

    check(
        arb_contains_si_value(
            t_right,
            5
        ),
        "left child ends at common endpoint 5"
    );

    /*
     * right child = [5,6].
     */
    check(
        interval_leaf_endpoints(
            t_left,
            t_right,
            &domain,
            &right,
            prec
        ),
        "right child endpoints can be evaluated"
    );

    check(
        arb_contains_si_value(
            t_left,
            5
        ),
        "right child begins at common endpoint 5"
    );

    check(
        arb_contains_si_value(
            t_right,
            6
        ),
        "right child ends at 6"
    );


    /*
     * ========================================================
     * Root partition
     * ========================================================
     */

    printf(
        "\n"
        "=== Root partition ===\n"
    );

    check(
        interval_leaf_bisect(
            &left,
            &right,
            &root
        ),
        "root bisects successfully"
    );

    /*
     * On [2,10]:
     *
     *     left  = [2,6],
     *     right = [6,10].
     */
    check(
        interval_leaf_endpoints(
            t_left,
            t_right,
            &domain,
            &left,
            prec
        ),
        "left half of root evaluates"
    );

    check(
        arb_contains_si_value(
            t_left,
            2
        ),
        "left half begins at root left endpoint"
    );

    check(
        arb_contains_si_value(
            t_right,
            6
        ),
        "left half ends at exact midpoint"
    );

    check(
        interval_leaf_endpoints(
            t_left,
            t_right,
            &domain,
            &right,
            prec
        ),
        "right half of root evaluates"
    );

    check(
        arb_contains_si_value(
            t_left,
            6
        ),
        "right half begins at exact midpoint"
    );

    check(
        arb_contains_si_value(
            t_right,
            10
        ),
        "right half ends at root right endpoint"
    );


    /*
     * ========================================================
     * Invalid symbolic leaves
     * ========================================================
     */

    printf(
        "\n"
        "=== Invalid symbolic leaves ===\n"
    );

    /*
     * depth 2 permits exactly indices 0,1,2,3.
     */
    invalid_leaf.depth = 2;

    fmpz_set_ui(
        invalid_leaf.index,
        4
    );

    check(
        !interval_leaf_is_valid(
            &invalid_leaf
        ),
        "index 4 is invalid at depth 2"
    );

    fmpz_set_si(
        invalid_leaf.index,
        -1
    );

    check(
        !interval_leaf_is_valid(
            &invalid_leaf
        ),
        "negative leaf index is rejected"
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
            "ALL INTERVAL GEOMETRY TESTS PASSED\n"
        );
    }
    else
    {
        printf(
            "%d INTERVAL GEOMETRY TEST(S) FAILED\n",
            failures
        );
    }

    printf(
        "============================================================\n"
    );


    fmpz_clear(index);

    arf_clear(hi);
    arf_clear(lo);

    arf_clear(b);
    arf_clear(a);

    arb_clear(right_ball);
    arb_clear(left_ball);

    arb_clear(interval);
    arb_clear(t_right);
    arb_clear(t_left);

    interval_leaf_clear(
        &invalid_leaf
    );

    interval_leaf_clear(
        &right
    );

    interval_leaf_clear(
        &left
    );

    interval_leaf_clear(
        &leaf
    );

    interval_leaf_clear(
        &root
    );

    interval_domain_clear(
        &domain
    );

    return failures == 0 ? 0 : 1;
}