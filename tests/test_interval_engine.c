#include <stdio.h>

#include "arb.h"
#include "arf.h"

#include "interval_engine.h"


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

    interval_domain_t domain;
    interval_engine_stats_t stats;

    arf_t left;
    arf_t right;

    arb_t kappa;


    printf(
        "============================================================\n"
        "INTERVAL ENGINE TESTS\n"
        "Adaptive Hardy-only dyadic certification\n"
        "precision = %ld bits\n"
        "============================================================\n",
        prec
    );


    interval_domain_init(
        &domain
    );

    interval_engine_stats_init(
        &stats
    );

    arf_init(left);
    arf_init(right);

    arb_init(kappa);


    /*
     * ========================================================
     * Small certifiable domain
     * ========================================================
     */

    printf(
        "\n"
        "=== Small certifiable domain ===\n"
    );

    arf_set_si(
        left,
        50
    );

    arf_set_d(
        right,
        50.01
    );

    check(
        interval_domain_set_arf(
            &domain,
            left,
            right
        ),
        "domain [50,50.01] is valid"
    );

    arb_zero(
        kappa
    );

    check(
        interval_engine_certify_domain(
            &stats,
            &domain,
            kappa,
            prec,
            20,
            0,
            NULL,
            NULL
        ),
        "adaptive engine certifies [50,50.01]"
    );

    check(
        stats.status ==
            INTERVAL_ENGINE_STATUS_SUCCESS,
        "successful run reports SUCCESS"
    );

    check(
        stats.certified_leaves > 0,
        "at least one terminal leaf was certified"
    );

    check(
        stats.has_minimum,
        "successful run records a minimum lower bound"
    );

    check(
        arf_sgn(
            stats.minimum_certified_lower
        ) > 0,
        "minimum certified lower bound is positive"
    );

    /*
     * Every evaluated node in a successful binary tree is
     * either certified or bisected.
     */
    check(
        stats.evaluations ==
            stats.certified_leaves +
            stats.bisected_leaves,
        "every successful evaluation is classified"
    );

    /*
     * A finite full binary tree satisfies
     *
     *     leaves = internal_nodes + 1.
     */
    check(
        stats.certified_leaves ==
            stats.bisected_leaves + 1,
        "successful adaptive tree satisfies leaves = bisections + 1"
    );


    /*
     * ========================================================
     * Deliberate depth limit
     * ========================================================
     */

    printf(
        "\n"
        "=== Deliberate depth limit ===\n"
    );

    arf_set_si(
        left,
        100
    );

    arf_set_si(
        right,
        101
    );

    check(
        interval_domain_set_arf(
            &domain,
            left,
            right
        ),
        "domain [100,101] is valid"
    );

    /*
     * The direct Hardy interval evaluation is known to be too
     * wide on [100,101].  With max_depth = 0 the engine is not
     * allowed to subdivide.
     */
    check(
        !interval_engine_certify_domain(
            &stats,
            &domain,
            kappa,
            prec,
            0,
            0,
            NULL,
            NULL
        ),
        "wide domain fails when subdivision is forbidden"
    );

    check(
        stats.status ==
            INTERVAL_ENGINE_STATUS_DEPTH_LIMIT,
        "failure is correctly classified as DEPTH_LIMIT"
    );

    check(
        stats.has_unresolved_leaf,
        "depth-limit run records the unresolved leaf"
    );

    check(
        stats.unresolved_leaf.depth == 0,
        "unresolved leaf is the root"
    );


    /*
     * ========================================================
     * Invalid domain
     * ========================================================
     */

    printf(
        "\n"
        "=== Invalid domain ===\n"
    );

    arf_set_si(
        left,
        10
    );

    arf_set_si(
        right,
        10
    );

    /*
     * interval_domain_set_arf must reject this domain, leaving
     * us free to invalidate the stored domain explicitly.
     */
    arf_set(
        domain.left,
        left
    );

    arf_set(
        domain.right,
        right
    );

    check(
        !interval_engine_certify_domain(
            &stats,
            &domain,
            kappa,
            prec,
            10,
            0,
            NULL,
            NULL
        ),
        "zero-width frozen domain is rejected"
    );

    check(
        stats.status ==
            INTERVAL_ENGINE_STATUS_INVALID,
        "invalid domain leaves status INVALID"
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
            "ALL INTERVAL ENGINE TESTS PASSED\n"
        );
    }
    else
    {
        printf(
            "%d INTERVAL ENGINE TEST(S) FAILED\n",
            failures
        );
    }

    printf(
        "============================================================\n"
    );


    arb_clear(kappa);

    arf_clear(right);
    arf_clear(left);

    interval_engine_stats_clear(
        &stats
    );

    interval_domain_clear(
        &domain
    );

    return failures == 0 ? 0 : 1;
}