#include <stdio.h>

#include "arb.h"
#include "fmpz.h"

#include "bilinear_certifier.h"
#include "interval_geometry.h"


#define LARGE_T_N 2
#define KAPPA_VALUE -100


static int
build_finite_domain(
    interval_domain_t *domain,
    ulong large_t_N,
    slong prec
)
{
    arb_t pi;
    arb_t left;
    arb_t right;

    int success = 0;

    arb_init(pi);
    arb_init(left);
    arb_init(right);

    arb_const_pi(
        pi,
        prec
    );

    /*
     * left = 2*pi
     */
    arb_mul_2exp_si(
        left,
        pi,
        1
    );

    /*
     * right = 2*pi*N^2
     */
    arb_set(
        right,
        left
    );

    arb_mul_ui(
        right,
        right,
        large_t_N,
        prec
    );

    arb_mul_ui(
        right,
        right,
        large_t_N,
        prec
    );

    if (!interval_domain_set_arb(
            domain,
            left,
            right,
            prec))
    {
        goto cleanup;
    }

    success = 1;

cleanup:

    arb_clear(right);
    arb_clear(left);
    arb_clear(pi);

    return success;
}


/*
 * Recursively test whether one stored leaf can be certified
 * directly or after exact dyadic refinement.
 *
 * extra_depth = number of levels below the stored leaf.
 */
static int
certify_with_refinement(
    const interval_domain_t *domain,
    const interval_leaf_t *leaf,
    const arb_t kappa,
    slong prec,
    ulong extra_depth,
    ulong max_extra_depth,
    ulong *evaluations,
    ulong *certified_leaves,
    ulong *subdivisions,
    ulong *max_extra_used
)
{
    arb_t interval;

    bilinear_certificate_t certificate;

    interval_leaf_t left;
    interval_leaf_t right;

    int completed;
    int result = 0;


    arb_init(interval);

    bilinear_certificate_init(
        &certificate
    );

    interval_leaf_init(
        &left
    );

    interval_leaf_init(
        &right
    );


    (*evaluations)++;


    if (!interval_leaf_arb(
            interval,
            domain,
            leaf,
            prec))
    {
        goto subdivide;
    }


    completed =
        bilinear_certify(
            &certificate,
            interval,
            kappa,
            prec
        );


    if (completed &&
        certificate.status ==
            BILINEAR_STATUS_CERTIFIED)
    {
        (*certified_leaves)++;

        if (extra_depth >
            *max_extra_used)
        {
            *max_extra_used =
                extra_depth;
        }

        result = 1;

        goto cleanup;
    }


    /*
     * Any rigorous non-certified result, or an unusable Hardy
     * enclosure, is handled by exact dyadic refinement.
     */
subdivide:

    if (extra_depth >=
        max_extra_depth)
    {
        printf(
            "    unresolved at extra depth %lu"
            " : leaf=(d=%lu,k=",
            extra_depth,
            leaf->depth
        );

        fmpz_print(
            leaf->index
        );

        printf(
            ")\n"
        );

        goto cleanup;
    }


    if (!interval_leaf_bisect(
            &left,
            &right,
            leaf))
    {
        goto cleanup;
    }


    (*subdivisions)++;


    if (!certify_with_refinement(
            domain,
            &left,
            kappa,
            prec,
            extra_depth + 1,
            max_extra_depth,
            evaluations,
            certified_leaves,
            subdivisions,
            max_extra_used))
    {
        goto cleanup;
    }


    if (!certify_with_refinement(
            domain,
            &right,
            kappa,
            prec,
            extra_depth + 1,
            max_extra_depth,
            evaluations,
            certified_leaves,
            subdivisions,
            max_extra_used))
    {
        goto cleanup;
    }


    result = 1;


cleanup:

    interval_leaf_clear(
        &right
    );

    interval_leaf_clear(
        &left
    );

    bilinear_certificate_clear(
        &certificate
    );

    arb_clear(
        interval
    );

    return result;
}


static int
run_precision_test(
    slong prec,
    ulong max_extra_depth
)
{
    interval_domain_t domain;

    interval_leaf_t leaf;

    arb_t kappa;

    fmpz_t index;

    ulong evaluations = 0;
    ulong certified_leaves = 0;
    ulong subdivisions = 0;
    ulong max_extra_used = 0;

    int success;


    interval_domain_init(
        &domain
    );

    interval_leaf_init(
        &leaf
    );

    arb_init(
        kappa
    );

    fmpz_init(
        index
    );


    if (!build_finite_domain(
            &domain,
            LARGE_T_N,
            prec))
    {
        printf(
            "precision %ld: could not build domain\n",
            prec
        );

        success = 0;
        goto cleanup;
    }


    /*
     * Stored problematic leaf:
     *
     *     (d,k) = (11,0)
     */
    fmpz_zero(
        index
    );

    if (!interval_leaf_set(
            &leaf,
            11,
            index))
    {
        printf(
            "precision %ld: could not build leaf\n",
            prec
        );

        success = 0;
        goto cleanup;
    }


    arb_set_si(
        kappa,
        KAPPA_VALUE
    );


    printf(
        "\n"
        "------------------------------------------------------------\n"
        "precision       = %ld bits\n"
        "stored leaf     = (d=11,k=0)\n"
        "max extra depth = %lu\n"
        "------------------------------------------------------------\n",
        prec,
        max_extra_depth
    );


    success =
        certify_with_refinement(
            &domain,
            &leaf,
            kappa,
            prec,
            0,
            max_extra_depth,
            &evaluations,
            &certified_leaves,
            &subdivisions,
            &max_extra_used
        );


    printf(
        "result           = %s\n"
        "evaluations      = %lu\n"
        "certified leaves = %lu\n"
        "subdivisions     = %lu\n"
        "max extra used   = %lu\n",
        success ? "SUCCESS" : "FAIL",
        evaluations,
        certified_leaves,
        subdivisions,
        max_extra_used
    );


cleanup:

    fmpz_clear(
        index
    );

    arb_clear(
        kappa
    );

    interval_leaf_clear(
        &leaf
    );

    interval_domain_clear(
        &domain
    );

    return success;
}


int
main(void)
{
    const ulong max_extra_depth = 12;

    int ok256;
    int ok320;
    int ok384;
    int ok512;


    printf(
        "============================================================\n"
        "VERIFIER LOCAL REFINEMENT TEST\n"
        "Stored leaf (d=11,k=0), domain [2*pi,8*pi], kappa=-100\n"
        "============================================================\n"
    );


    ok256 =
        run_precision_test(
            256,
            max_extra_depth
        );

    ok320 =
        run_precision_test(
            320,
            max_extra_depth
        );

    ok384 =
        run_precision_test(
            384,
            max_extra_depth
        );

    ok512 =
        run_precision_test(
            512,
            max_extra_depth
        );


    printf(
        "\n"
        "============================================================\n"
        "SUMMARY\n"
        "256 bits : %s\n"
        "320 bits : %s\n"
        "384 bits : %s\n"
        "512 bits : %s\n"
        "============================================================\n",
        ok256 ? "PASS" : "FAIL",
        ok320 ? "PASS" : "FAIL",
        ok384 ? "PASS" : "FAIL",
        ok512 ? "PASS" : "FAIL"
    );


    return
        (ok256 &&
         ok320 &&
         ok384 &&
         ok512)
        ? 0
        : 1;
}