#include "interval_geometry.h"

#include <limits.h>


/*
 * Build the exact dyadic coefficient
 *
 *     x = j / 2^depth
 *
 * as an Arb point.
 *
 * The value itself is exactly dyadic.  Arb is used here only
 * as the arithmetic container needed by the subsequent affine
 * evaluation.
 */
static int
interval_build_dyadic_coefficient(
    arb_t x,
    ulong depth,
    const fmpz_t j
)
{
    if (depth > (ulong) LONG_MAX)
        return 0;

    arb_set_fmpz(
        x,
        j
    );

    arb_mul_2exp_si(
        x,
        x,
        -(slong) depth
    );

    return arb_is_finite(x);
}


/*
 * Rigorously enclose
 *
 *     t = a + (b-a) x,
 *
 * where a and b are the frozen exact ARF endpoints of the
 * root domain and x is an exact dyadic coefficient.
 */
static int
interval_build_affine_endpoint(
    arb_t t,
    const interval_domain_t *domain,
    const arb_t x,
    slong prec
)
{
    arb_t a;
    arb_t b;
    arb_t width;
    arb_t tmp;

    int success = 0;

    arb_init(a);
    arb_init(b);
    arb_init(width);
    arb_init(tmp);

    arb_set_arf(
        a,
        domain->left
    );

    arb_set_arf(
        b,
        domain->right
    );

    /*
     * width = b - a.
     */
    arb_sub(
        width,
        b,
        a,
        prec
    );

    /*
     * tmp = (b-a) x.
     */
    arb_mul(
        tmp,
        width,
        x,
        prec
    );

    /*
     * t = a + (b-a) x.
     */
    arb_add(
        t,
        a,
        tmp,
        prec
    );

    if (!arb_is_finite(t))
        goto cleanup;

    success = 1;

cleanup:

    arb_clear(tmp);
    arb_clear(width);
    arb_clear(b);
    arb_clear(a);

    return success;
}


void
interval_domain_init(
    interval_domain_t *domain
)
{
    arf_init(
        domain->left
    );

    arf_init(
        domain->right
    );

    arf_zero(
        domain->left
    );

    arf_zero(
        domain->right
    );
}


void
interval_domain_clear(
    interval_domain_t *domain
)
{
    arf_clear(
        domain->right
    );

    arf_clear(
        domain->left
    );
}


int
interval_domain_is_valid(
    const interval_domain_t *domain
)
{
    if (!arf_is_finite(domain->left) ||
        !arf_is_finite(domain->right))
    {
        return 0;
    }

    return
        arf_cmp(
            domain->left,
            domain->right
        ) < 0;
}


int
interval_domain_set_arf(
    interval_domain_t *domain,
    const arf_t left,
    const arf_t right
)
{
    if (!arf_is_finite(left) ||
        !arf_is_finite(right))
    {
        return 0;
    }

    if (arf_cmp(left, right) >= 0)
        return 0;

    arf_set(
        domain->left,
        left
    );

    arf_set(
        domain->right,
        right
    );

    return 1;
}


int
interval_domain_set_arb(
    interval_domain_t *domain,
    const arb_t left,
    const arb_t right,
    slong prec
)
{
    arf_t frozen_left;
    arf_t frozen_right;

    int success = 0;

    if (prec <= 0)
        return 0;

    if (!arb_is_finite(left) ||
        !arb_is_finite(right))
    {
        return 0;
    }

    arf_init(frozen_left);
    arf_init(frozen_right);

    /*
     * Freeze outward binary endpoints.
     */
    arb_get_lbound_arf(
        frozen_left,
        left,
        prec
    );

    arb_get_ubound_arf(
        frozen_right,
        right,
        prec
    );

    if (!arf_is_finite(frozen_left) ||
        !arf_is_finite(frozen_right))
    {
        goto cleanup;
    }

    if (arf_cmp(
            frozen_left,
            frozen_right) >= 0)
    {
        goto cleanup;
    }

    arf_set(
        domain->left,
        frozen_left
    );

    arf_set(
        domain->right,
        frozen_right
    );

    success = 1;

cleanup:

    arf_clear(frozen_right);
    arf_clear(frozen_left);

    return success;
}


void
interval_leaf_init(
    interval_leaf_t *leaf
)
{
    leaf->depth = 0;

    fmpz_init(
        leaf->index
    );

    fmpz_zero(
        leaf->index
    );
}


void
interval_leaf_clear(
    interval_leaf_t *leaf
)
{
    fmpz_clear(
        leaf->index
    );
}


void
interval_leaf_set_root(
    interval_leaf_t *leaf
)
{
    leaf->depth = 0;

    fmpz_zero(
        leaf->index
    );
}


int
interval_leaf_is_valid(
    const interval_leaf_t *leaf
)
{
    fmpz_t count;

    int valid;

    /*
     * index must be nonnegative.
     */
    if (fmpz_sgn(leaf->index) < 0)
        return 0;

    fmpz_init(count);

    /*
     * count = 2^depth.
     */
    fmpz_one(count);

    fmpz_mul_2exp(
        count,
        count,
        leaf->depth
    );

    valid =
        fmpz_cmp(
            leaf->index,
            count
        ) < 0;

    fmpz_clear(count);

    return valid;
}


int
interval_leaf_set(
    interval_leaf_t *leaf,
    ulong depth,
    const fmpz_t index
)
{
    interval_leaf_t candidate;

    int success = 0;

    interval_leaf_init(
        &candidate
    );

    candidate.depth =
        depth;

    fmpz_set(
        candidate.index,
        index
    );

    if (!interval_leaf_is_valid(
            &candidate))
    {
        goto cleanup;
    }

    leaf->depth =
        candidate.depth;

    fmpz_set(
        leaf->index,
        candidate.index
    );

    success = 1;

cleanup:

    interval_leaf_clear(
        &candidate
    );

    return success;
}


int
interval_leaf_bisect(
    interval_leaf_t *left,
    interval_leaf_t *right,
    const interval_leaf_t *parent
)
{
    if (!interval_leaf_is_valid(parent))
        return 0;

    if (parent->depth == ULONG_MAX)
        return 0;

    left->depth =
        parent->depth + 1;

    right->depth =
        parent->depth + 1;

    /*
     * left.index = 2 * parent.index.
     */
    fmpz_mul_2exp(
        left->index,
        parent->index,
        1
    );

    /*
     * right.index = 2 * parent.index + 1.
     */
    fmpz_set(
        right->index,
        left->index
    );

    fmpz_add_ui(
        right->index,
        right->index,
        1
    );

    return
        interval_leaf_is_valid(left) &&
        interval_leaf_is_valid(right);
}


int
interval_leaf_endpoints(
    arb_t t_left,
    arb_t t_right,
    const interval_domain_t *domain,
    const interval_leaf_t *leaf,
    slong prec
)
{
    fmpz_t right_index;

    arb_t x_left;
    arb_t x_right;

    int success = 0;

    if (!interval_domain_is_valid(domain))
        return 0;

    if (!interval_leaf_is_valid(leaf))
        return 0;

    if (prec <= 0)
        return 0;

    if (leaf->depth > (ulong) LONG_MAX)
        return 0;

    fmpz_init(right_index);

    arb_init(x_left);
    arb_init(x_right);

    /*
     * right_index = index + 1.
     */
    fmpz_add_ui(
        right_index,
        leaf->index,
        1
    );

    if (!interval_build_dyadic_coefficient(
            x_left,
            leaf->depth,
            leaf->index))
    {
        goto cleanup;
    }

    if (!interval_build_dyadic_coefficient(
            x_right,
            leaf->depth,
            right_index))
    {
        goto cleanup;
    }

    if (!interval_build_affine_endpoint(
            t_left,
            domain,
            x_left,
            prec))
    {
        goto cleanup;
    }

    if (!interval_build_affine_endpoint(
            t_right,
            domain,
            x_right,
            prec))
    {
        goto cleanup;
    }

    success = 1;

cleanup:

    arb_clear(x_right);
    arb_clear(x_left);

    fmpz_clear(right_index);

    return success;
}


int
interval_leaf_arb(
    arb_t interval,
    const interval_domain_t *domain,
    const interval_leaf_t *leaf,
    slong prec
)
{
    arb_t left;
    arb_t right;

    arf_t lower;
    arf_t upper;

    int success = 0;

    if (prec <= 0)
        return 0;

    arb_init(left);
    arb_init(right);

    arf_init(lower);
    arf_init(upper);

    if (!interval_leaf_endpoints(
            left,
            right,
            domain,
            leaf,
            prec))
    {
        goto cleanup;
    }

    /*
     * Extract a downward-rounded lower endpoint and an
     * upward-rounded upper endpoint.
     */
    arb_get_lbound_arf(
        lower,
        left,
        prec
    );

    arb_get_ubound_arf(
        upper,
        right,
        prec
    );

    if (!arf_is_finite(lower) ||
        !arf_is_finite(upper))
    {
        goto cleanup;
    }

    if (arf_cmp(lower, upper) > 0)
        goto cleanup;

    /*
     * Build a single rigorous Arb enclosure of the complete
     * mathematical leaf.
     */
    arb_set_interval_arf(
        interval,
        lower,
        upper,
        prec
    );

    if (!arb_is_finite(interval))
        goto cleanup;

    success = 1;

cleanup:

    arf_clear(upper);
    arf_clear(lower);

    arb_clear(right);
    arb_clear(left);

    return success;
}