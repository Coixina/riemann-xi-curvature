#ifndef INTERVAL_GEOMETRY_H
#define INTERVAL_GEOMETRY_H

#include "arb.h"
#include "arf.h"
#include "fmpz.h"


/*
 * Frozen root interval for adaptive certification.
 *
 * The logical certification domain is
 *
 *     [left, right],
 *
 * where left and right are finite ARF numbers with
 *
 *     left < right.
 *
 * ARF values are binary floating-point numbers and therefore
 * represent exact dyadic rationals.  Once the domain has been
 * initialized, these endpoints remain fixed throughout the
 * adaptive run.
 *
 * A mathematical target interval [a,b] may be enclosed by
 *
 *     left  <= a,
 *     b <= right.
 *
 * Certification of the frozen root interval then certifies
 * the original target interval as well.
 */
typedef struct
{
    arf_t left;
    arf_t right;

} interval_domain_t;


/*
 * Exact dyadic leaf inside one frozen root interval.
 *
 * The leaf
 *
 *     (depth, index)
 *
 * represents the exact affine subinterval
 *
 *     u in [
 *         index / 2^depth,
 *         (index + 1) / 2^depth
 *     ]
 *
 * under
 *
 *     t(u)
 *       =
 *     left + (right - left) u.
 *
 * Validity condition:
 *
 *     0 <= index < 2^depth.
 *
 * The logical geometry is defined by this symbolic dyadic
 * identifier, not by an Arb midpoint/radius representation.
 */
typedef struct
{
    ulong depth;
    fmpz_t index;

} interval_leaf_t;


/*
 * Initialize and clear a frozen interval domain.
 */
void interval_domain_init(
    interval_domain_t *domain
);

void interval_domain_clear(
    interval_domain_t *domain
);


/*
 * Return 1 exactly when
 *
 *     left < right
 *
 * and both frozen endpoints are finite.
 */
int interval_domain_is_valid(
    const interval_domain_t *domain
);


/*
 * Freeze a domain directly from exact finite ARF endpoints.
 *
 * Requires
 *
 *     left < right.
 *
 * Returns 1 on success.
 */
int interval_domain_set_arf(
    interval_domain_t *domain,
    const arf_t left,
    const arf_t right
);


/*
 * Freeze a domain from two rigorous Arb endpoint enclosures.
 *
 * The stored endpoints are
 *
 *     domain->left
 *         = downward lower bound of left,
 *
 *     domain->right
 *         = upward upper bound of right.
 *
 * Thus, if the Arb arguments rigorously contain mathematical
 * endpoints a and b, the frozen domain contains [a,b].
 *
 * prec controls extraction of the outward ARF bounds.
 *
 * Returns 1 on success.
 */
int interval_domain_set_arb(
    interval_domain_t *domain,
    const arb_t left,
    const arb_t right,
    slong prec
);


/*
 * Initialize and clear a symbolic leaf.
 */
void interval_leaf_init(
    interval_leaf_t *leaf
);

void interval_leaf_clear(
    interval_leaf_t *leaf
);


/*
 * Set the root leaf
 *
 *     (0,0).
 */
void interval_leaf_set_root(
    interval_leaf_t *leaf
);


/*
 * Set an arbitrary dyadic leaf.
 *
 * Validity condition:
 *
 *     0 <= index < 2^depth.
 *
 * Returns 1 on success.
 */
int interval_leaf_set(
    interval_leaf_t *leaf,
    ulong depth,
    const fmpz_t index
);


/*
 * Return 1 exactly when the symbolic dyadic identifier is
 * valid.
 */
int interval_leaf_is_valid(
    const interval_leaf_t *leaf
);


/*
 * Exact combinatorial bisection:
 *
 *     (d,k)
 *
 * becomes
 *
 *     left  = (d+1, 2k),
 *     right = (d+1, 2k+1).
 *
 * Returns 1 on success.
 */
int interval_leaf_bisect(
    interval_leaf_t *left,
    interval_leaf_t *right,
    const interval_leaf_t *parent
);


/*
 * Build rigorous Arb enclosures of the exact mathematical
 * endpoints represented by one leaf in one frozen domain.
 *
 * If
 *
 *     a = domain->left,
 *     b = domain->right,
 *
 * then the exact endpoints are
 *
 *     t_left
 *       =
 *     a + (b-a) index / 2^depth,
 *
 *     t_right
 *       =
 *     a + (b-a) (index+1) / 2^depth.
 *
 * Arb arithmetic is used only to enclose these exact affine
 * quantities rigorously.
 */
int interval_leaf_endpoints(
    arb_t t_left,
    arb_t t_right,
    const interval_domain_t *domain,
    const interval_leaf_t *leaf,
    slong prec
);


/*
 * Build one Arb ball containing the complete exact leaf.
 *
 * This ball is an evaluation enclosure only.  It is not the
 * logical definition of the leaf.
 */
int interval_leaf_arb(
    arb_t interval,
    const interval_domain_t *domain,
    const interval_leaf_t *leaf,
    slong prec
);


#endif