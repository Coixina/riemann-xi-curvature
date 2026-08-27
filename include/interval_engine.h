#ifndef INTERVAL_ENGINE_H
#define INTERVAL_ENGINE_H

#include "arb.h"
#include "arf.h"

#include "interval_geometry.h"


/*
 * Final status of one adaptive certification run.
 */
typedef enum
{
    INTERVAL_ENGINE_STATUS_INVALID = 0,

    /*
     * Every terminal leaf was rigorously certified.
     */
    INTERVAL_ENGINE_STATUS_SUCCESS,

    /*
     * At least one unresolved leaf reached max_depth.
     *
     * This is not a mathematical counterexample.  It means
     * only that the current subdivision depth was insufficient
     * to complete the certificate.
     */
    INTERVAL_ENGINE_STATUS_DEPTH_LIMIT,

    /*
     * An unexpected internal failure occurred.
     */
    INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE

} interval_engine_status_t;


/*
 * Statistics and diagnostic information for one adaptive run.
 *
 * The structure must be initialized with
 * interval_engine_stats_init() before use.
 */
typedef struct
{
    interval_engine_status_t status;

    /*
     * Number of attempted leaf certifications.
     */
    ulong evaluations;

    /*
     * Number of terminal certified leaves.
     */
    ulong certified_leaves;

    /*
    * Fraction of the root interval already covered by
    * rigorously certified terminal leaves.
    *
    * Diagnostic only; it has no role in the certificate.
    */
    long double certified_fraction;

    /*
     * Number of leaves that were bisected.
     */
    ulong bisected_leaves;

    /*
     * Largest dyadic depth actually evaluated.
     */
    ulong max_depth_reached;

    /*
     * Smallest rigorous certified lower bound among all
     * accepted leaves.
     *
     * Valid only when has_minimum != 0.
     */
    arf_t minimum_certified_lower;

    interval_leaf_t minimum_leaf;

    int has_minimum;

    /*
     * Leaf responsible for unsuccessful termination.
     *
     * Valid only when has_unresolved_leaf != 0.
     */
    interval_leaf_t unresolved_leaf;

    int has_unresolved_leaf;

} interval_engine_stats_t;

/*
 * Optional sink called whenever a terminal leaf has been
 * rigorously certified.
 *
 * The interval engine knows nothing about the consumer of the
 * leaf.  A caller may use this callback, for example, to write
 * certified leaves to a certificate file.
 *
 * Returning 1 means that the leaf was accepted by the sink.
 * Returning 0 aborts the adaptive run with
 * INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE.
 */
typedef int (*interval_engine_leaf_sink_t)(
    const interval_leaf_t *leaf,
    void *user_data
);

/*
 * Initialize and clear an interval-engine statistics object.
 */
void interval_engine_stats_init(
    interval_engine_stats_t *stats
);

void interval_engine_stats_clear(
    interval_engine_stats_t *stats
);


/*
 * Certify one exact dyadic subtree of a frozen domain.
 *
 * domain:
 *     frozen root interval.
 *
 * root:
 *     exact symbolic dyadic leaf defining the subtree.
 *
 * kappa:
 *     curvature parameter used by bilinear_certifier.
 *
 * max_depth:
 *     absolute dyadic depth.  It is not a number of
 *     additional subdivision levels.
 *
 * progress_every:
 *     print one diagnostic progress line after this many
 *     leaf evaluations.  Set to zero to disable progress
 *     output.
 *
 * The engine uses iterative depth-first search:
 *
 *     CERTIFIED
 *         -> accept the leaf.
 *
 *     INCONCLUSIVE
 *         -> bisect.
 *
 *     unusable Hardy enclosure
 *         -> bisect.
 *
 *     unresolved at max_depth
 *         -> terminate with DEPTH_LIMIT.
 *
 * Returns 1 exactly when the complete subtree is rigorously
 * certified.
 */
int interval_engine_certify_subtree(
    interval_engine_stats_t *stats,
    const interval_domain_t *domain,
    const interval_leaf_t *root,
    const arb_t kappa,
    slong prec,
    ulong max_depth,
    ulong progress_every,
    interval_engine_leaf_sink_t leaf_sink,
    void *leaf_sink_data
);


/*
 * Certify the complete frozen domain.
 *
 * This is equivalent to certifying the symbolic root leaf
 *
 *     (depth,index) = (0,0).
 *
 * leaf_sink:
 *     optional callback invoked for every certified terminal
 *     leaf, in left-to-right order.  May be NULL.
 *
 * leaf_sink_data:
 *     opaque pointer passed unchanged to leaf_sink.
 */
int interval_engine_certify_domain(
    interval_engine_stats_t *stats,
    const interval_domain_t *domain,
    const arb_t kappa,
    slong prec,
    ulong max_depth,
    ulong progress_every,
    interval_engine_leaf_sink_t leaf_sink,
    void *leaf_sink_data
);


#endif