#include "interval_engine.h"

#include "bilinear_certifier.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


/*
 * Copy one symbolic interval leaf.
 */
static void
interval_engine_leaf_set(
    interval_leaf_t *destination,
    const interval_leaf_t *source
)
{
    destination->depth =
        source->depth;

    fmpz_set(
        destination->index,
        source->index
    );
}


/*
 * Reset statistics while preserving initialized FLINT objects.
 */
static void
interval_engine_stats_reset(
    interval_engine_stats_t *stats
)
{
    stats->status =
        INTERVAL_ENGINE_STATUS_INVALID;

    stats->evaluations = 0;
    stats->certified_leaves = 0;
    stats->certified_fraction = 0.0L;
    stats->bisected_leaves = 0;

    stats->max_depth_reached = 0;

    arf_zero(
        stats->minimum_certified_lower
    );

    stats->has_minimum = 0;
    stats->has_unresolved_leaf = 0;
}


/*
 * Update the minimum rigorous certified lower bound.
 */
static void
interval_engine_update_minimum(
    interval_engine_stats_t *stats,
    const interval_leaf_t *leaf,
    const arf_t certified_lower
)
{
    if (!stats->has_minimum ||
        arf_cmp(
            certified_lower,
            stats->minimum_certified_lower
        ) < 0)
    {
        arf_set(
            stats->minimum_certified_lower,
            certified_lower
        );

        interval_engine_leaf_set(
            &stats->minimum_leaf,
            leaf
        );

        stats->has_minimum = 1;
    }
}


/*
 * Save the leaf responsible for unsuccessful termination.
 */
static void
interval_engine_set_unresolved_leaf(
    interval_engine_stats_t *stats,
    const interval_leaf_t *leaf
)
{
    interval_engine_leaf_set(
        &stats->unresolved_leaf,
        leaf
    );

    stats->has_unresolved_leaf = 1;
}


/*
 * Simple iterative DFS stack.
 *
 * In binary depth-first traversal there is at most one pending
 * right sibling per depth.  Therefore max_depth + 2 entries
 * are sufficient.
 */
typedef struct
{
    interval_leaf_t *items;

    size_t size;
    size_t capacity;

} interval_engine_stack_t;


static int
interval_engine_stack_init(
    interval_engine_stack_t *stack,
    ulong max_depth
)
{
    size_t i;
    size_t capacity;

    stack->items = NULL;
    stack->size = 0;
    stack->capacity = 0;

    if (max_depth >
        (ulong) (SIZE_MAX - 2))
    {
        return 0;
    }

    capacity =
        (size_t) max_depth + 2;

    stack->items =
        (interval_leaf_t *)
        malloc(
            capacity *
            sizeof(interval_leaf_t)
        );

    if (stack->items == NULL)
        return 0;

    for (i = 0;
         i < capacity;
         ++i)
    {
        interval_leaf_init(
            stack->items + i
        );
    }

    stack->capacity =
        capacity;

    return 1;
}


static void
interval_engine_stack_clear(
    interval_engine_stack_t *stack
)
{
    size_t i;

    if (stack->items != NULL)
    {
        for (i = 0;
             i < stack->capacity;
             ++i)
        {
            interval_leaf_clear(
                stack->items + i
            );
        }

        free(
            stack->items
        );
    }

    stack->items = NULL;
    stack->size = 0;
    stack->capacity = 0;
}


static int
interval_engine_stack_push(
    interval_engine_stack_t *stack,
    const interval_leaf_t *leaf
)
{
    if (stack->size >=
        stack->capacity)
    {
        return 0;
    }

    interval_engine_leaf_set(
        stack->items + stack->size,
        leaf
    );

    stack->size++;

    return 1;
}


static int
interval_engine_stack_pop(
    interval_leaf_t *leaf,
    interval_engine_stack_t *stack
)
{
    if (stack->size == 0)
        return 0;

    stack->size--;

    interval_engine_leaf_set(
        leaf,
        stack->items + stack->size
    );

    return 1;
}


void
interval_engine_stats_init(
    interval_engine_stats_t *stats
)
{
    arf_init(
        stats->minimum_certified_lower
    );

    interval_leaf_init(
        &stats->minimum_leaf
    );

    interval_leaf_init(
        &stats->unresolved_leaf
    );

    interval_engine_stats_reset(
        stats
    );
}


void
interval_engine_stats_clear(
    interval_engine_stats_t *stats
)
{
    interval_leaf_clear(
        &stats->unresolved_leaf
    );

    interval_leaf_clear(
        &stats->minimum_leaf
    );

    arf_clear(
        stats->minimum_certified_lower
    );
}


int
interval_engine_certify_subtree(
    interval_engine_stats_t *stats,
    const interval_domain_t *domain,
    const interval_leaf_t *root,
    const arb_t kappa,
    slong prec,
    ulong max_depth,
    ulong progress_every,
    interval_engine_leaf_sink_t leaf_sink,
    void *leaf_sink_data
)
{
    interval_engine_stack_t stack;

    interval_leaf_t current;
    interval_leaf_t left;
    interval_leaf_t right;

    bilinear_certificate_t certificate;

    arb_t interval;

    int stack_initialized = 0;
    int result = 0;

    clock_t start_clock;
    clock_t current_clock;

    double elapsed_seconds;
    double evaluations_per_second;


    interval_engine_stats_reset(
        stats
    );


    if (!interval_domain_is_valid(domain))
        return 0;

    if (!interval_leaf_is_valid(root))
        return 0;

    if (!arb_is_finite(kappa))
        return 0;

    if (prec <= 0)
        return 0;

    /*
     * max_depth is absolute.
     */
    if (root->depth > max_depth)
        return 0;


    interval_leaf_init(
        &current
    );

    interval_leaf_init(
        &left
    );

    interval_leaf_init(
        &right
    );

    bilinear_certificate_init(
        &certificate
    );

    arb_init(
        interval
    );


    if (!interval_engine_stack_init(
            &stack,
            max_depth))
    {
        goto cleanup;
    }

    stack_initialized = 1;


    /*
     * Start DFS from the exact symbolic root.
     */
    if (!interval_engine_stack_push(
            &stack,
            root))
    {
        goto cleanup;
    }


    /*
     * Diagnostic timer only.
     *
     * It has no role in any rigorous certification decision.
     */
    start_clock =
        clock();


    while (stack.size > 0)
    {
        int evaluation_completed;


        if (!interval_engine_stack_pop(
                &current,
                &stack))
        {
            stats->status =
                INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE;

            goto cleanup;
        }


        if (current.depth >
            stats->max_depth_reached)
        {
            stats->max_depth_reached =
                current.depth;
        }


        stats->evaluations++;


        /*
         * Optional progress diagnostics.
         */
        if (progress_every != 0 &&
            stats->evaluations %
                progress_every == 0)
        {
            current_clock =
                clock();

            elapsed_seconds =
                (double)
                (current_clock - start_clock)
                /
                (double) CLOCKS_PER_SEC;

            if (elapsed_seconds > 0.0)
            {
                evaluations_per_second =
                    (double) stats->evaluations /
                    elapsed_seconds;
            }
            else
            {
                evaluations_per_second =
                    0.0;
            }

            printf(
                "[progress]"
                " eval=%lu"
                " certified=%lu"
                " covered=%.2Lf%%"
                " bisected=%lu"
                " depth=%lu"
                " max_depth=%lu"
                " stack=%zu"
                " elapsed=%.1f s"
                " rate=%.1f eval/s"
                " leaf=(d=%lu,k=",
                stats->evaluations,
                stats->certified_leaves,
                100.0L * stats->certified_fraction,
                stats->bisected_leaves,
                current.depth,
                stats->max_depth_reached,
                stack.size,
                elapsed_seconds,
                evaluations_per_second,
                current.depth
            );

            fmpz_print(
                current.index
            );

            printf(")\n");

            fflush(stdout);
        }


        /*
         * Convert the exact symbolic leaf into one rigorous Arb
         * enclosure solely for function evaluation.
         */
        if (!interval_leaf_arb(
                interval,
                domain,
                &current,
                prec))
        {
            stats->status =
                INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE;

            interval_engine_set_unresolved_leaf(
                stats,
                &current
            );

            goto cleanup;
        }


        /*
         * Attempt the rigorous local certification.
         *
         * A return value of zero may simply mean that the direct
         * Hardy evaluation was not usable on an interval this
         * wide.  In that case subdivision is appropriate.
         */
        evaluation_completed =
            bilinear_certify(
                &certificate,
                interval,
                kappa,
                prec
            );


        if (evaluation_completed &&
            certificate.status ==
                BILINEAR_STATUS_CERTIFIED)
        {
            /*
            * The sink is called before the leaf is counted as
            * successfully accepted by the engine.
            *
            * This guarantees that a writer failure cannot leave the
            * run marked as successfully certified.
            */
            if (leaf_sink != NULL)
            {
                if (!leaf_sink(
                        &current,
                        leaf_sink_data))
                {
                    stats->status =
                        INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE;

                    interval_engine_set_unresolved_leaf(
                        stats,
                        &current
                    );

                    goto cleanup;
                }
            }

            stats->certified_leaves++;

            stats->certified_fraction +=
                ldexpl(
                    1.0L,
                    -(int) (current.depth - root->depth)
                );


            interval_engine_update_minimum(
                stats,
                &current,
                certificate.certified_lower
            );


            continue;
        }


        /*
         * If bilinear_certify completed rigorously, the only
         * nonterminal result allowed here is INCONCLUSIVE.
         */
        if (evaluation_completed &&
            certificate.status !=
                BILINEAR_STATUS_INCONCLUSIVE)
        {
            stats->status =
                INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE;

            interval_engine_set_unresolved_leaf(
                stats,
                &current
            );

            goto cleanup;
        }


        /*
         * Either:
         *
         *   - the rigorous enclosure was inconclusive, or
         *   - the direct Hardy evaluation was unusable at this
         *     interval width.
         *
         * Neither is a mathematical counterexample.
         */
        if (current.depth >=
            max_depth)
        {
            stats->status =
                INTERVAL_ENGINE_STATUS_DEPTH_LIMIT;

            interval_engine_set_unresolved_leaf(
                stats,
                &current
            );

            goto cleanup;
        }


        /*
         * Exact symbolic bisection.
         */
        if (!interval_leaf_bisect(
                &left,
                &right,
                &current))
        {
            stats->status =
                INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE;

            interval_engine_set_unresolved_leaf(
                stats,
                &current
            );

            goto cleanup;
        }


        stats->bisected_leaves++;


        /*
         * Push right first and left second.
         *
         * Since the stack is LIFO, leaves are processed from
         * left to right.
         */
        if (!interval_engine_stack_push(
                &stack,
                &right))
        {
            stats->status =
                INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE;

            goto cleanup;
        }

        if (!interval_engine_stack_push(
                &stack,
                &left))
        {
            stats->status =
                INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE;

            goto cleanup;
        }
    }


    /*
     * The stack becomes empty only after every terminal branch
     * has been rigorously certified.
     */
    if (stats->certified_leaves == 0 ||
        !stats->has_minimum)
    {
        stats->status =
            INTERVAL_ENGINE_STATUS_EVALUATION_FAILURE;

        goto cleanup;
    }


    stats->status =
        INTERVAL_ENGINE_STATUS_SUCCESS;

    result = 1;


cleanup:

    if (stack_initialized)
    {
        interval_engine_stack_clear(
            &stack
        );
    }

    arb_clear(
        interval
    );

    bilinear_certificate_clear(
        &certificate
    );

    interval_leaf_clear(
        &right
    );

    interval_leaf_clear(
        &left
    );

    interval_leaf_clear(
        &current
    );

    return result;
}

int
interval_engine_certify_domain(
    interval_engine_stats_t *stats,
    const interval_domain_t *domain,
    const arb_t kappa,
    slong prec,
    ulong max_depth,
    ulong progress_every,
    interval_engine_leaf_sink_t leaf_sink,
    void *leaf_sink_data
)
{
    interval_leaf_t root;

    int result;


    interval_leaf_init(
        &root
    );

    interval_leaf_set_root(
        &root
    );


    result =
        interval_engine_certify_subtree(
            stats,
            domain,
            &root,
            kappa,
            prec,
            max_depth,
            progress_every,
            leaf_sink,
            leaf_sink_data
        );


    interval_leaf_clear(
        &root
    );

    return result;
}