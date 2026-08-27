#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arb.h"
#include "arf.h"
#include "fmpz.h"

#include "bilinear_certifier.h"
#include "interval_geometry.h"


#define VERIFY_DEFAULT_PRECISION 256
#define VERIFY_MINIMUM_PRECISION 256
#define VERIFY_DEFAULT_MAX_EXTRA_DEPTH 12

#define TOKEN_SIZE 4096
#define MAGIC_SIZE 64

#define RH_CERTIFICATE_FORMAT_VERSION 1


typedef struct
{
    ulong large_t_N;

    fmpz_t kappa_num;
    fmpz_t kappa_den;

    slong generator_precision;

    ulong leaf_count;

} certificate_header_t;


static void
certificate_header_init(
    certificate_header_t *header
)
{
    header->large_t_N = 0;

    fmpz_init(
        header->kappa_num
    );

    fmpz_init(
        header->kappa_den
    );

    header->generator_precision = 0;
    header->leaf_count = 0;
}


static void
certificate_header_clear(
    certificate_header_t *header
)
{
    fmpz_clear(
        header->kappa_den
    );

    fmpz_clear(
        header->kappa_num
    );
}


static int
parse_ulong_token(
    ulong *value,
    const char *text
)
{
    char *end;

    unsigned long parsed;


    errno = 0;

    parsed =
        strtoul(
            text,
            &end,
            10
        );

    if (errno != 0 ||
        end == text ||
        *end != '\0')
    {
        return 0;
    }

    *value =
        (ulong) parsed;

    return 1;
}


static int
parse_slong_token(
    slong *value,
    const char *text
)
{
    char *end;

    long parsed;


    errno = 0;

    parsed =
        strtol(
            text,
            &end,
            10
        );

    if (errno != 0 ||
        end == text ||
        *end != '\0')
    {
        return 0;
    }

    *value =
        (slong) parsed;

    return 1;
}


static int
read_named_token(
    FILE *file,
    const char *expected_name,
    char *value
)
{
    char name[MAGIC_SIZE];


    if (fscanf(
            file,
            "%63s %4095s",
            name,
            value
        ) != 2)
    {
        return 0;
    }

    return
        strcmp(
            name,
            expected_name
        ) == 0;
}


static int
read_certificate_header(
    FILE *file,
    certificate_header_t *header
)
{
    char magic[MAGIC_SIZE];
    char token[TOKEN_SIZE];

    ulong version;


    if (fscanf(
            file,
            "%63s %lu",
            magic,
            &version
        ) != 2)
    {
        return 0;
    }


    if (strcmp(
            magic,
            "RH_FINITE_CERTIFICATE") != 0)
    {
        return 0;
    }


    if (version !=
        RH_CERTIFICATE_FORMAT_VERSION)
    {
        return 0;
    }


    if (!read_named_token(
            file,
            "large_t_N",
            token))
    {
        return 0;
    }

    if (!parse_ulong_token(
            &header->large_t_N,
            token))
    {
        return 0;
    }

    if (header->large_t_N < 2)
        return 0;


    if (!read_named_token(
            file,
            "kappa_num",
            token))
    {
        return 0;
    }

    if (fmpz_set_str(
            header->kappa_num,
            token,
            10) != 0)
    {
        return 0;
    }


    if (!read_named_token(
            file,
            "kappa_den",
            token))
    {
        return 0;
    }

    if (fmpz_set_str(
            header->kappa_den,
            token,
            10) != 0)
    {
        return 0;
    }

    if (fmpz_sgn(
            header->kappa_den) <= 0)
    {
        return 0;
    }


    if (!read_named_token(
            file,
            "generator_precision",
            token))
    {
        return 0;
    }

    if (!parse_slong_token(
            &header->generator_precision,
            token))
    {
        return 0;
    }

    if (header->generator_precision <= 0)
        return 0;


    if (!read_named_token(
            file,
            "leaf_count",
            token))
    {
        return 0;
    }

    if (!parse_ulong_token(
            &header->leaf_count,
            token))
    {
        return 0;
    }

    if (header->leaf_count == 0)
        return 0;


    if (fscanf(
            file,
            "%4095s",
            token
        ) != 1)
    {
        return 0;
    }

    if (strcmp(
            token,
            "BEGIN_LEAVES") != 0)
    {
        return 0;
    }


    return 1;
}


static int
read_leaf(
    FILE *file,
    interval_leaf_t *leaf
)
{
    char depth_text[TOKEN_SIZE];
    char index_text[TOKEN_SIZE];

    ulong depth;

    fmpz_t index;

    int success = 0;


    if (fscanf(
            file,
            "%4095s %4095s",
            depth_text,
            index_text
        ) != 2)
    {
        return 0;
    }


    if (!parse_ulong_token(
            &depth,
            depth_text))
    {
        return 0;
    }


    fmpz_init(
        index
    );


    if (fmpz_set_str(
            index,
            index_text,
            10) != 0)
    {
        goto cleanup;
    }


    if (!interval_leaf_set(
            leaf,
            depth,
            index))
    {
        goto cleanup;
    }


    success = 1;


cleanup:

    fmpz_clear(
        index
    );

    return success;
}


/*
 * Exact integer adjacency test.
 *
 * Previous:
 *
 *     (d,k)
 *
 * Current:
 *
 *     (e,l)
 *
 * Adjacency means
 *
 *     (k+1)/2^d = l/2^e
 *
 * and is checked exactly as
 *
 *     (k+1) 2^e = l 2^d.
 */
static int
leaves_are_adjacent(
    const interval_leaf_t *previous,
    const interval_leaf_t *current
)
{
    fmpz_t lhs;
    fmpz_t rhs;
    fmpz_t tmp;

    int result;


    fmpz_init(lhs);
    fmpz_init(rhs);
    fmpz_init(tmp);


    fmpz_add_ui(
        tmp,
        previous->index,
        1
    );

    fmpz_mul_2exp(
        lhs,
        tmp,
        current->depth
    );


    fmpz_mul_2exp(
        rhs,
        current->index,
        previous->depth
    );


    result =
        fmpz_equal(
            lhs,
            rhs
        );


    fmpz_clear(tmp);
    fmpz_clear(rhs);
    fmpz_clear(lhs);


    return result;
}


static int
leaf_ends_at_one(
    const interval_leaf_t *leaf
)
{
    fmpz_t lhs;
    fmpz_t rhs;

    int result;


    fmpz_init(lhs);
    fmpz_init(rhs);


    fmpz_add_ui(
        lhs,
        leaf->index,
        1
    );


    fmpz_one(
        rhs
    );

    fmpz_mul_2exp(
        rhs,
        rhs,
        leaf->depth
    );


    result =
        fmpz_equal(
            lhs,
            rhs
        );


    fmpz_clear(rhs);
    fmpz_clear(lhs);


    return result;
}


/*
 * First pass:
 *
 * validate format and prove exact dyadic partition using only
 * exact integer arithmetic.
 */
static int
verify_partition(
    const char *path,
    certificate_header_t *header
)
{
    FILE *file;

    interval_leaf_t previous;
    interval_leaf_t current;

    char token[TOKEN_SIZE];

    ulong j;

    int have_previous = 0;
    int success = 0;


    file =
        fopen(
            path,
            "r"
        );

    if (file == NULL)
        return 0;


    interval_leaf_init(
        &previous
    );

    interval_leaf_init(
        &current
    );


    if (!read_certificate_header(
            file,
            header))
    {
        goto cleanup;
    }


    for (j = 0;
         j < header->leaf_count;
         ++j)
    {
        if (!read_leaf(
                file,
                &current))
        {
            goto cleanup;
        }


        if (!have_previous)
        {
            /*
             * First leaf must begin at zero.
             */
            if (!fmpz_is_zero(
                    current.index))
            {
                goto cleanup;
            }

            have_previous = 1;
        }
        else
        {
            if (!leaves_are_adjacent(
                    &previous,
                    &current))
            {
                goto cleanup;
            }
        }


        previous.depth =
            current.depth;

        fmpz_set(
            previous.index,
            current.index
        );
    }


    if (!have_previous)
        goto cleanup;


    /*
     * Last leaf must end exactly at one.
     */
    if (!leaf_ends_at_one(
            &previous))
    {
        goto cleanup;
    }


    if (fscanf(
            file,
            "%4095s",
            token
        ) != 1)
    {
        goto cleanup;
    }

    if (strcmp(
            token,
            "END_LEAVES") != 0)
    {
        goto cleanup;
    }


    /*
     * No non-whitespace tokens are allowed after END_LEAVES.
     */
    if (fscanf(
            file,
            "%4095s",
            token
        ) == 1)
    {
        goto cleanup;
    }


    success = 1;


cleanup:

    interval_leaf_clear(
        &current
    );

    interval_leaf_clear(
        &previous
    );

    fclose(
        file
    );


    return success;
}


/*
 * Construct a rigorous outer frozen domain containing exactly
 * the mathematical target
 *
 *     [2*pi, 2*pi*N^2].
 */
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
     * left = 2*pi.
     */
    arb_mul_2exp_si(
        left,
        pi,
        1
    );


    /*
     * right = 2*pi*N^2.
     *
     * Multiply by N twice to avoid native integer overflow in
     * forming N^2.
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
 * Convert exact rational kappa_num / kappa_den to one rigorous
 * Arb value.
 */
static int
build_kappa(
    arb_t kappa,
    const certificate_header_t *header,
    slong prec
)
{
    arb_t denominator;


    arb_init(
        denominator
    );


    arb_set_fmpz(
        kappa,
        header->kappa_num
    );

    arb_set_fmpz(
        denominator,
        header->kappa_den
    );


    arb_div(
        kappa,
        kappa,
        denominator,
        prec
    );


    arb_clear(
        denominator
    );


    return
        arb_is_finite(kappa);
}

/*
 * Diagnostics for the analytic verification pass.
 *
 * stored_leaves_verified:
 *     number of leaves read from the certificate whose complete
 *     regions have been rigorously verified.
 *
 * analytic_evaluations:
 *     total calls to bilinear_certify(), including descendants
 *     created by local verification refinement.
 *
 * analytic_certified_leaves:
 *     number of terminal analytic leaves actually accepted.
 *
 * additional_subdivisions:
 *     number of local dyadic subdivisions performed by the
 *     verifier beyond the stored certificate partition.
 *
 * max_extra_depth_used:
 *     deepest refinement level below any stored leaf.
 */
typedef struct
{
    ulong stored_leaves_verified;

    ulong analytic_evaluations;
    ulong analytic_certified_leaves;
    ulong additional_subdivisions;

    ulong max_extra_depth_used;

} verification_stats_t;

/*
 * Verify one stored certificate leaf, refining it locally when
 * the requested verification precision does not yield a usable
 * or decisive enclosure on the stored leaf itself.
 *
 * The refinement changes only the analytic verification.
 * The certificate partition has already been checked exactly
 * and is not modified.
 */
static int
verify_leaf_with_refinement(
    const interval_domain_t *domain,
    const interval_leaf_t *leaf,
    const arb_t kappa,
    slong prec,
    ulong extra_depth,
    ulong max_extra_depth,
    verification_stats_t *stats,
    arf_t minimum_lower,
    interval_leaf_t *minimum_leaf,
    int *has_minimum
)
{
    arb_t interval;

    bilinear_certificate_t certificate;

    interval_leaf_t left;
    interval_leaf_t right;

    int completed;
    int success = 0;


    arb_init(
        interval
    );

    bilinear_certificate_init(
        &certificate
    );

    interval_leaf_init(
        &left
    );

    interval_leaf_init(
        &right
    );


    /*
     * Build the rigorous enclosure of the exact symbolic leaf.
     *
     * Failure here is unexpected geometry/evaluation failure,
     * not merely an inconclusive Hardy evaluation.
     */
    if (!interval_leaf_arb(
            interval,
            domain,
            leaf,
            prec))
    {
        goto cleanup;
    }


    stats->analytic_evaluations++;


    completed =
        bilinear_certify(
            &certificate,
            interval,
            kappa,
            prec
        );


    /*
     * Terminal successful analytic leaf.
     */
    if (completed &&
        certificate.status ==
            BILINEAR_STATUS_CERTIFIED)
    {
        stats->analytic_certified_leaves++;


        if (extra_depth >
            stats->max_extra_depth_used)
        {
            stats->max_extra_depth_used =
                extra_depth;
        }


        if (!*has_minimum ||
            arf_cmp(
                certificate.certified_lower,
                minimum_lower
            ) < 0)
        {
            arf_set(
                minimum_lower,
                certificate.certified_lower
            );

            minimum_leaf->depth =
                leaf->depth;

            fmpz_set(
                minimum_leaf->index,
                leaf->index
            );

            *has_minimum = 1;
        }


        success = 1;

        goto cleanup;
    }


    /*
     * A completed evaluation may legitimately be inconclusive.
     * Any other completed status is an internal failure.
     */
    if (completed &&
        certificate.status !=
            BILINEAR_STATUS_INCONCLUSIVE)
    {
        goto cleanup;
    }


    /*
     * Either the interval evaluation was unusable, or the
     * rigorous result was inconclusive.
     *
     * Refine locally if allowed.
     */
    if (extra_depth >=
        max_extra_depth)
    {
        fprintf(
            stderr,
            "\n[verify] maximum extra depth reached"
            " at leaf (d=%lu,k=",
            leaf->depth
        );

        fmpz_fprint(
            stderr,
            leaf->index
        );

        fprintf(
            stderr,
            "), extra_depth=%lu\n",
            extra_depth
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


    stats->additional_subdivisions++;


    /*
     * Both children must verify.
     */
    if (!verify_leaf_with_refinement(
            domain,
            &left,
            kappa,
            prec,
            extra_depth + 1,
            max_extra_depth,
            stats,
            minimum_lower,
            minimum_leaf,
            has_minimum))
    {
        goto cleanup;
    }


    if (!verify_leaf_with_refinement(
            domain,
            &right,
            kappa,
            prec,
            extra_depth + 1,
            max_extra_depth,
            stats,
            minimum_lower,
            minimum_leaf,
            has_minimum))
    {
        goto cleanup;
    }


    success = 1;


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


    return success;
}

/*
 * Second pass:
 *
 * independently recompute the rigorous analytic certificate on
 * every stored leaf.
 *
 * A stored leaf may be refined locally if the requested
 * verification precision does not certify it directly.
 */
static int
verify_analytic_leaves(
    const char *path,
    const certificate_header_t *expected_header,
    slong prec,
    ulong max_extra_depth,
    verification_stats_t *stats,
    arf_t minimum_lower,
    interval_leaf_t *minimum_leaf,
    int *has_minimum
)
{
    FILE *file;

    certificate_header_t header;

    interval_domain_t domain;

    interval_leaf_t leaf;

    arb_t kappa;

    char token[TOKEN_SIZE];

    ulong j;

    int success = 0;


    stats->stored_leaves_verified = 0;
    stats->analytic_evaluations = 0;
    stats->analytic_certified_leaves = 0;
    stats->additional_subdivisions = 0;
    stats->max_extra_depth_used = 0;

    *has_minimum = 0;

    arf_zero(
        minimum_lower
    );


    certificate_header_init(
        &header
    );

    interval_domain_init(
        &domain
    );

    interval_leaf_init(
        &leaf
    );

    arb_init(
        kappa
    );


    file =
        fopen(
            path,
            "r"
        );

    if (file == NULL)
        goto cleanup;


    if (!read_certificate_header(
            file,
            &header))
    {
        goto close_file;
    }


    /*
     * The file must not change between verification passes.
     */
    if (header.large_t_N !=
            expected_header->large_t_N ||
        header.generator_precision !=
            expected_header->generator_precision ||
        header.leaf_count !=
            expected_header->leaf_count ||
        !fmpz_equal(
            header.kappa_num,
            expected_header->kappa_num) ||
        !fmpz_equal(
            header.kappa_den,
            expected_header->kappa_den))
    {
        goto close_file;
    }


    if (!build_finite_domain(
            &domain,
            header.large_t_N,
            prec))
    {
        goto close_file;
    }


    if (!build_kappa(
            kappa,
            &header,
            prec))
    {
        goto close_file;
    }


    for (j = 0;
         j < header.leaf_count;
         ++j)
    {
        if (!read_leaf(
                file,
                &leaf))
        {
            goto close_file;
        }


        if (!verify_leaf_with_refinement(
                &domain,
                &leaf,
                kappa,
                prec,
                0,
                max_extra_depth,
                stats,
                minimum_lower,
                minimum_leaf,
                has_minimum))
        {
            fprintf(
                stderr,
                "\n[verify] analytic verification failed"
                " for stored leaf %lu / %lu: (d=%lu,k=",
                j + 1,
                header.leaf_count,
                leaf.depth
            );

            fmpz_fprint(
                stderr,
                leaf.index
            );

            fprintf(
                stderr,
                ")\n"
            );

            goto close_file;
        }


        stats->stored_leaves_verified++;


        if (stats->stored_leaves_verified %
            10000 == 0)
        {
            printf(
                "[verify] stored leaves=%lu / %lu"
                " analytic eval=%lu"
                " extra split=%lu"
                " max_extra=%lu\n",
                stats->stored_leaves_verified,
                header.leaf_count,
                stats->analytic_evaluations,
                stats->additional_subdivisions,
                stats->max_extra_depth_used
            );

            fflush(stdout);
        }
    }


    if (fscanf(
            file,
            "%4095s",
            token
        ) != 1)
    {
        goto close_file;
    }


    if (strcmp(
            token,
            "END_LEAVES") != 0)
    {
        goto close_file;
    }


    if (stats->stored_leaves_verified !=
        header.leaf_count)
    {
        goto close_file;
    }


    success = 1;


close_file:

    fclose(
        file
    );


cleanup:

    arb_clear(
        kappa
    );

    interval_leaf_clear(
        &leaf
    );

    interval_domain_clear(
        &domain
    );

    certificate_header_clear(
        &header
    );


    return success;
}

int
main(
    int argc,
    char **argv
)
{
    certificate_header_t header;

    interval_leaf_t minimum_leaf;

    arf_t minimum_lower;

    slong verification_precision =
        VERIFY_DEFAULT_PRECISION;

    verification_stats_t verification_stats;

    ulong max_extra_depth =
    VERIFY_DEFAULT_MAX_EXTRA_DEPTH;

    int has_minimum = 0;

    int partition_ok;
    int analytic_ok;


    if (argc < 2 ||
        argc > 4)
    {
        fprintf(
            stderr,
            "Usage: %s certificate.txt [precision] [max_extra_depth]\n",
            argv[0]
        );

        return 2;
    }

    if (argc >= 3)
    {
        if (!parse_slong_token(
                &verification_precision,
                argv[2]))
        {
            fprintf(
                stderr,
                "Invalid precision: %s\n",
                argv[2]
            );

            return 2;
        }
    }


    if (argc == 4)
    {
        if (!parse_ulong_token(
                &max_extra_depth,
                argv[3]))
        {
            fprintf(
                stderr,
                "Invalid max_extra_depth: %s\n",
                argv[3]
            );

            return 2;
        }
    }

    if (verification_precision <
        VERIFY_MINIMUM_PRECISION)
    {
        fprintf(
            stderr,
            "Verification precision must be at least %d bits.\n",
            VERIFY_MINIMUM_PRECISION
        );

        return 2;
    }


    certificate_header_init(
        &header
    );

    interval_leaf_init(
        &minimum_leaf
    );

    arf_init(
        minimum_lower
    );


    printf(
        "============================================================\n"
        "RH FINITE CERTIFICATE VERIFIER\n"
        "============================================================\n"
    );


    partition_ok =
        verify_partition(
            argv[1],
            &header
        );


    if (!partition_ok)
    {
        printf(
            "PARTITION CHECK        : FAIL\n"
            "CERTIFICATE            : INVALID\n"
        );

        arf_clear(minimum_lower);
        interval_leaf_clear(&minimum_leaf);
        certificate_header_clear(&header);

        return 1;
    }

    printf(
        "Certificate format      : RH-FINITE-1\n"
        "large_t_N               : %lu\n"
        "Generator precision     : %ld bits\n"
        "Verification precision  : %ld bits\n"
        "Max extra depth         : %lu\n"
        "Terminal leaves         : %lu\n"
        "PARTITION CHECK         : PASS\n",
        header.large_t_N,
        header.generator_precision,
        verification_precision,
        max_extra_depth,
        header.leaf_count
    );


    analytic_ok =
        verify_analytic_leaves(
            argv[1],
            &header,
            verification_precision,
            max_extra_depth,
            &verification_stats,
            minimum_lower,
            &minimum_leaf,
            &has_minimum
        );


    if (!analytic_ok ||
        !has_minimum)
    {
        printf(
            "RIGOROUS LEAF CHECK    : FAIL\n"
            "CERTIFICATE            : INVALID\n"
        );

        arf_clear(minimum_lower);
        interval_leaf_clear(&minimum_leaf);
        certificate_header_clear(&header);

        return 1;
    }

    printf(
        "RIGOROUS LEAF CHECK     : PASS\n"
        "Stored leaves verified  : %lu\n"
        "Analytic evaluations    : %lu\n"
        "Analytic terminal leaves: %lu\n"
        "Additional subdivisions : %lu\n"
        "Max extra depth used    : %lu\n"
        "Minimum verified lower  : ",
        verification_stats.stored_leaves_verified,
        verification_stats.analytic_evaluations,
        verification_stats.analytic_certified_leaves,
        verification_stats.additional_subdivisions,
        verification_stats.max_extra_depth_used
    );

    arf_printd(
        minimum_lower,
        30
    );

    printf(
        "\n"
        "Minimum leaf            : (d=%lu,k=",
        minimum_leaf.depth
    );

    fmpz_print(
        minimum_leaf.index
    );

    printf(
        ")\n"
        "CERTIFICATE             : VALID\n"
        "============================================================\n"
    );


    arf_clear(
        minimum_lower
    );

    interval_leaf_clear(
        &minimum_leaf
    );

    certificate_header_clear(
        &header
    );


    return 0;
}