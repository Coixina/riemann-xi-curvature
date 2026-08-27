#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arb.h"
#include "arf.h"
#include "fmpz.h"

#include "certificate_writer.h"
#include "interval_engine.h"


#define RH_DEFAULT_PRECISION       256
#define RH_MINIMUM_PRECISION       256
#define RH_DEFAULT_MAX_DEPTH       40
#define RH_DEFAULT_PROGRESS_EVERY  10000


/*
 * Parse an unsigned long exactly from decimal text.
 */
static int
parse_ulong_token(
    ulong *value,
    const char *text
)
{
    char *end;

    unsigned long parsed;


    if (value == NULL ||
        text == NULL ||
        *text == '\0')
    {
        return 0;
    }


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


/*
 * Parse a signed long exactly from decimal text.
 */
static int
parse_slong_token(
    slong *value,
    const char *text
)
{
    char *end;

    long parsed;


    if (value == NULL ||
        text == NULL ||
        *text == '\0')
    {
        return 0;
    }


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


/*
 * Canonicalize
 *
 *     numerator / denominator
 *
 * so that
 *
 *     gcd(|numerator|, denominator) = 1,
 *     denominator > 0.
 */
static int
canonicalize_rational(
    fmpz_t numerator,
    fmpz_t denominator
)
{
    fmpz_t gcd;
    fmpz_t abs_numerator;


    if (fmpz_is_zero(
            denominator))
    {
        return 0;
    }


    /*
     * Move the sign to the numerator.
     */
    if (fmpz_sgn(
            denominator) < 0)
    {
        fmpz_neg(
            numerator,
            numerator
        );

        fmpz_neg(
            denominator,
            denominator
        );
    }


    /*
     * Canonical zero.
     */
    if (fmpz_is_zero(
            numerator))
    {
        fmpz_one(
            denominator
        );

        return 1;
    }


    fmpz_init(gcd);
    fmpz_init(abs_numerator);


    fmpz_abs(
        abs_numerator,
        numerator
    );


    fmpz_gcd(
        gcd,
        abs_numerator,
        denominator
    );


    fmpz_divexact(
        numerator,
        numerator,
        gcd
    );

    fmpz_divexact(
        denominator,
        denominator,
        gcd
    );


    fmpz_clear(abs_numerator);
    fmpz_clear(gcd);


    return 1;
}


/*
 * Parse an exact rational written as
 *
 *     numerator/denominator.
 */
static int
parse_fraction_kappa(
    fmpz_t numerator,
    fmpz_t denominator,
    const char *text
)
{
    const char *slash;

    char *num_text;
    char *den_text;

    size_t num_length;
    size_t den_length;

    int success = 0;


    slash =
        strchr(
            text,
            '/'
        );


    if (slash == NULL)
        return 0;


    /*
     * Exactly one slash.
     */
    if (strchr(
            slash + 1,
            '/') != NULL)
    {
        return 0;
    }


    num_length =
        (size_t) (slash - text);

    den_length =
        strlen(
            slash + 1
        );


    if (num_length == 0 ||
        den_length == 0)
    {
        return 0;
    }


    num_text =
        (char *)
        malloc(
            num_length + 1
        );

    den_text =
        (char *)
        malloc(
            den_length + 1
        );


    if (num_text == NULL ||
        den_text == NULL)
    {
        free(num_text);
        free(den_text);

        return 0;
    }


    memcpy(
        num_text,
        text,
        num_length
    );

    num_text[num_length] =
        '\0';


    memcpy(
        den_text,
        slash + 1,
        den_length + 1
    );


    if (fmpz_set_str(
            numerator,
            num_text,
            10) != 0)
    {
        goto cleanup;
    }


    if (fmpz_set_str(
            denominator,
            den_text,
            10) != 0)
    {
        goto cleanup;
    }


    if (!canonicalize_rational(
            numerator,
            denominator))
    {
        goto cleanup;
    }


    success = 1;


cleanup:

    free(den_text);
    free(num_text);


    return success;
}


/*
 * Parse an exact decimal.
 *
 * Examples:
 *
 *     0.092
 *       -> 92 / 1000
 *       -> 23 / 250
 *
 *     -0.09250
 *       -> -9250 / 100000
 *       -> -37 / 400
 *
 * No binary floating-point conversion is used.
 *
 * Scientific notation is deliberately not accepted.
 */
static int
parse_decimal_kappa(
    fmpz_t numerator,
    fmpz_t denominator,
    const char *text
)
{
    const char *p;

    char *digits;

    size_t length;
    size_t digits_used = 0;

    ulong fractional_digits = 0;

    int negative = 0;
    int seen_dot = 0;
    int seen_digit = 0;

    int success = 0;


    if (text == NULL ||
        *text == '\0')
    {
        return 0;
    }


    length =
        strlen(text);


    digits =
        (char *)
        malloc(
            length + 2
        );


    if (digits == NULL)
        return 0;


    p = text;


    if (*p == '+' ||
        *p == '-')
    {
        negative =
            (*p == '-');

        p++;

        if (*p == '\0')
            goto cleanup;
    }


    while (*p != '\0')
    {
        if (*p == '.')
        {
            if (seen_dot)
                goto cleanup;

            seen_dot = 1;

            p++;
            continue;
        }


        if (!isdigit(
                (unsigned char) *p))
        {
            goto cleanup;
        }


        seen_digit = 1;

        digits[digits_used++] =
            *p;


        if (seen_dot)
            fractional_digits++;


        p++;
    }


    if (!seen_digit)
        goto cleanup;


    /*
     * fmpz_set_str needs at least one digit.
     */
    digits[digits_used] =
        '\0';


    if (fmpz_set_str(
            numerator,
            digits,
            10) != 0)
    {
        goto cleanup;
    }


    if (negative)
    {
        fmpz_neg(
            numerator,
            numerator
        );
    }


    /*
     * denominator = 10^fractional_digits.
     *
     * Build it exactly with fmpz arithmetic.
     */
    fmpz_one(
        denominator
    );


    while (fractional_digits > 0)
    {
        fmpz_mul_ui(
            denominator,
            denominator,
            10
        );

        fractional_digits--;
    }


    if (!canonicalize_rational(
            numerator,
            denominator))
    {
        goto cleanup;
    }


    success = 1;


cleanup:

    free(digits);


    return success;
}


/*
 * Accepted kappa syntax:
 *
 *     integer
 *     decimal
 *     numerator/denominator
 *
 * Examples:
 *
 *     0
 *     0.092
 *     -0.125
 *     23/250
 *     -1/8
 *
 * Scientific notation is intentionally rejected.
 */
static int
parse_kappa(
    fmpz_t numerator,
    fmpz_t denominator,
    const char *text
)
{
    if (strchr(
            text,
            '/') != NULL)
    {
        return
            parse_fraction_kappa(
                numerator,
                denominator,
                text
            );
    }


    return
        parse_decimal_kappa(
            numerator,
            denominator,
            text
        );
}


/*
 * Convert exact rational kappa to one rigorous Arb value.
 */
static int
build_kappa_arb(
    arb_t kappa,
    const fmpz_t numerator,
    const fmpz_t denominator,
    slong prec
)
{
    arb_t den;


    arb_init(
        den
    );


    arb_set_fmpz(
        kappa,
        numerator
    );

    arb_set_fmpz(
        den,
        denominator
    );


    arb_div(
        kappa,
        kappa,
        den,
        prec
    );


    arb_clear(
        den
    );


    return
        arb_is_finite(
            kappa
        );
}


/*
 * Construct a rigorous frozen domain containing the exact
 * mathematical interval
 *
 *     [2*pi, 2*pi*N^2].
 *
 * interval_domain_set_arb() freezes outward ARF bounds.
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
     * Multiply by N twice rather than forming N^2 in a native
     * integer, avoiding native integer overflow.
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
 * Adapter between interval_engine's generic certified-leaf
 * callback and certificate_writer.
 */
static int
certificate_leaf_sink(
    const interval_leaf_t *leaf,
    void *user_data
)
{
    certificate_writer_t *writer;


    if (user_data == NULL)
        return 0;


    writer =
        (certificate_writer_t *)
        user_data;


    return
        certificate_writer_write_leaf(
            writer,
            leaf
        );
}


static void
print_usage(
    const char *program
)
{
    fprintf(
        stderr,
        "Usage:\n"
        "  %s large_t_N kappa output_certificate"
        " [precision] [max_depth] [progress_every]\n"
        "\n"
        "kappa syntax:\n"
        "  exact fraction : 23/250\n"
        "  exact decimal  : 0.092\n"
        "  integer        : 0\n"
        "\n"
        "Defaults:\n"
        "  precision      = %d bits\n"
        "  max_depth      = %d\n"
        "  progress_every = %d\n",
        program,
        RH_DEFAULT_PRECISION,
        RH_DEFAULT_MAX_DEPTH,
        RH_DEFAULT_PROGRESS_EVERY
    );
}


int
main(
    int argc,
    char **argv
)
{
    ulong large_t_N;

    slong precision =
        RH_DEFAULT_PRECISION;

    ulong max_depth =
        RH_DEFAULT_MAX_DEPTH;

    ulong progress_every =
        RH_DEFAULT_PROGRESS_EVERY;


    const char *kappa_text;
    const char *output_path;


    fmpz_t kappa_num;
    fmpz_t kappa_den;

    arb_t kappa;

    interval_domain_t domain;

    interval_engine_stats_t stats;

    certificate_writer_t writer;


    int writer_open = 0;
    int engine_success;

    int exit_code = 1;


    /*
     * Required:
     *
     *   large_t_N
     *   kappa
     *   output_certificate
     *
     * Optional:
     *
     *   precision
     *   max_depth
     *   progress_every
     */
    if (argc < 4 ||
        argc > 7)
    {
        print_usage(
            argv[0]
        );

        return 2;
    }


    /*
     * large_t_N.
     */
    if (!parse_ulong_token(
            &large_t_N,
            argv[1]) ||
        large_t_N < 2)
    {
        fprintf(
            stderr,
            "Invalid large_t_N: %s\n",
            argv[1]
        );

        return 2;
    }


    kappa_text =
        argv[2];

    output_path =
        argv[3];


    /*
     * Optional precision.
     */
    if (argc >= 5)
    {
        if (!parse_slong_token(
                &precision,
                argv[4]))
        {
            fprintf(
                stderr,
                "Invalid precision: %s\n",
                argv[4]
            );

            return 2;
        }
    }


    if (precision <
        RH_MINIMUM_PRECISION)
    {
        fprintf(
            stderr,
            "Production precision must be at least %d bits.\n",
            RH_MINIMUM_PRECISION
        );

        return 2;
    }


    /*
     * Optional maximum depth.
     */
    if (argc >= 6)
    {
        if (!parse_ulong_token(
                &max_depth,
                argv[5]))
        {
            fprintf(
                stderr,
                "Invalid max_depth: %s\n",
                argv[5]
            );

            return 2;
        }
    }


    /*
     * Optional progress cadence.
     *
     * Zero disables progress output.
     */
    if (argc >= 7)
    {
        if (!parse_ulong_token(
                &progress_every,
                argv[6]))
        {
            fprintf(
                stderr,
                "Invalid progress_every: %s\n",
                argv[6]
            );

            return 2;
        }
    }


    fmpz_init(
        kappa_num
    );

    fmpz_init(
        kappa_den
    );

    arb_init(
        kappa
    );

    interval_domain_init(
        &domain
    );

    interval_engine_stats_init(
        &stats
    );

    certificate_writer_init(
        &writer
    );


    /*
     * Parse kappa exactly from user text.
     */
    if (!parse_kappa(
            kappa_num,
            kappa_den,
            kappa_text))
    {
        fprintf(
            stderr,
            "Invalid kappa: %s\n"
            "Use an integer, exact decimal, or fraction p/q.\n",
            kappa_text
        );

        exit_code = 2;
        goto cleanup;
    }


    /*
     * Build the Arb representation used in the rigorous
     * calculation.
     */
    if (!build_kappa_arb(
            kappa,
            kappa_num,
            kappa_den,
            precision))
    {
        fprintf(
            stderr,
            "Could not construct rigorous kappa.\n"
        );

        goto cleanup;
    }


    /*
     * Build the complete finite target interval
     *
     *     [2*pi, 2*pi*large_t_N^2].
     */
    if (!build_finite_domain(
            &domain,
            large_t_N,
            precision))
    {
        fprintf(
            stderr,
            "Could not construct finite certification domain.\n"
        );

        goto cleanup;
    }


    /*
     * Print the complete run configuration before starting a
     * potentially long computation.
     */
    printf(
        "============================================================\n"
        "RH HARDY FINITE CERTIFIER\n"
        "============================================================\n"
        "large_t_N       = %lu\n"
        "domain          = [2*pi, 2*pi*%lu^2]\n"
        "input kappa     = %s\n"
        "exact kappa     = ",
        large_t_N,
        large_t_N,
        kappa_text
    );

    fmpz_print(
        kappa_num
    );

    printf(" / ");

    fmpz_print(
        kappa_den
    );

    printf(
        "\n"
        "precision       = %ld bits\n"
        "max depth       = %lu\n"
        "progress every  = %lu\n"
        "certificate     = %s\n"
        "============================================================\n",
        precision,
        max_depth,
        progress_every,
        output_path
    );

    fflush(stdout);


    /*
     * Open the certificate before running the engine.
     *
     * Leaves will be streamed directly to disk as they are
     * certified from left to right.
     */
    if (!certificate_writer_open(
            &writer,
            output_path,
            large_t_N,
            kappa_num,
            kappa_den,
            precision))
    {
        fprintf(
            stderr,
            "Could not open certificate file: %s\n",
            output_path
        );

        goto cleanup;
    }


    writer_open = 1;


    /*
     * Run the complete adaptive finite certification.
     */
    engine_success =
        interval_engine_certify_domain(
            &stats,
            &domain,
            kappa,
            precision,
            max_depth,
            progress_every,
            certificate_leaf_sink,
            &writer
        );


    if (!engine_success ||
        stats.status !=
            INTERVAL_ENGINE_STATUS_SUCCESS)
    {
        fprintf(
            stderr,
            "\nCertification did not complete successfully.\n"
            "status            = %d\n"
            "evaluations       = %lu\n"
            "certified leaves  = %lu\n"
            "bisected leaves   = %lu\n"
            "max depth reached = %lu\n",
            (int) stats.status,
            stats.evaluations,
            stats.certified_leaves,
            stats.bisected_leaves,
            stats.max_depth_reached
        );


        if (stats.has_unresolved_leaf)
        {
            fprintf(
                stderr,
                "unresolved leaf   = (d=%lu,k=",
                stats.unresolved_leaf.depth
            );

            fmpz_fprint(
                stderr,
                stats.unresolved_leaf.index
            );

            fprintf(
                stderr,
                ")\n"
            );
        }


        /*
         * Never leave an unfinished file looking like a valid
         * production certificate.
         */
        certificate_writer_abort(
            &writer
        );

        writer_open = 0;


        /*
         * Remove the incomplete output file if possible.
         */
        remove(
            output_path
        );


        goto cleanup;
    }


    /*
     * Internal consistency:
     *
     * every accepted engine leaf must have been written.
     */
    if (writer.leaf_count !=
        stats.certified_leaves)
    {
        fprintf(
            stderr,
            "Internal error: writer leaf count does not match "
            "engine leaf count.\n"
        );


        certificate_writer_abort(
            &writer
        );

        writer_open = 0;

        remove(
            output_path
        );


        goto cleanup;
    }


    /*
     * Only a completely successful certification receives the
     * END_LEAVES marker and final leaf_count.
     */
    if (!certificate_writer_finalize(
            &writer))
    {
        fprintf(
            stderr,
            "Could not finalize certificate file.\n"
        );

        writer_open = 0;

        remove(
            output_path
        );

        goto cleanup;
    }


    writer_open = 0;


    /*
     * Final report.
     */
    printf(
        "\n"
        "============================================================\n"
        "CERTIFICATION SUCCESSFUL\n"
        "============================================================\n"
        "evaluations       = %lu\n"
        "certified leaves  = %lu\n"
        "interval covered  = %.6Lf %%\n"
        "bisected leaves   = %lu\n"
        "max depth used    = %lu\n",
        stats.evaluations,
        stats.certified_leaves,
        100.0L * stats.certified_fraction,
        stats.bisected_leaves,
        stats.max_depth_reached
    );


    if (stats.has_minimum)
    {
        printf(
            "minimum lower    = "
        );

        arf_printd(
            stats.minimum_certified_lower,
            30
        );

        printf(
            "\n"
            "minimum leaf     = (d=%lu,k=",
            stats.minimum_leaf.depth
        );

        fmpz_print(
            stats.minimum_leaf.index
        );

        printf(
            ")\n"
        );
    }


    printf(
        "certificate      = %s\n"
        "============================================================\n",
        output_path
    );


    exit_code = 0;


cleanup:

    if (writer_open)
    {
        certificate_writer_abort(
            &writer
        );
    }


    interval_engine_stats_clear(
        &stats
    );

    interval_domain_clear(
        &domain
    );

    arb_clear(
        kappa
    );

    fmpz_clear(
        kappa_den
    );

    fmpz_clear(
        kappa_num
    );


    return exit_code;
}