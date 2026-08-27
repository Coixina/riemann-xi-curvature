#include <stdio.h>
#include <string.h>

#include "arb.h"
#include "fmpz.h"

#include "bilinear_certifier.h"
#include "certificate_writer.h"
#include "interval_geometry.h"


#define FIXTURE_PRECISION 256
#define FIXTURE_LARGE_T_N 2
#define FIXTURE_MAX_DEPTH 20


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
 * Small adaptive generator used only for the integration fixture.
 *
 * It deliberately does not use interval_engine: the purpose here is
 * to test certificate_writer + verify_certificate independently.
 */
static int
certify_and_write(
    certificate_writer_t *writer,
    const interval_domain_t *domain,
    const interval_leaf_t *leaf,
    const arb_t kappa,
    slong prec,
    ulong max_depth
)
{
    arb_t interval;

    bilinear_certificate_t certificate;

    interval_leaf_t left;
    interval_leaf_t right;

    int completed;
    int success = 0;

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


    if (!interval_leaf_arb(
            interval,
            domain,
            leaf,
            prec))
    {
        goto cleanup;
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
        success =
            certificate_writer_write_leaf(
                writer,
                leaf
            );

        goto cleanup;
    }


    if (completed &&
        certificate.status !=
            BILINEAR_STATUS_INCONCLUSIVE)
    {
        goto cleanup;
    }


    if (leaf->depth >=
        max_depth)
    {
        goto cleanup;
    }


    if (!interval_leaf_bisect(
            &left,
            &right,
            leaf))
    {
        goto cleanup;
    }


    /*
     * Left first, then right, so the certificate is emitted
     * exactly in geometric left-to-right order.
     */
    if (!certify_and_write(
            writer,
            domain,
            &left,
            kappa,
            prec,
            max_depth))
    {
        goto cleanup;
    }


    if (!certify_and_write(
            writer,
            domain,
            &right,
            kappa,
            prec,
            max_depth))
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

    arb_clear(interval);

    return success;
}


static int
write_bad_gap(
    const char *path
)
{
    FILE *file =
        fopen(path, "w");

    if (file == NULL)
        return 0;

    fprintf(
        file,
        "RH_FINITE_CERTIFICATE 1\n"
        "large_t_N 2\n"
        "kappa_num -100\n"
        "kappa_den 1\n"
        "generator_precision 256\n"
        "leaf_count 2\n"
        "BEGIN_LEAVES\n"
        "2 0\n"
        "1 1\n"
        "END_LEAVES\n"
    );

    return fclose(file) == 0;
}


static int
write_bad_overlap(
    const char *path
)
{
    FILE *file =
        fopen(path, "w");

    if (file == NULL)
        return 0;

    fprintf(
        file,
        "RH_FINITE_CERTIFICATE 1\n"
        "large_t_N 2\n"
        "kappa_num -100\n"
        "kappa_den 1\n"
        "generator_precision 256\n"
        "leaf_count 3\n"
        "BEGIN_LEAVES\n"
        "1 0\n"
        "2 1\n"
        "1 1\n"
        "END_LEAVES\n"
    );

    return fclose(file) == 0;
}


static int
write_bad_duplicate(
    const char *path
)
{
    FILE *file =
        fopen(path, "w");

    if (file == NULL)
        return 0;

    fprintf(
        file,
        "RH_FINITE_CERTIFICATE 1\n"
        "large_t_N 2\n"
        "kappa_num -100\n"
        "kappa_den 1\n"
        "generator_precision 256\n"
        "leaf_count 3\n"
        "BEGIN_LEAVES\n"
        "1 0\n"
        "1 0\n"
        "1 1\n"
        "END_LEAVES\n"
    );

    return fclose(file) == 0;
}


static int
write_bad_count(
    const char *path
)
{
    FILE *file =
        fopen(path, "w");

    if (file == NULL)
        return 0;

    fprintf(
        file,
        "RH_FINITE_CERTIFICATE 1\n"
        "large_t_N 2\n"
        "kappa_num -100\n"
        "kappa_den 1\n"
        "generator_precision 256\n"
        "leaf_count 2\n"
        "BEGIN_LEAVES\n"
        "0 0\n"
        "END_LEAVES\n"
    );

    return fclose(file) == 0;
}


int
main(
    int argc,
    char **argv
)
{
    char valid_path[4096];
    char gap_path[4096];
    char overlap_path[4096];
    char duplicate_path[4096];
    char count_path[4096];

    interval_domain_t domain;
    interval_leaf_t root;

    certificate_writer_t writer;

    arb_t kappa;

    fmpz_t kappa_num;
    fmpz_t kappa_den;

    int success = 0;


    if (argc != 2)
    {
        fprintf(
            stderr,
            "Usage: %s output_directory\n",
            argv[0]
        );

        return 2;
    }


    if (snprintf(
            valid_path,
            sizeof(valid_path),
            "%s/valid_certificate.txt",
            argv[1]
        ) >= (int) sizeof(valid_path) ||
        snprintf(
            gap_path,
            sizeof(gap_path),
            "%s/bad_gap.txt",
            argv[1]
        ) >= (int) sizeof(gap_path) ||
        snprintf(
            overlap_path,
            sizeof(overlap_path),
            "%s/bad_overlap.txt",
            argv[1]
        ) >= (int) sizeof(overlap_path) ||
        snprintf(
            duplicate_path,
            sizeof(duplicate_path),
            "%s/bad_duplicate.txt",
            argv[1]
        ) >= (int) sizeof(duplicate_path) ||
        snprintf(
            count_path,
            sizeof(count_path),
            "%s/bad_count.txt",
            argv[1]
        ) >= (int) sizeof(count_path))
    {
        fprintf(
            stderr,
            "Fixture path is too long.\n"
        );

        return 2;
    }


    interval_domain_init(
        &domain
    );

    interval_leaf_init(
        &root
    );

    certificate_writer_init(
        &writer
    );

    arb_init(kappa);

    fmpz_init(kappa_num);
    fmpz_init(kappa_den);


    /*
     * Exact fixture kappa = -100.
     *
     * This value is chosen solely to make the integration
     * certificate quick to generate.  It is not a theorem
     * parameter.
     */
    fmpz_set_si(
        kappa_num,
        -100
    );

    fmpz_one(
        kappa_den
    );

    arb_set_si(
        kappa,
        -100
    );


    if (!build_finite_domain(
            &domain,
            FIXTURE_LARGE_T_N,
            FIXTURE_PRECISION))
    {
        fprintf(
            stderr,
            "Could not build fixture domain.\n"
        );

        goto cleanup;
    }


    interval_leaf_set_root(
        &root
    );


    if (!certificate_writer_open(
            &writer,
            valid_path,
            FIXTURE_LARGE_T_N,
            kappa_num,
            kappa_den,
            FIXTURE_PRECISION))
    {
        fprintf(
            stderr,
            "Could not open valid fixture certificate.\n"
        );

        goto cleanup;
    }


    if (!certify_and_write(
            &writer,
            &domain,
            &root,
            kappa,
            FIXTURE_PRECISION,
            FIXTURE_MAX_DEPTH))
    {
        fprintf(
            stderr,
            "Could not generate valid fixture certificate.\n"
        );

        certificate_writer_abort(
            &writer
        );

        goto cleanup;
    }


    if (!certificate_writer_finalize(
            &writer))
    {
        fprintf(
            stderr,
            "Could not finalize valid fixture certificate.\n"
        );

        goto cleanup;
    }


    if (!write_bad_gap(gap_path) ||
        !write_bad_overlap(overlap_path) ||
        !write_bad_duplicate(duplicate_path) ||
        !write_bad_count(count_path))
    {
        fprintf(
            stderr,
            "Could not create malformed fixtures.\n"
        );

        goto cleanup;
    }


    printf(
        "Certificate fixtures created successfully.\n"
        "valid      : %s\n"
        "gap        : %s\n"
        "overlap    : %s\n"
        "duplicate  : %s\n"
        "bad count  : %s\n",
        valid_path,
        gap_path,
        overlap_path,
        duplicate_path,
        count_path
    );


    success = 1;


cleanup:

    if (writer.is_open)
    {
        certificate_writer_abort(
            &writer
        );
    }

    fmpz_clear(kappa_den);
    fmpz_clear(kappa_num);

    arb_clear(kappa);

    interval_leaf_clear(
        &root
    );

    interval_domain_clear(
        &domain
    );

    return success ? 0 : 1;
}