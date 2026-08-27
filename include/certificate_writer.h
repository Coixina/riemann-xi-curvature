#ifndef CERTIFICATE_WRITER_H
#define CERTIFICATE_WRITER_H

#include <stdio.h>

#include "flint.h"
#include "fmpz.h"

#include "interval_geometry.h"


#define RH_CERTIFICATE_FORMAT_VERSION 1
#define RH_CERTIFICATE_COUNT_WIDTH 20


typedef struct
{
    FILE *file;

    long leaf_count_offset;

    ulong leaf_count;

    int is_open;
    int is_finalized;

} certificate_writer_t;


/*
 * Initialize an inactive writer.
 */
void certificate_writer_init(
    certificate_writer_t *writer
);


/*
 * Begin a new finite certificate.
 *
 * kappa_num / kappa_den represent the exact rational kappa.
 * Requires
 *
 *     large_t_N >= 2,
 *     kappa_den > 0,
 *     generator_precision > 0.
 *
 * The output format is
 *
 *     RH_FINITE_CERTIFICATE 1
 *     large_t_N ...
 *     kappa_num ...
 *     kappa_den ...
 *     generator_precision ...
 *     leaf_count ...
 *     BEGIN_LEAVES
 *     ...
 *     END_LEAVES
 */
int certificate_writer_open(
    certificate_writer_t *writer,
    const char *path,
    ulong large_t_N,
    const fmpz_t kappa_num,
    const fmpz_t kappa_den,
    slong generator_precision
);


/*
 * Append one certified terminal dyadic leaf.
 *
 * Leaves are expected to be supplied in left-to-right order.
 * The independent verifier checks that ordering and the exact
 * partition structure; the writer does not trust or reconstruct
 * geometry here.
 */
int certificate_writer_write_leaf(
    certificate_writer_t *writer,
    const interval_leaf_t *leaf
);


/*
 * Finish the certificate, write END_LEAVES and patch the final
 * leaf count into the header.
 */
int certificate_writer_finalize(
    certificate_writer_t *writer
);


/*
 * Close an unfinished writer.
 *
 * An unfinished file is intentionally not a valid certificate.
 */
void certificate_writer_abort(
    certificate_writer_t *writer
);


#endif