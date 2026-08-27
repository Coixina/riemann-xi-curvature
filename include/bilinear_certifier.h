#ifndef BILINEAR_CERTIFIER_H
#define BILINEAR_CERTIFIER_H

#include "arb.h"
#include "arf.h"


/*
 * Logical result of one interval certification attempt.
 *
 * INVALID:
 *   The rigorous computation could not be completed.
 *
 * CERTIFIED:
 *   The exact Hardy bilinear functional is rigorously
 *   positive on the complete input interval.
 *
 * INCONCLUSIVE:
 *   The computation was rigorous, but the current interval
 *   was too wide to prove positivity.
 *
 * INCONCLUSIVE never means that the functional is negative.
 */
typedef enum
{
    BILINEAR_STATUS_INVALID = 0,
    BILINEAR_STATUS_CERTIFIED,
    BILINEAR_STATUS_INCONCLUSIVE

} bilinear_status_t;


/*
 * Rigorous certificate for
 *
 *     B_kappa[Z]
 *       =
 *     Z'^2 - Z Z'' + lambda_kappa Z^2
 *
 * on one complete Arb interval.
 *
 * value is a rigorous Arb enclosure of the functional.
 *
 * certified_lower satisfies
 *
 *     certified_lower
 *       <=
 *     inf_I B_kappa[Z].
 */
typedef struct
{
    arb_t value;
    arb_t lambda;

    arf_t certified_lower;

    bilinear_status_t status;

} bilinear_certificate_t;


/*
 * Initialize and clear a certificate object.
 */
void bilinear_certificate_init(
    bilinear_certificate_t *certificate
);

void bilinear_certificate_clear(
    bilinear_certificate_t *certificate
);


/*
 * Attempt to certify
 *
 *     B_kappa[Z](t) > 0
 *
 * uniformly on the complete real Arb interval t.
 *
 * The Hardy-Z jet is evaluated directly and rigorously by
 * hardy_jet.
 *
 * kappa is an Arb value. In particular, kappa = 0 can be
 * supplied exactly.
 *
 * Returns 1 if the rigorous computation completed.
 * Returns 0 on invalid input or if the direct Hardy interval
 * evaluation did not produce a usable finite enclosure.
 *
 * If the function returns 1, inspect certificate->status:
 *
 *     BILINEAR_STATUS_CERTIFIED
 *         certified_lower > 0.
 *
 *     BILINEAR_STATUS_INCONCLUSIVE
 *         a rigorous enclosure was obtained, but positivity
 *         was not proved on this interval.
 */
int bilinear_certify(
    bilinear_certificate_t *certificate,
    const arb_t t,
    const arb_t kappa,
    slong prec
);


#endif