#ifndef LAMBDA_MODEL_H
#define LAMBDA_MODEL_H

#include "arb.h"

/*
 * Hardy prefactor logarithmic-curvature model.
 *
 * With
 *
 *   Xi(t) = A(t) Z(t),
 *
 *   A(t)
 *     = c_A (t^2 + 1/4)
 *       |Gamma(1/4 + i t/2)|,
 *
 * define
 *
 *   lambda_0(t) = -(log A(t))'',
 *
 * and
 *
 *   lambda_kappa(t)
 *     = lambda_0(t) - kappa/2.
 *
 * Exact formula used for rigorous evaluation:
 *
 *   lambda_0(t)
 *     = 2 (t^2 - 1/4) / (t^2 + 1/4)^2
 *       + (1/4) Re psi^(1)(1/4 + i t/2),
 *
 * where psi^(1) is the trigamma function.
 */


/*
 * Compute lambda_0(t) rigorously for an Arb ball t.
 *
 * The input may be either a point ball or a genuine interval.
 */
void lambda_0(
    arb_t result,
    const arb_t t,
    slong prec
);


/*
 * Compute
 *
 *   lambda_kappa(t) = lambda_0(t) - kappa/2
 *
 * rigorously for Arb balls t and kappa.
 */
void lambda_kappa(
    arb_t result,
    const arb_t t,
    const arb_t kappa,
    slong prec
);


#endif