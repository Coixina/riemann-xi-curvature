#ifndef HARDY_JET_H
#define HARDY_JET_H

#include "arb.h"


/*
 * Rigorous real Taylor jet of the Hardy Z-function:
 *
 *     value = Z(t),
 *     dt    = Z'(t),
 *     d2t   = Z''(t).
 *
 * Each component is an Arb enclosure valid on the complete
 * input Arb ball t.
 */
typedef struct
{
    arb_t value;
    arb_t dt;
    arb_t d2t;

} hardy_jet_t;


/*
 * Initialize and clear a Hardy-Z jet.
 */
void hardy_jet_init(
    hardy_jet_t *jet
);

void hardy_jet_clear(
    hardy_jet_t *jet
);


/*
 * Rigorously evaluate
 *
 *     Z(t),
 *     Z'(t),
 *     Z''(t)
 *
 * on the complete real Arb ball t.
 *
 * FLINT returns Taylor coefficients
 *
 *     Z(t+x)
 *       =
 *     c_0 + c_1 x + c_2 x^2 + ...
 *
 * so that
 *
 *     Z   = c_0,
 *     Z'  = c_1,
 *     Z'' = 2 c_2.
 *
 * The evaluation is performed internally with Acb.  Since
 * Hardy Z is real on the real axis, every imaginary enclosure
 * must contain zero before its real component is accepted.
 *
 * There is no additional approximation or remainder budget:
 * the returned Arb balls directly enclose the exact Hardy-Z
 * jet.
 *
 Returns 1 when a finite rigorous enclosure of the complete
 Hardy-Z jet is obtained.
 Returns 0 on invalid input or when the direct interval
 evaluation does not produce a usable finite enclosure.
 A return value of 0 is not a mathematical counterexample;
 an adaptive caller may subdivide the input interval and retry.
 */
int hardy_jet_build(
    hardy_jet_t *jet,
    const arb_t t,
    slong prec
);


#endif