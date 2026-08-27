#include <stdio.h>

#include "arb.h"
#include "lambda_model.h"

#define TEST_PREC 256

static int failures = 0;


static void
print_ball(
    const char *name,
    const arb_t x
)
{
    printf("  %-36s = ", name);
    arb_printn(x, 30, 0);
    printf("\n");
}


static void
check_true(
    const char *name,
    int condition
)
{
    printf("[%-56s] ", name);

    if (condition)
    {
        printf("PASS\n");
    }
    else
    {
        printf("FAIL\n");
        failures++;
    }
}


static void
check_overlap(
    const char *name,
    const arb_t x,
    const arb_t y
)
{
    printf("[%-56s] ", name);

    if (arb_overlaps(x, y))
    {
        printf("PASS\n");
    }
    else
    {
        printf("FAIL\n");

        print_ball("x", x);
        print_ball("y", y);

        failures++;
    }
}


static void
test_lambda0_at_zero(void)
{
    arb_t t;
    arb_t value;
    arb_t expected;

    arb_t pi;
    arb_t catalan;
    arb_t tmp;

    arb_init(t);
    arb_init(value);
    arb_init(expected);

    arb_init(pi);
    arb_init(catalan);
    arb_init(tmp);

    printf("\n=== lambda_0 at t = 0 ===\n");

    arb_zero(t);

    lambda_0(
        value,
        t,
        TEST_PREC
    );

    /*
     * The special value
     *
     *   psi^(1)(1/4) = pi^2 + 8 G
     *
     * gives
     *
     *   lambda_0(0)
     *     = -8 + pi^2/4 + 2 G,
     *
     * where G is Catalan's constant.
     */

    arb_const_pi(pi, TEST_PREC);
    arb_const_catalan(catalan, TEST_PREC);

    arb_mul(expected, pi, pi, TEST_PREC);
    arb_div_ui(expected, expected, 4, TEST_PREC);

    arb_mul_ui(tmp, catalan, 2, TEST_PREC);
    arb_add(expected, expected, tmp, TEST_PREC);

    arb_sub_ui(expected, expected, 8, TEST_PREC);

    print_ball("lambda_0(0)", value);
    print_ball("special-value expression", expected);

    check_overlap(
        "lambda_0(0) matches special trigamma identity",
        value,
        expected
    );

    arb_clear(tmp);
    arb_clear(catalan);
    arb_clear(pi);

    arb_clear(expected);
    arb_clear(value);
    arb_clear(t);
}


static void
test_evenness(void)
{
    arb_t t_pos;
    arb_t t_neg;

    arb_t value_pos;
    arb_t value_neg;

    arb_init(t_pos);
    arb_init(t_neg);

    arb_init(value_pos);
    arb_init(value_neg);

    printf("\n=== Evenness of lambda_0 ===\n");

    arb_set_ui(t_pos, 37);
    arb_neg(t_neg, t_pos);

    lambda_0(
        value_pos,
        t_pos,
        TEST_PREC
    );

    lambda_0(
        value_neg,
        t_neg,
        TEST_PREC
    );

    print_ball("lambda_0(37)", value_pos);
    print_ball("lambda_0(-37)", value_neg);

    check_overlap(
        "lambda_0(-t) = lambda_0(t)",
        value_pos,
        value_neg
    );

    arb_clear(value_neg);
    arb_clear(value_pos);

    arb_clear(t_neg);
    arb_clear(t_pos);
}


static void
test_kappa_shift(void)
{
    arb_t t;
    arb_t kappa;

    arb_t lambda0_value;
    arb_t lambda_kappa_value;
    arb_t expected;
    arb_t tmp;

    arb_init(t);
    arb_init(kappa);

    arb_init(lambda0_value);
    arb_init(lambda_kappa_value);
    arb_init(expected);
    arb_init(tmp);

    printf("\n=== kappa shift ===\n");

    arb_set_ui(t, 50);

    /*
     * kappa = 1/8.
     */
    arb_one(kappa);
    arb_div_ui(kappa, kappa, 8, TEST_PREC);

    lambda_0(
        lambda0_value,
        t,
        TEST_PREC
    );

    lambda_kappa(
        lambda_kappa_value,
        t,
        kappa,
        TEST_PREC
    );

    /*
     * expected = lambda_0 - kappa/2.
     */
    arb_div_ui(tmp, kappa, 2, TEST_PREC);
    arb_sub(
        expected,
        lambda0_value,
        tmp,
        TEST_PREC
    );

    print_ball("lambda_0(50)", lambda0_value);
    print_ball("lambda_kappa(50)", lambda_kappa_value);

    check_overlap(
        "lambda_kappa = lambda_0 - kappa/2",
        lambda_kappa_value,
        expected
    );

    arb_clear(tmp);
    arb_clear(expected);
    arb_clear(lambda_kappa_value);
    arb_clear(lambda0_value);

    arb_clear(kappa);
    arb_clear(t);
}


static void
test_large_t_bound(void)
{
    arb_t t;
    arb_t t2;

    arb_t value;
    arb_t lower;
    arb_t upper;

    arb_init(t);
    arb_init(t2);

    arb_init(value);
    arb_init(lower);
    arb_init(upper);

    printf("\n=== Audited large-t bound at t = 200 ===\n");

    arb_set_ui(t, 200);

    lambda_0(
        value,
        t,
        TEST_PREC
    );

    arb_mul(t2, t, t, TEST_PREC);

    /*
     * lower = 0.99 / t^2.
     */
    arb_set_ui(lower, 99);
    arb_div_ui(lower, lower, 100, TEST_PREC);
    arb_div(lower, lower, t2, TEST_PREC);

    /*
     * upper = 2.01 / t^2.
     */
    arb_set_ui(upper, 201);
    arb_div_ui(upper, upper, 100, TEST_PREC);
    arb_div(upper, upper, t2, TEST_PREC);

    print_ball("0.99/t^2", lower);
    print_ball("lambda_0(200)", value);
    print_ball("2.01/t^2", upper);

    check_true(
        "0.99/t^2 < lambda_0(200)",
        arb_lt(lower, value)
    );

    check_true(
        "lambda_0(200) < 2.01/t^2",
        arb_lt(value, upper)
    );

    check_true(
        "lambda_0(200) is positive",
        arb_is_positive(value)
    );

    arb_clear(upper);
    arb_clear(lower);
    arb_clear(value);

    arb_clear(t2);
    arb_clear(t);
}


static void
test_interval_evaluation(void)
{
    arb_t t;
    arb_t value;

    arb_t left_point;
    arb_t right_point;

    arb_t left_value;
    arb_t right_value;

    arb_init(t);
    arb_init(value);

    arb_init(left_point);
    arb_init(right_point);

    arb_init(left_value);
    arb_init(right_value);

    printf("\n=== Interval evaluation ===\n");

    /*
     * t = 200 +/- 2^-10.
     */
    arb_set_ui(t, 200);
    arb_add_error_2exp_si(t, -10);

    lambda_0(
        value,
        t,
        TEST_PREC
    );

    arb_set_ui(left_point, 200);
    arb_sub_ui(left_point, left_point, 1, TEST_PREC);
    arb_add_ui(left_point, left_point, 1, TEST_PREC);

    /*
     * Construct the two rational test points
     *
     *   200 - 2^-10,
     *   200 + 2^-10.
     */
    arb_set_ui(left_point, 204800 - 1);
    arb_div_ui(
        left_point,
        left_point,
        1024,
        TEST_PREC
    );

    arb_set_ui(right_point, 204800 + 1);
    arb_div_ui(
        right_point,
        right_point,
        1024,
        TEST_PREC
    );

    lambda_0(
        left_value,
        left_point,
        TEST_PREC
    );

    lambda_0(
        right_value,
        right_point,
        TEST_PREC
    );

    print_ball("lambda_0(interval)", value);

    check_true(
        "interval result contains left-endpoint value",
        arb_contains(value, left_value)
    );

    check_true(
        "interval result contains right-endpoint value",
        arb_contains(value, right_value)
    );

    arb_clear(right_value);
    arb_clear(left_value);

    arb_clear(right_point);
    arb_clear(left_point);

    arb_clear(value);
    arb_clear(t);
}


int
main(void)
{
    printf("============================================================\n");
    printf("LAMBDA MODEL TESTS\n");
    printf("Arb/FLINT rigorous ball arithmetic\n");
    printf("precision = %d bits\n", TEST_PREC);
    printf("============================================================\n");

    test_lambda0_at_zero();
    test_evenness();
    test_kappa_shift();
    test_large_t_bound();
    test_interval_evaluation();

    printf("\n============================================================\n");

    if (failures == 0)
    {
        printf("ALL LAMBDA MODEL TESTS PASSED\n");
        printf("============================================================\n");
        return 0;
    }

    printf(
        "%d LAMBDA MODEL TEST(S) FAILED\n",
        failures
    );

    printf("============================================================\n");

    return 1;
}