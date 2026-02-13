from testing import assert_true, TestSuite
from src.modular import find_suitable_q
from src.polynomial import Polynomial
from collections import List
from random import random_si64


fn test_multidim_polynomial_multiplication() raises:
    print("Running test_multidim_polynomial_multiplication...")

    # Parameters
    var n_dim = 4
    var p_dim = 9  # Phi_9 degree = 6
    var target_bits = 20

    var q_modulus_u64 = find_suitable_q(n_dim, p_dim, target_bits)
    var q_modulus = UInt32(q_modulus_u64)
    print("Found suitable q:", q_modulus)

    # Create two polynomials
    var poly_first = Polynomial(n_dim, p_dim, q_modulus)
    var poly_second = Polynomial(n_dim, p_dim, q_modulus)

    # Initialize with some data (just random coefficients)
    for index in range(poly_first.total_length):
        var rand_coeff_1 = Int(random_si64(0, Int(q_modulus) - 1))
        var rand_coeff_2 = Int(random_si64(0, Int(q_modulus) - 1))
        poly_first.set_coefficient(index, rand_coeff_1)
        poly_second.set_coefficient(index, rand_coeff_2)

    # Perform multiplication (lazy NTT triggers here)
    var poly_product = poly_first.multiply(poly_second)

    assert_true(poly_product.is_in_ntt_i, "Product should be in NTT form (i)")
    assert_true(poly_product.is_in_ntt_x, "Product should be in NTT form (X)")
    assert_true(poly_product.is_in_ntt_w, "Product should be in NTT form (W)")

    print("test_multidim_polynomial_multiplication passed basic execution.")


fn test_polynomial_add_and_multiply_correctness() raises:
    print("Running test_polynomial_add_and_multiply_correctness...")

    var n_dim = 4
    var p_dim = 9
    var target_bits = 20
    var q_modulus_u64 = find_suitable_q(n_dim, p_dim, target_bits)
    var q_modulus = UInt32(q_modulus_u64)

    var polynomial_left = Polynomial(n_dim, p_dim, q_modulus)
    var polynomial_right = Polynomial(n_dim, p_dim, q_modulus)
    for coefficient_index in range(polynomial_left.total_length):
        var left_seed = (coefficient_index * 37 + Int(q_modulus) - 5) % Int(
            q_modulus
        )
        var right_seed = (coefficient_index * 53 + Int(q_modulus) - 7) % Int(
            q_modulus
        )
        polynomial_left.set_coefficient(coefficient_index, left_seed)
        polynomial_right.set_coefficient(coefficient_index, right_seed)

    var sum_poly = polynomial_left.add(polynomial_right)
    for coefficient_index in range(sum_poly.total_length):
        var left_value = polynomial_left.coefficient_values[coefficient_index]
        var right_value = polynomial_right.coefficient_values[coefficient_index]
        var expected_sum = left_value + right_value
        if expected_sum >= q_modulus:
            expected_sum -= q_modulus
        assert_true(
            sum_poly.coefficient_values[coefficient_index] == expected_sum,
            "add() should match conditional-subtraction modulo reduction",
        )

    var inplace_sum_left = Polynomial(n_dim, p_dim, q_modulus)
    var inplace_sum_right = Polynomial(n_dim, p_dim, q_modulus)
    for coefficient_index in range(inplace_sum_left.total_length):
        var left_seed = (coefficient_index * 37 + Int(q_modulus) - 5) % Int(
            q_modulus
        )
        var right_seed = (coefficient_index * 53 + Int(q_modulus) - 7) % Int(
            q_modulus
        )
        inplace_sum_left.set_coefficient(coefficient_index, left_seed)
        inplace_sum_right.set_coefficient(coefficient_index, right_seed)
    inplace_sum_left.add_inplace(inplace_sum_right)
    for coefficient_index in range(inplace_sum_left.total_length):
        assert_true(
            inplace_sum_left.coefficient_values[coefficient_index]
            == sum_poly.coefficient_values[coefficient_index],
            "add_inplace() should match add()",
        )

    var multiply_left = Polynomial(n_dim, p_dim, q_modulus)
    var multiply_right = Polynomial(n_dim, p_dim, q_modulus)
    for coefficient_index in range(multiply_left.total_length):
        var left_seed = (coefficient_index * 37 + Int(q_modulus) - 5) % Int(
            q_modulus
        )
        var right_seed = (coefficient_index * 53 + Int(q_modulus) - 7) % Int(
            q_modulus
        )
        multiply_left.set_coefficient(coefficient_index, left_seed)
        multiply_right.set_coefficient(coefficient_index, right_seed)

    var expected_left = Polynomial(n_dim, p_dim, q_modulus)
    var expected_right = Polynomial(n_dim, p_dim, q_modulus)
    for coefficient_index in range(expected_left.total_length):
        var left_seed = (coefficient_index * 37 + Int(q_modulus) - 5) % Int(
            q_modulus
        )
        var right_seed = (coefficient_index * 53 + Int(q_modulus) - 7) % Int(
            q_modulus
        )
        expected_left.set_coefficient(coefficient_index, left_seed)
        expected_right.set_coefficient(coefficient_index, right_seed)
    expected_left.transform_to_full_ntt()
    expected_right.transform_to_full_ntt()

    var product_poly = multiply_left.multiply(multiply_right)
    for coefficient_index in range(product_poly.total_length):
        var expected_product = UInt32(
            (
                UInt64(expected_left.coefficient_values[coefficient_index])
                * UInt64(expected_right.coefficient_values[coefficient_index])
            )
            % UInt64(q_modulus)
        )
        assert_true(
            product_poly.coefficient_values[coefficient_index]
            == expected_product,
            "multiply() should equal pointwise product in transformed domain",
        )

    print("test_polynomial_add_and_multiply_correctness passed.")


fn test_polynomial_multiplication_full_trials() raises:
    print("Running test_polynomial_multiplication_full_trials...")

    var n_dim = 4
    var p_dim = 9
    var target_bits = 20
    var q_modulus_u64 = find_suitable_q(n_dim, p_dim, target_bits)
    var q_modulus = UInt32(q_modulus_u64)
    var number_of_trials = 6

    for trial_index in range(number_of_trials):
        _ = trial_index
        var poly_left = Polynomial(n_dim, p_dim, q_modulus)
        var poly_right = Polynomial(n_dim, p_dim, q_modulus)
        var poly_left_reference = Polynomial(n_dim, p_dim, q_modulus)
        var poly_right_reference = Polynomial(n_dim, p_dim, q_modulus)
        var left_input_values = List[Int]()
        var right_input_values = List[Int]()

        for coefficient_index in range(poly_left.total_length):
            var left_random_value = Int(random_si64(0, Int(q_modulus) - 1))
            var right_random_value = Int(random_si64(0, Int(q_modulus) - 1))
            left_input_values.append(left_random_value)
            right_input_values.append(right_random_value)
            poly_left.set_coefficient(coefficient_index, left_random_value)
            poly_right.set_coefficient(coefficient_index, right_random_value)
            poly_left_reference.set_coefficient(
                coefficient_index, left_random_value
            )
            poly_right_reference.set_coefficient(
                coefficient_index, right_random_value
            )

        poly_left_reference.transform_to_full_ntt()
        poly_right_reference.transform_to_full_ntt()

        var poly_product_left_right = poly_left.multiply(poly_right)
        for coefficient_index in range(poly_product_left_right.total_length):
            var expected_product = UInt32(
                (
                    UInt64(
                        poly_left_reference.coefficient_values[
                            coefficient_index
                        ]
                    )
                    * UInt64(
                        poly_right_reference.coefficient_values[
                            coefficient_index
                        ]
                    )
                )
                % UInt64(q_modulus)
            )
            assert_true(
                poly_product_left_right.coefficient_values[coefficient_index]
                == expected_product,
                "multiply() must match explicit transformed-domain reference",
            )

        var poly_left_commute = Polynomial(n_dim, p_dim, q_modulus)
        var poly_right_commute = Polynomial(n_dim, p_dim, q_modulus)
        for coefficient_index in range(poly_left_commute.total_length):
            poly_left_commute.set_coefficient(
                coefficient_index, left_input_values[coefficient_index]
            )
            poly_right_commute.set_coefficient(
                coefficient_index, right_input_values[coefficient_index]
            )

        var poly_product_right_left = poly_right_commute.multiply(
            poly_left_commute
        )
        for coefficient_index in range(poly_product_left_right.total_length):
            assert_true(
                poly_product_left_right.coefficient_values[coefficient_index]
                == poly_product_right_left.coefficient_values[
                    coefficient_index
                ],
                "multiplication should be commutative in NTT pointwise form",
            )

        var poly_left_for_zero = Polynomial(n_dim, p_dim, q_modulus)
        var poly_zero_for_test = Polynomial(n_dim, p_dim, q_modulus)
        for coefficient_index in range(poly_left_for_zero.total_length):
            poly_left_for_zero.set_coefficient(
                coefficient_index, left_input_values[coefficient_index]
            )
            poly_zero_for_test.set_coefficient(coefficient_index, 0)
        var poly_product_with_zero = poly_left_for_zero.multiply(
            poly_zero_for_test
        )
        for coefficient_index in range(poly_product_with_zero.total_length):
            assert_true(
                poly_product_with_zero.coefficient_values[coefficient_index]
                == 0,
                "polynomial multiplied by zero should remain zero",
            )

    print("test_polynomial_multiplication_full_trials passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
