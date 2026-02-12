from testing import assert_true, TestSuite
from src.modular import find_suitable_q
from src.polynomial import Polynomial
from random import random_si64


fn test_multidim_polynomial_multiplication() raises:
    print("Running test_multidim_polynomial_multiplication...")

    # Parameters
    var n_dim = 4
    var p_dim = 9  # Phi_9 degree = 6
    var target_bits = 20

    var q_modulus = find_suitable_q(n_dim, p_dim, target_bits)
    print("Found suitable q:", q_modulus)

    # Create two polynomials
    var poly_first = Polynomial(n_dim, p_dim, q_modulus)
    var poly_second = Polynomial(n_dim, p_dim, q_modulus)

    # Initialize with some data (just random coefficients)
    for index in range(poly_first.total_length):
        poly_first.set_coefficient(index, Int(random_si64(0, q_modulus - 1)))
        poly_second.set_coefficient(index, Int(random_si64(0, q_modulus - 1)))

    # Perform multiplication (lazy NTT triggers here)
    var poly_product = poly_first.multiply(poly_second)

    # In a full test, we would compare against naive multiplication.
    # For this prototype, we verify that the transform steps didn't crash
    # and the metadata is set correctly.
    assert_true(poly_product.is_in_ntt_i, "Product should be in NTT form (i)")
    assert_true(poly_product.is_in_ntt_x, "Product should be in NTT form (X)")
    assert_true(poly_product.is_in_ntt_w, "Product should be in NTT form (W)")

    print("test_multidim_polynomial_multiplication passed basic execution.")


fn test_polynomial_add_and_multiply_correctness() raises:
    print("Running test_polynomial_add_and_multiply_correctness...")

    var n_dim = 4
    var p_dim = 9
    var target_bits = 20
    var q_modulus = find_suitable_q(n_dim, p_dim, target_bits)

    var polynomial_left = Polynomial(n_dim, p_dim, q_modulus)
    var polynomial_right = Polynomial(n_dim, p_dim, q_modulus)
    for coefficient_index in range(polynomial_left.total_length):
        var left_seed = (coefficient_index * 37 + q_modulus - 5) % q_modulus
        var right_seed = (coefficient_index * 53 + q_modulus - 7) % q_modulus
        polynomial_left.set_coefficient(coefficient_index, left_seed)
        polynomial_right.set_coefficient(coefficient_index, right_seed)

    var sum_poly = polynomial_left.add(polynomial_right)
    for coefficient_index in range(sum_poly.total_length):
        var left_value = Int(
            polynomial_left.coefficient_values[coefficient_index]
        )
        var right_value = Int(
            polynomial_right.coefficient_values[coefficient_index]
        )
        var expected_sum = left_value + right_value
        if expected_sum >= q_modulus:
            expected_sum -= q_modulus
        assert_true(
            Int(sum_poly.coefficient_values[coefficient_index]) == expected_sum,
            "add() should match conditional-subtraction modulo reduction",
        )

    var inplace_sum_left = Polynomial(n_dim, p_dim, q_modulus)
    var inplace_sum_right = Polynomial(n_dim, p_dim, q_modulus)
    for coefficient_index in range(inplace_sum_left.total_length):
        var left_seed = (coefficient_index * 37 + q_modulus - 5) % q_modulus
        var right_seed = (coefficient_index * 53 + q_modulus - 7) % q_modulus
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
        var left_seed = (coefficient_index * 37 + q_modulus - 5) % q_modulus
        var right_seed = (coefficient_index * 53 + q_modulus - 7) % q_modulus
        multiply_left.set_coefficient(coefficient_index, left_seed)
        multiply_right.set_coefficient(coefficient_index, right_seed)

    var expected_left = Polynomial(n_dim, p_dim, q_modulus)
    var expected_right = Polynomial(n_dim, p_dim, q_modulus)
    for coefficient_index in range(expected_left.total_length):
        var left_seed = (coefficient_index * 37 + q_modulus - 5) % q_modulus
        var right_seed = (coefficient_index * 53 + q_modulus - 7) % q_modulus
        expected_left.set_coefficient(coefficient_index, left_seed)
        expected_right.set_coefficient(coefficient_index, right_seed)
    expected_left.transform_to_full_ntt()
    expected_right.transform_to_full_ntt()

    var product_poly = multiply_left.multiply(multiply_right)
    for coefficient_index in range(product_poly.total_length):
        var expected_product = Int(
            (
                Int64(Int(expected_left.coefficient_values[coefficient_index]))
                * Int64(
                    Int(expected_right.coefficient_values[coefficient_index])
                )
            )
            % Int64(q_modulus)
        )
        assert_true(
            Int(product_poly.coefficient_values[coefficient_index])
            == expected_product,
            "multiply() should equal pointwise product in transformed domain",
        )

    var inplace_multiply_left = Polynomial(n_dim, p_dim, q_modulus)
    var inplace_multiply_right = Polynomial(n_dim, p_dim, q_modulus)
    for coefficient_index in range(inplace_multiply_left.total_length):
        var left_seed = (coefficient_index * 37 + q_modulus - 5) % q_modulus
        var right_seed = (coefficient_index * 53 + q_modulus - 7) % q_modulus
        inplace_multiply_left.set_coefficient(coefficient_index, left_seed)
        inplace_multiply_right.set_coefficient(coefficient_index, right_seed)
    inplace_multiply_left.multiply_inplace(inplace_multiply_right)
    for coefficient_index in range(inplace_multiply_left.total_length):
        assert_true(
            inplace_multiply_left.coefficient_values[coefficient_index]
            == product_poly.coefficient_values[coefficient_index],
            "multiply_inplace() should match multiply()",
        )

    print("test_polynomial_add_and_multiply_correctness passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
