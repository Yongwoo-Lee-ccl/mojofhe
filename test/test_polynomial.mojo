from testing import assert_true, TestSuite
from src.modular import find_suitable_q, find_primitive_root
from src.polynomial import Polynomial
from collections import List
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


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
