from testing import assert_true, TestSuite
from src.modular import find_suitable_q
from src.polynomial import Polynomial


fn test_sample_ternary() raises:
    print("Running test_sample_ternary...")
    var n_dim = 4
    var p_dim = 9
    var target_bits = 20
    var q_modulus = find_suitable_q(n_dim, p_dim, target_bits)

    var poly = Polynomial(n_dim, p_dim, q_modulus)
    poly.sample_ternary()

    for i in range(poly.total_length):
        var val = Int(poly.coefficient_values[i])
        assert_true(
            val == 0 or val == 1 or val == q_modulus - 1,
            "Ternary sample should be 0, 1, or q-1",
        )


fn test_sample_gaussian() raises:
    print("Running test_sample_gaussian...")
    var n_dim = 4
    var p_dim = 9
    var target_bits = 20
    var q_modulus = find_suitable_q(n_dim, p_dim, target_bits)

    var poly = Polynomial(n_dim, p_dim, q_modulus)
    poly.sample_gaussian()

    # Gaussian sample should be within a reasonable range
    # sum of 42 bits (bit1 - bit2) has max 42 and min -42.
    for i in range(poly.total_length):
        var val = Int(poly.coefficient_values[i])
        if val > q_modulus // 2:
            val -= q_modulus
        assert_true(
            val >= -42 and val <= 42,
            "Gaussian sample should be within [-42, 42]",
        )


fn test_sample_sparse() raises:
    print("Running test_sample_sparse...")
    var n_dim = 4
    var p_dim = 9
    var target_bits = 20
    var q_modulus = find_suitable_q(n_dim, p_dim, target_bits)

    var poly = Polynomial(n_dim, p_dim, q_modulus)
    var hamming_weight = 5
    poly.sample_sparse(hamming_weight)

    var count = 0
    for i in range(poly.total_length):
        var val = Int(poly.coefficient_values[i])
        if val != 0:
            assert_true(
                val == 1 or val == q_modulus - 1,
                "Sparse sample non-zero should be 1 or q-1",
            )
            count += 1

    assert_true(
        count == hamming_weight, "Hamming weight should match input weight"
    )


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
