from testing import assert_true, TestSuite
from src.modular import find_suitable_q
from src.ntt import apply_xyw_quotient_ntt_inplace
from collections import List
from random import random_si64


fn test_xyw_ring_inplace_ntt() raises:
    print("Running test_xyw_ring_inplace_ntt...")

    var n_power_of_2 = 4
    var p_power_of_3 = 9
    var target_bit_length = 20
    var q_modulus = find_suitable_q(
        n_power_of_2 * n_power_of_2, p_power_of_3, target_bit_length
    )

    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var total_length = 2 * n_power_of_2 * n_power_of_2 * phi_p_degree
    var coefficient_values = List[Int32]()
    for coefficient_index in range(total_length):
        coefficient_values.append(Int32(Int(random_si64(0, q_modulus - 1))))

    var transform_succeeded = apply_xyw_quotient_ntt_inplace(
        coefficient_values, n_power_of_2, p_power_of_3, q_modulus
    )
    assert_true(transform_succeeded, "In-place XYW quotient NTT should succeed")

    for coefficient_index in range(total_length):
        var transformed_value = Int(coefficient_values[coefficient_index])
        assert_true(
            transformed_value >= 0 and transformed_value < q_modulus,
            "Transformed coefficient must remain reduced modulo q",
        )

    print("test_xyw_ring_inplace_ntt passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
