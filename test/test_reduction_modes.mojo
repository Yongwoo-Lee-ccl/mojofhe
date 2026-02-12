from testing import assert_true, TestSuite
from src.modular import find_suitable_q
from src.ntt import apply_ixw_quotient_ntt_inplace
from collections import List
from random import random_si64


fn test_barrett_and_montgomery_match_on_ntt_workload() raises:
    print("Running test_barrett_and_montgomery_match_on_ntt_workload...")

    var n_power_of_2 = 1 << 7
    var p_power_of_3 = 81  # 3^4
    var target_bits = 25
    var q_modulus = find_suitable_q(n_power_of_2, p_power_of_3, target_bits)

    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var total_length = 2 * n_power_of_2 * phi_p_degree

    var input_coefficients = List[Int32]()
    var barrett_output = List[Int32]()
    var montgomery_output = List[Int32]()
    for coefficient_index in range(total_length):
        _ = coefficient_index
        var random_value = Int32(Int(random_si64(0, q_modulus - 1)))
        input_coefficients.append(random_value)
        barrett_output.append(random_value)
        montgomery_output.append(random_value)

    var barrett_success = apply_ixw_quotient_ntt_inplace(
        barrett_output, n_power_of_2, p_power_of_3, q_modulus, False
    )
    var montgomery_success = apply_ixw_quotient_ntt_inplace(
        montgomery_output, n_power_of_2, p_power_of_3, q_modulus, True
    )
    assert_true(barrett_success, "Barrett NTT should succeed")
    assert_true(montgomery_success, "Montgomery NTT should succeed")

    for coefficient_index in range(total_length):
        assert_true(
            barrett_output[coefficient_index]
            == montgomery_output[coefficient_index],
            "Barrett and Montgomery NTT outputs should match",
        )

    print("test_barrett_and_montgomery_match_on_ntt_workload passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
