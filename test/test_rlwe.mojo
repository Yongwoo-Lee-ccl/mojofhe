from collections import List
from random import random_si64
from testing import TestSuite, assert_true
from src.modular import find_suitable_q
from src.rlwe import (
    rlwe_keygen,
    rlwe_encrypt_binary,
    rlwe_decrypt_ntt,
    rlwe_ring_multiply_ntt,
    encode_binary_plaintext_ntt,
)


fn _centered_coefficient(value: UInt32, modulus_value: UInt32) -> Int:
    var centered_value = Int(value)
    var half_modulus = Int(modulus_value) // 2
    if centered_value > half_modulus:
        centered_value -= Int(modulus_value)
    return centered_value


fn _poly_add_modulus(
    left_values: List[UInt32],
    right_values: List[UInt32],
    modulus_value: UInt32,
) -> List[UInt32]:
    var result_values = List[UInt32](length=len(left_values), fill=0)
    for coefficient_index in range(len(left_values)):
        var summed_value = (
            left_values[coefficient_index] + right_values[coefficient_index]
        )
        if summed_value >= modulus_value:
            summed_value -= modulus_value
        result_values[coefficient_index] = summed_value
    return result_values^


fn _poly_equal(left_values: List[UInt32], right_values: List[UInt32]) -> Bool:
    if len(left_values) != len(right_values):
        return False
    for coefficient_index in range(len(left_values)):
        if left_values[coefficient_index] != right_values[coefficient_index]:
            return False
    return True


fn _sample_binary_message(polynomial_length: Int) -> List[UInt32]:
    var message_bits = List[UInt32](length=polynomial_length, fill=0)
    for coefficient_index in range(polynomial_length):
        message_bits[coefficient_index] = UInt32(Int(random_si64(0, 1)))
    return message_bits^


fn _poly_mul_ntt(
    left_values: List[UInt32],
    right_values: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    q_modulus: UInt32,
) -> List[UInt32]:
    return rlwe_ring_multiply_ntt(
        left_values,
        right_values,
        n_power_of_2,
        p_power_of_3,
        q_modulus,
    )


fn _expected_decrypt_rhs_ntt(
    message_bits: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    q_modulus: UInt32,
    public_a_ntt_values: List[UInt32],
    public_b_ntt_values: List[UInt32],
    secret_key_ntt_values: List[UInt32],
    ephemeral_u_ntt_values: List[UInt32],
    error_1_ntt_values: List[UInt32],
    error_2_ntt_values: List[UInt32],
) -> List[UInt32]:
    var message_ntt_values = encode_binary_plaintext_ntt(
        message_bits, n_power_of_2, p_power_of_3, q_modulus
    )
    var a_times_secret_ntt_values = rlwe_ring_multiply_ntt(
        public_a_ntt_values,
        secret_key_ntt_values,
        n_power_of_2,
        p_power_of_3,
        q_modulus,
    )
    var key_noise_ntt_values = _poly_add_modulus(
        public_b_ntt_values, a_times_secret_ntt_values, q_modulus
    )
    var key_noise_times_u_ntt_values = rlwe_ring_multiply_ntt(
        key_noise_ntt_values,
        ephemeral_u_ntt_values,
        n_power_of_2,
        p_power_of_3,
        q_modulus,
    )
    var error_2_times_secret_ntt_values = rlwe_ring_multiply_ntt(
        error_2_ntt_values,
        secret_key_ntt_values,
        n_power_of_2,
        p_power_of_3,
        q_modulus,
    )
    return _poly_add_modulus(
        _poly_add_modulus(
            message_ntt_values, key_noise_times_u_ntt_values, q_modulus
        ),
        _poly_add_modulus(
            error_1_ntt_values, error_2_times_secret_ntt_values, q_modulus
        ),
        q_modulus,
    )


fn test_rlwe_secret_key_is_ternary() raises:
    print("Running test_rlwe_secret_key_is_ternary...")

    var n_power_of_2 = 4
    var p_power_of_3 = 9
    var q_modulus = UInt32(find_suitable_q(n_power_of_2, p_power_of_3, 20))
    var key_pair = rlwe_keygen(n_power_of_2, p_power_of_3, q_modulus)

    for coefficient_index in range(len(key_pair.secret_key_coeff_values)):
        var centered_value = _centered_coefficient(
            key_pair.secret_key_coeff_values[coefficient_index], q_modulus
        )
        assert_true(
            centered_value >= -1 and centered_value <= 1,
            "secret key coefficient must be ternary (-1, 0, 1)",
        )

    print("test_rlwe_secret_key_is_ternary passed.")


fn test_rlwe_decrypt_equation_matches_ntt_reference() raises:
    print("Running test_rlwe_decrypt_equation_matches_ntt_reference...")

    var n_power_of_2 = 4
    var p_power_of_3 = 9
    var q_modulus = UInt32(find_suitable_q(n_power_of_2, p_power_of_3, 20))
    var key_pair = rlwe_keygen(n_power_of_2, p_power_of_3, q_modulus)
    var message_bits = _sample_binary_message(
        len(key_pair.secret_key_coeff_values)
    )

    var encryption_result = rlwe_encrypt_binary(message_bits, key_pair)
    var decrypted_ntt_values = rlwe_decrypt_ntt(
        encryption_result.c0_ntt_values,
        encryption_result.c1_ntt_values,
        key_pair,
    )
    var expected_rhs_ntt_values = _expected_decrypt_rhs_ntt(
        message_bits,
        n_power_of_2,
        p_power_of_3,
        q_modulus,
        key_pair.public_a_ntt_values,
        key_pair.public_b_ntt_values,
        key_pair.secret_key_ntt_values,
        encryption_result.ephemeral_u_ntt_values,
        encryption_result.error_1_ntt_values,
        encryption_result.error_2_ntt_values,
    )

    assert_true(
        _poly_equal(decrypted_ntt_values, expected_rhs_ntt_values),
        "decryption in NTT domain must match RLWE equation",
    )

    print("test_rlwe_decrypt_equation_matches_ntt_reference passed.")


fn test_rlwe_decrypt_equation_multiple_trials() raises:
    print("Running test_rlwe_decrypt_equation_multiple_trials...")

    var n_power_of_2 = 4
    var p_power_of_3 = 9
    var q_modulus = UInt32(find_suitable_q(n_power_of_2, p_power_of_3, 20))
    var number_of_trials = 6

    for trial_index in range(number_of_trials):
        _ = trial_index
        var key_pair = rlwe_keygen(n_power_of_2, p_power_of_3, q_modulus)
        var message_bits = _sample_binary_message(
            len(key_pair.secret_key_coeff_values)
        )

        var encryption_result = rlwe_encrypt_binary(message_bits, key_pair)
        var decrypted_ntt_values = rlwe_decrypt_ntt(
            encryption_result.c0_ntt_values,
            encryption_result.c1_ntt_values,
            key_pair,
        )
        var expected_rhs_ntt_values = _expected_decrypt_rhs_ntt(
            message_bits,
            n_power_of_2,
            p_power_of_3,
            q_modulus,
            key_pair.public_a_ntt_values,
            key_pair.public_b_ntt_values,
            key_pair.secret_key_ntt_values,
            encryption_result.ephemeral_u_ntt_values,
            encryption_result.error_1_ntt_values,
            encryption_result.error_2_ntt_values,
        )
        assert_true(
            _poly_equal(decrypted_ntt_values, expected_rhs_ntt_values),
            "all trials must satisfy RLWE decryption equation in NTT domain",
        )

    print("test_rlwe_decrypt_equation_multiple_trials passed.")


fn test_rlwe_ciphertext_addition_homomorphic() raises:
    print("Running test_rlwe_ciphertext_addition_homomorphic...")

    var n_power_of_2 = 4
    var p_power_of_3 = 9
    var q_modulus = UInt32(find_suitable_q(n_power_of_2, p_power_of_3, 20))
    var key_pair = rlwe_keygen(n_power_of_2, p_power_of_3, q_modulus)
    var message_0 = _sample_binary_message(
        len(key_pair.secret_key_coeff_values)
    )
    var message_1 = _sample_binary_message(
        len(key_pair.secret_key_coeff_values)
    )

    var encryption_0 = rlwe_encrypt_binary(message_0, key_pair)
    var encryption_1 = rlwe_encrypt_binary(message_1, key_pair)

    var decrypted_0 = rlwe_decrypt_ntt(
        encryption_0.c0_ntt_values,
        encryption_0.c1_ntt_values,
        key_pair,
    )
    var decrypted_1 = rlwe_decrypt_ntt(
        encryption_1.c0_ntt_values,
        encryption_1.c1_ntt_values,
        key_pair,
    )

    var added_c0 = _poly_add_modulus(
        encryption_0.c0_ntt_values, encryption_1.c0_ntt_values, q_modulus
    )
    var added_c1 = _poly_add_modulus(
        encryption_0.c1_ntt_values, encryption_1.c1_ntt_values, q_modulus
    )
    var decrypted_added = rlwe_decrypt_ntt(added_c0, added_c1, key_pair)
    var expected_added = _poly_add_modulus(decrypted_0, decrypted_1, q_modulus)

    assert_true(
        _poly_equal(decrypted_added, expected_added),
        "ciphertext addition should match addition after decryption",
    )

    print("test_rlwe_ciphertext_addition_homomorphic passed.")


fn test_rlwe_ciphertext_hadamard_multiply_three_component_decrypt() raises:
    print(
        "Running"
        " test_rlwe_ciphertext_hadamard_multiply_three_component_decrypt..."
    )

    var n_power_of_2 = 4
    var p_power_of_3 = 9
    var q_modulus = UInt32(find_suitable_q(n_power_of_2, p_power_of_3, 20))
    var key_pair = rlwe_keygen(n_power_of_2, p_power_of_3, q_modulus)
    var message_0 = _sample_binary_message(
        len(key_pair.secret_key_coeff_values)
    )
    var message_1 = _sample_binary_message(
        len(key_pair.secret_key_coeff_values)
    )

    var encryption_0 = rlwe_encrypt_binary(message_0, key_pair)
    var encryption_1 = rlwe_encrypt_binary(message_1, key_pair)

    var b0 = encryption_0.c0_ntt_values.copy()
    var a0 = encryption_0.c1_ntt_values.copy()
    var b1 = encryption_1.c0_ntt_values.copy()
    var a1 = encryption_1.c1_ntt_values.copy()

    var d0 = _poly_mul_ntt(b0, b1, n_power_of_2, p_power_of_3, q_modulus)
    var a0b1 = _poly_mul_ntt(a0, b1, n_power_of_2, p_power_of_3, q_modulus)
    var a1b0 = _poly_mul_ntt(a1, b0, n_power_of_2, p_power_of_3, q_modulus)
    var d1 = _poly_add_modulus(a0b1, a1b0, q_modulus)
    var d2 = _poly_mul_ntt(a0, a1, n_power_of_2, p_power_of_3, q_modulus)

    var s_ntt = key_pair.secret_key_ntt_values.copy()
    var s_square_ntt = _poly_mul_ntt(
        s_ntt, s_ntt, n_power_of_2, p_power_of_3, q_modulus
    )
    var d1_times_s = _poly_mul_ntt(
        d1, s_ntt, n_power_of_2, p_power_of_3, q_modulus
    )
    var d2_times_s_square = _poly_mul_ntt(
        d2, s_square_ntt, n_power_of_2, p_power_of_3, q_modulus
    )
    var three_component_decrypt = _poly_add_modulus(
        _poly_add_modulus(d0, d1_times_s, q_modulus),
        d2_times_s_square,
        q_modulus,
    )

    var decrypted_0 = rlwe_decrypt_ntt(
        encryption_0.c0_ntt_values,
        encryption_0.c1_ntt_values,
        key_pair,
    )
    var decrypted_1 = rlwe_decrypt_ntt(
        encryption_1.c0_ntt_values,
        encryption_1.c1_ntt_values,
        key_pair,
    )
    var expected_product = _poly_mul_ntt(
        decrypted_0,
        decrypted_1,
        n_power_of_2,
        p_power_of_3,
        q_modulus,
    )

    assert_true(
        _poly_equal(three_component_decrypt, expected_product),
        "d0 + d1*s + d2*s^2 should match product of decrypted polynomials",
    )

    print(
        "test_rlwe_ciphertext_hadamard_multiply_three_component_decrypt passed."
    )


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
