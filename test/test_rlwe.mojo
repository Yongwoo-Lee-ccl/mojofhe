from collections import List
from testing import TestSuite, assert_true
from random import random_si64
from src.rlwe import rlwe_keygen, rlwe_encrypt_binary, rlwe_decrypt_binary


fn _centered_coefficient(value: UInt32, modulus_value: UInt32) -> Int:
    var centered_value = Int(value)
    var half_modulus = Int(modulus_value) // 2
    if centered_value > half_modulus:
        centered_value -= Int(modulus_value)
    return centered_value


fn test_rlwe_secret_key_is_ternary() raises:
    print("Running test_rlwe_secret_key_is_ternary...")
    var polynomial_degree = 16
    var modulus_value: UInt32 = 40961
    var key_pair = rlwe_keygen(polynomial_degree, modulus_value)

    for coefficient_index in range(polynomial_degree):
        var centered_value = _centered_coefficient(
            key_pair.secret_key_values[coefficient_index], modulus_value
        )
        assert_true(
            centered_value >= -1 and centered_value <= 1,
            "secret key coefficient must be ternary (-1, 0, 1)",
        )

    print("test_rlwe_secret_key_is_ternary passed.")


fn test_rlwe_encrypt_decrypt_binary_message() raises:
    print("Running test_rlwe_encrypt_decrypt_binary_message...")
    var polynomial_degree = 16
    var modulus_value: UInt32 = 40961
    var key_pair = rlwe_keygen(polynomial_degree, modulus_value)

    var message_values = List[UInt32](length=polynomial_degree, fill=0)
    var message_reference = List[UInt32](length=polynomial_degree, fill=0)
    for coefficient_index in range(polynomial_degree):
        var sampled_bit = UInt32(Int(random_si64(0, 1)))
        message_values[coefficient_index] = sampled_bit
        message_reference[coefficient_index] = sampled_bit

    var ciphertext_value = rlwe_encrypt_binary(
        message_values,
        key_pair.public_a_values,
        key_pair.public_b_values,
        modulus_value,
    )
    var decrypted_values = rlwe_decrypt_binary(
        ciphertext_value,
        key_pair.secret_key_values,
        modulus_value,
    )

    for coefficient_index in range(polynomial_degree):
        assert_true(
            decrypted_values[coefficient_index]
            == message_reference[coefficient_index],
            "decrypted coefficient should match original bit",
        )

    print("test_rlwe_encrypt_decrypt_binary_message passed.")


fn test_rlwe_encrypt_decrypt_multiple_trials() raises:
    print("Running test_rlwe_encrypt_decrypt_multiple_trials...")
    var polynomial_degree = 16
    var modulus_value: UInt32 = 40961
    var number_of_trials = 10

    for trial_index in range(number_of_trials):
        _ = trial_index
        var key_pair = rlwe_keygen(polynomial_degree, modulus_value)
        var message_values = List[UInt32](length=polynomial_degree, fill=0)
        var message_reference = List[UInt32](length=polynomial_degree, fill=0)
        for coefficient_index in range(polynomial_degree):
            var sampled_bit = UInt32(Int(random_si64(0, 1)))
            message_values[coefficient_index] = sampled_bit
            message_reference[coefficient_index] = sampled_bit

        var ciphertext_value = rlwe_encrypt_binary(
            message_values,
            key_pair.public_a_values,
            key_pair.public_b_values,
            modulus_value,
        )
        var decrypted_values = rlwe_decrypt_binary(
            ciphertext_value,
            key_pair.secret_key_values,
            modulus_value,
        )
        for coefficient_index in range(polynomial_degree):
            assert_true(
                decrypted_values[coefficient_index]
                == message_reference[coefficient_index],
                "all trials must decrypt correctly",
            )

    print("test_rlwe_encrypt_decrypt_multiple_trials passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
