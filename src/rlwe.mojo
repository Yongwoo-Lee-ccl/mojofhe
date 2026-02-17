from collections import List
from src.polynomial import Polynomial


fn _normalize_modulus(value: Int64, modulus_value: UInt32) -> UInt32:
    var modulus_int64 = Int64(modulus_value)
    var reduced_value = value % modulus_int64
    if reduced_value < 0:
        reduced_value += modulus_int64
    return UInt32(reduced_value)


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


fn _poly_subtract_modulus(
    left_values: List[UInt32],
    right_values: List[UInt32],
    modulus_value: UInt32,
) -> List[UInt32]:
    var result_values = List[UInt32](length=len(left_values), fill=0)
    for coefficient_index in range(len(left_values)):
        var difference_value = Int64(left_values[coefficient_index]) - Int64(
            right_values[coefficient_index]
        )
        result_values[coefficient_index] = _normalize_modulus(
            difference_value, modulus_value
        )
    return result_values^


fn _negacyclic_multiply_modulus(
    left_values: List[UInt32],
    right_values: List[UInt32],
    modulus_value: UInt32,
) -> List[UInt32]:
    var polynomial_degree = len(left_values)
    var result_values = List[UInt32](length=polynomial_degree, fill=0)
    for left_index in range(polynomial_degree):
        for right_index in range(polynomial_degree):
            var target_index = left_index + right_index
            var signed_product = Int64(left_values[left_index]) * Int64(
                right_values[right_index]
            )
            if target_index >= polynomial_degree:
                target_index -= polynomial_degree
                signed_product = -signed_product
            var accumulated_value = (
                Int64(result_values[target_index]) + signed_product
            )
            result_values[target_index] = _normalize_modulus(
                accumulated_value, modulus_value
            )
    return result_values^


fn _sample_via_existing_polynomial(
    polynomial_degree: Int, modulus_value: UInt32
) -> List[UInt32]:
    if polynomial_degree != 16:
        return List[UInt32](length=polynomial_degree, fill=0)
    var sampled_poly = Polynomial(4, 3, modulus_value)
    return sampled_poly.coefficient_values.copy()


fn _sample_uniform_polynomial(
    polynomial_degree: Int, modulus_value: UInt32
) -> List[UInt32]:
    var sampled_poly = Polynomial(4, 3, modulus_value)
    if polynomial_degree != sampled_poly.total_length:
        return _sample_via_existing_polynomial(polynomial_degree, modulus_value)
    sampled_poly.sample_uniform()
    return sampled_poly.coefficient_values.copy()


fn _sample_ternary_polynomial(
    polynomial_degree: Int, modulus_value: UInt32
) -> List[UInt32]:
    var sampled_poly = Polynomial(4, 3, modulus_value)
    if polynomial_degree != sampled_poly.total_length:
        return _sample_via_existing_polynomial(polynomial_degree, modulus_value)
    sampled_poly.sample_ternary()
    return sampled_poly.coefficient_values.copy()


fn _sample_gaussian_polynomial(
    polynomial_degree: Int, modulus_value: UInt32
) -> List[UInt32]:
    var sampled_poly = Polynomial(4, 3, modulus_value)
    if polynomial_degree != sampled_poly.total_length:
        return _sample_via_existing_polynomial(polynomial_degree, modulus_value)
    sampled_poly.sample_gaussian()
    return sampled_poly.coefficient_values.copy()


fn _encode_plaintext_binary(
    plaintext_bits: List[UInt32], modulus_value: UInt32
) -> List[UInt32]:
    var encoded_values = List[UInt32](length=len(plaintext_bits), fill=0)
    var scaling_value = modulus_value // 2
    for coefficient_index in range(len(plaintext_bits)):
        if plaintext_bits[coefficient_index] % 2 == 1:
            encoded_values[coefficient_index] = scaling_value
    return encoded_values^


fn _decode_plaintext_binary(
    encoded_values: List[UInt32], modulus_value: UInt32
) -> List[UInt32]:
    var decoded_values = List[UInt32](length=len(encoded_values), fill=0)
    var scaling_value = Int64(modulus_value // 2)
    var modulus_int64 = Int64(modulus_value)
    for coefficient_index in range(len(encoded_values)):
        var value = Int64(encoded_values[coefficient_index])
        var distance_to_zero = value
        if modulus_int64 - value < distance_to_zero:
            distance_to_zero = modulus_int64 - value

        var centered_delta_distance = value - scaling_value
        if centered_delta_distance < 0:
            centered_delta_distance = -centered_delta_distance
        if modulus_int64 - centered_delta_distance < centered_delta_distance:
            centered_delta_distance = modulus_int64 - centered_delta_distance

        if centered_delta_distance < distance_to_zero:
            decoded_values[coefficient_index] = 1
    return decoded_values^


struct RLWEKeyPair(Movable):
    var secret_key_values: List[UInt32]
    var public_a_values: List[UInt32]
    var public_b_values: List[UInt32]

    fn __init__(
        out self,
        secret_key_values: List[UInt32],
        public_a_values: List[UInt32],
        public_b_values: List[UInt32],
    ):
        self.secret_key_values = secret_key_values.copy()
        self.public_a_values = public_a_values.copy()
        self.public_b_values = public_b_values.copy()


struct RLWECiphertext(Movable):
    var c0_values: List[UInt32]
    var c1_values: List[UInt32]

    fn __init__(out self, c0_values: List[UInt32], c1_values: List[UInt32]):
        self.c0_values = c0_values.copy()
        self.c1_values = c1_values.copy()


fn rlwe_keygen(polynomial_degree: Int, modulus_value: UInt32) -> RLWEKeyPair:
    var secret_values = _sample_ternary_polynomial(
        polynomial_degree, modulus_value
    )
    var public_a_values = _sample_uniform_polynomial(
        polynomial_degree, modulus_value
    )
    var error_values = _sample_gaussian_polynomial(
        polynomial_degree, modulus_value
    )

    var a_times_secret = _negacyclic_multiply_modulus(
        public_a_values, secret_values, modulus_value
    )
    var error_minus_a_times_secret = _poly_subtract_modulus(
        error_values, a_times_secret, modulus_value
    )

    return RLWEKeyPair(
        secret_values, public_a_values, error_minus_a_times_secret
    )


fn rlwe_encrypt_binary(
    plaintext_bits: List[UInt32],
    public_a_values: List[UInt32],
    public_b_values: List[UInt32],
    modulus_value: UInt32,
) -> RLWECiphertext:
    var polynomial_degree = len(plaintext_bits)
    var ephemeral_values = _sample_ternary_polynomial(
        polynomial_degree, modulus_value
    )
    var error_1_values = _sample_gaussian_polynomial(
        polynomial_degree, modulus_value
    )
    var error_2_values = _sample_gaussian_polynomial(
        polynomial_degree, modulus_value
    )

    var public_b_times_ephemeral = _negacyclic_multiply_modulus(
        public_b_values, ephemeral_values, modulus_value
    )
    var public_a_times_ephemeral = _negacyclic_multiply_modulus(
        public_a_values, ephemeral_values, modulus_value
    )
    var encoded_plaintext = _encode_plaintext_binary(
        plaintext_bits, modulus_value
    )
    var c0_with_error = _poly_add_modulus(
        public_b_times_ephemeral, error_1_values, modulus_value
    )
    var c0_values = _poly_add_modulus(
        c0_with_error, encoded_plaintext, modulus_value
    )
    var c1_values = _poly_add_modulus(
        public_a_times_ephemeral, error_2_values, modulus_value
    )

    return RLWECiphertext(c0_values, c1_values)


fn rlwe_decrypt_binary(
    ciphertext_value: RLWECiphertext,
    secret_key_values: List[UInt32],
    modulus_value: UInt32,
) -> List[UInt32]:
    var c1_times_secret = _negacyclic_multiply_modulus(
        ciphertext_value.c1_values, secret_key_values, modulus_value
    )
    var decrypted_scaled_values = _poly_add_modulus(
        ciphertext_value.c0_values, c1_times_secret, modulus_value
    )
    return _decode_plaintext_binary(decrypted_scaled_values, modulus_value)
