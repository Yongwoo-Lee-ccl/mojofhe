from collections import List
from src.modular import compute_barrett_ratio, multiply_mod_barrett
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


fn _make_polynomial_from_coefficients(
    coefficient_values: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> Polynomial:
    var polynomial_value = Polynomial(n_power_of_2, p_power_of_3, modulus_value)
    for coefficient_index in range(polynomial_value.total_length):
        polynomial_value.coefficient_values[
            coefficient_index
        ] = coefficient_values[coefficient_index]
    return polynomial_value^


fn _sample_existing_uniform(
    n_power_of_2: Int, p_power_of_3: Int, modulus_value: UInt32
) -> List[UInt32]:
    var sampled_polynomial = Polynomial(
        n_power_of_2, p_power_of_3, modulus_value
    )
    sampled_polynomial.sample_uniform()
    return sampled_polynomial.coefficient_values.copy()


fn _sample_existing_ternary(
    n_power_of_2: Int, p_power_of_3: Int, modulus_value: UInt32
) -> List[UInt32]:
    var sampled_polynomial = Polynomial(
        n_power_of_2, p_power_of_3, modulus_value
    )
    sampled_polynomial.sample_ternary()
    return sampled_polynomial.coefficient_values.copy()


fn _sample_existing_gaussian(
    n_power_of_2: Int, p_power_of_3: Int, modulus_value: UInt32
) -> List[UInt32]:
    var sampled_polynomial = Polynomial(
        n_power_of_2, p_power_of_3, modulus_value
    )
    sampled_polynomial.sample_gaussian()
    return sampled_polynomial.coefficient_values.copy()


fn _to_ntt_with_existing_transform(
    coefficient_values: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    var polynomial_value = _make_polynomial_from_coefficients(
        coefficient_values, n_power_of_2, p_power_of_3, modulus_value
    )
    polynomial_value.transform_to_full_ntt()
    return polynomial_value.coefficient_values.copy()


fn _multiply_with_existing_ntt(
    left_coefficient_values: List[UInt32],
    right_coefficient_values: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    # Uses existing Polynomial.multiply(), which internally applies the
    # project's multivariate NTT pipeline for this quotient ring.
    var left_poly = _make_polynomial_from_coefficients(
        left_coefficient_values, n_power_of_2, p_power_of_3, modulus_value
    )
    var right_poly = _make_polynomial_from_coefficients(
        right_coefficient_values, n_power_of_2, p_power_of_3, modulus_value
    )
    var product_poly = left_poly.multiply(right_poly)
    return product_poly.coefficient_values.copy()


fn _pointwise_multiply_ntt_modulus(
    left_ntt_values: List[UInt32],
    right_ntt_values: List[UInt32],
    modulus_value: UInt32,
) -> List[UInt32]:
    var result_values = List[UInt32](length=len(left_ntt_values), fill=0)
    var barrett_ratio = compute_barrett_ratio(modulus_value)
    for coefficient_index in range(len(left_ntt_values)):
        result_values[coefficient_index] = multiply_mod_barrett(
            left_ntt_values[coefficient_index],
            right_ntt_values[coefficient_index],
            modulus_value,
            barrett_ratio,
        )
    return result_values^


fn _encode_binary_plaintext_coefficients(
    plaintext_bits: List[UInt32], modulus_value: UInt32
) -> List[UInt32]:
    var encoded_values = List[UInt32](length=len(plaintext_bits), fill=0)
    var scaling_value = modulus_value // 2
    for coefficient_index in range(len(plaintext_bits)):
        if plaintext_bits[coefficient_index] % 2 == 1:
            encoded_values[coefficient_index] = scaling_value
    return encoded_values^


struct RLWEKeyPair(Movable):
    var n_power_of_2: Int
    var p_power_of_3: Int
    var q_modulus: UInt32
    var secret_key_coeff_values: List[UInt32]
    var secret_key_ntt_values: List[UInt32]
    var public_a_ntt_values: List[UInt32]
    var public_b_ntt_values: List[UInt32]

    fn __init__(
        out self,
        n_power_of_2: Int,
        p_power_of_3: Int,
        q_modulus: UInt32,
        secret_key_coeff_values: List[UInt32],
        secret_key_ntt_values: List[UInt32],
        public_a_ntt_values: List[UInt32],
        public_b_ntt_values: List[UInt32],
    ):
        self.n_power_of_2 = n_power_of_2
        self.p_power_of_3 = p_power_of_3
        self.q_modulus = q_modulus
        self.secret_key_coeff_values = secret_key_coeff_values.copy()
        self.secret_key_ntt_values = secret_key_ntt_values.copy()
        self.public_a_ntt_values = public_a_ntt_values.copy()
        self.public_b_ntt_values = public_b_ntt_values.copy()


struct RLWEEncryptionResult(Movable):
    var c0_ntt_values: List[UInt32]
    var c1_ntt_values: List[UInt32]
    var ephemeral_u_ntt_values: List[UInt32]
    var error_1_ntt_values: List[UInt32]
    var error_2_ntt_values: List[UInt32]

    fn __init__(
        out self,
        c0_ntt_values: List[UInt32],
        c1_ntt_values: List[UInt32],
        ephemeral_u_ntt_values: List[UInt32],
        error_1_ntt_values: List[UInt32],
        error_2_ntt_values: List[UInt32],
    ):
        self.c0_ntt_values = c0_ntt_values.copy()
        self.c1_ntt_values = c1_ntt_values.copy()
        self.ephemeral_u_ntt_values = ephemeral_u_ntt_values.copy()
        self.error_1_ntt_values = error_1_ntt_values.copy()
        self.error_2_ntt_values = error_2_ntt_values.copy()


fn rlwe_keygen(
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> RLWEKeyPair:
    var secret_key_coeff_values = _sample_existing_ternary(
        n_power_of_2, p_power_of_3, modulus_value
    )
    var public_a_coeff_values = _sample_existing_uniform(
        n_power_of_2, p_power_of_3, modulus_value
    )
    var error_coeff_values = _sample_existing_gaussian(
        n_power_of_2, p_power_of_3, modulus_value
    )

    var secret_key_ntt_values = _to_ntt_with_existing_transform(
        secret_key_coeff_values,
        n_power_of_2,
        p_power_of_3,
        modulus_value,
    )
    var public_a_ntt_values = _to_ntt_with_existing_transform(
        public_a_coeff_values,
        n_power_of_2,
        p_power_of_3,
        modulus_value,
    )
    var error_ntt_values = _to_ntt_with_existing_transform(
        error_coeff_values,
        n_power_of_2,
        p_power_of_3,
        modulus_value,
    )

    var a_times_secret_ntt_values = _pointwise_multiply_ntt_modulus(
        public_a_ntt_values, secret_key_ntt_values, modulus_value
    )
    var public_b_ntt_values = _poly_subtract_modulus(
        error_ntt_values, a_times_secret_ntt_values, modulus_value
    )

    return RLWEKeyPair(
        n_power_of_2,
        p_power_of_3,
        modulus_value,
        secret_key_coeff_values,
        secret_key_ntt_values,
        public_a_ntt_values,
        public_b_ntt_values,
    )


fn rlwe_encrypt_binary(
    plaintext_bits: List[UInt32], key_pair: RLWEKeyPair
) -> RLWEEncryptionResult:
    var encoded_plaintext_coeff_values = _encode_binary_plaintext_coefficients(
        plaintext_bits, key_pair.q_modulus
    )
    var encoded_plaintext_ntt_values = _to_ntt_with_existing_transform(
        encoded_plaintext_coeff_values,
        key_pair.n_power_of_2,
        key_pair.p_power_of_3,
        key_pair.q_modulus,
    )

    var ephemeral_u_coeff_values = _sample_existing_ternary(
        key_pair.n_power_of_2, key_pair.p_power_of_3, key_pair.q_modulus
    )
    var error_1_coeff_values = _sample_existing_gaussian(
        key_pair.n_power_of_2, key_pair.p_power_of_3, key_pair.q_modulus
    )
    var error_2_coeff_values = _sample_existing_gaussian(
        key_pair.n_power_of_2, key_pair.p_power_of_3, key_pair.q_modulus
    )

    var ephemeral_u_ntt_values = _to_ntt_with_existing_transform(
        ephemeral_u_coeff_values,
        key_pair.n_power_of_2,
        key_pair.p_power_of_3,
        key_pair.q_modulus,
    )
    var error_1_ntt_values = _to_ntt_with_existing_transform(
        error_1_coeff_values,
        key_pair.n_power_of_2,
        key_pair.p_power_of_3,
        key_pair.q_modulus,
    )
    var error_2_ntt_values = _to_ntt_with_existing_transform(
        error_2_coeff_values,
        key_pair.n_power_of_2,
        key_pair.p_power_of_3,
        key_pair.q_modulus,
    )

    var public_b_times_u_ntt_values = _pointwise_multiply_ntt_modulus(
        key_pair.public_b_ntt_values,
        ephemeral_u_ntt_values,
        key_pair.q_modulus,
    )
    var public_a_times_u_ntt_values = _pointwise_multiply_ntt_modulus(
        key_pair.public_a_ntt_values,
        ephemeral_u_ntt_values,
        key_pair.q_modulus,
    )

    var c0_with_error_values = _poly_add_modulus(
        public_b_times_u_ntt_values, error_1_ntt_values, key_pair.q_modulus
    )
    var c0_ntt_values = _poly_add_modulus(
        c0_with_error_values, encoded_plaintext_ntt_values, key_pair.q_modulus
    )
    var c1_ntt_values = _poly_add_modulus(
        public_a_times_u_ntt_values, error_2_ntt_values, key_pair.q_modulus
    )

    return RLWEEncryptionResult(
        c0_ntt_values,
        c1_ntt_values,
        ephemeral_u_ntt_values,
        error_1_ntt_values,
        error_2_ntt_values,
    )


fn rlwe_decrypt_ntt(
    c0_ntt_values: List[UInt32],
    c1_ntt_values: List[UInt32],
    key_pair: RLWEKeyPair,
) -> List[UInt32]:
    var c1_times_secret_ntt_values = _pointwise_multiply_ntt_modulus(
        c1_ntt_values,
        key_pair.secret_key_ntt_values,
        key_pair.q_modulus,
    )
    return _poly_add_modulus(
        c0_ntt_values,
        c1_times_secret_ntt_values,
        key_pair.q_modulus,
    )


fn rlwe_ring_multiply_ntt(
    left_ntt_values: List[UInt32],
    right_ntt_values: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    _ = n_power_of_2
    _ = p_power_of_3
    return _pointwise_multiply_ntt_modulus(
        left_ntt_values, right_ntt_values, modulus_value
    )


fn encode_binary_plaintext_ntt(
    plaintext_bits: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    return _to_ntt_with_existing_transform(
        _encode_binary_plaintext_coefficients(plaintext_bits, modulus_value),
        n_power_of_2,
        p_power_of_3,
        modulus_value,
    )
