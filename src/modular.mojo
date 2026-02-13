from src.arithmetic import is_prime, power_modular


fn mod_pow(
    base_value: UInt64, exponent_value: UInt64, modulus_value: UInt64
) -> UInt64:
    return power_modular(base_value, exponent_value, modulus_value)


fn compute_barrett_ratio(modulus_value: UInt32) -> UInt64:
    # k = 31 for q < 2^30. Ratio uses floor(2^(2k) / q).
    return (UInt64(1) << 62) // UInt64(modulus_value)


fn multiply_mod_barrett(
    left_value: UInt32,
    right_value: UInt32,
    modulus_value: UInt32,
    barrett_ratio: UInt64,
) -> UInt32:
    # Barrett reduction specialized for int32-lane friendly moduli.
    var wide_product = UInt64(left_value) * UInt64(right_value)
    var q1_estimate = wide_product >> 30
    var q2_estimate = q1_estimate * barrett_ratio
    var q3_estimate = q2_estimate >> 32
    var reduced_value = Int64(wide_product) - Int64(
        q3_estimate * UInt64(modulus_value)
    )

    while reduced_value >= Int64(modulus_value):
        reduced_value -= Int64(modulus_value)
    while reduced_value < 0:
        reduced_value += Int64(modulus_value)

    return UInt32(reduced_value)


fn compute_montgomery_neg_inv(modulus_value: UInt32) -> UInt64:
    # Compute -q^{-1} mod 2^32 using Newton iteration in the 32-bit ring.
    var lower_word_mask = (UInt64(1) << 32) - 1
    var inverse_estimate: UInt64 = 1
    var modulus_word = UInt64(modulus_value) & lower_word_mask
    for iteration_index in range(6):
        _ = iteration_index
        inverse_estimate = (
            inverse_estimate * (2 - modulus_word * inverse_estimate)
        ) & lower_word_mask
    return (-inverse_estimate) & lower_word_mask


fn compute_montgomery_r2(modulus_value: UInt32) -> UInt32:
    var r_mod_q = UInt32((UInt64(1) << 32) % UInt64(modulus_value))
    var barrett_ratio = compute_barrett_ratio(modulus_value)
    return multiply_mod_barrett(r_mod_q, r_mod_q, modulus_value, barrett_ratio)


fn _montgomery_reduce(
    wide_value: UInt64, modulus_value: UInt32, montgomery_neg_inv: UInt64
) -> UInt32:
    var lower_word_mask = (UInt64(1) << 32) - 1
    var correction_factor = (
        (wide_value & lower_word_mask) * montgomery_neg_inv
    ) & lower_word_mask
    var reduced_value = (
        wide_value + correction_factor * UInt64(modulus_value)
    ) >> 32
    if reduced_value >= UInt64(modulus_value):
        reduced_value -= UInt64(modulus_value)
    return UInt32(reduced_value)


fn montgomery_multiply_raw(
    left_value: UInt32,
    right_value: UInt32,
    modulus_value: UInt32,
    montgomery_neg_inv: UInt64,
) -> UInt32:
    return _montgomery_reduce(
        UInt64(left_value) * UInt64(right_value),
        modulus_value,
        montgomery_neg_inv,
    )


fn to_montgomery_domain(
    value: UInt32,
    modulus_value: UInt32,
    montgomery_neg_inv: UInt64,
    montgomery_r2: UInt32,
) -> UInt32:
    return montgomery_multiply_raw(
        value, montgomery_r2, modulus_value, montgomery_neg_inv
    )


fn multiply_mod_montgomery_with_rhs_mont(
    left_value: UInt32,
    right_value_montgomery: UInt32,
    modulus_value: UInt32,
    montgomery_neg_inv: UInt64,
) -> UInt32:
    # If right operand is in Montgomery domain, output keeps left domain:
    # normal x mont -> normal, mont x mont -> mont.
    return montgomery_multiply_raw(
        left_value,
        right_value_montgomery,
        modulus_value,
        montgomery_neg_inv,
    )


fn multiply_mod_montgomery(
    left_value: UInt32,
    right_value: UInt32,
    modulus_value: UInt32,
    montgomery_neg_inv: UInt64,
    montgomery_r2: UInt32,
) -> UInt32:
    # Convert both operands to Montgomery domain and reduce back to normal.
    var left_montgomery = to_montgomery_domain(
        left_value, modulus_value, montgomery_neg_inv, montgomery_r2
    )
    var right_montgomery = to_montgomery_domain(
        right_value, modulus_value, montgomery_neg_inv, montgomery_r2
    )
    var product_montgomery = montgomery_multiply_raw(
        left_montgomery,
        right_montgomery,
        modulus_value,
        montgomery_neg_inv,
    )
    return _montgomery_reduce(
        UInt64(product_montgomery), modulus_value, montgomery_neg_inv
    )


fn find_suitable_q(n_dim: Int, p_dim: Int, target_bit_length: Int) -> UInt64:
    """
    Finds a prime q such that q = 1 mod (4 * n_dim * p_dim) and q has
    roughly target_bit_length bits.
    """
    var target_magnitude = UInt64(1) << target_bit_length
    var cyclic_order = UInt64(4 * n_dim * p_dim)

    var search_start_multiplier = target_magnitude // cyclic_order
    var search_offset = UInt64(0)
    while True:
        var candidate_up = (
            search_start_multiplier + search_offset
        ) * cyclic_order + 1
        if candidate_up > 0 and is_prime(candidate_up):
            return candidate_up

        if search_offset > 0:
            var candidate_down = (
                search_start_multiplier - search_offset
            ) * cyclic_order + 1
            if candidate_down > 0 and is_prime(candidate_down):
                return candidate_down

        search_offset += 1


fn find_primitive_root(modulus_value: UInt64, order_value: UInt64) -> UInt64:
    """
    Finds a primitive n-th root of unity modulo q.
    Returns 0 if no primitive root is found.
    """
    if (modulus_value - 1) % order_value != 0:
        return 0

    var search_candidate: UInt64 = 2
    while search_candidate < modulus_value:
        var exponent_for_root = (modulus_value - 1) // order_value
        var root_candidate = mod_pow(
            search_candidate, exponent_for_root, modulus_value
        )

        if root_candidate != 1:
            var is_actually_primitive = True

            # Check prime factors of order_value to ensure it's primitive
            # We assume order_value consists of factors 2, 3, 5 based
            # on problem constraints
            if order_value % 2 == 0:
                if (
                    mod_pow(root_candidate, order_value // 2, modulus_value)
                    == 1
                ):
                    is_actually_primitive = False

            if is_actually_primitive and order_value % 3 == 0:
                if (
                    mod_pow(root_candidate, order_value // 3, modulus_value)
                    == 1
                ):
                    is_actually_primitive = False

            if is_actually_primitive and order_value % 5 == 0:
                if (
                    mod_pow(root_candidate, order_value // 5, modulus_value)
                    == 1
                ):
                    is_actually_primitive = False

            if is_actually_primitive:
                return root_candidate

        search_candidate += 1

    return 0
