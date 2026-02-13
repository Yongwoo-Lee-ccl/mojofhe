from collections import List
from src.modular import (
    mod_pow,
    find_primitive_root,
    compute_barrett_ratio,
    multiply_mod_barrett,
    compute_montgomery_neg_inv,
    compute_montgomery_r2,
    to_montgomery_domain,
    multiply_mod_montgomery_with_rhs_mont,
)


fn _multiply_modulo_int32(
    left_value: Int,
    right_value: Int,
    modulus_int: Int,
    barrett_ratio: Int64,
    use_montgomery: Bool,
    montgomery_neg_inv: Int64,
    montgomery_r2: Int,
) -> Int:
    _ = montgomery_r2
    if barrett_ratio < 0:
        return Int(
            (Int64(left_value) * Int64(right_value)) % Int64(modulus_int)
        )
    if use_montgomery:
        return multiply_mod_montgomery_with_rhs_mont(
            left_value,
            right_value,
            modulus_int,
            montgomery_neg_inv,
        )
    return multiply_mod_barrett(
        left_value, right_value, modulus_int, barrett_ratio
    )


fn _add_modulo(left_value: Int, right_value: Int, modulus_int: Int) -> Int:
    var sum_value = left_value + right_value
    if sum_value >= modulus_int:
        sum_value -= modulus_int
    return sum_value


fn _sub_modulo(left_value: Int, right_value: Int, modulus_int: Int) -> Int:
    var diff_value = left_value - right_value
    if diff_value < 0:
        diff_value += modulus_int
    return diff_value


fn _add_three_modulo(
    first_value: Int, second_value: Int, third_value: Int, modulus_int: Int
) -> Int:
    return _add_modulo(
        _add_modulo(first_value, second_value, modulus_int),
        third_value,
        modulus_int,
    )


fn apply_i_axis_transform(
    mut coefficient_values: List[Int32],
    total_length: Int,
    root_imaginary_unit: Int32,
    modulus_value: Int32,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
    if use_naive_modulo:
        barrett_ratio = -1
    var montgomery_neg_inv: Int64 = 0
    var montgomery_r2: Int = 0
    var root_multiplier = Int(root_imaginary_unit)
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_int)
        montgomery_r2 = compute_montgomery_r2(modulus_int)
        root_multiplier = to_montgomery_domain(
            root_multiplier,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
    var half_stride = total_length // 2
    for offset_index in range(half_stride):
        var real_part = Int(coefficient_values[offset_index])
        var imag_part = Int(coefficient_values[offset_index + half_stride])
        var weighted_imag = _multiply_modulo_int32(
            imag_part,
            root_multiplier,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )
        coefficient_values[offset_index] = Int32(
            (real_part + weighted_imag) % modulus_int
        )
        coefficient_values[offset_index + half_stride] = Int32(
            (real_part - weighted_imag + modulus_int) % modulus_int
        )


fn _mod_inverse(value: Int, modulus_value: Int) -> Int:
    var normalized_value = value % modulus_value
    if normalized_value < 0:
        normalized_value += modulus_value

    var previous_remainder = modulus_value
    var current_remainder = normalized_value
    var previous_coefficient = 0
    var current_coefficient = 1
    while current_remainder != 0:
        var quotient_value = previous_remainder // current_remainder
        var next_remainder = (
            previous_remainder - quotient_value * current_remainder
        )
        previous_remainder = current_remainder
        current_remainder = next_remainder

        var next_coefficient = (
            previous_coefficient - quotient_value * current_coefficient
        )
        previous_coefficient = current_coefficient
        current_coefficient = next_coefficient

    if previous_remainder != 1:
        return -1

    var normalized_coefficient = previous_coefficient % modulus_value
    if normalized_coefficient < 0:
        normalized_coefficient += modulus_value
    return normalized_coefficient


fn apply_inverse_i_axis_transform(
    mut coefficient_values: List[Int32],
    total_length: Int,
    root_imaginary_unit: Int32,
    modulus_value: Int32,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
    if use_naive_modulo:
        barrett_ratio = -1
    var montgomery_neg_inv: Int64 = 0
    var montgomery_r2: Int = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_int)
        montgomery_r2 = compute_montgomery_r2(modulus_int)

    var inverse_two = _mod_inverse(2, modulus_int)
    var inverse_root = _mod_inverse(Int(root_imaginary_unit), modulus_int)
    if inverse_two < 0 or inverse_root < 0:
        return

    var inverse_two_multiplier = inverse_two
    var inverse_root_multiplier = inverse_root
    if use_montgomery:
        inverse_two_multiplier = to_montgomery_domain(
            inverse_two,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        inverse_root_multiplier = to_montgomery_domain(
            inverse_root,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )

    var half_stride = total_length // 2
    for offset_index in range(half_stride):
        var first_transformed = Int(coefficient_values[offset_index])
        var second_transformed = Int(
            coefficient_values[offset_index + half_stride]
        )
        var summed_value = (
            first_transformed + second_transformed
        ) % modulus_int
        var difference_value = (
            first_transformed - second_transformed + modulus_int
        ) % modulus_int
        var real_part = _multiply_modulo_int32(
            summed_value,
            inverse_two_multiplier,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )
        var imag_scaled = _multiply_modulo_int32(
            difference_value,
            inverse_two_multiplier,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )
        var imag_part = _multiply_modulo_int32(
            imag_scaled,
            inverse_root_multiplier,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )
        coefficient_values[offset_index] = Int32(real_part)
        coefficient_values[offset_index + half_stride] = Int32(imag_part)


fn apply_radix2_dif_ntt(
    mut coefficient_values: List[Int32],
    transform_length: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
    if use_naive_modulo:
        barrett_ratio = -1
    var montgomery_neg_inv: Int64 = 0
    var montgomery_r2: Int = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_int)
        montgomery_r2 = compute_montgomery_r2(modulus_int)
    var root_int = Int(root_of_unity)
    var current_step_size = transform_length
    var stage_twiddle = root_int
    var montgomery_one = 1
    if use_montgomery:
        montgomery_one = to_montgomery_domain(
            1,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        stage_twiddle = to_montgomery_domain(
            root_int,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
    while current_step_size > 1:
        var half_step = current_step_size // 2
        for group_start in range(0, transform_length, current_step_size):
            var current_twiddle = montgomery_one
            for butterfly_index in range(half_step):
                var upper_index_base = (
                    base_offset + (group_start + butterfly_index) * block_size
                )
                var lower_index_base = (
                    base_offset
                    + (group_start + butterfly_index + half_step) * block_size
                )
                # TODO: SIMD this contiguous butterfly block update loop.
                for inner_offset in range(block_size):
                    var upper_val = Int(
                        coefficient_values[upper_index_base + inner_offset]
                    )
                    var lower_val = Int(
                        coefficient_values[lower_index_base + inner_offset]
                    )
                    var sum_val = _add_modulo(upper_val, lower_val, modulus_int)
                    var diff_val = _sub_modulo(
                        upper_val, lower_val, modulus_int
                    )
                    coefficient_values[upper_index_base + inner_offset] = Int32(
                        sum_val
                    )
                    coefficient_values[lower_index_base + inner_offset] = Int32(
                        _multiply_modulo_int32(
                            diff_val,
                            current_twiddle,
                            modulus_int,
                            barrett_ratio,
                            use_montgomery,
                            montgomery_neg_inv,
                            montgomery_r2,
                        )
                    )
                current_twiddle = _multiply_modulo_int32(
                    current_twiddle,
                    stage_twiddle,
                    modulus_int,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
                    montgomery_r2,
                )
        stage_twiddle = _multiply_modulo_int32(
            stage_twiddle,
            stage_twiddle,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )
        current_step_size = half_step


fn apply_radix2_dif_intt(
    mut coefficient_values: List[Int32],
    transform_length: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
    if use_naive_modulo:
        barrett_ratio = -1
    var montgomery_neg_inv: Int64 = 0
    var montgomery_r2: Int = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_int)
        montgomery_r2 = compute_montgomery_r2(modulus_int)

    var inverse_two = _mod_inverse(2, modulus_int)
    var root_int_standard = Int(root_of_unity)
    var inverse_root_standard = _mod_inverse(root_int_standard, modulus_int)
    if inverse_two < 0 or inverse_root_standard < 0:
        return

    var inverse_two_multiplier = inverse_two
    var montgomery_one = 1
    if use_montgomery:
        inverse_two_multiplier = to_montgomery_domain(
            inverse_two,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        montgomery_one = to_montgomery_domain(
            1,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )

    var current_step_size = 2
    while current_step_size <= transform_length:
        var stage_exponent = transform_length // current_step_size
        var inverse_stage_twiddle_standard = mod_pow(
            inverse_root_standard, stage_exponent, modulus_int
        )
        var inverse_stage_twiddle = inverse_stage_twiddle_standard
        if use_montgomery:
            inverse_stage_twiddle = to_montgomery_domain(
                inverse_stage_twiddle_standard,
                modulus_int,
                montgomery_neg_inv,
                montgomery_r2,
            )

        var half_step = current_step_size // 2
        for group_start in range(0, transform_length, current_step_size):
            var current_inverse_twiddle = 1
            if use_montgomery:
                current_inverse_twiddle = montgomery_one
            for butterfly_index in range(half_step):
                var upper_index_base = (
                    base_offset + (group_start + butterfly_index) * block_size
                )
                var lower_index_base = (
                    base_offset
                    + (group_start + butterfly_index + half_step) * block_size
                )
                for inner_offset in range(block_size):
                    var summed_value = Int(
                        coefficient_values[upper_index_base + inner_offset]
                    )
                    var difference_twiddled = Int(
                        coefficient_values[lower_index_base + inner_offset]
                    )
                    var difference_value = _multiply_modulo_int32(
                        difference_twiddled,
                        current_inverse_twiddle,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var upper_recovered = _multiply_modulo_int32(
                        _add_modulo(
                            summed_value, difference_value, modulus_int
                        ),
                        inverse_two_multiplier,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var lower_recovered = _multiply_modulo_int32(
                        _sub_modulo(
                            summed_value, difference_value, modulus_int
                        ),
                        inverse_two_multiplier,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    coefficient_values[upper_index_base + inner_offset] = Int32(
                        upper_recovered
                    )
                    coefficient_values[lower_index_base + inner_offset] = Int32(
                        lower_recovered
                    )
                current_inverse_twiddle = _multiply_modulo_int32(
                    current_inverse_twiddle,
                    inverse_stage_twiddle,
                    modulus_int,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
                    montgomery_r2,
                )
        current_step_size *= 2


fn apply_radix3_dif_ntt(
    mut coefficient_values: List[Int32],
    transform_length: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
    if use_naive_modulo:
        barrett_ratio = -1
    var montgomery_neg_inv: Int64 = 0
    var montgomery_r2: Int = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_int)
        montgomery_r2 = compute_montgomery_r2(modulus_int)
    var root_int = Int(root_of_unity)
    var current_step_size = transform_length
    var root_order_3 = mod_pow(root_int, transform_length // 3, modulus_int)
    var montgomery_one = 1
    if use_montgomery:
        montgomery_one = to_montgomery_domain(
            1,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        root_int = to_montgomery_domain(
            root_int,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        root_order_3 = to_montgomery_domain(
            root_order_3,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
    var root_order_3_sq = _multiply_modulo_int32(
        root_order_3,
        root_order_3,
        modulus_int,
        barrett_ratio,
        use_montgomery,
        montgomery_neg_inv,
        montgomery_r2,
    )
    var stage_twiddle = root_int
    while current_step_size >= 3:
        var third_step = current_step_size // 3
        for group_start in range(0, transform_length, current_step_size):
            var twiddle_1 = montgomery_one
            for butterfly_index in range(third_step):
                var twiddle_2 = _multiply_modulo_int32(
                    twiddle_1,
                    twiddle_1,
                    modulus_int,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
                    montgomery_r2,
                )
                var first_index_base = (
                    base_offset + (group_start + butterfly_index) * block_size
                )
                var second_index_base = (
                    base_offset
                    + (group_start + butterfly_index + third_step) * block_size
                )
                var third_index_base = (
                    base_offset
                    + (group_start + butterfly_index + 2 * third_step)
                    * block_size
                )
                for inner_offset in range(block_size):
                    var first_value = Int(
                        coefficient_values[first_index_base + inner_offset]
                    )
                    var second_value = Int(
                        coefficient_values[second_index_base + inner_offset]
                    )
                    var third_value = Int(
                        coefficient_values[third_index_base + inner_offset]
                    )
                    var sum_all = _add_three_modulo(
                        first_value, second_value, third_value, modulus_int
                    )
                    var second_weighted_zeta = _multiply_modulo_int32(
                        second_value,
                        root_order_3,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var third_weighted_zeta_sq = _multiply_modulo_int32(
                        third_value,
                        root_order_3_sq,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var sum_zeta = _add_three_modulo(
                        first_value,
                        second_weighted_zeta,
                        third_weighted_zeta_sq,
                        modulus_int,
                    )
                    var second_weighted_zeta_sq = _multiply_modulo_int32(
                        second_value,
                        root_order_3_sq,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var third_weighted_zeta = _multiply_modulo_int32(
                        third_value,
                        root_order_3,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var sum_zeta_sq = _add_three_modulo(
                        first_value,
                        second_weighted_zeta_sq,
                        third_weighted_zeta,
                        modulus_int,
                    )
                    coefficient_values[first_index_base + inner_offset] = Int32(
                        sum_all
                    )
                    coefficient_values[
                        second_index_base + inner_offset
                    ] = Int32(
                        _multiply_modulo_int32(
                            sum_zeta,
                            twiddle_1,
                            modulus_int,
                            barrett_ratio,
                            use_montgomery,
                            montgomery_neg_inv,
                            montgomery_r2,
                        )
                    )
                    coefficient_values[third_index_base + inner_offset] = Int32(
                        _multiply_modulo_int32(
                            sum_zeta_sq,
                            twiddle_2,
                            modulus_int,
                            barrett_ratio,
                            use_montgomery,
                            montgomery_neg_inv,
                            montgomery_r2,
                        )
                    )
                twiddle_1 = _multiply_modulo_int32(
                    twiddle_1,
                    stage_twiddle,
                    modulus_int,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
                    montgomery_r2,
                )
        var stage_twiddle_squared = _multiply_modulo_int32(
            stage_twiddle,
            stage_twiddle,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )
        stage_twiddle = _multiply_modulo_int32(
            stage_twiddle_squared,
            stage_twiddle,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )
        current_step_size = third_step


fn apply_radix3_dif_intt(
    mut coefficient_values: List[Int32],
    transform_length: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
    if use_naive_modulo:
        barrett_ratio = -1
    var montgomery_neg_inv: Int64 = 0
    var montgomery_r2: Int = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_int)
        montgomery_r2 = compute_montgomery_r2(modulus_int)

    var inverse_three = _mod_inverse(3, modulus_int)
    var root_int_standard = Int(root_of_unity)
    var inverse_root_standard = _mod_inverse(root_int_standard, modulus_int)
    if inverse_three < 0 or inverse_root_standard < 0:
        return

    var root_order_3_standard = mod_pow(
        root_int_standard, transform_length // 3, modulus_int
    )
    var root_order_3_sq_standard = Int(
        (Int64(root_order_3_standard) * Int64(root_order_3_standard))
        % Int64(modulus_int)
    )

    var root_order_3 = root_order_3_standard
    var root_order_3_sq = root_order_3_sq_standard
    var inverse_three_multiplier = inverse_three
    var montgomery_one = 1
    if use_montgomery:
        root_order_3 = to_montgomery_domain(
            root_order_3_standard,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        root_order_3_sq = to_montgomery_domain(
            root_order_3_sq_standard,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        inverse_three_multiplier = to_montgomery_domain(
            inverse_three,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        montgomery_one = to_montgomery_domain(
            1,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )

    var current_step_size = 3
    while current_step_size <= transform_length:
        var stage_exponent = transform_length // current_step_size
        var inverse_stage_twiddle_standard = mod_pow(
            inverse_root_standard, stage_exponent, modulus_int
        )
        var inverse_stage_twiddle = inverse_stage_twiddle_standard
        if use_montgomery:
            inverse_stage_twiddle = to_montgomery_domain(
                inverse_stage_twiddle_standard,
                modulus_int,
                montgomery_neg_inv,
                montgomery_r2,
            )

        var third_step = current_step_size // 3
        for group_start in range(0, transform_length, current_step_size):
            var current_inverse_twiddle = 1
            if use_montgomery:
                current_inverse_twiddle = montgomery_one
            for butterfly_index in range(third_step):
                var current_inverse_twiddle_sq = _multiply_modulo_int32(
                    current_inverse_twiddle,
                    current_inverse_twiddle,
                    modulus_int,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
                    montgomery_r2,
                )
                var first_index_base = (
                    base_offset + (group_start + butterfly_index) * block_size
                )
                var second_index_base = (
                    base_offset
                    + (group_start + butterfly_index + third_step) * block_size
                )
                var third_index_base = (
                    base_offset
                    + (group_start + butterfly_index + 2 * third_step)
                    * block_size
                )
                for inner_offset in range(block_size):
                    var first_transformed = Int(
                        coefficient_values[first_index_base + inner_offset]
                    )
                    var second_transformed = Int(
                        coefficient_values[second_index_base + inner_offset]
                    )
                    var third_transformed = Int(
                        coefficient_values[third_index_base + inner_offset]
                    )
                    var second_unweighted = _multiply_modulo_int32(
                        second_transformed,
                        current_inverse_twiddle,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var third_unweighted = _multiply_modulo_int32(
                        third_transformed,
                        current_inverse_twiddle_sq,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )

                    var recovered_first = _multiply_modulo_int32(
                        _add_three_modulo(
                            first_transformed,
                            second_unweighted,
                            third_unweighted,
                            modulus_int,
                        ),
                        inverse_three_multiplier,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )

                    var second_with_zeta_sq = _multiply_modulo_int32(
                        second_unweighted,
                        root_order_3_sq,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var third_with_zeta = _multiply_modulo_int32(
                        third_unweighted,
                        root_order_3,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var recovered_second = _multiply_modulo_int32(
                        _add_three_modulo(
                            first_transformed,
                            second_with_zeta_sq,
                            third_with_zeta,
                            modulus_int,
                        ),
                        inverse_three_multiplier,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )

                    var second_with_zeta = _multiply_modulo_int32(
                        second_unweighted,
                        root_order_3,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var third_with_zeta_sq = _multiply_modulo_int32(
                        third_unweighted,
                        root_order_3_sq,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )
                    var recovered_third = _multiply_modulo_int32(
                        _add_three_modulo(
                            first_transformed,
                            second_with_zeta,
                            third_with_zeta_sq,
                            modulus_int,
                        ),
                        inverse_three_multiplier,
                        modulus_int,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                        montgomery_r2,
                    )

                    coefficient_values[first_index_base + inner_offset] = Int32(
                        recovered_first
                    )
                    coefficient_values[
                        second_index_base + inner_offset
                    ] = Int32(recovered_second)
                    coefficient_values[third_index_base + inner_offset] = Int32(
                        recovered_third
                    )

                current_inverse_twiddle = _multiply_modulo_int32(
                    current_inverse_twiddle,
                    inverse_stage_twiddle,
                    modulus_int,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
                    montgomery_r2,
                )
        current_step_size *= 3


fn apply_cyclotomic_pruned_ntt(
    mut coefficient_values: List[Int32],
    m_parameter: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
    if use_naive_modulo:
        barrett_ratio = -1
    var montgomery_neg_inv: Int64 = 0
    var montgomery_r2: Int = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_int)
        montgomery_r2 = compute_montgomery_r2(modulus_int)
    var root_int_standard = Int(root_of_unity)
    var root_int = root_int_standard
    var root_order_3 = mod_pow(root_int_standard, m_parameter, modulus_int)
    var montgomery_one = 1
    if use_montgomery:
        montgomery_one = to_montgomery_domain(
            1,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        root_int = to_montgomery_domain(
            root_int,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
        root_order_3 = to_montgomery_domain(
            root_order_3,
            modulus_int,
            montgomery_neg_inv,
            montgomery_r2,
        )
    var root_order_3_sq = _multiply_modulo_int32(
        root_order_3,
        root_order_3,
        modulus_int,
        barrett_ratio,
        use_montgomery,
        montgomery_neg_inv,
        montgomery_r2,
    )
    var current_twiddle = montgomery_one
    for offset_index in range(m_parameter):
        var twiddle_sq = _multiply_modulo_int32(
            current_twiddle,
            current_twiddle,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )
        var first_index_base = base_offset + offset_index * block_size
        var second_index_base = (
            base_offset + (offset_index + m_parameter) * block_size
        )
        for inner_offset in range(block_size):
            var first_value = Int(
                coefficient_values[first_index_base + inner_offset]
            )
            var second_value = Int(
                coefficient_values[second_index_base + inner_offset]
            )
            var second_weighted_zeta = _multiply_modulo_int32(
                second_value,
                root_order_3,
                modulus_int,
                barrett_ratio,
                use_montgomery,
                montgomery_neg_inv,
                montgomery_r2,
            )
            var first_branch_value = _multiply_modulo_int32(
                _add_modulo(first_value, second_weighted_zeta, modulus_int),
                current_twiddle,
                modulus_int,
                barrett_ratio,
                use_montgomery,
                montgomery_neg_inv,
                montgomery_r2,
            )
            var second_weighted_zeta_sq = _multiply_modulo_int32(
                second_value,
                root_order_3_sq,
                modulus_int,
                barrett_ratio,
                use_montgomery,
                montgomery_neg_inv,
                montgomery_r2,
            )
            var second_branch_value = _multiply_modulo_int32(
                _add_modulo(first_value, second_weighted_zeta_sq, modulus_int),
                twiddle_sq,
                modulus_int,
                barrett_ratio,
                use_montgomery,
                montgomery_neg_inv,
                montgomery_r2,
            )
            coefficient_values[first_index_base + inner_offset] = Int32(
                first_branch_value
            )
            coefficient_values[second_index_base + inner_offset] = Int32(
                second_branch_value
            )
        current_twiddle = _multiply_modulo_int32(
            current_twiddle,
            root_int,
            modulus_int,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
            montgomery_r2,
        )

    var root_for_m = Int32(mod_pow(root_int_standard, 3, modulus_int))
    apply_radix3_dif_ntt(
        coefficient_values,
        m_parameter,
        root_for_m,
        modulus_value,
        block_size,
        base_offset,
        use_montgomery,
        use_naive_modulo,
    )


fn apply_cyclotomic_pruned_intt(
    mut coefficient_values: List[Int32],
    m_parameter: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var inverse_matrix = _build_cyclotomic_inverse_matrix(
        m_parameter,
        Int(root_of_unity),
        modulus_int,
        use_montgomery,
        use_naive_modulo,
    )
    var phi_degree = 2 * m_parameter
    if len(inverse_matrix) != phi_degree * phi_degree:
        return
    _apply_cyclotomic_inverse_matrix(
        coefficient_values,
        inverse_matrix,
        phi_degree,
        modulus_int,
        block_size,
        base_offset,
    )


fn _build_cyclotomic_inverse_matrix(
    m_parameter: Int,
    root_of_unity: Int,
    modulus_int: Int,
    use_montgomery: Bool,
    use_naive_modulo: Bool = False,
) -> List[Int]:
    var phi_degree = 2 * m_parameter
    var forward_matrix = List[Int](length=phi_degree * phi_degree, fill=0)

    for basis_index in range(phi_degree):
        var basis_vector = List[Int32](length=phi_degree, fill=0)
        basis_vector[basis_index] = 1
        apply_cyclotomic_pruned_ntt(
            basis_vector,
            m_parameter,
            Int32(root_of_unity),
            Int32(modulus_int),
            1,
            0,
            use_montgomery,
            use_naive_modulo,
        )
        for row_index in range(phi_degree):
            forward_matrix[row_index * phi_degree + basis_index] = Int(
                basis_vector[row_index]
            )

    var augmented_width = 2 * phi_degree
    var augmented_matrix = List[Int](
        length=phi_degree * augmented_width, fill=0
    )
    for row_index in range(phi_degree):
        for col_index in range(phi_degree):
            augmented_matrix[
                row_index * augmented_width + col_index
            ] = forward_matrix[row_index * phi_degree + col_index]
        augmented_matrix[
            row_index * augmented_width + phi_degree + row_index
        ] = 1

    for pivot_col in range(phi_degree):
        var pivot_row = pivot_col
        while (
            pivot_row < phi_degree
            and augmented_matrix[pivot_row * augmented_width + pivot_col] == 0
        ):
            pivot_row += 1
        if pivot_row == phi_degree:
            var empty_matrix = List[Int]()
            return empty_matrix^

        if pivot_row != pivot_col:
            for col_index in range(augmented_width):
                var top_index = pivot_col * augmented_width + col_index
                var swap_index = pivot_row * augmented_width + col_index
                var swap_value = augmented_matrix[top_index]
                augmented_matrix[top_index] = augmented_matrix[swap_index]
                augmented_matrix[swap_index] = swap_value

        var pivot_value = augmented_matrix[
            pivot_col * augmented_width + pivot_col
        ]
        var inverse_pivot = _mod_inverse(pivot_value, modulus_int)
        if inverse_pivot < 0:
            var empty_matrix = List[Int]()
            return empty_matrix^

        for col_index in range(augmented_width):
            augmented_matrix[pivot_col * augmented_width + col_index] = Int(
                (
                    Int64(
                        augmented_matrix[
                            pivot_col * augmented_width + col_index
                        ]
                    )
                    * Int64(inverse_pivot)
                )
                % Int64(modulus_int)
            )

        for row_index in range(phi_degree):
            if row_index == pivot_col:
                continue
            var elimination_factor = augmented_matrix[
                row_index * augmented_width + pivot_col
            ]
            if elimination_factor == 0:
                continue
            for col_index in range(augmented_width):
                var updated_value = (
                    augmented_matrix[row_index * augmented_width + col_index]
                    - Int(
                        (
                            Int64(elimination_factor)
                            * Int64(
                                augmented_matrix[
                                    pivot_col * augmented_width + col_index
                                ]
                            )
                        )
                        % Int64(modulus_int)
                    )
                ) % modulus_int
                if updated_value < 0:
                    updated_value += modulus_int
                augmented_matrix[
                    row_index * augmented_width + col_index
                ] = updated_value

    var inverse_matrix = List[Int](length=phi_degree * phi_degree, fill=0)
    for row_index in range(phi_degree):
        for col_index in range(phi_degree):
            inverse_matrix[
                row_index * phi_degree + col_index
            ] = augmented_matrix[
                row_index * augmented_width + phi_degree + col_index
            ]
    return inverse_matrix^


fn _apply_cyclotomic_inverse_matrix(
    mut coefficient_values: List[Int32],
    inverse_matrix: List[Int],
    phi_degree: Int,
    modulus_int: Int,
    block_size: Int,
    base_offset: Int,
):
    for inner_offset in range(block_size):
        var input_values = List[Int](length=phi_degree, fill=0)
        for col_index in range(phi_degree):
            input_values[col_index] = Int(
                coefficient_values[
                    base_offset + col_index * block_size + inner_offset
                ]
            )

        var output_values = List[Int](length=phi_degree, fill=0)
        for row_index in range(phi_degree):
            var accumulated_value = 0
            for col_index in range(phi_degree):
                accumulated_value = (
                    accumulated_value
                    + Int(
                        (
                            Int64(
                                inverse_matrix[
                                    row_index * phi_degree + col_index
                                ]
                            )
                            * Int64(input_values[col_index])
                        )
                        % Int64(modulus_int)
                    )
                ) % modulus_int
            output_values[row_index] = accumulated_value

        for row_index in range(phi_degree):
            coefficient_values[
                base_offset + row_index * block_size + inner_offset
            ] = Int32(output_values[row_index])


fn apply_ixw_quotient_ntt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: Int,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
) -> Bool:
    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var expected_total_length = 2 * n_power_of_2 * phi_p_degree
    if len(coefficient_values) != expected_total_length:
        return False

    var root_j = find_primitive_root(modulus_value, 4)
    var root_4n = find_primitive_root(modulus_value, 4 * n_power_of_2)
    var root_p = find_primitive_root(modulus_value, p_power_of_3)
    if root_j < 0 or root_4n < 0 or root_p < 0:
        return False

    apply_i_axis_transform(
        coefficient_values,
        expected_total_length,
        Int32(root_j),
        Int32(modulus_value),
        use_montgomery,
        use_naive_modulo,
    )

    var component_stride = expected_total_length // 2
    apply_radix2_dif_ntt(
        coefficient_values,
        n_power_of_2,
        Int32(root_4n),
        Int32(modulus_value),
        phi_p_degree,
        0,
        use_montgomery,
        use_naive_modulo,
    )
    var root_for_minus_j = mod_pow(root_4n, 3, modulus_value)
    apply_radix2_dif_ntt(
        coefficient_values,
        n_power_of_2,
        Int32(root_for_minus_j),
        Int32(modulus_value),
        phi_p_degree,
        component_stride,
        use_montgomery,
        use_naive_modulo,
    )

    var m_parameter = p_power_of_3 // 3
    var number_of_slices = 2 * n_power_of_2
    for slice_index in range(number_of_slices):
        var slice_base_offset = slice_index * phi_p_degree
        apply_cyclotomic_pruned_ntt(
            coefficient_values,
            m_parameter,
            Int32(root_p),
            Int32(modulus_value),
            1,
            slice_base_offset,
            use_montgomery,
            use_naive_modulo,
        )

    return True


fn apply_ixw_quotient_intt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: Int,
    use_montgomery: Bool = False,
    use_naive_modulo: Bool = False,
) -> Bool:
    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var expected_total_length = 2 * n_power_of_2 * phi_p_degree
    if len(coefficient_values) != expected_total_length:
        return False

    var root_j = find_primitive_root(modulus_value, 4)
    var root_4n = find_primitive_root(modulus_value, 4 * n_power_of_2)
    var root_p = find_primitive_root(modulus_value, p_power_of_3)
    if root_j < 0 or root_4n < 0 or root_p < 0:
        return False

    var m_parameter = p_power_of_3 // 3
    var inverse_cyclotomic_matrix = _build_cyclotomic_inverse_matrix(
        m_parameter,
        root_p,
        modulus_value,
        use_montgomery,
        use_naive_modulo,
    )
    if len(inverse_cyclotomic_matrix) != phi_p_degree * phi_p_degree:
        return False

    var number_of_slices = 2 * n_power_of_2
    for slice_index in range(number_of_slices):
        var slice_base_offset = slice_index * phi_p_degree
        _apply_cyclotomic_inverse_matrix(
            coefficient_values,
            inverse_cyclotomic_matrix,
            phi_p_degree,
            modulus_value,
            1,
            slice_base_offset,
        )

    var component_stride = expected_total_length // 2
    apply_radix2_dif_intt(
        coefficient_values,
        n_power_of_2,
        Int32(root_4n),
        Int32(modulus_value),
        phi_p_degree,
        0,
        use_montgomery,
        use_naive_modulo,
    )
    var root_for_minus_j = mod_pow(root_4n, 3, modulus_value)
    apply_radix2_dif_intt(
        coefficient_values,
        n_power_of_2,
        Int32(root_for_minus_j),
        Int32(modulus_value),
        phi_p_degree,
        component_stride,
        use_montgomery,
        use_naive_modulo,
    )

    apply_inverse_i_axis_transform(
        coefficient_values,
        expected_total_length,
        Int32(root_j),
        Int32(modulus_value),
        use_montgomery,
        use_naive_modulo,
    )
    return True


fn apply_xyw_quotient_ntt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_int: Int,
) -> Bool:
    # Layout: [i_component][y_axis][x_axis][w_axis]
    # Total length is 2 * n * n * phi_p where phi_p = 2 * (p / 3).
    if p_power_of_3 % 3 != 0:
        return False

    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var expected_total_length = 2 * n_power_of_2 * n_power_of_2 * phi_p_degree
    if len(coefficient_values) != expected_total_length:
        return False

    var root_imaginary_unit = find_primitive_root(modulus_int, 4)
    var root_4n = find_primitive_root(modulus_int, 4 * n_power_of_2)
    var root_for_minus_imaginary = mod_pow(root_4n, 3, modulus_int)
    var root_p = find_primitive_root(modulus_int, p_power_of_3)
    if (
        root_imaginary_unit < 0
        or root_4n < 0
        or root_for_minus_imaginary < 0
        or root_p < 0
    ):
        return False

    # i-axis first.
    apply_i_axis_transform(
        coefficient_values,
        expected_total_length,
        Int32(root_imaginary_unit),
        Int32(modulus_int),
    )

    # X-axis per i-component and Y slice.
    var component_length = n_power_of_2 * n_power_of_2 * phi_p_degree
    var y_slice_length = n_power_of_2 * phi_p_degree
    for component_index in range(2):
        var component_base_offset = component_index * component_length
        var x_root_for_component = Int32(root_4n)
        if component_index == 1:
            x_root_for_component = Int32(root_for_minus_imaginary)
        for y_axis_index in range(n_power_of_2):
            var y_slice_offset = (
                component_base_offset + y_axis_index * y_slice_length
            )
            apply_radix2_dif_ntt(
                coefficient_values,
                n_power_of_2,
                x_root_for_component,
                Int32(modulus_int),
                phi_p_degree,
                y_slice_offset,
            )

    # Y-axis per i-component, treating each Y point as a contiguous XW vector.
    for component_index in range(2):
        var component_base_offset = component_index * component_length
        var y_root_for_component = Int32(root_4n)
        if component_index == 1:
            y_root_for_component = Int32(root_for_minus_imaginary)
        apply_radix2_dif_ntt(
            coefficient_values,
            n_power_of_2,
            y_root_for_component,
            Int32(modulus_int),
            y_slice_length,
            component_base_offset,
        )

    # W-axis cyclotomic transform for every [i, y, x] lane.
    var cyclotomic_m_parameter = p_power_of_3 // 3
    var number_of_xy_slices = 2 * n_power_of_2 * n_power_of_2
    for xy_slice_index in range(number_of_xy_slices):
        var xy_slice_offset = xy_slice_index * phi_p_degree
        apply_cyclotomic_pruned_ntt(
            coefficient_values,
            cyclotomic_m_parameter,
            Int32(root_p),
            Int32(modulus_int),
            1,
            xy_slice_offset,
        )

    return True


fn apply_xyw_quotient_intt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_int: Int,
) -> Bool:
    if p_power_of_3 % 3 != 0:
        return False

    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var expected_total_length = 2 * n_power_of_2 * n_power_of_2 * phi_p_degree
    if len(coefficient_values) != expected_total_length:
        return False

    var root_imaginary_unit = find_primitive_root(modulus_int, 4)
    var root_4n = find_primitive_root(modulus_int, 4 * n_power_of_2)
    var root_for_minus_imaginary = mod_pow(root_4n, 3, modulus_int)
    var root_p = find_primitive_root(modulus_int, p_power_of_3)
    if (
        root_imaginary_unit < 0
        or root_4n < 0
        or root_for_minus_imaginary < 0
        or root_p < 0
    ):
        return False

    var component_length = n_power_of_2 * n_power_of_2 * phi_p_degree
    var y_slice_length = n_power_of_2 * phi_p_degree

    var cyclotomic_m_parameter = p_power_of_3 // 3
    var inverse_cyclotomic_matrix = _build_cyclotomic_inverse_matrix(
        cyclotomic_m_parameter, root_p, modulus_int, False
    )
    if len(inverse_cyclotomic_matrix) != phi_p_degree * phi_p_degree:
        return False

    var number_of_xy_slices = 2 * n_power_of_2 * n_power_of_2
    for xy_slice_index in range(number_of_xy_slices):
        var xy_slice_offset = xy_slice_index * phi_p_degree
        _apply_cyclotomic_inverse_matrix(
            coefficient_values,
            inverse_cyclotomic_matrix,
            phi_p_degree,
            modulus_int,
            1,
            xy_slice_offset,
        )

    for component_index in range(2):
        var component_base_offset = component_index * component_length
        var y_root_for_component = Int32(root_4n)
        if component_index == 1:
            y_root_for_component = Int32(root_for_minus_imaginary)
        apply_radix2_dif_intt(
            coefficient_values,
            n_power_of_2,
            y_root_for_component,
            Int32(modulus_int),
            y_slice_length,
            component_base_offset,
        )

    for component_index in range(2):
        var component_base_offset = component_index * component_length
        var x_root_for_component = Int32(root_4n)
        if component_index == 1:
            x_root_for_component = Int32(root_for_minus_imaginary)
        for y_axis_index in range(n_power_of_2):
            var y_slice_offset = (
                component_base_offset + y_axis_index * y_slice_length
            )
            apply_radix2_dif_intt(
                coefficient_values,
                n_power_of_2,
                x_root_for_component,
                Int32(modulus_int),
                phi_p_degree,
                y_slice_offset,
            )

    apply_inverse_i_axis_transform(
        coefficient_values,
        expected_total_length,
        Int32(root_imaginary_unit),
        Int32(modulus_int),
    )
    return True
