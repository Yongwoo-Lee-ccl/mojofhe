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


fn apply_i_axis_transform(
    mut coefficient_values: List[Int32],
    total_length: Int,
    root_imaginary_unit: Int32,
    modulus_value: Int32,
    use_montgomery: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
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


fn apply_radix2_dif_ntt(
    mut coefficient_values: List[Int32],
    transform_length: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
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
                    var sum_val = (upper_val + lower_val) % modulus_int
                    var diff_val = (
                        upper_val - lower_val + modulus_int
                    ) % modulus_int
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


fn apply_radix3_dif_ntt(
    mut coefficient_values: List[Int32],
    transform_length: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
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
                    var sum_all = (
                        first_value + second_value + third_value
                    ) % modulus_int
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
                    var sum_zeta = (
                        first_value
                        + second_weighted_zeta
                        + third_weighted_zeta_sq
                    ) % modulus_int
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
                    var sum_zeta_sq = (
                        first_value
                        + second_weighted_zeta_sq
                        + third_weighted_zeta
                    ) % modulus_int
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


fn apply_cyclotomic_pruned_ntt(
    mut coefficient_values: List[Int32],
    m_parameter: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
):
    var modulus_int = Int(modulus_value)
    var barrett_ratio = compute_barrett_ratio(modulus_int)
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
                first_value + second_weighted_zeta,
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
                first_value + second_weighted_zeta_sq,
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
    )
    apply_radix3_dif_ntt(
        coefficient_values,
        m_parameter,
        root_for_m,
        modulus_value,
        block_size,
        base_offset + (m_parameter * block_size),
        use_montgomery,
    )


fn apply_ixw_quotient_ntt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: Int,
    use_montgomery: Bool = False,
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
        )

    return True
