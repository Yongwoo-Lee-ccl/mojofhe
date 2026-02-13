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


fn _multiply_modulo_uint32(
    left_value: UInt32,
    right_value: UInt32,
    modulus_value: UInt32,
    barrett_ratio: UInt64,
    use_montgomery: Bool,
    montgomery_neg_inv: UInt64,
) -> UInt32:
    if use_montgomery:
        return multiply_mod_montgomery_with_rhs_mont(
            left_value,
            right_value,
            modulus_value,
            montgomery_neg_inv,
        )
    return multiply_mod_barrett(
        left_value, right_value, modulus_value, barrett_ratio
    )


fn apply_i_axis_transform(
    mut coefficient_values: List[UInt32],
    total_length: Int,
    root_imaginary_unit: UInt32,
    modulus_value: UInt32,
    use_montgomery: Bool = False,
):
    var barrett_ratio = compute_barrett_ratio(modulus_value)
    var montgomery_neg_inv: UInt64 = 0
    var montgomery_r2: UInt32 = 0
    var root_multiplier = root_imaginary_unit
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_value)
        montgomery_r2 = compute_montgomery_r2(modulus_value)
        root_multiplier = to_montgomery_domain(
            root_multiplier,
            modulus_value,
            montgomery_neg_inv,
            montgomery_r2,
        )
    else:
        _ = montgomery_r2

    var half_stride = total_length // 2
    for offset_index in range(half_stride):
        var real_part = coefficient_values[offset_index]
        var imag_part = coefficient_values[offset_index + half_stride]
        var weighted_imag = _multiply_modulo_uint32(
            imag_part,
            root_multiplier,
            modulus_value,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
        )
        coefficient_values[offset_index] = (
            real_part + weighted_imag
        ) % modulus_value
        coefficient_values[offset_index + half_stride] = (
            real_part + modulus_value - weighted_imag
        ) % modulus_value


fn apply_radix2_dif_ntt(
    mut coefficient_values: List[UInt32],
    transform_length: Int,
    root_of_unity: UInt32,
    modulus_value: UInt32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
):
    var barrett_ratio = compute_barrett_ratio(modulus_value)
    var montgomery_neg_inv: UInt64 = 0
    var montgomery_r2: UInt32 = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_value)
        montgomery_r2 = compute_montgomery_r2(modulus_value)
    else:
        _ = montgomery_r2

    var current_step_size = transform_length
    var stage_twiddle = root_of_unity
    var montgomery_one: UInt32 = 1
    if use_montgomery:
        montgomery_one = to_montgomery_domain(
            1,
            modulus_value,
            montgomery_neg_inv,
            montgomery_r2,
        )
        stage_twiddle = to_montgomery_domain(
            stage_twiddle,
            modulus_value,
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
                for inner_offset in range(block_size):
                    var upper_val = coefficient_values[
                        upper_index_base + inner_offset
                    ]
                    var lower_val = coefficient_values[
                        lower_index_base + inner_offset
                    ]
                    var sum_val = (upper_val + lower_val) % modulus_value
                    var diff_val = (
                        upper_val + modulus_value - lower_val
                    ) % modulus_value
                    coefficient_values[
                        upper_index_base + inner_offset
                    ] = sum_val
                    coefficient_values[
                        lower_index_base + inner_offset
                    ] = _multiply_modulo_uint32(
                        diff_val,
                        current_twiddle,
                        modulus_value,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                    )
                current_twiddle = _multiply_modulo_uint32(
                    current_twiddle,
                    stage_twiddle,
                    modulus_value,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
                )
        stage_twiddle = _multiply_modulo_uint32(
            stage_twiddle,
            stage_twiddle,
            modulus_value,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
        )
        current_step_size = half_step


fn apply_radix3_dif_ntt(
    mut coefficient_values: List[UInt32],
    transform_length: Int,
    root_of_unity: UInt32,
    modulus_value: UInt32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
):
    var barrett_ratio = compute_barrett_ratio(modulus_value)
    var montgomery_neg_inv: UInt64 = 0
    var montgomery_r2: UInt32 = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_value)
        montgomery_r2 = compute_montgomery_r2(modulus_value)
    else:
        _ = montgomery_r2

    var root_standard = root_of_unity
    var root_order_3 = UInt32(
        mod_pow(
            UInt64(root_standard),
            UInt64(transform_length // 3),
            UInt64(modulus_value),
        )
    )
    var montgomery_one: UInt32 = 1
    var stage_twiddle = root_standard

    if use_montgomery:
        montgomery_one = to_montgomery_domain(
            1,
            modulus_value,
            montgomery_neg_inv,
            montgomery_r2,
        )
        stage_twiddle = to_montgomery_domain(
            stage_twiddle,
            modulus_value,
            montgomery_neg_inv,
            montgomery_r2,
        )
        root_order_3 = to_montgomery_domain(
            root_order_3,
            modulus_value,
            montgomery_neg_inv,
            montgomery_r2,
        )

    var root_order_3_sq = _multiply_modulo_uint32(
        root_order_3,
        root_order_3,
        modulus_value,
        barrett_ratio,
        use_montgomery,
        montgomery_neg_inv,
    )

    var current_step_size = transform_length
    while current_step_size >= 3:
        var third_step = current_step_size // 3
        for group_start in range(0, transform_length, current_step_size):
            var twiddle_1 = montgomery_one
            for butterfly_index in range(third_step):
                var twiddle_2 = _multiply_modulo_uint32(
                    twiddle_1,
                    twiddle_1,
                    modulus_value,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
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
                    var first_value = coefficient_values[
                        first_index_base + inner_offset
                    ]
                    var second_value = coefficient_values[
                        second_index_base + inner_offset
                    ]
                    var third_value = coefficient_values[
                        third_index_base + inner_offset
                    ]
                    var sum_all = (
                        first_value + second_value + third_value
                    ) % modulus_value

                    var second_weighted_zeta = _multiply_modulo_uint32(
                        second_value,
                        root_order_3,
                        modulus_value,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                    )
                    var third_weighted_zeta_sq = _multiply_modulo_uint32(
                        third_value,
                        root_order_3_sq,
                        modulus_value,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                    )
                    var sum_zeta = (
                        first_value
                        + second_weighted_zeta
                        + third_weighted_zeta_sq
                    ) % modulus_value

                    var second_weighted_zeta_sq = _multiply_modulo_uint32(
                        second_value,
                        root_order_3_sq,
                        modulus_value,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                    )
                    var third_weighted_zeta = _multiply_modulo_uint32(
                        third_value,
                        root_order_3,
                        modulus_value,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                    )
                    var sum_zeta_sq = (
                        first_value
                        + second_weighted_zeta_sq
                        + third_weighted_zeta
                    ) % modulus_value

                    coefficient_values[
                        first_index_base + inner_offset
                    ] = sum_all
                    coefficient_values[
                        second_index_base + inner_offset
                    ] = _multiply_modulo_uint32(
                        sum_zeta,
                        twiddle_1,
                        modulus_value,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                    )
                    coefficient_values[
                        third_index_base + inner_offset
                    ] = _multiply_modulo_uint32(
                        sum_zeta_sq,
                        twiddle_2,
                        modulus_value,
                        barrett_ratio,
                        use_montgomery,
                        montgomery_neg_inv,
                    )
                twiddle_1 = _multiply_modulo_uint32(
                    twiddle_1,
                    stage_twiddle,
                    modulus_value,
                    barrett_ratio,
                    use_montgomery,
                    montgomery_neg_inv,
                )
        var stage_twiddle_squared = _multiply_modulo_uint32(
            stage_twiddle,
            stage_twiddle,
            modulus_value,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
        )
        stage_twiddle = _multiply_modulo_uint32(
            stage_twiddle_squared,
            stage_twiddle,
            modulus_value,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
        )
        current_step_size = third_step


fn apply_cyclotomic_pruned_ntt(
    mut coefficient_values: List[UInt32],
    m_parameter: Int,
    root_of_unity: UInt32,
    modulus_value: UInt32,
    block_size: Int,
    base_offset: Int = 0,
    use_montgomery: Bool = False,
):
    var barrett_ratio = compute_barrett_ratio(modulus_value)
    var montgomery_neg_inv: UInt64 = 0
    var montgomery_r2: UInt32 = 0
    if use_montgomery:
        montgomery_neg_inv = compute_montgomery_neg_inv(modulus_value)
        montgomery_r2 = compute_montgomery_r2(modulus_value)
    else:
        _ = montgomery_r2

    var root_standard = root_of_unity
    var root_order_3 = UInt32(
        mod_pow(
            UInt64(root_standard), UInt64(m_parameter), UInt64(modulus_value)
        )
    )
    var montgomery_one: UInt32 = 1
    var stage_twiddle = root_standard

    if use_montgomery:
        montgomery_one = to_montgomery_domain(
            1,
            modulus_value,
            montgomery_neg_inv,
            montgomery_r2,
        )
        stage_twiddle = to_montgomery_domain(
            stage_twiddle,
            modulus_value,
            montgomery_neg_inv,
            montgomery_r2,
        )
        root_order_3 = to_montgomery_domain(
            root_order_3,
            modulus_value,
            montgomery_neg_inv,
            montgomery_r2,
        )

    var root_order_3_sq = _multiply_modulo_uint32(
        root_order_3,
        root_order_3,
        modulus_value,
        barrett_ratio,
        use_montgomery,
        montgomery_neg_inv,
    )

    var current_twiddle = montgomery_one
    for offset_index in range(m_parameter):
        var twiddle_sq = _multiply_modulo_uint32(
            current_twiddle,
            current_twiddle,
            modulus_value,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
        )
        var first_index_base = base_offset + offset_index * block_size
        var second_index_base = (
            base_offset + (offset_index + m_parameter) * block_size
        )
        for inner_offset in range(block_size):
            var first_value = coefficient_values[
                first_index_base + inner_offset
            ]
            var second_value = coefficient_values[
                second_index_base + inner_offset
            ]
            var second_weighted_zeta = _multiply_modulo_uint32(
                second_value,
                root_order_3,
                modulus_value,
                barrett_ratio,
                use_montgomery,
                montgomery_neg_inv,
            )
            var first_branch_value = _multiply_modulo_uint32(
                (first_value + second_weighted_zeta) % modulus_value,
                current_twiddle,
                modulus_value,
                barrett_ratio,
                use_montgomery,
                montgomery_neg_inv,
            )
            var second_weighted_zeta_sq = _multiply_modulo_uint32(
                second_value,
                root_order_3_sq,
                modulus_value,
                barrett_ratio,
                use_montgomery,
                montgomery_neg_inv,
            )
            var second_branch_value = _multiply_modulo_uint32(
                (first_value + second_weighted_zeta_sq) % modulus_value,
                twiddle_sq,
                modulus_value,
                barrett_ratio,
                use_montgomery,
                montgomery_neg_inv,
            )
            coefficient_values[
                first_index_base + inner_offset
            ] = first_branch_value
            coefficient_values[
                second_index_base + inner_offset
            ] = second_branch_value
        current_twiddle = _multiply_modulo_uint32(
            current_twiddle,
            stage_twiddle,
            modulus_value,
            barrett_ratio,
            use_montgomery,
            montgomery_neg_inv,
        )

    var root_for_m = UInt32(
        mod_pow(UInt64(root_standard), UInt64(3), UInt64(modulus_value))
    )
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
    mut coefficient_values: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: UInt32,
    use_montgomery: Bool = False,
) -> Bool:
    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var expected_total_length = 2 * n_power_of_2 * phi_p_degree
    if len(coefficient_values) != expected_total_length:
        return False

    var root_j = find_primitive_root(UInt64(modulus_value), 4)
    var root_4n = find_primitive_root(
        UInt64(modulus_value), UInt64(4 * n_power_of_2)
    )
    var root_p = find_primitive_root(
        UInt64(modulus_value), UInt64(p_power_of_3)
    )
    if root_j == 0 or root_4n == 0 or root_p == 0:
        return False

    apply_i_axis_transform(
        coefficient_values,
        expected_total_length,
        UInt32(root_j),
        modulus_value,
        use_montgomery,
    )

    var component_stride = expected_total_length // 2
    apply_radix2_dif_ntt(
        coefficient_values,
        n_power_of_2,
        UInt32(root_4n),
        modulus_value,
        phi_p_degree,
        0,
        use_montgomery,
    )
    var root_for_minus_j = mod_pow(root_4n, 3, UInt64(modulus_value))
    apply_radix2_dif_ntt(
        coefficient_values,
        n_power_of_2,
        UInt32(root_for_minus_j),
        modulus_value,
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
            UInt32(root_p),
            modulus_value,
            1,
            slice_base_offset,
            use_montgomery,
        )

    return True


fn apply_xyw_quotient_ntt_inplace(
    mut coefficient_values: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: UInt32,
    use_montgomery: Bool = False,
) -> Bool:
    # Layout: [i_component][y_axis][x_axis][w_axis]
    # Total length is 2 * n * n * phi_p where phi_p = 2 * (p / 3).
    if p_power_of_3 % 3 != 0:
        return False

    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var expected_total_length = 2 * n_power_of_2 * n_power_of_2 * phi_p_degree
    if len(coefficient_values) != expected_total_length:
        return False

    var root_imaginary_unit = find_primitive_root(UInt64(modulus_value), 4)
    var root_4n = find_primitive_root(
        UInt64(modulus_value), UInt64(4 * n_power_of_2)
    )
    var root_for_minus_imaginary = mod_pow(
        root_4n, 3, UInt64(modulus_value)
    )
    var root_p = find_primitive_root(
        UInt64(modulus_value), UInt64(p_power_of_3)
    )
    if (
        root_imaginary_unit == 0
        or root_4n == 0
        or root_for_minus_imaginary == 0
        or root_p == 0
    ):
        return False

    # i-axis first.
    apply_i_axis_transform(
        coefficient_values,
        expected_total_length,
        UInt32(root_imaginary_unit),
        modulus_value,
        use_montgomery,
    )

    # X-axis per i-component and Y slice.
    var component_length = n_power_of_2 * n_power_of_2 * phi_p_degree
    var y_slice_length = n_power_of_2 * phi_p_degree
    for component_index in range(2):
        var component_base_offset = component_index * component_length
        var x_root_for_component = UInt32(root_4n)
        if component_index == 1:
            x_root_for_component = UInt32(root_for_minus_imaginary)
        for y_axis_index in range(n_power_of_2):
            var y_slice_offset = (
                component_base_offset + y_axis_index * y_slice_length
            )
            apply_radix2_dif_ntt(
                coefficient_values,
                n_power_of_2,
                x_root_for_component,
                modulus_value,
                phi_p_degree,
                y_slice_offset,
                use_montgomery,
            )

    # Y-axis per i-component, treating each Y point as a contiguous XW vector.
    for component_index in range(2):
        var component_base_offset = component_index * component_length
        var y_root_for_component = UInt32(root_4n)
        if component_index == 1:
            y_root_for_component = UInt32(root_for_minus_imaginary)
        apply_radix2_dif_ntt(
            coefficient_values,
            n_power_of_2,
            y_root_for_component,
            modulus_value,
            y_slice_length,
            component_base_offset,
            use_montgomery,
        )

    # W-axis cyclotomic transform for every [i, y, x] lane.
    var cyclotomic_m_parameter = p_power_of_3 // 3
    var number_of_xy_slices = 2 * n_power_of_2 * n_power_of_2
    for xy_slice_index in range(number_of_xy_slices):
        var xy_slice_offset = xy_slice_index * phi_p_degree
        apply_cyclotomic_pruned_ntt(
            coefficient_values,
            cyclotomic_m_parameter,
            UInt32(root_p),
            modulus_value,
            1,
            xy_slice_offset,
            use_montgomery,
        )

    return True
