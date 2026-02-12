from collections import List
from src.modular import find_primitive_root, mod_pow


fn _multiply_modulo_int32(
    left_value: Int, right_value: Int, modulus_int: Int
) -> Int:
    return Int((Int64(left_value) * Int64(right_value)) % Int64(modulus_int))


fn apply_i_axis_transform(
    mut coefficient_values: List[Int32],
    total_length: Int,
    root_imaginary_unit: Int32,
    modulus_value: Int32,
):
    var modulus_int = Int(modulus_value)
    var half_stride = total_length // 2
    for offset_index in range(half_stride):
        var real_part = Int(coefficient_values[offset_index])
        var imag_part = Int(coefficient_values[offset_index + half_stride])
        var weighted_imag = _multiply_modulo_int32(
            imag_part, Int(root_imaginary_unit), modulus_int
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
):
    var modulus_int = Int(modulus_value)
    var root_int = Int(root_of_unity)
    var current_step_size = transform_length
    var stage_twiddle = root_int
    while current_step_size > 1:
        var half_step = current_step_size // 2
        for group_start in range(0, transform_length, current_step_size):
            var current_twiddle = 1
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
                            diff_val, current_twiddle, modulus_int
                        )
                    )
                current_twiddle = _multiply_modulo_int32(
                    current_twiddle, stage_twiddle, modulus_int
                )
        stage_twiddle = _multiply_modulo_int32(
            stage_twiddle, stage_twiddle, modulus_int
        )
        current_step_size = half_step


fn apply_radix3_dif_ntt(
    mut coefficient_values: List[Int32],
    transform_length: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
):
    var modulus_int = Int(modulus_value)
    var root_int = Int(root_of_unity)
    var current_step_size = transform_length
    var root_order_3 = mod_pow(root_int, transform_length // 3, modulus_int)
    var root_order_3_sq = _multiply_modulo_int32(
        root_order_3, root_order_3, modulus_int
    )
    var stage_twiddle = root_int
    while current_step_size >= 3:
        var third_step = current_step_size // 3
        for group_start in range(0, transform_length, current_step_size):
            var twiddle_1 = 1
            for butterfly_index in range(third_step):
                var twiddle_2 = _multiply_modulo_int32(
                    twiddle_1, twiddle_1, modulus_int
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
                        second_value, root_order_3, modulus_int
                    )
                    var third_weighted_zeta_sq = _multiply_modulo_int32(
                        third_value, root_order_3_sq, modulus_int
                    )
                    var sum_zeta = (
                        first_value
                        + second_weighted_zeta
                        + third_weighted_zeta_sq
                    ) % modulus_int
                    var second_weighted_zeta_sq = _multiply_modulo_int32(
                        second_value, root_order_3_sq, modulus_int
                    )
                    var third_weighted_zeta = _multiply_modulo_int32(
                        third_value, root_order_3, modulus_int
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
                        _multiply_modulo_int32(sum_zeta, twiddle_1, modulus_int)
                    )
                    coefficient_values[third_index_base + inner_offset] = Int32(
                        _multiply_modulo_int32(
                            sum_zeta_sq, twiddle_2, modulus_int
                        )
                    )
                twiddle_1 = _multiply_modulo_int32(
                    twiddle_1, stage_twiddle, modulus_int
                )
        var stage_twiddle_squared = _multiply_modulo_int32(
            stage_twiddle, stage_twiddle, modulus_int
        )
        stage_twiddle = _multiply_modulo_int32(
            stage_twiddle_squared, stage_twiddle, modulus_int
        )
        current_step_size = third_step


fn apply_cyclotomic_pruned_ntt(
    mut coefficient_values: List[Int32],
    m_parameter: Int,
    root_of_unity: Int32,
    modulus_value: Int32,
    block_size: Int,
    base_offset: Int = 0,
):
    var modulus_int = Int(modulus_value)
    var root_int = Int(root_of_unity)
    var root_order_3 = mod_pow(root_int, m_parameter, modulus_int)
    var root_order_3_sq = _multiply_modulo_int32(
        root_order_3, root_order_3, modulus_int
    )
    var current_twiddle = 1
    for offset_index in range(m_parameter):
        var twiddle_sq = _multiply_modulo_int32(
            current_twiddle, current_twiddle, modulus_int
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
                second_value, root_order_3, modulus_int
            )
            var first_branch_value = _multiply_modulo_int32(
                first_value + second_weighted_zeta,
                current_twiddle,
                modulus_int,
            )
            var second_weighted_zeta_sq = _multiply_modulo_int32(
                second_value, root_order_3_sq, modulus_int
            )
            var second_branch_value = _multiply_modulo_int32(
                first_value + second_weighted_zeta_sq, twiddle_sq, modulus_int
            )
            coefficient_values[first_index_base + inner_offset] = Int32(
                first_branch_value
            )
            coefficient_values[second_index_base + inner_offset] = Int32(
                second_branch_value
            )
        current_twiddle = _multiply_modulo_int32(
            current_twiddle, root_int, modulus_int
        )

    var root_for_m = Int32(mod_pow(root_int, 3, modulus_int))
    apply_radix3_dif_ntt(
        coefficient_values,
        m_parameter,
        root_for_m,
        modulus_value,
        block_size,
        base_offset,
    )
    apply_radix3_dif_ntt(
        coefficient_values,
        m_parameter,
        root_for_m,
        modulus_value,
        block_size,
        base_offset + (m_parameter * block_size),
    )


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
