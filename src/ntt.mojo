from collections import List
from src.modular import mod_pow


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
        var weighted_imag = Int(
            (Int128(imag_part) * Int128(Int(root_imaginary_unit)))
            % Int128(modulus_int)
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
    while current_step_size > 1:
        var half_step = current_step_size // 2
        var step_twiddle = mod_pow(
            root_int, transform_length // current_step_size, modulus_int
        )
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
                        (Int128(diff_val) * Int128(current_twiddle))
                        % Int128(modulus_int)
                    )
                current_twiddle = Int(
                    (Int128(current_twiddle) * Int128(step_twiddle))
                    % Int128(modulus_int)
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
    var root_order_3_sq = Int(
        (Int128(root_order_3) * Int128(root_order_3)) % Int128(modulus_int)
    )
    while current_step_size >= 3:
        var third_step = current_step_size // 3
        var step_twiddle = mod_pow(
            root_int, transform_length // current_step_size, modulus_int
        )
        for group_start in range(0, transform_length, current_step_size):
            var twiddle_1 = 1
            for butterfly_index in range(third_step):
                var twiddle_2 = Int(
                    (Int128(twiddle_1) * Int128(twiddle_1))
                    % Int128(modulus_int)
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
                    var second_weighted_zeta = Int(
                        (Int128(second_value) * Int128(root_order_3))
                        % Int128(modulus_int)
                    )
                    var third_weighted_zeta_sq = Int(
                        (Int128(third_value) * Int128(root_order_3_sq))
                        % Int128(modulus_int)
                    )
                    var sum_zeta = (
                        first_value
                        + second_weighted_zeta
                        + third_weighted_zeta_sq
                    ) % modulus_int
                    var second_weighted_zeta_sq = Int(
                        (Int128(second_value) * Int128(root_order_3_sq))
                        % Int128(modulus_int)
                    )
                    var third_weighted_zeta = Int(
                        (Int128(third_value) * Int128(root_order_3))
                        % Int128(modulus_int)
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
                        (Int128(sum_zeta) * Int128(twiddle_1))
                        % Int128(modulus_int)
                    )
                    coefficient_values[third_index_base + inner_offset] = Int32(
                        (Int128(sum_zeta_sq) * Int128(twiddle_2))
                        % Int128(modulus_int)
                    )
                twiddle_1 = Int(
                    (Int128(twiddle_1) * Int128(step_twiddle))
                    % Int128(modulus_int)
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
    var root_order_3_sq = Int(
        (Int128(root_order_3) * Int128(root_order_3)) % Int128(modulus_int)
    )
    var current_twiddle = 1
    for offset_index in range(m_parameter):
        var twiddle_sq = Int(
            (Int128(current_twiddle) * Int128(current_twiddle))
            % Int128(modulus_int)
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
            var second_weighted_zeta = Int(
                (Int128(second_value) * Int128(root_order_3))
                % Int128(modulus_int)
            )
            var first_branch_value = Int(
                (
                    Int128(first_value + second_weighted_zeta)
                    * Int128(current_twiddle)
                )
                % Int128(modulus_int)
            )
            var second_weighted_zeta_sq = Int(
                (Int128(second_value) * Int128(root_order_3_sq))
                % Int128(modulus_int)
            )
            var second_branch_value = Int(
                (
                    Int128(first_value + second_weighted_zeta_sq)
                    * Int128(twiddle_sq)
                )
                % Int128(modulus_int)
            )
            coefficient_values[first_index_base + inner_offset] = Int32(
                first_branch_value
            )
            coefficient_values[second_index_base + inner_offset] = Int32(
                second_branch_value
            )
        current_twiddle = Int(
            (Int128(current_twiddle) * Int128(root_int)) % Int128(modulus_int)
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
