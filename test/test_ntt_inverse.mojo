from collections import List
from testing import assert_true, TestSuite
from src.modular import find_suitable_q
from src.ntt import (
    apply_ixw_quotient_ntt_inplace,
    apply_ixw_quotient_intt_inplace,
    apply_xyw_quotient_ntt_inplace,
    apply_xyw_quotient_intt_inplace,
)


fn _normalize_mod(value: Int, modulus_value: Int) -> Int:
    var reduced_value = value % modulus_value
    if reduced_value < 0:
        reduced_value += modulus_value
    return reduced_value


fn _ixw_index(
    component_index: Int,
    x_index: Int,
    w_index: Int,
    n_value: Int,
    phi_degree: Int,
) -> Int:
    return (
        component_index * (n_value * phi_degree)
        + x_index * phi_degree
        + w_index
    )


fn _xyw_index(
    component_index: Int,
    y_index: Int,
    x_index: Int,
    w_index: Int,
    n_value: Int,
    phi_degree: Int,
) -> Int:
    var component_stride = n_value * n_value * phi_degree
    return (
        component_index * component_stride
        + y_index * (n_value * phi_degree)
        + x_index * phi_degree
        + w_index
    )


fn _copy_list_u32(values: List[UInt32]) -> List[UInt32]:
    var copied_values = List[UInt32]()
    for value in values:
        copied_values.append(value)
    return copied_values^


fn _build_ixw_values(
    n_value: Int, p_value: Int, modulus_value: UInt32
) -> List[UInt32]:
    var phi_degree = 2 * (p_value // 3)
    var total_length = 2 * n_value * phi_degree
    var values = List[UInt32](length=total_length, fill=0)
    for component_idx in range(2):
        for x_idx in range(n_value):
            for w_idx in range(phi_degree):
                var seed_value = (
                    component_idx * 131
                    + x_idx * 47
                    + w_idx * 19
                    + n_value * 23
                    + p_value * 17
                )
                values[
                    _ixw_index(component_idx, x_idx, w_idx, n_value, phi_degree)
                ] = UInt32(_normalize_mod(seed_value, Int(modulus_value)))
    return values^


fn _build_xyw_values(
    n_value: Int, p_value: Int, modulus_value: UInt32
) -> List[UInt32]:
    var phi_degree = 2 * (p_value // 3)
    var total_length = 2 * n_value * n_value * phi_degree
    var values = List[UInt32](length=total_length, fill=0)
    for component_idx in range(2):
        for y_idx in range(n_value):
            for x_idx in range(n_value):
                for w_idx in range(phi_degree):
                    var seed_value = (
                        component_idx * 149
                        + y_idx * 41
                        + x_idx * 29
                        + w_idx * 13
                        + n_value * 31
                        + p_value * 11
                    )
                    values[
                        _xyw_index(
                            component_idx,
                            y_idx,
                            x_idx,
                            w_idx,
                            n_value,
                            phi_degree,
                        )
                    ] = UInt32(_normalize_mod(seed_value, Int(modulus_value)))
    return values^


fn _naive_ixw_multiply_mod_quotient(
    left_values: List[UInt32],
    right_values: List[UInt32],
    n_value: Int,
    p_value: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    var m_param = p_value // 3
    var phi_deg = 2 * m_param
    var w_temp_deg = 4 * m_param
    var x_temp_deg = 2 * n_value

    var temp_comp_0 = List[Int](length=x_temp_deg * w_temp_deg, fill=0)
    var temp_comp_1 = List[Int](length=x_temp_deg * w_temp_deg, fill=0)

    for left_x in range(n_value):
        for left_w in range(phi_deg):
            var left_0 = Int(
                left_values[_ixw_index(0, left_x, left_w, n_value, phi_deg)]
            )
            var left_1 = Int(
                left_values[_ixw_index(1, left_x, left_w, n_value, phi_deg)]
            )
            for right_x in range(n_value):
                for right_w in range(phi_deg):
                    var right_0 = Int(
                        right_values[
                            _ixw_index(0, right_x, right_w, n_value, phi_deg)
                        ]
                    )
                    var right_1 = Int(
                        right_values[
                            _ixw_index(1, right_x, right_w, n_value, phi_deg)
                        ]
                    )
                    var prod_0 = left_0 * right_0 - left_1 * right_1
                    var prod_1 = left_0 * right_1 + left_1 * right_0
                    var x_deg = left_x + right_x
                    var w_deg = left_w + right_w
                    var temp_idx = x_deg * w_temp_deg + w_deg
                    temp_comp_0[temp_idx] = _normalize_mod(
                        temp_comp_0[temp_idx] + prod_0,
                        Int(modulus_value),
                    )
                    temp_comp_1[temp_idx] = _normalize_mod(
                        temp_comp_1[temp_idx] + prod_1,
                        Int(modulus_value),
                    )

    var red_comp_0 = List[Int](length=n_value * w_temp_deg, fill=0)
    var red_comp_1 = List[Int](length=n_value * w_temp_deg, fill=0)
    for x_deg in range(2 * n_value - 1):
        var red_x = x_deg
        if x_deg >= n_value:
            red_x = x_deg - n_value
        for w_deg in range(w_temp_deg):
            var src_idx = x_deg * w_temp_deg + w_deg
            var target_idx = red_x * w_temp_deg + w_deg
            var comp_0_val = temp_comp_0[src_idx]
            var comp_1_val = temp_comp_1[src_idx]
            red_comp_0[target_idx] = _normalize_mod(
                red_comp_0[target_idx] + comp_0_val,
                Int(modulus_value),
            )
            red_comp_1[target_idx] = _normalize_mod(
                red_comp_1[target_idx] + comp_1_val,
                Int(modulus_value),
            )

    var expected_vals = List[UInt32](length=2 * n_value * phi_deg, fill=0)
    for x_idx in range(n_value):
        var w_deg = w_temp_deg - 1
        while w_deg >= phi_deg:
            var lane_idx = x_idx * w_temp_deg + w_deg
            var lane_0 = red_comp_0[lane_idx]
            var lane_1 = red_comp_1[lane_idx]
            if lane_0 != 0 or lane_1 != 0:
                var idx_minus_m = x_idx * w_temp_deg + (w_deg - m_param)
                var idx_minus_2m = x_idx * w_temp_deg + (w_deg - phi_deg)
                red_comp_0[idx_minus_m] = _normalize_mod(
                    red_comp_0[idx_minus_m] - lane_0,
                    Int(modulus_value),
                )
                red_comp_0[idx_minus_2m] = _normalize_mod(
                    red_comp_0[idx_minus_2m] - lane_0,
                    Int(modulus_value),
                )
                red_comp_1[idx_minus_m] = _normalize_mod(
                    red_comp_1[idx_minus_m] - lane_1,
                    Int(modulus_value),
                )
                red_comp_1[idx_minus_2m] = _normalize_mod(
                    red_comp_1[idx_minus_2m] - lane_1,
                    Int(modulus_value),
                )
            red_comp_0[lane_idx] = 0
            red_comp_1[lane_idx] = 0
            w_deg -= 1

        for w_idx in range(phi_deg):
            var red_idx = x_idx * w_temp_deg + w_idx
            expected_vals[
                _ixw_index(0, x_idx, w_idx, n_value, phi_deg)
            ] = UInt32(red_comp_0[red_idx])
            expected_vals[
                _ixw_index(1, x_idx, w_idx, n_value, phi_deg)
            ] = UInt32(red_comp_1[red_idx])

    return expected_vals^


fn _naive_xyw_multiply_mod_quotient(
    left_values: List[UInt32],
    right_values: List[UInt32],
    n_value: Int,
    p_value: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    var m_param = p_value // 3
    var phi_deg = 2 * m_param
    var w_temp_deg = 4 * m_param
    var x_temp_deg = 2 * n_value
    var y_temp_deg = 2 * n_value
    var temp_plane_sz = x_temp_deg * w_temp_deg
    var temp_total_sz = y_temp_deg * temp_plane_sz

    var temp_comp_0 = List[Int](length=temp_total_sz, fill=0)
    var temp_comp_1 = List[Int](length=temp_total_sz, fill=0)

    for left_y in range(n_value):
        for left_x in range(n_value):
            for left_w in range(phi_deg):
                var left_0 = Int(
                    left_values[
                        _xyw_index(0, left_y, left_x, left_w, n_value, phi_deg)
                    ]
                )
                var left_1 = Int(
                    left_values[
                        _xyw_index(1, left_y, left_x, left_w, n_value, phi_deg)
                    ]
                )
                for right_y in range(n_value):
                    for right_x in range(n_value):
                        for right_w in range(phi_deg):
                            var right_0 = Int(
                                right_values[
                                    _xyw_index(
                                        0,
                                        right_y,
                                        right_x,
                                        right_w,
                                        n_value,
                                        phi_deg,
                                    )
                                ]
                            )
                            var right_1 = Int(
                                right_values[
                                    _xyw_index(
                                        1,
                                        right_y,
                                        right_x,
                                        right_w,
                                        n_value,
                                        phi_deg,
                                    )
                                ]
                            )
                            var prod_0 = left_0 * right_0 - left_1 * right_1
                            var prod_1 = left_0 * right_1 + left_1 * right_0
                            var y_deg = left_y + right_y
                            var x_deg = left_x + right_x
                            var w_deg = left_w + right_w
                            var temp_idx = (
                                y_deg * temp_plane_sz
                                + x_deg * w_temp_deg
                                + w_deg
                            )
                            temp_comp_0[temp_idx] = _normalize_mod(
                                temp_comp_0[temp_idx] + prod_0,
                                Int(modulus_value),
                            )
                            temp_comp_1[temp_idx] = _normalize_mod(
                                temp_comp_1[temp_idx] + prod_1,
                                Int(modulus_value),
                            )

    var red_plane_sz = n_value * w_temp_deg
    var red_total_sz = n_value * red_plane_sz
    var red_comp_0 = List[Int](length=red_total_sz, fill=0)
    var red_comp_1 = List[Int](length=red_total_sz, fill=0)
    for y_deg in range(2 * n_value - 1):
        var red_y = y_deg
        if y_deg >= n_value:
            red_y = y_deg - n_value
        for x_deg in range(2 * n_value - 1):
            var red_x = x_deg
            if x_deg >= n_value:
                red_x = x_deg - n_value
            for w_deg in range(w_temp_deg):
                var src_idx = y_deg * temp_plane_sz + x_deg * w_temp_deg + w_deg
                var target_idx = (
                    red_y * red_plane_sz + red_x * w_temp_deg + w_deg
                )
                var comp_0_val = temp_comp_0[src_idx]
                var comp_1_val = temp_comp_1[src_idx]
                red_comp_0[target_idx] = _normalize_mod(
                    red_comp_0[target_idx] + comp_0_val,
                    Int(modulus_value),
                )
                red_comp_1[target_idx] = _normalize_mod(
                    red_comp_1[target_idx] + comp_1_val,
                    Int(modulus_value),
                )

    var expected_vals = List[UInt32](
        length=2 * n_value * n_value * phi_deg, fill=0
    )
    for y_idx in range(n_value):
        for x_idx in range(n_value):
            var w_deg = w_temp_deg - 1
            while w_deg >= phi_deg:
                var lane_idx = y_idx * red_plane_sz + x_idx * w_temp_deg + w_deg
                var lane_0 = red_comp_0[lane_idx]
                var lane_1 = red_comp_1[lane_idx]
                if lane_0 != 0 or lane_1 != 0:
                    var idx_minus_m = (
                        y_idx * red_plane_sz
                        + x_idx * w_temp_deg
                        + (w_deg - m_param)
                    )
                    var idx_minus_2m = (
                        y_idx * red_plane_sz
                        + x_idx * w_temp_deg
                        + (w_deg - phi_deg)
                    )
                    red_comp_0[idx_minus_m] = _normalize_mod(
                        red_comp_0[idx_minus_m] - lane_0,
                        Int(modulus_value),
                    )
                    red_comp_0[idx_minus_2m] = _normalize_mod(
                        red_comp_0[idx_minus_2m] - lane_0,
                        Int(modulus_value),
                    )
                    red_comp_1[idx_minus_m] = _normalize_mod(
                        red_comp_1[idx_minus_m] - lane_1,
                        Int(modulus_value),
                    )
                    red_comp_1[idx_minus_2m] = _normalize_mod(
                        red_comp_1[idx_minus_2m] - lane_1,
                        Int(modulus_value),
                    )
                red_comp_0[lane_idx] = 0
                red_comp_1[lane_idx] = 0
                w_deg -= 1

            for w_idx in range(phi_deg):
                var red_idx = y_idx * red_plane_sz + x_idx * w_temp_deg + w_idx
                expected_vals[
                    _xyw_index(0, y_idx, x_idx, w_idx, n_value, phi_deg)
                ] = UInt32(red_comp_0[red_idx])
                expected_vals[
                    _xyw_index(1, y_idx, x_idx, w_idx, n_value, phi_deg)
                ] = UInt32(red_comp_1[red_idx])

    return expected_vals^


fn _assert_ixw_roundtrip(
    n_value: Int, p_value: Int, modulus_value: UInt32, use_montgomery: Bool
) raises:
    var original_values = _build_ixw_values(n_value, p_value, modulus_value)
    var transformed_values = _copy_list_u32(original_values)
    var ntt_success = apply_ixw_quotient_ntt_inplace(
        transformed_values,
        n_value,
        p_value,
        modulus_value,
        use_montgomery,
    )
    assert_true(ntt_success, "IXW NTT should succeed")
    var intt_success = apply_ixw_quotient_intt_inplace(
        transformed_values,
        n_value,
        p_value,
        modulus_value,
        use_montgomery,
    )
    assert_true(intt_success, "IXW inverse NTT should succeed")
    for coeff_idx in range(len(original_values)):
        assert_true(
            transformed_values[coeff_idx] == original_values[coeff_idx],
            "IXW NTT->INTT must recover original coefficients",
        )


fn _assert_xyw_roundtrip(
    n_value: Int, p_value: Int, modulus_value: UInt32
) raises:
    var original_values = _build_xyw_values(n_value, p_value, modulus_value)
    var transformed_values = _copy_list_u32(original_values)
    var ntt_success = apply_xyw_quotient_ntt_inplace(
        transformed_values, n_value, p_value, modulus_value
    )
    assert_true(ntt_success, "XYW NTT should succeed")
    var intt_success = apply_xyw_quotient_intt_inplace(
        transformed_values, n_value, p_value, modulus_value
    )
    assert_true(intt_success, "XYW inverse NTT should succeed")
    for coeff_idx in range(len(original_values)):
        assert_true(
            transformed_values[coeff_idx] == original_values[coeff_idx],
            "XYW NTT->INTT must recover original coefficients",
        )


fn _assert_ixw_multiplication_matches_naive(
    n_value: Int, p_value: Int, modulus_value: UInt32
) raises:
    var left_values = _build_ixw_values(n_value, p_value, modulus_value)
    var right_values = _build_ixw_values(n_value, p_value, modulus_value)
    for coeff_idx in range(len(right_values)):
        right_values[coeff_idx] = UInt32(
            _normalize_mod(
                Int(right_values[coeff_idx]) + 7 * coeff_idx + 11,
                Int(modulus_value),
            )
        )
    var expected_values = _naive_ixw_multiply_mod_quotient(
        left_values, right_values, n_value, p_value, modulus_value
    )

    var left_transformed = _copy_list_u32(left_values)
    var right_transformed = _copy_list_u32(right_values)
    assert_true(
        apply_ixw_quotient_ntt_inplace(
            left_transformed, n_value, p_value, modulus_value
        ),
        "IXW forward NTT for left input should succeed",
    )
    assert_true(
        apply_ixw_quotient_ntt_inplace(
            right_transformed, n_value, p_value, modulus_value
        ),
        "IXW forward NTT for right input should succeed",
    )

    var pointwise_product = List[UInt32](length=len(left_transformed), fill=0)
    for coeff_idx in range(len(pointwise_product)):
        var multiplied_value = (
            UInt64(left_transformed[coeff_idx])
            * UInt64(right_transformed[coeff_idx])
        ) % UInt64(modulus_value)
        pointwise_product[coeff_idx] = UInt32(multiplied_value)

    assert_true(
        apply_ixw_quotient_intt_inplace(
            pointwise_product, n_value, p_value, modulus_value
        ),
        "IXW inverse NTT after pointwise multiply should succeed",
    )
    for coeff_idx in range(len(pointwise_product)):
        assert_true(
            pointwise_product[coeff_idx] == expected_values[coeff_idx],
            "IXW NTT multiplication must match naive quotient multiplication",
        )


fn _assert_xyw_multiplication_matches_naive(
    n_value: Int, p_value: Int, modulus_value: UInt32
) raises:
    var left_values = _build_xyw_values(n_value, p_value, modulus_value)
    var right_values = _build_xyw_values(n_value, p_value, modulus_value)
    for coeff_idx in range(len(right_values)):
        right_values[coeff_idx] = UInt32(
            _normalize_mod(
                Int(right_values[coeff_idx]) + 5 * coeff_idx + 17,
                Int(modulus_value),
            )
        )
    var expected_values = _naive_xyw_multiply_mod_quotient(
        left_values, right_values, n_value, p_value, modulus_value
    )

    var left_transformed = _copy_list_u32(left_values)
    var right_transformed = _copy_list_u32(right_values)
    assert_true(
        apply_xyw_quotient_ntt_inplace(
            left_transformed, n_value, p_value, modulus_value
        ),
        "XYW forward NTT for left input should succeed",
    )
    assert_true(
        apply_xyw_quotient_ntt_inplace(
            right_transformed, n_value, p_value, modulus_value
        ),
        "XYW forward NTT for right input should succeed",
    )

    var pointwise_product = List[UInt32](length=len(left_transformed), fill=0)
    for coeff_idx in range(len(pointwise_product)):
        var multiplied_value = (
            UInt64(left_transformed[coeff_idx])
            * UInt64(right_transformed[coeff_idx])
        ) % UInt64(modulus_value)
        pointwise_product[coeff_idx] = UInt32(multiplied_value)

    assert_true(
        apply_xyw_quotient_intt_inplace(
            pointwise_product, n_value, p_value, modulus_value
        ),
        "XYW inverse NTT after pointwise multiply should succeed",
    )
    for coeff_idx in range(len(pointwise_product)):
        assert_true(
            pointwise_product[coeff_idx] == expected_values[coeff_idx],
            "XYW NTT multiplication must match naive quotient multiplication",
        )


fn test_ixw_and_xyw_roundtrip_small_grids() raises:
    print("Running test_ixw_and_xyw_roundtrip_small_grids...")

    var n_values = List[Int]()
    n_values.append(4)
    var p_values = List[Int]()
    p_values.append(3)

    for n_value in n_values:
        for p_value in p_values:
            var mod_u64 = find_suitable_q(n_value, p_value, 17)
            var modulus_value = UInt32(mod_u64)
            _assert_ixw_roundtrip(n_value, p_value, modulus_value, False)
            _assert_ixw_roundtrip(n_value, p_value, modulus_value, True)
            _assert_xyw_roundtrip(n_value, p_value, modulus_value)

    print("test_ixw_and_xyw_roundtrip_small_grids passed.")


fn manual_check_ntt_multiplication_matches_naive_2d_3d_4d() raises:
    print("Running manual_check_ntt_multiplication_matches_naive_2d_3d_4d...")

    var p_values = List[Int]()
    p_values.append(3)
    p_values.append(9)
    var n_values_xyw = List[Int]()
    n_values_xyw.append(4)  # i/x/y/w path with n>=4 production setting

    for p_value in p_values:
        for n_value in n_values_xyw:
            var mod_u64 = find_suitable_q(n_value, p_value, 17)
            var modulus_value = UInt32(mod_u64)
            _assert_xyw_multiplication_matches_naive(
                n_value, p_value, modulus_value
            )

    print("manual_check_ntt_multiplication_matches_naive_2d_3d_4d passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
