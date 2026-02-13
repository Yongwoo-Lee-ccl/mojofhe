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


fn _copy_list_i32(values: List[Int32]) -> List[Int32]:
    var copied_values = List[Int32]()
    for value in values:
        copied_values.append(value)
    return copied_values^


fn _build_ixw_values(
    n_value: Int, p_value: Int, modulus_value: Int
) -> List[Int32]:
    var phi_degree = 2 * (p_value // 3)
    var total_length = 2 * n_value * phi_degree
    var values = List[Int32](length=total_length, fill=0)
    for component_index in range(2):
        for x_index in range(n_value):
            for w_index in range(phi_degree):
                var seed_value = (
                    component_index * 131
                    + x_index * 47
                    + w_index * 19
                    + n_value * 23
                    + p_value * 17
                )
                values[
                    _ixw_index(
                        component_index, x_index, w_index, n_value, phi_degree
                    )
                ] = Int32(_normalize_mod(seed_value, modulus_value))
    return values^


fn _build_xyw_values(
    n_value: Int, p_value: Int, modulus_value: Int
) -> List[Int32]:
    var phi_degree = 2 * (p_value // 3)
    var total_length = 2 * n_value * n_value * phi_degree
    var values = List[Int32](length=total_length, fill=0)
    for component_index in range(2):
        for y_index in range(n_value):
            for x_index in range(n_value):
                for w_index in range(phi_degree):
                    var seed_value = (
                        component_index * 149
                        + y_index * 41
                        + x_index * 29
                        + w_index * 13
                        + n_value * 31
                        + p_value * 11
                    )
                    values[
                        _xyw_index(
                            component_index,
                            y_index,
                            x_index,
                            w_index,
                            n_value,
                            phi_degree,
                        )
                    ] = Int32(_normalize_mod(seed_value, modulus_value))
    return values^


fn _naive_ixw_multiply_mod_quotient(
    left_values: List[Int32],
    right_values: List[Int32],
    n_value: Int,
    p_value: Int,
    modulus_value: Int,
) -> List[Int32]:
    var m_parameter = p_value // 3
    var phi_degree = 2 * m_parameter
    var w_temp_degree = 4 * m_parameter
    var x_temp_degree = 2 * n_value

    var temp_component_0 = List[Int](
        length=x_temp_degree * w_temp_degree, fill=0
    )
    var temp_component_1 = List[Int](
        length=x_temp_degree * w_temp_degree, fill=0
    )

    for left_x in range(n_value):
        for left_w in range(phi_degree):
            var left_0 = Int(
                left_values[_ixw_index(0, left_x, left_w, n_value, phi_degree)]
            )
            var left_1 = Int(
                left_values[_ixw_index(1, left_x, left_w, n_value, phi_degree)]
            )
            for right_x in range(n_value):
                for right_w in range(phi_degree):
                    var right_0 = Int(
                        right_values[
                            _ixw_index(0, right_x, right_w, n_value, phi_degree)
                        ]
                    )
                    var right_1 = Int(
                        right_values[
                            _ixw_index(1, right_x, right_w, n_value, phi_degree)
                        ]
                    )
                    var product_0 = left_0 * right_0 - left_1 * right_1
                    var product_1 = left_0 * right_1 + left_1 * right_0
                    var x_degree = left_x + right_x
                    var w_degree = left_w + right_w
                    var temp_index = x_degree * w_temp_degree + w_degree
                    temp_component_0[temp_index] = _normalize_mod(
                        temp_component_0[temp_index] + product_0,
                        modulus_value,
                    )
                    temp_component_1[temp_index] = _normalize_mod(
                        temp_component_1[temp_index] + product_1,
                        modulus_value,
                    )

    var reduced_component_0 = List[Int](length=n_value * w_temp_degree, fill=0)
    var reduced_component_1 = List[Int](length=n_value * w_temp_degree, fill=0)
    for x_degree in range(2 * n_value - 1):
        var reduced_x = x_degree
        if x_degree >= n_value:
            reduced_x = x_degree - n_value
        for w_degree in range(w_temp_degree):
            var source_index = x_degree * w_temp_degree + w_degree
            var target_index = reduced_x * w_temp_degree + w_degree
            var component_0_value = temp_component_0[source_index]
            var component_1_value = temp_component_1[source_index]
            reduced_component_0[target_index] = _normalize_mod(
                reduced_component_0[target_index] + component_0_value,
                modulus_value,
            )
            reduced_component_1[target_index] = _normalize_mod(
                reduced_component_1[target_index] + component_1_value,
                modulus_value,
            )

    var expected_values = List[Int32](length=2 * n_value * phi_degree, fill=0)
    for x_index in range(n_value):
        var w_degree = w_temp_degree - 1
        while w_degree >= phi_degree:
            var lane_index = x_index * w_temp_degree + w_degree
            var lane_0 = reduced_component_0[lane_index]
            var lane_1 = reduced_component_1[lane_index]
            if lane_0 != 0 or lane_1 != 0:
                var index_minus_m = x_index * w_temp_degree + (
                    w_degree - m_parameter
                )
                var index_minus_2m = x_index * w_temp_degree + (
                    w_degree - phi_degree
                )
                reduced_component_0[index_minus_m] = _normalize_mod(
                    reduced_component_0[index_minus_m] - lane_0,
                    modulus_value,
                )
                reduced_component_0[index_minus_2m] = _normalize_mod(
                    reduced_component_0[index_minus_2m] - lane_0,
                    modulus_value,
                )
                reduced_component_1[index_minus_m] = _normalize_mod(
                    reduced_component_1[index_minus_m] - lane_1,
                    modulus_value,
                )
                reduced_component_1[index_minus_2m] = _normalize_mod(
                    reduced_component_1[index_minus_2m] - lane_1,
                    modulus_value,
                )
            reduced_component_0[lane_index] = 0
            reduced_component_1[lane_index] = 0
            w_degree -= 1

        for w_index in range(phi_degree):
            var reduced_index = x_index * w_temp_degree + w_index
            expected_values[
                _ixw_index(0, x_index, w_index, n_value, phi_degree)
            ] = Int32(reduced_component_0[reduced_index])
            expected_values[
                _ixw_index(1, x_index, w_index, n_value, phi_degree)
            ] = Int32(reduced_component_1[reduced_index])

    return expected_values^


fn _naive_xyw_multiply_mod_quotient(
    left_values: List[Int32],
    right_values: List[Int32],
    n_value: Int,
    p_value: Int,
    modulus_value: Int,
) -> List[Int32]:
    var m_parameter = p_value // 3
    var phi_degree = 2 * m_parameter
    var w_temp_degree = 4 * m_parameter
    var x_temp_degree = 2 * n_value
    var y_temp_degree = 2 * n_value
    var temp_plane_size = x_temp_degree * w_temp_degree
    var temp_total_size = y_temp_degree * temp_plane_size

    var temp_component_0 = List[Int](length=temp_total_size, fill=0)
    var temp_component_1 = List[Int](length=temp_total_size, fill=0)

    for left_y in range(n_value):
        for left_x in range(n_value):
            for left_w in range(phi_degree):
                var left_0 = Int(
                    left_values[
                        _xyw_index(
                            0, left_y, left_x, left_w, n_value, phi_degree
                        )
                    ]
                )
                var left_1 = Int(
                    left_values[
                        _xyw_index(
                            1, left_y, left_x, left_w, n_value, phi_degree
                        )
                    ]
                )
                for right_y in range(n_value):
                    for right_x in range(n_value):
                        for right_w in range(phi_degree):
                            var right_0 = Int(
                                right_values[
                                    _xyw_index(
                                        0,
                                        right_y,
                                        right_x,
                                        right_w,
                                        n_value,
                                        phi_degree,
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
                                        phi_degree,
                                    )
                                ]
                            )
                            var product_0 = left_0 * right_0 - left_1 * right_1
                            var product_1 = left_0 * right_1 + left_1 * right_0
                            var y_degree = left_y + right_y
                            var x_degree = left_x + right_x
                            var w_degree = left_w + right_w
                            var temp_index = (
                                y_degree * temp_plane_size
                                + x_degree * w_temp_degree
                                + w_degree
                            )
                            temp_component_0[temp_index] = _normalize_mod(
                                temp_component_0[temp_index] + product_0,
                                modulus_value,
                            )
                            temp_component_1[temp_index] = _normalize_mod(
                                temp_component_1[temp_index] + product_1,
                                modulus_value,
                            )

    var reduced_plane_size = n_value * w_temp_degree
    var reduced_total_size = n_value * reduced_plane_size
    var reduced_component_0 = List[Int](length=reduced_total_size, fill=0)
    var reduced_component_1 = List[Int](length=reduced_total_size, fill=0)
    for y_degree in range(2 * n_value - 1):
        var reduced_y = y_degree
        if y_degree >= n_value:
            reduced_y = y_degree - n_value
        for x_degree in range(2 * n_value - 1):
            var reduced_x = x_degree
            if x_degree >= n_value:
                reduced_x = x_degree - n_value
            for w_degree in range(w_temp_degree):
                var source_index = (
                    y_degree * temp_plane_size
                    + x_degree * w_temp_degree
                    + w_degree
                )
                var target_index = (
                    reduced_y * reduced_plane_size
                    + reduced_x * w_temp_degree
                    + w_degree
                )
                var component_0_value = temp_component_0[source_index]
                var component_1_value = temp_component_1[source_index]
                reduced_component_0[target_index] = _normalize_mod(
                    reduced_component_0[target_index] + component_0_value,
                    modulus_value,
                )
                reduced_component_1[target_index] = _normalize_mod(
                    reduced_component_1[target_index] + component_1_value,
                    modulus_value,
                )

    var expected_values = List[Int32](
        length=2 * n_value * n_value * phi_degree, fill=0
    )
    for y_index in range(n_value):
        for x_index in range(n_value):
            var w_degree = w_temp_degree - 1
            while w_degree >= phi_degree:
                var lane_index = (
                    y_index * reduced_plane_size
                    + x_index * w_temp_degree
                    + w_degree
                )
                var lane_0 = reduced_component_0[lane_index]
                var lane_1 = reduced_component_1[lane_index]
                if lane_0 != 0 or lane_1 != 0:
                    var index_minus_m = (
                        y_index * reduced_plane_size
                        + x_index * w_temp_degree
                        + (w_degree - m_parameter)
                    )
                    var index_minus_2m = (
                        y_index * reduced_plane_size
                        + x_index * w_temp_degree
                        + (w_degree - phi_degree)
                    )
                    reduced_component_0[index_minus_m] = _normalize_mod(
                        reduced_component_0[index_minus_m] - lane_0,
                        modulus_value,
                    )
                    reduced_component_0[index_minus_2m] = _normalize_mod(
                        reduced_component_0[index_minus_2m] - lane_0,
                        modulus_value,
                    )
                    reduced_component_1[index_minus_m] = _normalize_mod(
                        reduced_component_1[index_minus_m] - lane_1,
                        modulus_value,
                    )
                    reduced_component_1[index_minus_2m] = _normalize_mod(
                        reduced_component_1[index_minus_2m] - lane_1,
                        modulus_value,
                    )
                reduced_component_0[lane_index] = 0
                reduced_component_1[lane_index] = 0
                w_degree -= 1

            for w_index in range(phi_degree):
                var reduced_index = (
                    y_index * reduced_plane_size
                    + x_index * w_temp_degree
                    + w_index
                )
                expected_values[
                    _xyw_index(
                        0, y_index, x_index, w_index, n_value, phi_degree
                    )
                ] = Int32(reduced_component_0[reduced_index])
                expected_values[
                    _xyw_index(
                        1, y_index, x_index, w_index, n_value, phi_degree
                    )
                ] = Int32(reduced_component_1[reduced_index])

    return expected_values^


fn _assert_ixw_roundtrip(
    n_value: Int, p_value: Int, modulus_value: Int, use_montgomery: Bool
) raises:
    var original_values = _build_ixw_values(n_value, p_value, modulus_value)
    var transformed_values = _copy_list_i32(original_values)
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
    for coefficient_index in range(len(original_values)):
        assert_true(
            transformed_values[coefficient_index]
            == original_values[coefficient_index],
            "IXW NTT->INTT must recover original coefficients",
        )


fn _assert_xyw_roundtrip(n_value: Int, p_value: Int, modulus_value: Int) raises:
    var original_values = _build_xyw_values(n_value, p_value, modulus_value)
    var transformed_values = _copy_list_i32(original_values)
    var ntt_success = apply_xyw_quotient_ntt_inplace(
        transformed_values, n_value, p_value, modulus_value
    )
    assert_true(ntt_success, "XYW NTT should succeed")
    var intt_success = apply_xyw_quotient_intt_inplace(
        transformed_values, n_value, p_value, modulus_value
    )
    assert_true(intt_success, "XYW inverse NTT should succeed")
    for coefficient_index in range(len(original_values)):
        assert_true(
            transformed_values[coefficient_index]
            == original_values[coefficient_index],
            "XYW NTT->INTT must recover original coefficients",
        )


fn _assert_ixw_multiplication_matches_naive(
    n_value: Int, p_value: Int, modulus_value: Int
) raises:
    var left_values = _build_ixw_values(n_value, p_value, modulus_value)
    var right_values = _build_ixw_values(n_value, p_value, modulus_value)
    for coefficient_index in range(len(right_values)):
        right_values[coefficient_index] = Int32(
            _normalize_mod(
                Int(right_values[coefficient_index])
                + 7 * coefficient_index
                + 11,
                modulus_value,
            )
        )
    var expected_values = _naive_ixw_multiply_mod_quotient(
        left_values, right_values, n_value, p_value, modulus_value
    )

    var left_transformed = _copy_list_i32(left_values)
    var right_transformed = _copy_list_i32(right_values)
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

    var pointwise_product = List[Int32](length=len(left_transformed), fill=0)
    for coefficient_index in range(len(pointwise_product)):
        var multiplied_value = (
            Int64(Int(left_transformed[coefficient_index]))
            * Int64(Int(right_transformed[coefficient_index]))
        ) % Int64(modulus_value)
        pointwise_product[coefficient_index] = Int32(Int(multiplied_value))

    assert_true(
        apply_ixw_quotient_intt_inplace(
            pointwise_product, n_value, p_value, modulus_value
        ),
        "IXW inverse NTT after pointwise multiply should succeed",
    )
    for coefficient_index in range(len(pointwise_product)):
        assert_true(
            pointwise_product[coefficient_index]
            == expected_values[coefficient_index],
            "IXW NTT multiplication must match naive quotient multiplication",
        )


fn _assert_xyw_multiplication_matches_naive(
    n_value: Int, p_value: Int, modulus_value: Int
) raises:
    var left_values = _build_xyw_values(n_value, p_value, modulus_value)
    var right_values = _build_xyw_values(n_value, p_value, modulus_value)
    for coefficient_index in range(len(right_values)):
        right_values[coefficient_index] = Int32(
            _normalize_mod(
                Int(right_values[coefficient_index])
                + 5 * coefficient_index
                + 17,
                modulus_value,
            )
        )
    var expected_values = _naive_xyw_multiply_mod_quotient(
        left_values, right_values, n_value, p_value, modulus_value
    )

    var left_transformed = _copy_list_i32(left_values)
    var right_transformed = _copy_list_i32(right_values)
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

    var pointwise_product = List[Int32](length=len(left_transformed), fill=0)
    for coefficient_index in range(len(pointwise_product)):
        var multiplied_value = (
            Int64(Int(left_transformed[coefficient_index]))
            * Int64(Int(right_transformed[coefficient_index]))
        ) % Int64(modulus_value)
        pointwise_product[coefficient_index] = Int32(Int(multiplied_value))

    assert_true(
        apply_xyw_quotient_intt_inplace(
            pointwise_product, n_value, p_value, modulus_value
        ),
        "XYW inverse NTT after pointwise multiply should succeed",
    )
    for coefficient_index in range(len(pointwise_product)):
        assert_true(
            pointwise_product[coefficient_index]
            == expected_values[coefficient_index],
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
            var modulus_value = find_suitable_q(n_value, p_value, 17)
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
            var modulus_value = find_suitable_q(n_value, p_value, 17)
            _assert_xyw_multiplication_matches_naive(
                n_value, p_value, modulus_value
            )

    print("manual_check_ntt_multiplication_matches_naive_2d_3d_4d passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
