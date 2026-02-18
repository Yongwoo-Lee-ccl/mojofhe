from collections import List
from testing import TestSuite, assert_true
from src.modular import (
    compute_barrett_ratio,
    find_suitable_q,
    multiply_mod_barrett,
)
from src.encoded_matrix import EncodedPolynomial


fn _add_modulus(
    left_value: UInt32, right_value: UInt32, modulus_value: UInt32
) -> UInt32:
    var summed_value = left_value + right_value
    if summed_value >= modulus_value:
        summed_value -= modulus_value
    return summed_value


fn _subtract_modulus(
    left_value: UInt32, right_value: UInt32, modulus_value: UInt32
) -> UInt32:
    var difference_value = Int64(left_value) - Int64(right_value)
    var reduced_value = difference_value % Int64(modulus_value)
    if reduced_value < 0:
        reduced_value += Int64(modulus_value)
    return UInt32(reduced_value)


fn _w_inverse_reference(
    coefficient_values: List[UInt32],
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    var phi_degree = 2 * (p_power_of_3 // 3)
    var reduction_shift = p_power_of_3 // 3

    var exponent_space_values = List[UInt32](length=p_power_of_3, fill=0)
    for coefficient_index in range(phi_degree):
        var coefficient_value = coefficient_values[coefficient_index]
        var mapped_degree = 0
        if coefficient_index != 0:
            mapped_degree = p_power_of_3 - coefficient_index
        exponent_space_values[mapped_degree] = _add_modulus(
            exponent_space_values[mapped_degree],
            coefficient_value,
            modulus_value,
        )

    var reduction_index = p_power_of_3 - 1
    while reduction_index >= phi_degree:
        var top_coefficient = exponent_space_values[reduction_index]
        if top_coefficient != 0:
            exponent_space_values[reduction_index] = 0
            var first_reduction_index = reduction_index - reduction_shift
            var second_reduction_index = reduction_index - phi_degree
            exponent_space_values[first_reduction_index] = _subtract_modulus(
                exponent_space_values[first_reduction_index],
                top_coefficient,
                modulus_value,
            )
            exponent_space_values[second_reduction_index] = _subtract_modulus(
                exponent_space_values[second_reduction_index],
                top_coefficient,
                modulus_value,
            )
        reduction_index -= 1

    var reduced_values = List[UInt32](length=phi_degree, fill=0)
    for coefficient_index in range(phi_degree):
        reduced_values[coefficient_index] = exponent_space_values[
            coefficient_index
        ]
    return reduced_values^


fn _w_multiply_reference(
    left_values: List[UInt32],
    right_values: List[UInt32],
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    var phi_degree = 2 * (p_power_of_3 // 3)
    var reduction_shift = p_power_of_3 // 3
    var barrett_ratio = compute_barrett_ratio(modulus_value)

    var product_degree_bound = 2 * phi_degree - 1
    var unreduced_values = List[UInt32](length=product_degree_bound, fill=0)
    for left_index in range(phi_degree):
        for right_index in range(phi_degree):
            var product_value = multiply_mod_barrett(
                left_values[left_index],
                right_values[right_index],
                modulus_value,
                barrett_ratio,
            )
            var product_index = left_index + right_index
            unreduced_values[product_index] = _add_modulus(
                unreduced_values[product_index],
                product_value,
                modulus_value,
            )

    var reduction_index = product_degree_bound - 1
    while reduction_index >= phi_degree:
        var top_coefficient = unreduced_values[reduction_index]
        if top_coefficient != 0:
            unreduced_values[reduction_index] = 0
            var first_reduction_index = reduction_index - reduction_shift
            var second_reduction_index = reduction_index - phi_degree
            unreduced_values[first_reduction_index] = _subtract_modulus(
                unreduced_values[first_reduction_index],
                top_coefficient,
                modulus_value,
            )
            unreduced_values[second_reduction_index] = _subtract_modulus(
                unreduced_values[second_reduction_index],
                top_coefficient,
                modulus_value,
            )
        reduction_index -= 1

    var reduced_values = List[UInt32](length=phi_degree, fill=0)
    for coefficient_index in range(phi_degree):
        reduced_values[coefficient_index] = unreduced_values[coefficient_index]
    return reduced_values^


fn _poly_equal(
    left_value: EncodedPolynomial, right_value: EncodedPolynomial
) -> Bool:
    if len(left_value.coefficient_values) != len(
        right_value.coefficient_values
    ):
        return False
    for coefficient_index in range(len(left_value.coefficient_values)):
        if (
            left_value.coefficient_values[coefficient_index]
            != right_value.coefficient_values[coefficient_index]
        ):
            return False
    return True


fn test_rhs_automorphism_involution() raises:
    print("Running test_rhs_automorphism_involution...")

    var n_dim = 4
    var p_power_of_3 = 9
    var q_modulus = UInt32(find_suitable_q(n_dim, p_power_of_3, 20))
    var poly_value = EncodedPolynomial(n_dim, p_power_of_3, q_modulus)

    var coefficient_seed = Int(q_modulus) - 19
    for x_index in range(n_dim):
        for y_index in range(n_dim):
            for w_index in range(poly_value.phi_p_degree):
                coefficient_seed = (
                    coefficient_seed * 73 + 17 + x_index + 3 * y_index
                ) % Int(q_modulus)
                poly_value.set_coefficient(
                    x_index,
                    y_index,
                    w_index,
                    coefficient_seed,
                )

    var transformed_once = poly_value.apply_rhs_automorphism()
    var transformed_twice = transformed_once.apply_rhs_automorphism()

    assert_true(
        _poly_equal(poly_value, transformed_twice),
        "b(X,Y,W) should equal b((X^{-1})^{-1},Y,(W^{-1})^{-1})",
    )

    print("test_rhs_automorphism_involution passed.")


fn test_batched_trace_multiply_matches_direct_formula() raises:
    print("Running test_batched_trace_multiply_matches_direct_formula...")

    var n_dim = 4
    var p_power_of_3 = 9
    var q_modulus = UInt32(find_suitable_q(n_dim, p_power_of_3, 20))

    var left_value = EncodedPolynomial(n_dim, p_power_of_3, q_modulus)
    var right_value = EncodedPolynomial(n_dim, p_power_of_3, q_modulus)
    var left_reference = EncodedPolynomial(n_dim, p_power_of_3, q_modulus)
    var right_reference = EncodedPolynomial(n_dim, p_power_of_3, q_modulus)

    var left_seed = Int(q_modulus) - 31
    var right_seed = Int(q_modulus) - 47
    for x_index in range(n_dim):
        for y_index in range(n_dim):
            for w_index in range(left_value.phi_p_degree):
                left_seed = (left_seed * 37 + 11 + 5 * x_index + y_index) % Int(
                    q_modulus
                )
                right_seed = (
                    right_seed * 41 + 13 + x_index + 7 * y_index
                ) % Int(q_modulus)
                left_value.set_coefficient(
                    x_index,
                    y_index,
                    w_index,
                    left_seed,
                )
                left_reference.set_coefficient(
                    x_index,
                    y_index,
                    w_index,
                    left_seed,
                )
                right_value.set_coefficient(
                    x_index,
                    y_index,
                    w_index,
                    right_seed,
                )
                right_reference.set_coefficient(
                    x_index,
                    y_index,
                    w_index,
                    right_seed,
                )

    var computed_value = left_value.trace_multiply(right_value)
    assert_true(
        computed_value.is_w_ntt,
        "trace_multiply() should output W-NTT-domain polynomial",
    )
    var expected_value = EncodedPolynomial(n_dim, p_power_of_3, q_modulus)

    for row_index in range(n_dim):
        for column_index in range(n_dim):
            var accumulated_values = List[UInt32](
                length=left_value.phi_p_degree,
                fill=0,
            )
            for inner_index in range(n_dim):
                var left_coefficients = left_reference.get_w_polynomial(
                    row_index,
                    inner_index,
                )

                # Right side uses Y^{-1} and Z^{-1}; in coefficient form
                # this maps the first index by negation modulo n.
                var source_x_index = (n_dim - column_index) % n_dim
                var right_source_coefficients = (
                    right_reference.get_w_polynomial(
                        source_x_index,
                        inner_index,
                    )
                )
                var right_mapped_coefficients = _w_inverse_reference(
                    right_source_coefficients,
                    p_power_of_3,
                    q_modulus,
                )

                var product_values = _w_multiply_reference(
                    left_coefficients,
                    right_mapped_coefficients,
                    p_power_of_3,
                    q_modulus,
                )
                for coefficient_index in range(left_value.phi_p_degree):
                    accumulated_values[coefficient_index] = _add_modulus(
                        accumulated_values[coefficient_index],
                        product_values[coefficient_index],
                        q_modulus,
                    )
            expected_value.set_w_polynomial(
                row_index,
                column_index,
                accumulated_values,
            )

    assert_true(
        expected_value.is_w_ntt == False,
        "reference should remain in coefficient domain before transform",
    )
    expected_value.transform_w_to_ntt()

    assert_true(
        expected_value.is_w_ntt,
        "reference should be transformed to W-NTT before comparison",
    )

    assert_true(
        _poly_equal(computed_value, expected_value),
        "trace_multiply() should match theorem formula in W-NTT domain",
    )

    print("test_batched_trace_multiply_matches_direct_formula passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
