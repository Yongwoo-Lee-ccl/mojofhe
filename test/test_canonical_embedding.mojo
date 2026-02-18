from testing import TestSuite, assert_true
from src.canonical_embedding import (
    ComplexTensor3D,
    CanonicalEmbeddingContext,
    encode_canonical_complex,
    decode_canonical_complex,
)

struct ComplexPair(Copyable, Movable):
    var real_part: Float64
    var imag_part: Float64

    fn __init__(out self, real_part: Float64 = 0.0, imag_part: Float64 = 0.0):
        self.real_part = real_part
        self.imag_part = imag_part


struct ComplexPolynomial(Movable):
    var real_values: List[Float64]
    var imag_values: List[Float64]

    fn __init__(out self, degree: Int):
        self.real_values = List[Float64](length=degree, fill=0.0)
        self.imag_values = List[Float64](length=degree, fill=0.0)


fn _abs_f64(value: Float64) -> Float64:
    if value < 0.0:
        return -value
    return value


fn _complex_close(
    left_real: Float64,
    left_imag: Float64,
    right_real: Float64,
    right_imag: Float64,
    tolerance: Float64,
) -> Bool:
    return (
        _abs_f64(left_real - right_real) <= tolerance
        and _abs_f64(left_imag - right_imag) <= tolerance
    )


fn _make_input_tensor(n_dim: Int, phi_degree: Int) -> ComplexTensor3D:
    var tensor_value = ComplexTensor3D(n_dim, n_dim, phi_degree)
    for row_index in range(n_dim):
        for column_index in range(n_dim):
            for batch_index in range(phi_degree):
                var real_seed = Float64(
                    row_index * 11 + column_index * 7 + batch_index * 3
                )
                var imag_seed = Float64(
                    row_index * 5 - column_index * 2 + batch_index * 13
                )
                tensor_value.set(
                    row_index,
                    column_index,
                    batch_index,
                    (real_seed % 17.0) / 8.0,
                    (imag_seed % 19.0) / 9.0,
                )
    return tensor_value^


fn _make_input_tensor_with_offset(
    n_dim: Int,
    phi_degree: Int,
    offset: Int,
) -> ComplexTensor3D:
    var tensor_value = ComplexTensor3D(n_dim, n_dim, phi_degree)
    for row_index in range(n_dim):
        for column_index in range(n_dim):
            for batch_index in range(phi_degree):
                var real_seed = Float64(
                    row_index * 11 + column_index * 7 + batch_index * 3 + offset
                )
                var imag_seed = Float64(
                    row_index * 5
                    - column_index * 2
                    + batch_index * 13
                    + 2 * offset
                )
                tensor_value.set(
                    row_index,
                    column_index,
                    batch_index,
                    (real_seed % 23.0) / 10.0,
                    (imag_seed % 29.0) / 11.0,
                )
    return tensor_value^


fn _complex_multiply(
    left_real: Float64,
    left_imag: Float64,
    right_real: Float64,
    right_imag: Float64,
) -> ComplexPair:
    return ComplexPair(
        left_real * right_real - left_imag * right_imag,
        left_real * right_imag + left_imag * right_real,
    )


fn _complex_conjugate(real_part: Float64, imag_part: Float64) -> ComplexPair:
    return ComplexPair(real_part, -imag_part)


fn _w_inverse_polynomial_complex(
    coefficient_real_values: List[Float64],
    coefficient_imag_values: List[Float64],
    p_value: Int,
) -> ComplexPolynomial:
    var phi_degree = len(coefficient_real_values)
    var reduction_shift = p_value // 3

    var exponent_real_values = List[Float64](length=p_value, fill=0.0)
    var exponent_imag_values = List[Float64](length=p_value, fill=0.0)

    for coefficient_index in range(phi_degree):
        var mapped_degree = 0
        if coefficient_index != 0:
            mapped_degree = p_value - coefficient_index
        exponent_real_values[mapped_degree] += coefficient_real_values[
            coefficient_index
        ]
        exponent_imag_values[mapped_degree] += coefficient_imag_values[
            coefficient_index
        ]

    var reduction_index = p_value - 1
    while reduction_index >= phi_degree:
        var top_real = exponent_real_values[reduction_index]
        var top_imag = exponent_imag_values[reduction_index]
        if _abs_f64(top_real) > 0.0 or _abs_f64(top_imag) > 0.0:
            exponent_real_values[reduction_index] = 0.0
            exponent_imag_values[reduction_index] = 0.0
            var first_reduction_index = reduction_index - reduction_shift
            var second_reduction_index = reduction_index - phi_degree

            exponent_real_values[first_reduction_index] -= top_real
            exponent_imag_values[first_reduction_index] -= top_imag
            exponent_real_values[second_reduction_index] -= top_real
            exponent_imag_values[second_reduction_index] -= top_imag
        reduction_index -= 1

    var output_values = ComplexPolynomial(phi_degree)
    for coefficient_index in range(phi_degree):
        output_values.real_values[coefficient_index] = exponent_real_values[
            coefficient_index
        ]
        output_values.imag_values[coefficient_index] = exponent_imag_values[
            coefficient_index
        ]
    return output_values^


fn _w_multiply_mod_cyclotomic_complex(
    left_real_values: List[Float64],
    left_imag_values: List[Float64],
    right_real_values: List[Float64],
    right_imag_values: List[Float64],
    p_value: Int,
) -> ComplexPolynomial:
    var phi_degree = len(left_real_values)
    var reduction_shift = p_value // 3
    var product_degree_bound = 2 * phi_degree - 1

    var unreduced_real_values = List[Float64](
        length=product_degree_bound,
        fill=0.0,
    )
    var unreduced_imag_values = List[Float64](
        length=product_degree_bound,
        fill=0.0,
    )

    for left_index in range(phi_degree):
        for right_index in range(phi_degree):
            var product_value = _complex_multiply(
                left_real_values[left_index],
                left_imag_values[left_index],
                right_real_values[right_index],
                right_imag_values[right_index],
            )
            var product_index = left_index + right_index
            unreduced_real_values[product_index] += product_value.real_part
            unreduced_imag_values[product_index] += product_value.imag_part

    var reduction_index = product_degree_bound - 1
    while reduction_index >= phi_degree:
        var top_real = unreduced_real_values[reduction_index]
        var top_imag = unreduced_imag_values[reduction_index]
        if _abs_f64(top_real) > 0.0 or _abs_f64(top_imag) > 0.0:
            unreduced_real_values[reduction_index] = 0.0
            unreduced_imag_values[reduction_index] = 0.0
            var first_reduction_index = reduction_index - reduction_shift
            var second_reduction_index = reduction_index - phi_degree

            unreduced_real_values[first_reduction_index] -= top_real
            unreduced_imag_values[first_reduction_index] -= top_imag
            unreduced_real_values[second_reduction_index] -= top_real
            unreduced_imag_values[second_reduction_index] -= top_imag
        reduction_index -= 1

    var output_values = ComplexPolynomial(phi_degree)
    for coefficient_index in range(phi_degree):
        output_values.real_values[coefficient_index] = unreduced_real_values[
            coefficient_index
        ]
        output_values.imag_values[coefficient_index] = unreduced_imag_values[
            coefficient_index
        ]
    return output_values^


fn _rhs_automorphism_for_trace(
    encoded_tensor: ComplexTensor3D,
    n_dim: Int,
    p_value: Int,
    phi_degree: Int,
) -> ComplexTensor3D:
    var output_tensor = ComplexTensor3D(n_dim, n_dim, phi_degree)

    for x_index in range(n_dim):
        var mapped_x_index = (n_dim - x_index) % n_dim
        for y_index in range(n_dim):
            var temp_real_values = List[Float64](length=phi_degree, fill=0.0)
            var temp_imag_values = List[Float64](length=phi_degree, fill=0.0)
            for w_index in range(phi_degree):
                var conjugated_value = _complex_conjugate(
                    encoded_tensor.get_real(x_index, y_index, w_index),
                    encoded_tensor.get_imag(x_index, y_index, w_index),
                )
                var mapped_real = conjugated_value.real_part
                var mapped_imag = conjugated_value.imag_part
                if x_index != 0:
                    # Multiply by -i for X^{-x} under X^n = i.
                    var multiplied = _complex_multiply(
                        mapped_real,
                        mapped_imag,
                        0.0,
                        -1.0,
                    )
                    mapped_real = multiplied.real_part
                    mapped_imag = multiplied.imag_part
                temp_real_values[w_index] = mapped_real
                temp_imag_values[w_index] = mapped_imag

            var mapped_w_values = _w_inverse_polynomial_complex(
                temp_real_values,
                temp_imag_values,
                p_value,
            )
            for w_index in range(phi_degree):
                output_tensor.set(
                    mapped_x_index,
                    y_index,
                    w_index,
                    mapped_w_values.real_values[w_index],
                    mapped_w_values.imag_values[w_index],
                )

    return output_tensor^


fn _trace_multiply_encoded_polynomials(
    encoded_left: ComplexTensor3D,
    encoded_right: ComplexTensor3D,
    n_dim: Int,
    p_value: Int,
    phi_degree: Int,
) -> ComplexTensor3D:
    var rhs_prime = _rhs_automorphism_for_trace(
        encoded_right,
        n_dim,
        p_value,
        phi_degree,
    )
    var output_tensor = ComplexTensor3D(n_dim, n_dim, phi_degree)

    for row_index in range(n_dim):
        for col_index in range(n_dim):
            var acc_real_values = List[Float64](length=phi_degree, fill=0.0)
            var acc_imag_values = List[Float64](length=phi_degree, fill=0.0)
            for inner_index in range(n_dim):
                var left_real_values = List[Float64](length=phi_degree, fill=0.0)
                var left_imag_values = List[Float64](length=phi_degree, fill=0.0)
                var right_real_values = List[Float64](length=phi_degree, fill=0.0)
                var right_imag_values = List[Float64](length=phi_degree, fill=0.0)
                for w_index in range(phi_degree):
                    left_real_values[w_index] = encoded_left.get_real(
                        row_index,
                        inner_index,
                        w_index,
                    )
                    left_imag_values[w_index] = encoded_left.get_imag(
                        row_index,
                        inner_index,
                        w_index,
                    )
                    right_real_values[w_index] = rhs_prime.get_real(
                        col_index,
                        inner_index,
                        w_index,
                    )
                    right_imag_values[w_index] = rhs_prime.get_imag(
                        col_index,
                        inner_index,
                        w_index,
                    )

                var product_values = _w_multiply_mod_cyclotomic_complex(
                    left_real_values,
                    left_imag_values,
                    right_real_values,
                    right_imag_values,
                    p_value,
                )
                for w_index in range(phi_degree):
                    acc_real_values[w_index] += product_values.real_values[
                        w_index
                    ]
                    acc_imag_values[w_index] += product_values.imag_values[
                        w_index
                    ]

            for w_index in range(phi_degree):
                output_tensor.set(
                    row_index,
                    col_index,
                    w_index,
                    acc_real_values[w_index],
                    acc_imag_values[w_index],
                )
    return output_tensor^


fn test_canonical_embedding_roundtrip_without_rounding() raises:
    print("Running test_canonical_embedding_roundtrip_without_rounding...")

    var n_dim = 4
    var p_value = 9
    var context = CanonicalEmbeddingContext(n_dim, p_value)
    var input_tensor = _make_input_tensor(n_dim, context.phi_p_degree)

    var encoded_tensor = encode_canonical_complex(
        input_tensor,
        context,
        1.0,
        round_to_gaussian_integer=False,
    )
    var decoded_tensor = decode_canonical_complex(encoded_tensor, context, 1.0)

    var tolerance = 1e-7
    for row_index in range(n_dim):
        for column_index in range(n_dim):
            for batch_index in range(context.phi_p_degree):
                assert_true(
                    _complex_close(
                        input_tensor.get_real(
                            row_index, column_index, batch_index
                        ),
                        input_tensor.get_imag(
                            row_index, column_index, batch_index
                        ),
                        decoded_tensor.get_real(
                            row_index, column_index, batch_index
                        ),
                        decoded_tensor.get_imag(
                            row_index, column_index, batch_index
                        ),
                        tolerance,
                    ),
                    "decode(encode(input)) should recover input without rounding",
                )

    print("test_canonical_embedding_roundtrip_without_rounding passed.")


fn test_canonical_embedding_roundtrip_with_scaling_no_rounding() raises:
    print("Running test_canonical_embedding_roundtrip_with_scaling_no_rounding...")

    var n_dim = 4
    var p_value = 9
    var scaling_factor = 64.0
    var context = CanonicalEmbeddingContext(n_dim, p_value)
    var input_tensor = _make_input_tensor(n_dim, context.phi_p_degree)

    var encoded_tensor = encode_canonical_complex(
        input_tensor,
        context,
        scaling_factor,
        round_to_gaussian_integer=False,
    )
    var decoded_tensor = decode_canonical_complex(
        encoded_tensor,
        context,
        scaling_factor,
    )

    var tolerance = 1e-6
    for row_index in range(n_dim):
        for column_index in range(n_dim):
            for batch_index in range(context.phi_p_degree):
                assert_true(
                    _complex_close(
                        input_tensor.get_real(
                            row_index, column_index, batch_index
                        ),
                        input_tensor.get_imag(
                            row_index, column_index, batch_index
                        ),
                        decoded_tensor.get_real(
                            row_index, column_index, batch_index
                        ),
                        decoded_tensor.get_imag(
                            row_index, column_index, batch_index
                        ),
                        tolerance,
                    ),
                    "scaling and inverse scaling should preserve values",
                )

    print("test_canonical_embedding_roundtrip_with_scaling_no_rounding passed.")


fn test_trace_matmul_matches_naive_complex_matmul() raises:
    print("Running test_trace_matmul_matches_naive_complex_matmul...")

    var n_dim = 4
    var p_value = 9
    var context = CanonicalEmbeddingContext(n_dim, p_value)

    var matrix_a = _make_input_tensor_with_offset(
        n_dim,
        context.phi_p_degree,
        2,
    )
    var matrix_b = _make_input_tensor_with_offset(
        n_dim,
        context.phi_p_degree,
        9,
    )

    var encoded_a = encode_canonical_complex(
        matrix_a,
        context,
        1.0,
        round_to_gaussian_integer=False,
    )
    var encoded_b = encode_canonical_complex(
        matrix_b,
        context,
        1.0,
        round_to_gaussian_integer=False,
    )

    # encode -> trace-matmul -> decode matches naive A * B* in C.
    var trace_product_encoded = _trace_multiply_encoded_polynomials(
        encoded_a,
        encoded_b,
        n_dim,
        p_value,
        context.phi_p_degree,
    )
    for coefficient_index in range(len(trace_product_encoded.real_values)):
        trace_product_encoded.real_values[coefficient_index] *= Float64(n_dim)
        trace_product_encoded.imag_values[coefficient_index] *= Float64(n_dim)

    var decoded_trace_product = decode_canonical_complex(
        trace_product_encoded,
        context,
        1.0,
    )

    var product_tolerance = 2e-5
    for row_index in range(n_dim):
        for col_index in range(n_dim):
            for batch_index in range(context.phi_p_degree):
                var expected_real = 0.0
                var expected_imag = 0.0
                for inner_index in range(n_dim):
                    var right_conjugated = _complex_conjugate(
                        matrix_b.get_real(col_index, inner_index, batch_index),
                        matrix_b.get_imag(col_index, inner_index, batch_index),
                    )
                    var product_value = _complex_multiply(
                        matrix_a.get_real(row_index, inner_index, batch_index),
                        matrix_a.get_imag(row_index, inner_index, batch_index),
                        right_conjugated.real_part,
                        right_conjugated.imag_part,
                    )
                    expected_real += product_value.real_part
                    expected_imag += product_value.imag_part

                assert_true(
                    _complex_close(
                        decoded_trace_product.get_real(
                            row_index,
                            col_index,
                            batch_index,
                        ),
                        decoded_trace_product.get_imag(
                            row_index,
                            col_index,
                            batch_index,
                        ),
                        expected_real,
                        expected_imag,
                        product_tolerance,
                    ),
                    "trace-based encoded matmul should match naive complex matmul",
                )

    print("test_trace_matmul_matches_naive_complex_matmul passed.")


fn test_encode_decode_roundtrip_for_matmul_inputs() raises:
    print("Running test_encode_decode_roundtrip_for_matmul_inputs...")

    var n_dim = 4
    var p_value = 9
    var context = CanonicalEmbeddingContext(n_dim, p_value)

    var matrix_a = _make_input_tensor_with_offset(
        n_dim,
        context.phi_p_degree,
        2,
    )
    var matrix_b = _make_input_tensor_with_offset(
        n_dim,
        context.phi_p_degree,
        9,
    )

    var encoded_a = encode_canonical_complex(
        matrix_a,
        context,
        1.0,
        round_to_gaussian_integer=False,
    )
    var encoded_b = encode_canonical_complex(
        matrix_b,
        context,
        1.0,
        round_to_gaussian_integer=False,
    )
    var decoded_a = decode_canonical_complex(encoded_a, context, 1.0)
    var decoded_b = decode_canonical_complex(encoded_b, context, 1.0)

    var roundtrip_tolerance = 1e-6
    for row_index in range(n_dim):
        for col_index in range(n_dim):
            for batch_index in range(context.phi_p_degree):
                assert_true(
                    _complex_close(
                        matrix_a.get_real(row_index, col_index, batch_index),
                        matrix_a.get_imag(row_index, col_index, batch_index),
                        decoded_a.get_real(row_index, col_index, batch_index),
                        decoded_a.get_imag(row_index, col_index, batch_index),
                        roundtrip_tolerance,
                    ),
                    "encode/decode should preserve matrix A",
                )
                assert_true(
                    _complex_close(
                        matrix_b.get_real(row_index, col_index, batch_index),
                        matrix_b.get_imag(row_index, col_index, batch_index),
                        decoded_b.get_real(row_index, col_index, batch_index),
                        decoded_b.get_imag(row_index, col_index, batch_index),
                        roundtrip_tolerance,
                    ),
                    "encode/decode should preserve matrix B",
                )

    print("test_encode_decode_roundtrip_for_matmul_inputs passed.")


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
