from testing import TestSuite, assert_true
from src.canonical_embedding import (
    ComplexTensor3D,
    CanonicalEmbeddingContext,
    encode_canonical_complex,
    decode_canonical_complex,
)


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


fn main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
