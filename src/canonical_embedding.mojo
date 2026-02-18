from collections import List
from math import cos, sin, pi
from src.complex_utils import ComplexPair, complex_multiply_pair, complex_divide_pair


fn _round_nearest(value: Float64) -> Int:
    if value >= 0.0:
        return Int(value + 0.5)
    return Int(value - 0.5)


fn _pow_mod(base_value: Int, exponent_value: Int, modulus_value: Int) -> Int:
    var result_value = 1
    var base_reduced = base_value % modulus_value
    var exponent_remaining = exponent_value

    while exponent_remaining > 0:
        if exponent_remaining % 2 == 1:
            result_value = (result_value * base_reduced) % modulus_value
        base_reduced = (base_reduced * base_reduced) % modulus_value
        exponent_remaining //= 2
    return result_value


fn _gcd(value_a: Int, value_b: Int) -> Int:
    var a = value_a
    var b = value_b
    while b != 0:
        var temp = a % b
        a = b
        b = temp
    if a < 0:
        return -a
    return a


fn _complex_pow(
    base_real: Float64,
    base_imag: Float64,
    exponent_value: Int,
) -> ComplexPair:
    var result_real = 1.0
    var result_imag = 0.0
    if exponent_value == 0:
        return ComplexPair(result_real, result_imag)

    for exponent_index in range(exponent_value):
        _ = exponent_index
        var next_value = complex_multiply_pair(
            result_real,
            result_imag,
            base_real,
            base_imag,
        )
        result_real = next_value.real_part
        result_imag = next_value.imag_part
    return ComplexPair(result_real, result_imag)


fn _primitive_root_real_imag(order_value: Int) -> ComplexPair:
    var angle = 2.0 * pi / Float64(order_value)
    return ComplexPair(cos(angle), sin(angle))


struct ComplexTensor3D(Movable):
    var dim_x: Int
    var dim_y: Int
    var dim_z: Int
    var real_values: List[Float64]
    var imag_values: List[Float64]

    fn __init__(out self, dim_x: Int, dim_y: Int, dim_z: Int):
        self.dim_x = dim_x
        self.dim_y = dim_y
        self.dim_z = dim_z
        self.real_values = List[Float64](
            length=dim_x * dim_y * dim_z,
            fill=0.0,
        )
        self.imag_values = List[Float64](
            length=dim_x * dim_y * dim_z,
            fill=0.0,
        )

    fn _index(self, x_index: Int, y_index: Int, z_index: Int) -> Int:
        return x_index * self.dim_y * self.dim_z + y_index * self.dim_z + z_index

    fn set(
        mut self,
        x_index: Int,
        y_index: Int,
        z_index: Int,
        real_part: Float64,
        imag_part: Float64,
    ):
        var index_value = self._index(x_index, y_index, z_index)
        self.real_values[index_value] = real_part
        self.imag_values[index_value] = imag_part

    fn get_real(self, x_index: Int, y_index: Int, z_index: Int) -> Float64:
        return self.real_values[self._index(x_index, y_index, z_index)]

    fn get_imag(self, x_index: Int, y_index: Int, z_index: Int) -> Float64:
        return self.imag_values[self._index(x_index, y_index, z_index)]


struct ComplexMatrix(Movable):
    var rows: Int
    var cols: Int
    var real_values: List[Float64]
    var imag_values: List[Float64]

    fn __init__(out self, rows: Int, cols: Int):
        self.rows = rows
        self.cols = cols
        self.real_values = List[Float64](length=rows * cols, fill=0.0)
        self.imag_values = List[Float64](length=rows * cols, fill=0.0)

    fn _index(self, row_index: Int, col_index: Int) -> Int:
        return row_index * self.cols + col_index

    fn set(
        mut self,
        row_index: Int,
        col_index: Int,
        real_part: Float64,
        imag_part: Float64,
    ):
        var index_value = self._index(row_index, col_index)
        self.real_values[index_value] = real_part
        self.imag_values[index_value] = imag_part

    fn get_real(self, row_index: Int, col_index: Int) -> Float64:
        return self.real_values[self._index(row_index, col_index)]

    fn get_imag(self, row_index: Int, col_index: Int) -> Float64:
        return self.imag_values[self._index(row_index, col_index)]


fn _build_vandermonde(
    root_real_values: List[Float64],
    root_imag_values: List[Float64],
    degree: Int,
) -> ComplexMatrix:
    var output_matrix = ComplexMatrix(degree, degree)
    for row_index in range(degree):
        var current_real = 1.0
        var current_imag = 0.0
        for column_index in range(degree):
            output_matrix.set(
                row_index,
                column_index,
                current_real,
                current_imag,
            )
            var next_value = complex_multiply_pair(
                current_real,
                current_imag,
                root_real_values[row_index],
                root_imag_values[row_index],
            )
            current_real = next_value.real_part
            current_imag = next_value.imag_part
    return output_matrix^


fn _invert_square_matrix(input_matrix: ComplexMatrix, dimension: Int) -> ComplexMatrix:
    var left_matrix = ComplexMatrix(dimension, dimension)
    var right_matrix = ComplexMatrix(dimension, dimension)

    for row_index in range(dimension):
        for col_index in range(dimension):
            left_matrix.set(
                row_index,
                col_index,
                input_matrix.get_real(row_index, col_index),
                input_matrix.get_imag(row_index, col_index),
            )
            if row_index == col_index:
                right_matrix.set(row_index, col_index, 1.0, 0.0)

    for pivot_index in range(dimension):
        var pivot_row_index = pivot_index
        var best_norm = (
            left_matrix.get_real(pivot_row_index, pivot_index)
            * left_matrix.get_real(pivot_row_index, pivot_index)
            + left_matrix.get_imag(pivot_row_index, pivot_index)
            * left_matrix.get_imag(pivot_row_index, pivot_index)
        )
        for candidate_row in range(pivot_index + 1, dimension):
            var candidate_norm = (
                left_matrix.get_real(candidate_row, pivot_index)
                * left_matrix.get_real(candidate_row, pivot_index)
                + left_matrix.get_imag(candidate_row, pivot_index)
                * left_matrix.get_imag(candidate_row, pivot_index)
            )
            if candidate_norm > best_norm:
                best_norm = candidate_norm
                pivot_row_index = candidate_row

        if pivot_row_index != pivot_index:
            for column_index in range(dimension):
                var temp_left_real = left_matrix.get_real(
                    pivot_index, column_index
                )
                var temp_left_imag = left_matrix.get_imag(
                    pivot_index, column_index
                )
                left_matrix.set(
                    pivot_index,
                    column_index,
                    left_matrix.get_real(pivot_row_index, column_index),
                    left_matrix.get_imag(pivot_row_index, column_index),
                )
                left_matrix.set(
                    pivot_row_index,
                    column_index,
                    temp_left_real,
                    temp_left_imag,
                )

                var temp_right_real = right_matrix.get_real(
                    pivot_index, column_index
                )
                var temp_right_imag = right_matrix.get_imag(
                    pivot_index, column_index
                )
                right_matrix.set(
                    pivot_index,
                    column_index,
                    right_matrix.get_real(pivot_row_index, column_index),
                    right_matrix.get_imag(pivot_row_index, column_index),
                )
                right_matrix.set(
                    pivot_row_index,
                    column_index,
                    temp_right_real,
                    temp_right_imag,
                )

        var pivot_real = left_matrix.get_real(pivot_index, pivot_index)
        var pivot_imag = left_matrix.get_imag(pivot_index, pivot_index)
        for column_index in range(dimension):
            var left_value = complex_divide_pair(
                left_matrix.get_real(pivot_index, column_index),
                left_matrix.get_imag(pivot_index, column_index),
                pivot_real,
                pivot_imag,
            )
            left_matrix.set(
                pivot_index,
                column_index,
                left_value.real_part,
                left_value.imag_part,
            )

            var right_value = complex_divide_pair(
                right_matrix.get_real(pivot_index, column_index),
                right_matrix.get_imag(pivot_index, column_index),
                pivot_real,
                pivot_imag,
            )
            right_matrix.set(
                pivot_index,
                column_index,
                right_value.real_part,
                right_value.imag_part,
            )

        for row_index in range(dimension):
            if row_index == pivot_index:
                continue
            var elimination_real = left_matrix.get_real(row_index, pivot_index)
            var elimination_imag = left_matrix.get_imag(row_index, pivot_index)
            for column_index in range(dimension):
                var left_product = complex_multiply_pair(
                    elimination_real,
                    elimination_imag,
                    left_matrix.get_real(pivot_index, column_index),
                    left_matrix.get_imag(pivot_index, column_index),
                )
                left_matrix.set(
                    row_index,
                    column_index,
                    left_matrix.get_real(row_index, column_index)
                    - left_product.real_part,
                    left_matrix.get_imag(row_index, column_index)
                    - left_product.imag_part,
                )

                var right_product = complex_multiply_pair(
                    elimination_real,
                    elimination_imag,
                    right_matrix.get_real(pivot_index, column_index),
                    right_matrix.get_imag(pivot_index, column_index),
                )
                right_matrix.set(
                    row_index,
                    column_index,
                    right_matrix.get_real(row_index, column_index)
                    - right_product.real_part,
                    right_matrix.get_imag(row_index, column_index)
                    - right_product.imag_part,
                )

    return right_matrix^


struct CanonicalEmbeddingContext(Movable):
    var n_dim: Int
    var p_value: Int
    var phi_p_degree: Int
    var eval_x: ComplexMatrix
    var eval_y: ComplexMatrix
    var eval_w: ComplexMatrix
    var inverse_x: ComplexMatrix
    var inverse_y: ComplexMatrix
    var inverse_w: ComplexMatrix

    fn __init__(out self, n_dim: Int, p_value: Int):
        self.n_dim = n_dim
        self.p_value = p_value

        var unit_exponents = List[Int]()
        for exponent_value in range(1, p_value):
            if _gcd(exponent_value, p_value) == 1:
                unit_exponents.append(exponent_value)
        self.phi_p_degree = len(unit_exponents)

        var primitive_4n_root = _primitive_root_real_imag(4 * n_dim)
        var zeta_real_values = List[Float64](length=n_dim, fill=0.0)
        var zeta_imag_values = List[Float64](length=n_dim, fill=0.0)
        for root_index in range(n_dim):
            var exponent_value = _pow_mod(5, root_index, 4 * n_dim)
            var zeta_power = _complex_pow(
                primitive_4n_root.real_part,
                primitive_4n_root.imag_part,
                exponent_value,
            )
            zeta_real_values[root_index] = zeta_power.real_part
            zeta_imag_values[root_index] = zeta_power.imag_part

        var primitive_p_root = _primitive_root_real_imag(p_value)
        var eta_real_values = List[Float64](
            length=self.phi_p_degree,
            fill=0.0,
        )
        var eta_imag_values = List[Float64](
            length=self.phi_p_degree,
            fill=0.0,
        )
        for exponent_index in range(self.phi_p_degree):
            var eta_power = _complex_pow(
                primitive_p_root.real_part,
                primitive_p_root.imag_part,
                unit_exponents[exponent_index],
            )
            eta_real_values[exponent_index] = eta_power.real_part
            eta_imag_values[exponent_index] = eta_power.imag_part

        self.eval_x = _build_vandermonde(
            zeta_real_values,
            zeta_imag_values,
            self.n_dim,
        )
        self.eval_y = _build_vandermonde(
            zeta_real_values,
            zeta_imag_values,
            self.n_dim,
        )
        self.eval_w = _build_vandermonde(
            eta_real_values,
            eta_imag_values,
            self.phi_p_degree,
        )

        self.inverse_x = _invert_square_matrix(self.eval_x, self.n_dim)
        self.inverse_y = _invert_square_matrix(self.eval_y, self.n_dim)
        self.inverse_w = _invert_square_matrix(self.eval_w, self.phi_p_degree)


fn _apply_separable_transform(
    input_tensor: ComplexTensor3D,
    matrix_x: ComplexMatrix,
    matrix_y: ComplexMatrix,
    matrix_w: ComplexMatrix,
) -> ComplexTensor3D:
    var output_tensor = ComplexTensor3D(
        input_tensor.dim_x,
        input_tensor.dim_y,
        input_tensor.dim_z,
    )

    var stage_x = ComplexTensor3D(
        input_tensor.dim_x,
        input_tensor.dim_y,
        input_tensor.dim_z,
    )
    for output_x in range(input_tensor.dim_x):
        for y_index in range(input_tensor.dim_y):
            for w_index in range(input_tensor.dim_z):
                var accum_real = 0.0
                var accum_imag = 0.0
                for input_x in range(input_tensor.dim_x):
                    var product_value = complex_multiply_pair(
                        matrix_x.get_real(output_x, input_x),
                        matrix_x.get_imag(output_x, input_x),
                        input_tensor.get_real(input_x, y_index, w_index),
                        input_tensor.get_imag(input_x, y_index, w_index),
                    )
                    accum_real += product_value.real_part
                    accum_imag += product_value.imag_part
                stage_x.set(output_x, y_index, w_index, accum_real, accum_imag)

    var stage_y = ComplexTensor3D(
        input_tensor.dim_x,
        input_tensor.dim_y,
        input_tensor.dim_z,
    )
    for x_index in range(input_tensor.dim_x):
        for output_y in range(input_tensor.dim_y):
            for w_index in range(input_tensor.dim_z):
                var accum_real = 0.0
                var accum_imag = 0.0
                for input_y in range(input_tensor.dim_y):
                    var product_value = complex_multiply_pair(
                        matrix_y.get_real(output_y, input_y),
                        matrix_y.get_imag(output_y, input_y),
                        stage_x.get_real(x_index, input_y, w_index),
                        stage_x.get_imag(x_index, input_y, w_index),
                    )
                    accum_real += product_value.real_part
                    accum_imag += product_value.imag_part
                stage_y.set(x_index, output_y, w_index, accum_real, accum_imag)

    for x_index in range(input_tensor.dim_x):
        for y_index in range(input_tensor.dim_y):
            for output_w in range(input_tensor.dim_z):
                var accum_real = 0.0
                var accum_imag = 0.0
                for input_w in range(input_tensor.dim_z):
                    var product_value = complex_multiply_pair(
                        matrix_w.get_real(output_w, input_w),
                        matrix_w.get_imag(output_w, input_w),
                        stage_y.get_real(x_index, y_index, input_w),
                        stage_y.get_imag(x_index, y_index, input_w),
                    )
                    accum_real += product_value.real_part
                    accum_imag += product_value.imag_part
                output_tensor.set(
                    x_index,
                    y_index,
                    output_w,
                    accum_real,
                    accum_imag,
                )

    return output_tensor^


fn encode_canonical_complex(
    input_matrices: ComplexTensor3D,
    context: CanonicalEmbeddingContext,
    scaling_factor: Float64,
    round_to_gaussian_integer: Bool = True,
) -> ComplexTensor3D:
    var coefficient_tensor = _apply_separable_transform(
        input_matrices,
        context.inverse_x,
        context.inverse_y,
        context.inverse_w,
    )

    for coefficient_index in range(len(coefficient_tensor.real_values)):
        var scaled_real = (
            coefficient_tensor.real_values[coefficient_index] * scaling_factor
        )
        var scaled_imag = (
            coefficient_tensor.imag_values[coefficient_index] * scaling_factor
        )
        if round_to_gaussian_integer:
            coefficient_tensor.real_values[coefficient_index] = Float64(
                _round_nearest(scaled_real)
            )
            coefficient_tensor.imag_values[coefficient_index] = Float64(
                _round_nearest(scaled_imag)
            )
        else:
            coefficient_tensor.real_values[coefficient_index] = scaled_real
            coefficient_tensor.imag_values[coefficient_index] = scaled_imag
    return coefficient_tensor^


fn decode_canonical_complex(
    encoded_coefficients: ComplexTensor3D,
    context: CanonicalEmbeddingContext,
    scaling_factor: Float64,
) -> ComplexTensor3D:
    var decoded_tensor = _apply_separable_transform(
        encoded_coefficients,
        context.eval_x,
        context.eval_y,
        context.eval_w,
    )

    var inverse_scaling = 1.0 / scaling_factor
    for coefficient_index in range(len(decoded_tensor.real_values)):
        decoded_tensor.real_values[coefficient_index] *= inverse_scaling
        decoded_tensor.imag_values[coefficient_index] *= inverse_scaling
    return decoded_tensor^
