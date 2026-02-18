from collections import List
from src.modular import (
    compute_barrett_ratio,
    multiply_mod_barrett,
    find_primitive_root,
)
from src.ntt import apply_cyclotomic_pruned_ntt
from src.canonical_embedding import ComplexTensor3D
from src.complex_utils import complex_multiply_pair, complex_conjugate_pair
from src.number_theory_utils import normalize_modulus_i64, mod_inverse_int64


fn _add_modulus(
    left_value: UInt32, right_value: UInt32, modulus_value: UInt32
) -> UInt32:
    var summed_value = left_value + right_value
    if summed_value >= modulus_value:
        summed_value -= modulus_value
    return summed_value


fn _apply_w_inverse_automorphism(
    coefficient_values: List[UInt32],
    p_power_of_3: Int,
    modulus_value: UInt32,
) -> List[UInt32]:
    var phi_degree = 2 * (p_power_of_3 // 3)
    var reduction_shift = p_power_of_3 // 3

    var exponent_space_values = List[UInt32](length=p_power_of_3, fill=0)
    for coefficient_index in range(phi_degree):
        var coefficient_value = coefficient_values[coefficient_index]
        if coefficient_value == 0:
            continue
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
            exponent_space_values[first_reduction_index] = normalize_modulus_i64(
                Int64(exponent_space_values[first_reduction_index])
                - Int64(top_coefficient),
                modulus_value,
            )
            exponent_space_values[second_reduction_index] = normalize_modulus_i64(
                Int64(exponent_space_values[second_reduction_index])
                - Int64(top_coefficient),
                modulus_value,
            )
        reduction_index -= 1

    var reduced_values = List[UInt32](length=phi_degree, fill=0)
    for coefficient_index in range(phi_degree):
        reduced_values[coefficient_index] = exponent_space_values[
            coefficient_index
        ]
    return reduced_values^


fn _pointwise_w_multiply_ntt(
    left_values: List[UInt32],
    right_values: List[UInt32],
    modulus_value: UInt32,
    barrett_ratio: UInt64,
) -> List[UInt32]:
    var product_values = List[UInt32](length=len(left_values), fill=0)
    for coefficient_index in range(len(left_values)):
        product_values[coefficient_index] = multiply_mod_barrett(
            left_values[coefficient_index],
            right_values[coefficient_index],
            modulus_value,
            barrett_ratio,
        )
    return product_values^


struct EncodedPolynomial(Movable):
    # Coefficients for a(X,Y,W)=sum_{x,y,w} a_{xyw} X^x Y^y W^w
    # with x,y in [0,n), w in [0,phi(p)).
    var coefficient_values: List[UInt32]
    var n_dim: Int
    var p_power_of_3: Int
    var phi_p_degree: Int
    var q_modulus: UInt32
    var barrett_ratio: UInt64
    var is_w_ntt: Bool

    fn __init__(out self, n_dim: Int, p_power_of_3: Int, q_modulus: UInt32):
        self.n_dim = n_dim
        self.p_power_of_3 = p_power_of_3
        self.phi_p_degree = 2 * (p_power_of_3 // 3)
        self.q_modulus = q_modulus
        self.barrett_ratio = compute_barrett_ratio(q_modulus)
        self.coefficient_values = List[UInt32](
            length=n_dim * n_dim * self.phi_p_degree,
            fill=0,
        )
        self.is_w_ntt = False

    fn _index(self, x_index: Int, y_index: Int, w_index: Int) -> Int:
        return (
            x_index * self.n_dim * self.phi_p_degree
            + y_index * self.phi_p_degree
            + w_index
        )

    fn set_coefficient(
        mut self, x_index: Int, y_index: Int, w_index: Int, value: Int
    ):
        self.coefficient_values[
            self._index(x_index, y_index, w_index)
        ] = normalize_modulus_i64(Int64(value), self.q_modulus)

    fn get_coefficient(
        self, x_index: Int, y_index: Int, w_index: Int
    ) -> UInt32:
        return self.coefficient_values[self._index(x_index, y_index, w_index)]

    fn get_w_polynomial(self, x_index: Int, y_index: Int) -> List[UInt32]:
        var result_values = List[UInt32](length=self.phi_p_degree, fill=0)
        var base_index = self._index(x_index, y_index, 0)
        for coefficient_index in range(self.phi_p_degree):
            result_values[coefficient_index] = self.coefficient_values[
                base_index + coefficient_index
            ]
        return result_values^

    fn set_w_polynomial(
        mut self,
        x_index: Int,
        y_index: Int,
        coefficient_values: List[UInt32],
    ):
        var base_index = self._index(x_index, y_index, 0)
        for coefficient_index in range(self.phi_p_degree):
            self.coefficient_values[
                base_index + coefficient_index
            ] = coefficient_values[coefficient_index]

    fn transform_w_to_ntt(mut self):
        if self.is_w_ntt:
            return

        var root_p = UInt32(
            find_primitive_root(
                UInt64(self.q_modulus), UInt64(self.p_power_of_3)
            )
        )
        if root_p == 0:
            return

        var m_parameter = self.p_power_of_3 // 3
        var number_of_slices = self.n_dim * self.n_dim
        for slice_index in range(number_of_slices):
            var slice_base_offset = slice_index * self.phi_p_degree
            apply_cyclotomic_pruned_ntt(
                self.coefficient_values,
                m_parameter,
                root_p,
                self.q_modulus,
                1,
                slice_base_offset,
            )
        self.is_w_ntt = True

    fn apply_rhs_automorphism(mut self) -> EncodedPolynomial:
        # Computes b'(X,Y,W)=b(X^{-1},Y,W^{-1}) in coefficient form.
        # This must be applied before W-NTT conversion.
        if self.is_w_ntt:
            var copy_value = EncodedPolynomial(
                self.n_dim, self.p_power_of_3, self.q_modulus
            )
            for coefficient_index in range(len(self.coefficient_values)):
                copy_value.coefficient_values[
                    coefficient_index
                ] = self.coefficient_values[coefficient_index]
            copy_value.is_w_ntt = True
            return copy_value^

        var transformed_value = EncodedPolynomial(
            self.n_dim, self.p_power_of_3, self.q_modulus
        )
        for x_index in range(self.n_dim):
            var mapped_x_index = (self.n_dim - x_index) % self.n_dim
            for y_index in range(self.n_dim):
                var w_coefficients = self.get_w_polynomial(x_index, y_index)
                var mapped_w_coefficients = _apply_w_inverse_automorphism(
                    w_coefficients,
                    self.p_power_of_3,
                    self.q_modulus,
                )
                transformed_value.set_w_polynomial(
                    mapped_x_index,
                    y_index,
                    mapped_w_coefficients,
                )
        return transformed_value^

    fn trace_multiply(
        mut self, mut rhs: EncodedPolynomial
    ) -> EncodedPolynomial:
        # Computes c = Tr_Z(a(X,Z,W) * b(Y^{-1},Z^{-1},W^{-1})).
        # Uses existing W-axis NTT pipeline and returns c in W-NTT domain.
        var rhs_prime = rhs.apply_rhs_automorphism()

        self.transform_w_to_ntt()
        rhs_prime.transform_w_to_ntt()

        var result_value = EncodedPolynomial(
            self.n_dim, self.p_power_of_3, self.q_modulus
        )
        var n_dim = self.n_dim
        var phi_degree = self.phi_p_degree

        for row_index in range(n_dim):
            var left_row_base = row_index * n_dim * phi_degree
            for column_index in range(n_dim):
                var right_row_base = column_index * n_dim * phi_degree
                var accumulated_values = List[UInt32](
                    length=phi_degree,
                    fill=0,
                )
                for inner_index in range(n_dim):
                    var left_lane_base = left_row_base + inner_index * phi_degree
                    var right_lane_base = (
                        right_row_base + inner_index * phi_degree
                    )
                    for coefficient_index in range(phi_degree):
                        var product_value = multiply_mod_barrett(
                            self.coefficient_values[
                                left_lane_base + coefficient_index
                            ],
                            rhs_prime.coefficient_values[
                                right_lane_base + coefficient_index
                            ],
                            self.q_modulus,
                            self.barrett_ratio,
                        )
                        accumulated_values[coefficient_index] = _add_modulus(
                            accumulated_values[coefficient_index],
                            product_value,
                            self.q_modulus,
                        )
                result_value.set_w_polynomial(
                    row_index,
                    column_index,
                    accumulated_values,
                )

        result_value.is_w_ntt = True
        return result_value^

    fn scale_by_inverse_n(mut self):
        var n_inverse_int64 = mod_inverse_int64(
            UInt32(self.n_dim), self.q_modulus
        )
        if n_inverse_int64 < 0:
            return
        var n_inverse = UInt32(UInt64(n_inverse_int64))
        for coefficient_index in range(len(self.coefficient_values)):
            self.coefficient_values[coefficient_index] = multiply_mod_barrett(
                self.coefficient_values[coefficient_index],
                n_inverse,
                self.q_modulus,
                self.barrett_ratio,
            )


struct _ComplexPolynomial(Movable):
    var real_values: List[Float64]
    var imag_values: List[Float64]

    fn __init__(out self, degree: Int):
        self.real_values = List[Float64](length=degree, fill=0.0)
        self.imag_values = List[Float64](length=degree, fill=0.0)


fn _w_inverse_polynomial_complex(
    coefficient_real_values: List[Float64],
    coefficient_imag_values: List[Float64],
    p_value: Int,
) -> _ComplexPolynomial:
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
        if top_real != 0.0 or top_imag != 0.0:
            exponent_real_values[reduction_index] = 0.0
            exponent_imag_values[reduction_index] = 0.0
            var first_reduction_index = reduction_index - reduction_shift
            var second_reduction_index = reduction_index - phi_degree

            exponent_real_values[first_reduction_index] -= top_real
            exponent_imag_values[first_reduction_index] -= top_imag
            exponent_real_values[second_reduction_index] -= top_real
            exponent_imag_values[second_reduction_index] -= top_imag
        reduction_index -= 1

    var output_values = _ComplexPolynomial(phi_degree)
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
) -> _ComplexPolynomial:
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
            var product_value = complex_multiply_pair(
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
        if top_real != 0.0 or top_imag != 0.0:
            unreduced_real_values[reduction_index] = 0.0
            unreduced_imag_values[reduction_index] = 0.0
            var first_reduction_index = reduction_index - reduction_shift
            var second_reduction_index = reduction_index - phi_degree

            unreduced_real_values[first_reduction_index] -= top_real
            unreduced_imag_values[first_reduction_index] -= top_imag
            unreduced_real_values[second_reduction_index] -= top_real
            unreduced_imag_values[second_reduction_index] -= top_imag
        reduction_index -= 1

    var output_values = _ComplexPolynomial(phi_degree)
    for coefficient_index in range(phi_degree):
        output_values.real_values[coefficient_index] = unreduced_real_values[
            coefficient_index
        ]
        output_values.imag_values[coefficient_index] = unreduced_imag_values[
            coefficient_index
        ]
    return output_values^


fn _rhs_automorphism_for_trace_complex(
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
                var conjugated_value = complex_conjugate_pair(
                    encoded_tensor.get_real(x_index, y_index, w_index),
                    encoded_tensor.get_imag(x_index, y_index, w_index),
                )
                var mapped_real = conjugated_value.real_part
                var mapped_imag = conjugated_value.imag_part
                if x_index != 0:
                    var multiplied = complex_multiply_pair(
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


fn trace_multiply_encoded_complex(
    encoded_left: ComplexTensor3D,
    encoded_right: ComplexTensor3D,
    n_dim: Int,
    p_value: Int,
    phi_degree: Int,
) -> ComplexTensor3D:
    # Computes c = Tr_Z(a(X,Z,W) * conjugate(b)(Y^{-1},Z^{-1},W^{-1}))
    # on encoded coefficient tensors in C.
    var rhs_prime = _rhs_automorphism_for_trace_complex(
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
