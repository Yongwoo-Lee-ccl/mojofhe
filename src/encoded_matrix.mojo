from collections import List
from src.modular import compute_barrett_ratio, multiply_mod_barrett


fn _normalize_modulus(value: Int64, modulus_value: UInt32) -> UInt32:
    var modulus_int64 = Int64(modulus_value)
    var reduced_value = value % modulus_int64
    if reduced_value < 0:
        reduced_value += modulus_int64
    return UInt32(reduced_value)


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
    return _normalize_modulus(
        Int64(left_value) - Int64(right_value), modulus_value
    )


fn _multiply_modulus(
    left_value: UInt32,
    right_value: UInt32,
    modulus_value: UInt32,
    barrett_ratio: UInt64,
) -> UInt32:
    return multiply_mod_barrett(
        left_value, right_value, modulus_value, barrett_ratio
    )


fn _mod_inverse(value: UInt32, modulus_value: UInt32) -> UInt32:
    var normalized_value = Int64(value % modulus_value)
    var modulus_int64 = Int64(modulus_value)

    var previous_remainder = modulus_int64
    var current_remainder = normalized_value
    var previous_coefficient: Int64 = 0
    var current_coefficient: Int64 = 1

    while current_remainder != 0:
        var quotient_value = previous_remainder // current_remainder
        var next_remainder = (
            previous_remainder - quotient_value * current_remainder
        )
        previous_remainder = current_remainder
        current_remainder = next_remainder

        var next_coefficient = (
            previous_coefficient - quotient_value * current_coefficient
        )
        previous_coefficient = current_coefficient
        current_coefficient = next_coefficient

    if previous_remainder != 1:
        return 0

    var normalized_coefficient = previous_coefficient % modulus_int64
    if normalized_coefficient < 0:
        normalized_coefficient += modulus_int64
    return UInt32(normalized_coefficient)


fn _multiply_w_polynomials_mod_cyclotomic(
    left_values: List[UInt32],
    right_values: List[UInt32],
    p_power_of_3: Int,
    modulus_value: UInt32,
    barrett_ratio: UInt64,
) -> List[UInt32]:
    var phi_degree = 2 * (p_power_of_3 // 3)
    var reduction_shift = p_power_of_3 // 3
    var product_degree_bound = 2 * phi_degree - 1

    var unreduced_values = List[UInt32](length=product_degree_bound, fill=0)
    for left_index in range(phi_degree):
        var left_coefficient = left_values[left_index]
        if left_coefficient == 0:
            continue
        for right_index in range(phi_degree):
            var right_coefficient = right_values[right_index]
            if right_coefficient == 0:
                continue
            var product_value = _multiply_modulus(
                left_coefficient,
                right_coefficient,
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


struct EncodedPolynomial(Movable):
    # Coefficients for a(X,Y,W)=sum_{x,y,w} a_{xyw} X^x Y^y W^w
    # with x,y in [0,n), w in [0,phi(p)).
    var coefficient_values: List[UInt32]
    var n_dim: Int
    var p_power_of_3: Int
    var phi_p_degree: Int
    var q_modulus: UInt32
    var barrett_ratio: UInt64

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
        ] = _normalize_modulus(Int64(value), self.q_modulus)

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

    fn apply_rhs_automorphism(mut self) -> EncodedPolynomial:
        # Computes b'(X,Y,W)=b(X^{-1},Y,W^{-1}) in coefficient form.
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
        # Computes c = Tr_Z(a(X,Z,W) * b(Y^{-1},Z^{-1},W^{-1}))
        # as matrix multiplication over R^(p): C = A * (B')^T.
        var rhs_prime = rhs.apply_rhs_automorphism()
        var result_value = EncodedPolynomial(
            self.n_dim, self.p_power_of_3, self.q_modulus
        )

        for row_index in range(self.n_dim):
            for column_index in range(self.n_dim):
                var accumulated_values = List[UInt32](
                    length=self.phi_p_degree,
                    fill=0,
                )
                for inner_index in range(self.n_dim):
                    var left_w_coefficients = self.get_w_polynomial(
                        row_index, inner_index
                    )
                    var right_w_coefficients = rhs_prime.get_w_polynomial(
                        column_index, inner_index
                    )
                    var product_w_coefficients = (
                        _multiply_w_polynomials_mod_cyclotomic(
                            left_w_coefficients,
                            right_w_coefficients,
                            self.p_power_of_3,
                            self.q_modulus,
                            self.barrett_ratio,
                        )
                    )
                    for coefficient_index in range(self.phi_p_degree):
                        accumulated_values[coefficient_index] = _add_modulus(
                            accumulated_values[coefficient_index],
                            product_w_coefficients[coefficient_index],
                            self.q_modulus,
                        )
                result_value.set_w_polynomial(
                    row_index,
                    column_index,
                    accumulated_values,
                )

        return result_value^

    fn scale_by_inverse_n(mut self):
        var n_inverse = _mod_inverse(UInt32(self.n_dim), self.q_modulus)
        if n_inverse == 0:
            return
        for coefficient_index in range(len(self.coefficient_values)):
            self.coefficient_values[coefficient_index] = _multiply_modulus(
                self.coefficient_values[coefficient_index],
                n_inverse,
                self.q_modulus,
                self.barrett_ratio,
            )
