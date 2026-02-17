from collections import List
from random import random_si64
from src.modular import (
    find_primitive_root,
    mod_pow,
    compute_barrett_ratio,
    multiply_mod_barrett,
)
from src.ntt import (
    apply_i_axis_transform,
    apply_radix2_dif_ntt,
    apply_cyclotomic_pruned_ntt,
)


fn _mod_inverse(value: UInt32, modulus_value: UInt32) -> Int64:
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
        return -1

    var normalized_coefficient = previous_coefficient % modulus_int64
    if normalized_coefficient < 0:
        normalized_coefficient += modulus_int64
    return normalized_coefficient


fn integer_to_rns(value: UInt128, moduli: List[UInt32]) -> List[UInt32]:
    var residue_values = List[UInt32]()
    for modulus_index in range(len(moduli)):
        var modulus_value = moduli[modulus_index]
        residue_values.append(UInt32(value % UInt128(modulus_value)))
    return residue_values^


fn integer_from_rns(
    residue_values: List[UInt32], moduli: List[UInt32]
) -> UInt128:
    var combined_modulus: UInt128 = 1
    for modulus_index in range(len(moduli)):
        combined_modulus *= UInt128(moduli[modulus_index])

    var reconstructed_value: UInt128 = 0
    for modulus_index in range(len(moduli)):
        var modulus_value = moduli[modulus_index]
        var partial_modulus = combined_modulus // UInt128(modulus_value)
        var partial_remainder = UInt32(partial_modulus % UInt128(modulus_value))
        var inverse_res = _mod_inverse(partial_remainder, modulus_value)
        if inverse_res < 0:
            continue
        var inverse_value = UInt128(UInt64(inverse_res))
        var crt_term = (
            UInt128(residue_values[modulus_index])
            * partial_modulus
            * inverse_value
        ) % combined_modulus
        reconstructed_value = (
            reconstructed_value + crt_term
        ) % combined_modulus
    return reconstructed_value


struct Polynomial(Movable):
    # Coefficients are stored as UInt32 lanes.
    # We assume q_modulus < 2^32.
    var coefficient_values: List[UInt32]
    var total_length: Int
    var q_modulus: UInt32
    var n_power_of_2: Int
    var p_power_of_3: Int
    var phi_p_degree: Int
    var barrett_ratio: UInt64

    var is_in_ntt_i: Bool
    var is_in_ntt_x: Bool
    var is_in_ntt_w: Bool

    fn __init__(out self, n_val: Int, p_val: Int, q_val: UInt32):
        self.n_power_of_2 = n_val
        self.p_power_of_3 = p_val
        self.q_modulus = q_val
        self.phi_p_degree = 2 * (p_val // 3)
        self.total_length = 2 * self.n_power_of_2 * self.phi_p_degree
        self.barrett_ratio = compute_barrett_ratio(self.q_modulus)
        self.coefficient_values = List[UInt32](length=self.total_length, fill=0)
        self.is_in_ntt_i = False
        self.is_in_ntt_x = False
        self.is_in_ntt_w = False

    fn set_coefficient(mut self, index: Int, value: Int):
        """
        Converter from signed input to unsigned coefficient modulo q.
        """
        var q_int = Int(self.q_modulus)
        var normalized_value = value % q_int
        if normalized_value < 0:
            normalized_value += q_int
        self.coefficient_values[index] = UInt32(normalized_value)

    fn _add_modulus_if_needed(
        self, left_value: UInt32, right_value: UInt32
    ) -> UInt32:
        var summed_value = left_value + right_value
        if summed_value >= self.q_modulus:
            summed_value -= self.q_modulus
        return summed_value

    fn _multiply_modulus_u32(
        self, left_value: UInt32, right_value: UInt32
    ) -> UInt32:
        return multiply_mod_barrett(
            left_value, right_value, self.q_modulus, self.barrett_ratio
        )

    fn transform_to_full_ntt(mut self):
        if not self.is_in_ntt_i:
            var root_j = UInt32(find_primitive_root(UInt64(self.q_modulus), 4))
            apply_i_axis_transform(
                self.coefficient_values,
                self.total_length,
                root_j,
                self.q_modulus,
            )
            self.is_in_ntt_i = True
        if not self.is_in_ntt_x:
            var root_4n = find_primitive_root(
                UInt64(self.q_modulus), UInt64(4 * self.n_power_of_2)
            )
            var component_stride = self.total_length // 2
            apply_radix2_dif_ntt(
                self.coefficient_values,
                self.n_power_of_2,
                UInt32(root_4n),
                self.q_modulus,
                self.phi_p_degree,
                0,
            )
            var root_for_minus_j = mod_pow(root_4n, 3, UInt64(self.q_modulus))
            apply_radix2_dif_ntt(
                self.coefficient_values,
                self.n_power_of_2,
                UInt32(root_for_minus_j),
                self.q_modulus,
                self.phi_p_degree,
                component_stride,
            )
            self.is_in_ntt_x = True
        if not self.is_in_ntt_w:
            var root_p = UInt32(
                find_primitive_root(
                    UInt64(self.q_modulus), UInt64(self.p_power_of_3)
                )
            )
            var m_parameter = self.p_power_of_3 // 3
            var number_of_slices = 2 * self.n_power_of_2
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
            self.is_in_ntt_w = True

    fn add(self, other: Polynomial) -> Polynomial:
        var result_poly = Polynomial(
            self.n_power_of_2, self.p_power_of_3, self.q_modulus
        )
        for offset in range(self.total_length):
            result_poly.coefficient_values[
                offset
            ] = self._add_modulus_if_needed(
                self.coefficient_values[offset],
                other.coefficient_values[offset],
            )
        result_poly.is_in_ntt_i = self.is_in_ntt_i
        result_poly.is_in_ntt_x = self.is_in_ntt_x
        result_poly.is_in_ntt_w = self.is_in_ntt_w
        return result_poly^

    fn add_inplace(mut self, other: Polynomial):
        for offset in range(self.total_length):
            self.coefficient_values[offset] = self._add_modulus_if_needed(
                self.coefficient_values[offset],
                other.coefficient_values[offset],
            )
        self.is_in_ntt_i = self.is_in_ntt_i and other.is_in_ntt_i
        self.is_in_ntt_x = self.is_in_ntt_x and other.is_in_ntt_x
        self.is_in_ntt_w = self.is_in_ntt_w and other.is_in_ntt_w

    fn multiply(mut self, mut other: Polynomial) -> Polynomial:
        self.transform_to_full_ntt()
        other.transform_to_full_ntt()
        var result_poly = Polynomial(
            self.n_power_of_2, self.p_power_of_3, self.q_modulus
        )
        for offset in range(self.total_length):
            var product_value = self._multiply_modulus_u32(
                self.coefficient_values[offset],
                other.coefficient_values[offset],
            )
            result_poly.coefficient_values[offset] = product_value
        result_poly.is_in_ntt_i = True
        result_poly.is_in_ntt_x = True
        result_poly.is_in_ntt_w = True
        return result_poly^

    fn multiply_inplace(mut self, mut other: Polynomial):
        self.transform_to_full_ntt()
        other.transform_to_full_ntt()
        for offset in range(self.total_length):
            self.coefficient_values[offset] = self._multiply_modulus_u32(
                self.coefficient_values[offset],
                other.coefficient_values[offset],
            )
        self.is_in_ntt_i = True
        self.is_in_ntt_x = True
        self.is_in_ntt_w = True

    fn sample_ternary(mut self):
        for coefficient_index in range(self.total_length):
            var random_value = Int(random_si64(-1, 1))
            self.set_coefficient(coefficient_index, random_value)

    fn sample_gaussian(mut self):
        # Approximate centered noise as sum(42 Bernoulli bits) - 21.
        for coefficient_index in range(self.total_length):
            var sum_value: Int = 0
            for bit_index in range(42):
                _ = bit_index
                sum_value += Int(random_si64(0, 1))
            self.set_coefficient(coefficient_index, sum_value - 21)

    fn sample_sparse(mut self, hamming_weight: Int):
        # First reset all coefficients to 0
        for coefficient_index in range(self.total_length):
            self.coefficient_values[coefficient_index] = 0

        var current_weight = 0
        while current_weight < hamming_weight:
            var random_index = Int(random_si64(0, self.total_length - 1))
            if self.coefficient_values[random_index] == 0:
                var sign = Int(random_si64(0, 1))
                if sign == 0:
                    self.set_coefficient(random_index, -1)
                else:
                    self.set_coefficient(random_index, 1)
                current_weight += 1

    fn sample_uniform(mut self):
        var q_int = Int(self.q_modulus)
        for coefficient_index in range(self.total_length):
            var random_value = Int(random_si64(0, q_int - 1))
            self.set_coefficient(coefficient_index, random_value)


struct RNSPolynomial(Movable):
    # Outer axis: RNS channels (q0, q1, q2, ...).
    # Inner axis: polynomial coefficients.
    var residue_channels: List[List[UInt32]]
    var rns_moduli: List[UInt32]
    var n_power_of_2: Int
    var p_power_of_3: Int
    var total_length: Int
    var combined_modulus: UInt128

    fn __init__(
        out self,
        n_val: Int,
        p_val: Int,
        moduli: List[UInt32],
    ):
        self.residue_channels = List[List[UInt32]]()
        self.rns_moduli = List[UInt32]()
        self.n_power_of_2 = n_val
        self.p_power_of_3 = p_val
        self.combined_modulus = 1

        var phi_p_degree = 2 * (self.p_power_of_3 // 3)
        self.total_length = 2 * self.n_power_of_2 * phi_p_degree

        for modulus_index in range(len(moduli)):
            var modulus_value = moduli[modulus_index]
            self.rns_moduli.append(modulus_value)
            self.residue_channels.append(
                List[UInt32](length=self.total_length, fill=0)
            )
            self.combined_modulus *= UInt128(modulus_value)

    fn set_coefficient_from_integer(mut self, index: Int, value: UInt128):
        var residue_values = integer_to_rns(value, self.rns_moduli)
        for modulus_index in range(len(self.rns_moduli)):
            self.residue_channels[modulus_index][index] = residue_values[
                modulus_index
            ]

    fn get_coefficient_rns(self, index: Int) -> List[UInt32]:
        var residue_values = List[UInt32]()
        for modulus_index in range(len(self.rns_moduli)):
            residue_values.append(self.residue_channels[modulus_index][index])
        return residue_values^

    fn get_coefficient_integer(self, index: Int) -> UInt128:
        return integer_from_rns(
            self.get_coefficient_rns(index), self.rns_moduli
        )

    fn add_residuewise(self, other: RNSPolynomial) -> RNSPolynomial:
        var result_poly = RNSPolynomial(
            self.n_power_of_2, self.p_power_of_3, self.rns_moduli
        )
        for modulus_index in range(len(self.rns_moduli)):
            var modulus_value = self.rns_moduli[modulus_index]
            for coefficient_index in range(self.total_length):
                var left_value = self.residue_channels[modulus_index][
                    coefficient_index
                ]
                var right_value = other.residue_channels[modulus_index][
                    coefficient_index
                ]
                var summed_value = left_value + right_value
                if summed_value >= modulus_value:
                    summed_value -= modulus_value
                result_poly.residue_channels[modulus_index][
                    coefficient_index
                ] = summed_value
        return result_poly^

    fn multiply_residuewise(self, other: RNSPolynomial) -> RNSPolynomial:
        var result_poly = RNSPolynomial(
            self.n_power_of_2, self.p_power_of_3, self.rns_moduli
        )
        for modulus_index in range(len(self.rns_moduli)):
            var modulus_value = self.rns_moduli[modulus_index]
            var barrett_ratio = compute_barrett_ratio(modulus_value)
            for coefficient_index in range(self.total_length):
                var left_value = self.residue_channels[modulus_index][
                    coefficient_index
                ]
                var right_value = other.residue_channels[modulus_index][
                    coefficient_index
                ]
                result_poly.residue_channels[modulus_index][
                    coefficient_index
                ] = multiply_mod_barrett(
                    left_value,
                    right_value,
                    modulus_value,
                    barrett_ratio,
                )
        return result_poly^
