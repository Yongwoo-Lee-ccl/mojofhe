from collections import List
from src.modular import find_primitive_root, mod_pow
from src.ntt import (
    apply_i_axis_transform,
    apply_radix2_dif_ntt,
    apply_cyclotomic_pruned_ntt,
)


struct Polynomial(Movable):
    # Coefficients are stored as Int32 lanes for SIMD/GPU-friendly arithmetic.
    # We assume q_modulus < 2^30 so adds/subs and twiddle products stay in safe ranges.
    var coefficient_values: List[Int32]
    var total_length: Int
    var q_modulus: Int
    var n_power_of_2: Int
    var p_power_of_3: Int
    var phi_p_degree: Int

    var is_in_ntt_i: Bool
    var is_in_ntt_x: Bool
    var is_in_ntt_w: Bool

    fn __init__(out self, n_val: Int, p_val: Int, q_val: Int):
        if q_val >= (1 << 30):
            print(
                "Warning: q_modulus should be < 2^30 for Int32 lane"
                " acceleration."
            )
        self.n_power_of_2 = n_val
        self.p_power_of_3 = p_val
        self.q_modulus = q_val
        self.phi_p_degree = 2 * (p_val // 3)
        self.total_length = 2 * self.n_power_of_2 * self.phi_p_degree
        self.coefficient_values = List[Int32](length=self.total_length, fill=0)
        self.is_in_ntt_i = False
        self.is_in_ntt_x = False
        self.is_in_ntt_w = False

    fn set_coefficient(mut self, index: Int, value: Int):
        var normalized_value = value % self.q_modulus
        if normalized_value < 0:
            normalized_value += self.q_modulus
        self.coefficient_values[index] = Int32(normalized_value)

    fn _add_modulus_if_needed(self, left_value: Int, right_value: Int) -> Int32:
        var summed_value = left_value + right_value
        if summed_value >= self.q_modulus:
            summed_value -= self.q_modulus
        return Int32(summed_value)

    fn _multiply_modulus_i32(self, left_value: Int, right_value: Int) -> Int32:
        return Int32(
            (Int64(left_value) * Int64(right_value)) % Int64(self.q_modulus)
        )

    fn transform_to_full_ntt(mut self):
        if not self.is_in_ntt_i:
            var root_j = Int32(find_primitive_root(self.q_modulus, 4))
            apply_i_axis_transform(
                self.coefficient_values,
                self.total_length,
                root_j,
                Int32(self.q_modulus),
            )
            self.is_in_ntt_i = True
        if not self.is_in_ntt_x:
            var root_4n = find_primitive_root(
                self.q_modulus, 4 * self.n_power_of_2
            )
            var component_stride = self.total_length // 2
            apply_radix2_dif_ntt(
                self.coefficient_values,
                self.n_power_of_2,
                Int32(root_4n),
                Int32(self.q_modulus),
                self.phi_p_degree,
                0,
            )
            var root_for_minus_j = mod_pow(root_4n, 3, self.q_modulus)
            apply_radix2_dif_ntt(
                self.coefficient_values,
                self.n_power_of_2,
                Int32(root_for_minus_j),
                Int32(self.q_modulus),
                self.phi_p_degree,
                component_stride,
            )
            self.is_in_ntt_x = True
        if not self.is_in_ntt_w:
            var root_p = Int32(
                find_primitive_root(self.q_modulus, self.p_power_of_3)
            )
            var m_parameter = self.p_power_of_3 // 3
            var number_of_slices = 2 * self.n_power_of_2
            for slice_index in range(number_of_slices):
                var slice_base_offset = slice_index * self.phi_p_degree
                apply_cyclotomic_pruned_ntt(
                    self.coefficient_values,
                    m_parameter,
                    root_p,
                    Int32(self.q_modulus),
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
                Int(self.coefficient_values[offset]),
                Int(other.coefficient_values[offset]),
            )
        result_poly.is_in_ntt_i = self.is_in_ntt_i
        result_poly.is_in_ntt_x = self.is_in_ntt_x
        result_poly.is_in_ntt_w = self.is_in_ntt_w
        return result_poly^

    fn add_inplace(mut self, other: Polynomial):
        for offset in range(self.total_length):
            self.coefficient_values[offset] = self._add_modulus_if_needed(
                Int(self.coefficient_values[offset]),
                Int(other.coefficient_values[offset]),
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
            var product_value = self._multiply_modulus_i32(
                Int(self.coefficient_values[offset]),
                Int(other.coefficient_values[offset]),
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
            self.coefficient_values[offset] = self._multiply_modulus_i32(
                Int(self.coefficient_values[offset]),
                Int(other.coefficient_values[offset]),
            )
        self.is_in_ntt_i = True
        self.is_in_ntt_x = True
        self.is_in_ntt_w = True
