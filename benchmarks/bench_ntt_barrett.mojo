from collections import List
from src.modular import find_suitable_q
from src.ntt import apply_ixw_quotient_ntt_inplace


fn _make_input(total_length: Int, q_modulus: Int) -> List[Int32]:
    var input_values = List[Int32]()
    for coefficient_index in range(total_length):
        var seeded_value = (coefficient_index * 8191 + 12345) % q_modulus
        input_values.append(Int32(seeded_value))
    return input_values^


fn _copy_values(source_values: List[Int32]) -> List[Int32]:
    var copied_values = List[Int32]()
    for coefficient_index in range(len(source_values)):
        copied_values.append(source_values[coefficient_index])
    return copied_values^


fn main() raises:
    var n_power_of_2 = 1 << 7
    var p_power_of_3 = 81  # 3^4
    var target_bits = 25
    var q_modulus = find_suitable_q(n_power_of_2, p_power_of_3, target_bits)

    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var total_length = 2 * n_power_of_2 * phi_p_degree
    var base_input = _make_input(total_length, q_modulus)

    var iterations = 6
    var checksum_accumulator: Int = 0
    for iteration_index in range(iterations):
        _ = iteration_index
        var working_values = _copy_values(base_input)
        var success_flag = apply_ixw_quotient_ntt_inplace(
            working_values, n_power_of_2, p_power_of_3, q_modulus, False
        )
        if not success_flag:
            print("Barrett benchmark setup failed.")
            return
        checksum_accumulator += Int(working_values[0])
        checksum_accumulator += Int(working_values[total_length - 1])

    print("barrett checksum:", checksum_accumulator)
