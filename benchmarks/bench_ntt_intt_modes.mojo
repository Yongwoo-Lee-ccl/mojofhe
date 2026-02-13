from collections import List
from src.modular import find_suitable_q
from src.ntt import (
    apply_ixw_quotient_ntt_inplace,
    apply_ixw_quotient_intt_inplace,
)


fn _make_input(total_length: Int, q_modulus: UInt32) -> List[UInt32]:
    var input_values = List[UInt32]()
    var q_int = Int(q_modulus)
    for coefficient_index in range(total_length):
        var seeded_value = (coefficient_index * 8191 + 12345) % q_int
        input_values.append(UInt32(seeded_value))
    return input_values^


fn _copy_values(source_values: List[UInt32]) -> List[UInt32]:
    var copied_values = List[UInt32]()
    for coefficient_index in range(len(source_values)):
        copied_values.append(source_values[coefficient_index])
    return copied_values^


fn _run_mode(
    mode_name: StringLiteral,
    base_input: List[UInt32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    q_modulus: UInt32,
    use_montgomery: Bool,
):
    var iterations = 4
    var checksum_accumulator: Int = 0
    for iteration_index in range(iterations):
        _ = iteration_index
        var working_values = _copy_values(base_input)
        var ntt_success = apply_ixw_quotient_ntt_inplace(
            working_values,
            n_power_of_2,
            p_power_of_3,
            q_modulus,
            use_montgomery,
        )
        if not ntt_success:
            print(mode_name, "NTT setup failed.")
            return
        var intt_success = apply_ixw_quotient_intt_inplace(
            working_values,
            n_power_of_2,
            p_power_of_3,
            q_modulus,
            use_montgomery,
        )
        if not intt_success:
            print(mode_name, "INTT setup failed.")
            return
        checksum_accumulator += Int(working_values[0])
        checksum_accumulator += Int(working_values[len(working_values) - 1])

    print(mode_name, "checksum:", checksum_accumulator)


fn main() raises:
    var n_power_of_2 = 1 << 7
    var p_power_of_3 = 81  # 3^4
    var target_bits = 25
    var q_modulus_u64 = find_suitable_q(n_power_of_2, p_power_of_3, target_bits)
    var q_modulus = UInt32(q_modulus_u64)

    var phi_p_degree = 2 * (p_power_of_3 // 3)
    var total_length = 2 * n_power_of_2 * phi_p_degree
    var base_input = _make_input(total_length, q_modulus)

    _run_mode(
        "barrett",
        base_input,
        n_power_of_2,
        p_power_of_3,
        q_modulus,
        False,
    )
    _run_mode(
        "montgomery",
        base_input,
        n_power_of_2,
        p_power_of_3,
        q_modulus,
        True,
    )
