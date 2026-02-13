from math import sqrt
from collections import List


fn _multiply_and_modulo(
    multiplicand: UInt64, multiplier: UInt64, modulus_value: UInt64
) -> UInt64:
    return UInt64(
        (UInt128(multiplicand) * UInt128(multiplier)) % UInt128(modulus_value)
    )


fn power_modular(
    base_value: UInt64, exponent_value: UInt64, modulus_value: UInt64
) -> UInt64:
    var result_accumulator: UInt64 = 1
    var base_shifted: UInt64 = base_value % modulus_value
    var exponent_remaining: UInt64 = exponent_value

    while exponent_remaining > 0:
        if exponent_remaining % 2 == 1:
            result_accumulator = _multiply_and_modulo(
                result_accumulator, base_shifted, modulus_value
            )
        base_shifted = _multiply_and_modulo(
            base_shifted, base_shifted, modulus_value
        )
        exponent_remaining //= 2
    return result_accumulator


fn is_prime(number_to_test: UInt64) -> Bool:
    if number_to_test <= 3:
        return number_to_test > 1
    if number_to_test % 2 == 0 or number_to_test % 3 == 0:
        return False

    var odd_component: UInt64 = number_to_test - 1
    var power_of_two_factor: UInt64 = 0
    while odd_component % 2 == 0:
        odd_component //= 2
        power_of_two_factor += 1

    var witness_bases = List[UInt64]()
    witness_bases.append(2)
    witness_bases.append(3)
    witness_bases.append(5)
    witness_bases.append(7)
    witness_bases.append(11)
    witness_bases.append(13)
    witness_bases.append(17)
    witness_bases.append(19)
    witness_bases.append(23)

    for base_index in range(len(witness_bases)):
        var current_base = witness_bases[base_index]
        if current_base >= number_to_test:
            break

        var miller_rabin_value = power_modular(
            current_base, odd_component, number_to_test
        )
        if miller_rabin_value == 1 or miller_rabin_value == number_to_test - 1:
            continue

        var is_probably_composite = True
        for repeat_index in range(Int(power_of_two_factor) - 1):
            _ = repeat_index
            miller_rabin_value = _multiply_and_modulo(
                miller_rabin_value, miller_rabin_value, number_to_test
            )
            if miller_rabin_value == number_to_test - 1:
                is_probably_composite = False
                break

        if is_probably_composite:
            return False

    return True
