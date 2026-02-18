fn normalize_modulus_i64(value: Int64, modulus_value: UInt32) -> UInt32:
    var modulus_int64 = Int64(modulus_value)
    var reduced_value = value % modulus_int64
    if reduced_value < 0:
        reduced_value += modulus_int64
    return UInt32(reduced_value)


fn mod_inverse_int64(value: UInt32, modulus_value: UInt32) -> Int64:
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
