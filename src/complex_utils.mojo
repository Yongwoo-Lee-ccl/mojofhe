from complex import ComplexFloat64


struct ComplexPair(Copyable, Movable):
    var real_part: Float64
    var imag_part: Float64

    fn __init__(out self, real_part: Float64 = 0.0, imag_part: Float64 = 0.0):
        self.real_part = real_part
        self.imag_part = imag_part


fn abs_f64(value: Float64) -> Float64:
    if value < 0.0:
        return -value
    return value


fn complex_multiply_pair(
    left_real: Float64,
    left_imag: Float64,
    right_real: Float64,
    right_imag: Float64,
) -> ComplexPair:
    var left_value = ComplexFloat64(left_real, left_imag)
    var right_value = ComplexFloat64(right_real, right_imag)
    var result_value = left_value * right_value
    return ComplexPair(result_value.re, result_value.im)


fn complex_divide_pair(
    numerator_real: Float64,
    numerator_imag: Float64,
    denominator_real: Float64,
    denominator_imag: Float64,
) -> ComplexPair:
    var denominator_abs_squared = denominator_real * denominator_real + (
        denominator_imag * denominator_imag
    )
    if denominator_abs_squared == 0.0:
        return ComplexPair()

    var numerator_value = ComplexFloat64(numerator_real, numerator_imag)
    var denominator_value = ComplexFloat64(denominator_real, denominator_imag)
    var result_value = numerator_value / denominator_value
    return ComplexPair(result_value.re, result_value.im)


fn complex_conjugate_pair(real_part: Float64, imag_part: Float64) -> ComplexPair:
    return ComplexPair(real_part, -imag_part)


fn complex_close_values(
    left_real: Float64,
    left_imag: Float64,
    right_real: Float64,
    right_imag: Float64,
    tolerance: Float64,
) -> Bool:
    return (
        abs_f64(left_real - right_real) <= tolerance
        and abs_f64(left_imag - right_imag) <= tolerance
    )
