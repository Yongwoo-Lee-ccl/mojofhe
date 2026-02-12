from src.arithmetic import is_prime


fn mod_pow(base_value: Int, exponent_value: Int, modulus_value: Int) -> Int:
    var result_value: Int = 1
    var base_copy: Int = base_value % modulus_value
    var exponent_copy: Int = exponent_value

    while exponent_copy > 0:
        if exponent_copy % 2 == 1:
            result_value = Int(
                (Int128(result_value) * Int128(base_copy))
                % Int128(modulus_value)
            )
        base_copy = Int(
            (Int128(base_copy) * Int128(base_copy)) % Int128(modulus_value)
        )
        exponent_copy //= 2
    return result_value


fn find_suitable_q(n_dim: Int, p_dim: Int, target_bit_length: Int) -> Int:
    """
    Finds a prime q such that q = 1 mod (4 * n_dim * p_dim) and q has roughly target_bit_length bits.
    """
    var target_magnitude = 1 << target_bit_length
    var cyclic_order = 4 * n_dim * p_dim

    var search_start_multiplier = target_magnitude // cyclic_order
    var search_offset = 0
    while True:
        var candidate_up = (
            search_start_multiplier + search_offset
        ) * cyclic_order + 1
        if candidate_up > 0 and is_prime(candidate_up):
            return candidate_up

        if search_offset > 0:
            var candidate_down = (
                search_start_multiplier - search_offset
            ) * cyclic_order + 1
            if candidate_down > 0 and is_prime(candidate_down):
                return candidate_down

        search_offset += 1


fn find_primitive_root(modulus_value: Int, order_value: Int) -> Int:
    """
    Finds a primitive n-th root of unity modulo q.
    """
    if (modulus_value - 1) % order_value != 0:
        return -1

    var search_candidate = 2
    while search_candidate < modulus_value:
        var exponent_for_root = (modulus_value - 1) // order_value
        var root_candidate = mod_pow(
            search_candidate, exponent_for_root, modulus_value
        )

        if root_candidate != 1:
            var is_actually_primitive = True

            # Check prime factors of order_value to ensure it's primitive
            # We assume order_value consists of factors 2, 3, 5 based on problem constraints
            if order_value % 2 == 0:
                if (
                    mod_pow(root_candidate, order_value // 2, modulus_value)
                    == 1
                ):
                    is_actually_primitive = False

            if is_actually_primitive and order_value % 3 == 0:
                if (
                    mod_pow(root_candidate, order_value // 3, modulus_value)
                    == 1
                ):
                    is_actually_primitive = False

            if is_actually_primitive and order_value % 5 == 0:
                if (
                    mod_pow(root_candidate, order_value // 5, modulus_value)
                    == 1
                ):
                    is_actually_primitive = False

            if is_actually_primitive:
                return root_candidate

        search_candidate += 1

    return -1
