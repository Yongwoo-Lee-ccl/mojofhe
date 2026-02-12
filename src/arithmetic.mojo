from math import sqrt
from collections import List


fn _mul_mod(a: Int, b: Int, m: Int) -> Int:
    return Int((Int128(a) * Int128(b)) % Int128(m))


fn _power(base: Int, exp: Int, mod: Int) -> Int:
    var res: Int = 1
    var b: Int = base % mod
    var e: Int = exp

    while e > 0:
        if e % 2 == 1:
            res = _mul_mod(res, b, mod)
        b = _mul_mod(b, b, mod)
        e //= 2
    return res


fn is_prime(num: Int) -> Bool:
    if num <= 3:
        return num > 1
    if num % 2 == 0 or num % 3 == 0:
        return False

    var d: Int = num - 1
    var s: Int = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    var bases = List[Int]()
    bases.append(2)
    bases.append(3)
    bases.append(5)
    bases.append(7)
    bases.append(11)
    bases.append(13)
    bases.append(17)
    bases.append(19)
    bases.append(23)

    for i in range(len(bases)):
        var a = bases[i]
        if a >= num:
            break

        var x = _power(a, d, num)
        if x == 1 or x == num - 1:
            continue

        var composite = True
        for _ in range(s - 1):
            x = _mul_mod(x, x, num)
            if x == num - 1:
                composite = False
                break

        if composite:
            return False

    return True
