from math import sqrt


fn is_prime(num: Int) -> Bool:
    """
    Checks if a number is prime.
    """
    if num <= 1:
        return False
    if num == 2:
        return True
    if num % 2 == 0:
        return False
    var limit = Int(sqrt(num)) + 1
    for i in range(3, limit, 2):
        if num % i == 0:
            return False
    return True
