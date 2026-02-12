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
    let limit = Int(sqrt(num)) + 1
    for i in range(3, limit, 2):
        if num % i == 0:
            return False
    return True

fn test_is_prime():
    """
    Tests the is_prime function.
    """
    print("Running test_is_prime...")
    # Test prime numbers
    assert is_prime(2), "2 should be prime"
    assert is_prime(3), "3 should be prime"
    assert is_prime(5), "5 should be prime"
    assert is_prime(7), "7 should be prime"
    assert is_prime(11), "11 should be prime"
    assert is_prime(13), "13 should be prime"

    # Test non-prime numbers
    assert not is_prime(1), "1 should not be prime"
    assert not is_prime(4), "4 should not be prime"
    assert not is_prime(6), "6 should not be prime"
    assert not is_prime(8), "8 should not be prime"
    assert not is_prime(9), "9 should not be prime"
    assert not is_prime(10), "10 should not be prime"
    assert not is_prime(12), "12 should not be prime"

    print("test_is_prime passed.")

fn main():
    test_is_prime()
