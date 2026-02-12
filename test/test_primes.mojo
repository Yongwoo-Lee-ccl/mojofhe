from testing import assert_true, assert_false, TestSuite
from src.arithmetic import is_prime


def test_is_prime():
    """
    Tests the is_prime function.
    """
    print("Running test_is_prime...")
    # Test small prime numbers
    assert_true(is_prime(2), "2 should be prime")
    assert_true(is_prime(3), "3 should be prime")
    assert_true(is_prime(5), "5 should be prime")
    assert_true(is_prime(7), "7 should be prime")
    assert_true(is_prime(11), "11 should be prime")
    assert_true(is_prime(13), "13 should be prime")

    # Test larger prime numbers
    assert_true(is_prime(65537), "2^16 + 1 (65537) should be prime")

    # Test non-prime numbers
    assert_false(is_prime(1), "1 should not be prime")
    assert_false(is_prime(4), "4 should not be prime")
    assert_false(is_prime(6), "6 should not be prime")
    assert_false(is_prime(8), "8 should not be prime")
    assert_false(is_prime(9), "9 should not be prime")
    assert_false(is_prime(10), "10 should not be prime")
    assert_false(is_prime(12), "12 should not be prime")
    assert_false(is_prime(65535), "2^16 - 1 (65535) should not be prime")

    print("test_is_prime passed.")


def main():
    TestSuite.discover_tests[__functions_in_module()]().run()
