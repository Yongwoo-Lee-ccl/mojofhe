# NTT Function Reference (Recent Additions)

This page documents recently added NTT entry points and options in
`src/ntt.mojo`.

## Public Entry Points

### `apply_ixw_quotient_ntt_inplace(...) -> Bool`

```mojo
fn apply_ixw_quotient_ntt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: Int,
    use_montgomery: Bool = False,
) -> Bool
```

Purpose:
- In-place NTT for IXW quotient-ring layout.
- Supports two modular multiplication paths:
  - Barrett (default)
  - Montgomery (`use_montgomery=True`)

Input expectations:
- `coefficient_values` length must be `2 * n_power_of_2 * (2 * (p_power_of_3 / 3))`.
- Primitive roots of orders `4`, `4 * n_power_of_2`, and `p_power_of_3`
  must exist under `modulus_value`.

Returns:
- `True` on success.
- `False` if size/root preconditions are not satisfied.

### `apply_xyw_quotient_ntt_inplace(...) -> Bool`

```mojo
fn apply_xyw_quotient_ntt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_int: Int,
) -> Bool
```

Purpose:
- In-place NTT for XYW quotient-ring layout.
- Uses the standard modular multiplication path (no `use_montgomery` flag).

Input expectations:
- `p_power_of_3` must be divisible by `3`.
- `coefficient_values` length must be
  `2 * n_power_of_2 * n_power_of_2 * (2 * (p_power_of_3 / 3))`.
- Primitive roots of orders `4`, `4 * n_power_of_2`, and `p_power_of_3`
  must exist under `modulus_int`.

Returns:
- `True` on success.
- `False` if size/root preconditions are not satisfied.

## Supporting Internal Functions

These functions are internal building blocks used by the public entry points:

- `_multiply_modulo_int32(...)`
  - Chooses Barrett or Montgomery multiplication.
- `apply_i_axis_transform(...)`
  - Performs the i-axis split/merge transform.
- `apply_radix2_dif_ntt(...)`
  - Radix-2 DIF stage engine with configurable block size and base offset.
- `apply_radix3_dif_ntt(...)`
  - Radix-3 DIF stage engine.
- `apply_cyclotomic_pruned_ntt(...)`
  - Pruned cyclotomic transform used in W-axis processing.

## Tests and Benchmarks

Relevant tests:
- `test/test_reduction_modes.mojo`
- `test/test_ntt_xyw.mojo`

Relevant benchmarks:
- `benchmarks/bench_ntt_barrett.mojo`
- `benchmarks/bench_ntt_montgomery.mojo`
