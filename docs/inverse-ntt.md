# Inverse NTT and Multiplication Validation

This project now includes inverse IXW/XYW quotient transforms and round-trip
validation tests.

## New APIs

### IXW inverse transform

```mojo
fn apply_ixw_quotient_intt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_value: Int,
    use_montgomery: Bool = False,
) -> Bool
```

- In-place inverse of `apply_ixw_quotient_ntt_inplace`.
- Returns `False` when roots or input-size preconditions are not satisfied.
- Supports two modular arithmetic modes:
  - Barrett: `use_montgomery=False`
  - Montgomery: `use_montgomery=True`

### XYW inverse transform

```mojo
fn apply_xyw_quotient_intt_inplace(
    mut coefficient_values: List[Int32],
    n_power_of_2: Int,
    p_power_of_3: Int,
    modulus_int: Int,
) -> Bool
```

- In-place inverse of `apply_xyw_quotient_ntt_inplace`.

## Tests

`test/test_ntt_inverse.mojo` adds:

- Round-trip checks (`NTT -> INTT`) over small grids:
  - IXW: `n in {1,2,4}`, `p in {3,9}`
  - XYW: `n in {2,4}`, `p in {3,9}`
- A manual helper for multiplication checks against a naive quotient model
  is included but not part of the default automated test run yet.

## Benchmark

`benchmarks/bench_ntt_intt_modes.mojo` benchmarks full `NTT + INTT` cycles for:

- Barrett reduction
- Montgomery reduction

Run via:

```bash
source .venv/bin/activate
bash benchmarks/run_reduction_bench.sh
```
