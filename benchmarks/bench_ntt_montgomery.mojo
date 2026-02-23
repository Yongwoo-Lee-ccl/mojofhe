from src.benchmark_utils import run_ixw_ntt_benchmark


fn main() raises:
    run_ixw_ntt_benchmark(
        True,
        "montgomery",
        "Montgomery benchmark setup failed.",
    )
