#!/bin/bash

set -e

if [ -z "$VIRTUAL_ENV" ]; then
    echo "Not in a virtual environment. Please activate it first."
    echo "source .venv/bin/activate"
    exit 1
fi

echo "Benchmarking NTT reduction modes (n=2^7, p=3^5)"
for run_index in 1 2 3; do
    echo "barrett run $run_index"
    /usr/bin/time -f "elapsed %e s" \
        mojo run -I . benchmarks/bench_ntt_barrett.mojo
done

for run_index in 1 2 3; do
    echo "montgomery run $run_index"
    /usr/bin/time -f "elapsed %e s" \
        mojo run -I . benchmarks/bench_ntt_montgomery.mojo
done
