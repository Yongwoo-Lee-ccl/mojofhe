#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Check if we are in the virtual environment
if [ -z "$VIRTUAL_ENV" ]; then
    echo "Not in a virtual environment. Please activate it first."
    echo "source .venv/bin/activate"
    exit 1
fi

echo "Running all tests in the 'test' folder..."

for test_file in test/*; do
    if [ -f "$test_file" ]; then
        echo "Running $test_file..."
        mojo run -I . "$test_file"
    fi
done

echo "All tests passed."
