#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Check if we are in the virtual environment
if [ -z "$VIRTUAL_ENV" ]; then
    echo "Not in a virtual environment. Please activate it first."
    echo "source .venv/bin/activate"
    exit 1
fi

echo "Running Mojo formatter..."

staged_mojo_files=()
while IFS= read -r file; do
    staged_mojo_files+=("$file")
done < <(git diff --cached --name-only --diff-filter=ACM | grep '\.mojo$' || true)

if [ ${#staged_mojo_files[@]} -gt 0 ]; then
    echo "Formatting staged .mojo files..."
    printf "%s\n" "${staged_mojo_files[@]}" | xargs mojo format
    echo "Adding formatted files to the commit..."
    printf "%s\n" "${staged_mojo_files[@]}" | xargs git add
fi

exit 0
