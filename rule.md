# Contribution Rules

This document outlines the rules for contributing to this project. Please read these rules carefully before contributing.

## 1. Linting
Always run the linter before committing your changes. This is enforced by a pre-commit hook, but you can also run it manually:
```bash
./lint.sh
```
This is especially important before creating a pull request.

The maximum line length is 80 characters (ruler length: 80). Ensure your code adheres to this limit for better readability.

## 2. Testing
Always run the tests and verify that all tests pass before creating a pull request.
```bash
./run_tests.sh
```

## 3. New Features
If you are adding a new feature, you must also add corresponding tests.

## 4. Variable Naming
All variable names should be self-explaining. Do not use single-letter variable names like `a`, `b`, or even `i` for indexing. Instead, use descriptive names like `slot_index`.

## 5. DRY (Do Not Repeat Yourself)
Avoid re-implementing existing logic. Before adding new helper functions, check whether equivalent functions already exist in `src/` and reuse them.

If duplicate logic is discovered during review, refactor to a shared function/module rather than keeping multiple copies.
