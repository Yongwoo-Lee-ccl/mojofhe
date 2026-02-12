# Contribution Rules

This document outlines the rules for contributing to this project. Please read these rules carefully before contributing.

## 1. Linting
Always run the linter before committing your changes. This is enforced by a pre-commit hook, but you can also run it manually:
```bash
./lint.sh
```
This is especially important before creating a pull request.

## 2. Testing
Always run the tests and verify that all tests pass before creating a pull request.
```bash
./run_tests.sh
```

## 3. New Features
If you are adding a new feature, you must also add corresponding tests.

## 4. Variable Naming
All variable names should be self-explaining. Do not use single-letter variable names like `a`, `b`, or even `i` for indexing. Instead, use descriptive names like `slot_index`.
