# mojofhe
Mojo Implementation of FHE, vibe-coded

## Project Structure
This project follows a standard layout for Mojo development:

- `src/`: Contains the core Mojo source code for the project.
- `test/`: Holds unit tests and integration tests for the Mojo codebase.
- `examples/`: Provides example usage of the project's Mojo modules.
- `benchmarks/`: Contains performance benchmarks for critical sections of the code.
- `.venv/`: Python virtual environment for managing dependencies.

## Getting Started

### Prerequisites
- Python 3.8+
- Git

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/mojofhe.git
    cd mojofhe
    ```

2.  **Set up a Python virtual environment:**
    It's recommended to use a virtual environment to manage project dependencies.
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```

3.  **Install Mojo:**
    With your virtual environment activated, install Mojo using pip:
    ```bash
    pip install mojo
    ```
    (Note: If you encounter issues, refer to the official Mojo installation guide: [https://docs.modular.com/mojo/manual/install](https://docs.modular.com/mojo/manual/install))

### Hello World Example

Here's a simple "Hello, World!" example in Mojo.

1.  Create a file named `src/hello.mojo`:
    ```mojo
    fn main():
        print("Hello, Mojo!")
    ```

2.  **Run the example:**
    Make sure your virtual environment is activated.
    ```bash
    mojo src/hello.mojo
    ```
    Expected output:
    ```
    Hello, Mojo!
    ```

## Linting and Formatting
This project uses `mojo format` to ensure consistent code style. A pre-commit hook is set up to automatically format your staged `.mojo` files before you commit.

To enable the pre-commit hook, run the following command from the root of the project:
```bash
ln -s ../../lint.sh .git/hooks/pre-commit
```

Now, every time you commit a `.mojo` file, it will be automatically formatted according to the project's style.

## Development
(Add more sections for development guidelines, testing, contributing, etc., as the project evolves.)

## Documentation
- `docs/inverse-ntt.md`: inverse NTT APIs, round-trip tests, and NTT/INTT
  benchmark modes.
