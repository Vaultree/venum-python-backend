# VENum Python backend

Welcome to the Python backend for VENum! VENum stands for Vaultree Encrypted Numbers and is Vaultree's homomorphic encryption library. This open-source project provides a straightforward Python implementation of homomorphic encryption, designed for ease of use.

## General Features and goals

- **Homomorphic Encryption**: Perform computations on encrypted data without decrypting it, ensuring data privacy and security.
- **Python Backend**: A simple and accessible Python-based implementation for developers and researchers.
- **Seamless Interchangeability**: The current code is a low-level API for homomorphic encryption and while it can be used directly, the intended use is through our VENumpy API which allows the Python backend to be swapped for Vaultree's closed-source Rust implementation without modifying your code.
- **Flexible Licensing**: The Python backend is open-source and free to use, while the Rust backend offers significant performance improvements for real-world scenarios and requires a license to access.

## Acknowledgments

This project is part of Vaultree's mission to make privacy-preserving computation accessible and practical for everyone.

## Homomorphic Encryption Features

### Current
- Addition of ciphers
- Subtraction of ciphers

### Upcoming
- Multiplication of ciphers

### Future Direction
We aim to provide all features available in the Rust backend. They are described in detail in the following Vaultree research papers:  
- [Efficient and Practical Homomorphic Encryption Framework](https://eprint.iacr.org/2024/1105.pdf)  
- [Advanced Cryptographic Techniques for Scalable Privacy](https://eprint.iacr.org/2024/1622.pdf)  

## Usage

### Dependencies

Before using the `justfile`, ensure you have the following installed:

- [Just command runner](https://github.com/casey/just)
- Python 3
- Container Engine: `docker`, `podman` or another compatible container engine.

Run `just setup` to automatically set up the required dependencies in a virtual environment.

### Recipes
  
- **`just setup`**: Sets up the Python virtual environment, installs the project locally in editable mode, and installs required dependencies for building packages.

- **`just test [path='./']`**: Runs tests in the specified directory (default: `tests/`).

- **`just build`**: Builds a Python wheel package for the project.

- **`just build-container`**: Builds a container image for the project using the container engine (`podman` by default). Extracts the build artifacts into the `dist_container/` directory.

- **`just clean`**: Cleans up build artifacts, including the `dist/`, `dist_container/`, and `build/` directories, and any `.egg-info` files.

## Documentation

For detailed usage examples and understanding of the library's functionality, please refer to the following resources:

- **Test Cases**: Explore the `tests/` directory for a variety of test cases that demonstrate how to use the library's features in practical scenarios. These tests provide concrete examples of encrypting, performing operations on encrypted data, and decrypting results.

- **Docstrings**: The codebase includes docstrings that document functions and classes. To view the docstrings interactively, use the Python REPL or your preferred IDE to inspect functions and classes.

The documentation is an ongoing effort and will continue to evolve over time. However, the codebase has been thoughtfully designed with clarity and usability in mind, ensuring that it is intuitive and straightforward for developers to work with, even in the absence of instructions.

We encourage users to contribute to improving both the documentation and test cases to enhance clarity and usability for the community.

## License

Please check the [License file](LICENSE.md).

## Support

For issues or questions, please open an issue on [GitHub](https://github.com/Vaultree/venum-python-backend) or contact us at [our support page](https://support.vaultree.com).
