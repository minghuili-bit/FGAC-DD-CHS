## Credits and Source

Portions of this project are derived from the [JEDI library](https://github.com/ucbrise/jedi-pairing), originally developed by Sam Kumar and collaborators. The core pairing implementation, WKD-IBE primitives, and C/Go wrappers are based on that library and used here under the BSD 3-Clause License.

For reference, the original JEDI project README is available [here](https://github.com/ucbrise/jedi-pairing#readme).

## Modifications Summary

### 1. Platform Compatibility

The original project was designed for **Linux** systems; we adapted the build configuration and assembly code to enable compilation and execution on **macOS**

### 2. Code Modifications

We extended the core WKD-IBE logic based on the original research implementation. The modifications were primarily made to the following files: /include/wkdibe/api.hpp, /src/wkdibe/api.cpp, and tests/test_integration.cpp.