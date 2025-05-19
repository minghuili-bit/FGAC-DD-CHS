## Credits and Source

Portions of this project are derived from the [JEDI library](https://github.com/ucbrise/jedi-pairing), originally developed by Sam Kumar and collaborators. The core pairing implementation, WKD-IBE primitives, and C/Go wrappers are based on that library and used here under the BSD 3-Clause License.

For reference, the original JEDI project README is available [here](https://github.com/ucbrise/jedi-pairing#readme).

## Modifications Summary

### 1. Platform Compatibility

- The original project was designed for **Linux** systems; we adapted the build configuration and assembly code to enable compilation and execution on **macOS**

- Modified `.s` files in `src/core/arch/x86_64/` to comply with macOS assembly syntax, such as replacing:

  ```asm
  .globl my_function         →  .globl _my_function
  my_function:               →  _my_function: