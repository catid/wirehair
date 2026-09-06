#ifndef WH2_THUE_MORSE_MATVEC_GFNI_R0_H
#define WH2_THUE_MORSE_MATVEC_GFNI_R0_H

#include <cstdint>

// Benchmark-only fixed-six arithmetic. Initialize the shared GF256 runtime
// before calling any entry point; no entry point initializes or allocates it.
namespace wh2_matvec_gfni_r0 {

// True only when this build has the target helper and the initialized runtime
// enables both GFNI (including its AVX-512/OS prerequisites) and AVX2.
bool Available();

// Replace vector with matrix * vector in Wirehair's GF(256), polynomial 0x14d.
// matrix is exactly 36 readable row-major bytes; vector is exactly six readable
// and writable bytes. Both may have arbitrary alignment and may overlap.
// All input reads complete before the final six-byte write. These fixed-size
// primitives do not validate pointers, dimensions, capacities or initialization.
void Apply(const std::uint8_t* matrix, std::uint8_t* vector);

// Exact original scalar operand order, also satisfying the overlap contract.
void ApplyScalar(const std::uint8_t* matrix, std::uint8_t* vector);

} // namespace wh2_matvec_gfni_r0

#endif
