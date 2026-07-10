# Shipped Heavy Matrix Contract

Wirehair's legacy wire format contains one fixed 6x18 GF(256) matrix. The
bytes are Cauchy-derived and must not change: `Generate_HeavyRows()` first
recreates them from seed `2318331135281` and compares every byte with the
matrix compiled into the codec.

The byte-reproduction check is the generator's compatibility guarantee. It is
not a reliability guarantee after binary Gaussian elimination mixes the left
columns into the right 6x6 corner. The old search claimed a perturbation-proof
"perfect" seed, but its rank check eliminated the wrong columns. With that bug
fixed, the shipped matrix becomes singular in roughly 1/256 random binary
perturbation trials, the expected floor for a random 6x6 GF(256) matrix.

`gen_tables --heavy-trials N` reports this singular rate as an empirical,
deterministic measurement. The result is informational and never certifies
all perturbations. Small samples may observe zero failures. Production keeps
the shipped bytes for wire compatibility and handles occasional rank loss by
accepting additional recovery rows.

The precode simulator's full-coverage MDS heavy patch is a separate modeling
assumption. Its results must not be described as a guarantee of the legacy
18-column production matrix.
