"""ADR-0035 Oracle worker shared package.

Holds the library-agnostic ``certified_round_f32`` routine and the
hand-derived tricky-case corpus. Each oracle worker (Arb, mpmath,
Maxima) imports from this package; the package itself has no
dependency on FLINT/Arb, mpmath, or Maxima.

The separation is the load-bearing verification handle: the rounding
routine and the corpus can be exercised without invoking any function
library, so a bug in the routine cannot be confused with a bug in the
ball arithmetic of any particular oracle.
"""
