# ADR-0032: Libm reciprocal and root kernels (cot, sec, csc, cbrt, hypot, rootn) ship as direct primary kernels, not derived aliases

- **Status**: accepted
- **Date**: 2026-05-22

## Context

The Phase 1 correctness sweep work breakdown (when filed under
`docs/decisions/plans/`) proposed a single pre-freeze task: add
`cot`, `sec`, `csc`, `cbrt`, `hypot`, `rootn` as "trivial aliases
over already-present kernels" so the frozen unary surface looks
complete to a user scanning it against libm or MPFR. The framing
was "a few lines each."

The framing is recalled-shortcut math. `cot(x) = 1 / tan(x)` is
mathematically true; a correctly rounded `tan(x)` followed by a
correctly rounded reciprocal is *not* a correctly rounded `cot(x)`
for hard-to-round inputs. The two roundings compose into up to
1 ULP error in the reciprocal direction. The standard libm
references that ship correctly rounded transcendentals (CRlibm
under INRIA, Sun's `fdlibm`, the IBM Accurate Portable Math
Library) implement each of `cot`, `sec`, `csc` as a direct primary
kernel for exactly this reason: their domain reduction, polynomial
approximant, and Ziv ladder are all designed for the function
itself, not for a wrapped reciprocal of a sibling function.

The same shape recurs across the rest of the alias list:

- `sec(x) = 1 / cos(x)` and `csc(x) = 1 / sin(x)`: identical
  double-rounding hazard, same direct-kernel resolution in the
  reference libms.
- `hypot(x, y) = sqrt(x*x + y*y)`: the naive composition overflows
  or underflows where a direct kernel (Moler's `hypot`, IEEE
  754-2019 §9.2.1 recommended) does not. Even setting the overflow
  question aside, the squaring loses half the input precision and
  the composition cannot recover it.
- `cbrt(x) = pow(x, 1.0 / 3.0)`: `1 / 3` cannot be exactly
  represented in floating point, so the naive composition cannot
  be correctly rounded over the reals for any input. The reference
  libms use a direct table-driven Newton iteration with a
  specifically derived initial approximation.
- `rootn(x, n) = pow(x, 1.0 / n)`: same trap as `cbrt`, plus the
  domain question for even `n` and negative `x`. IEEE 754-2019
  §9.2 specifies it directly and it earns the spec entry by being
  a real kernel, not a wrapper.

This is the same recalled-shortcut trap pfloat just spent slice 8a
relearning twice (the `beta_case4` reciprocal product as `O(m)` in
caller supplied `m`; the parser cap framed as derived from
`i64::MAX / log2(10)` when in fact the i32-bounded parsed exponent
makes that threshold unreachable). The pattern is constant: a
recalled identity ("`B(-n,m) = (-1)^m / (m * C(n,m))`",
"`exponent fits the binary exponent type`",
"`cot(x) = 1 / tan(x)`") looks like a one-line implementation and
hides a property the spec or input domain demands.

Adding the aliases pre-freeze, even as a cosmetic surface-completion
step, enables a specific overclaim downstream. The exhaustive `f32`
sweep is sparse against the double-rounding boundaries (the boundary
inputs are countable in a way the `f32` grid does not directly
sample); a derived alias passes the sweep, is recorded as
`correctly-rounded` in the per-function status table, and is still
wrong against an independent direct kernel or in an `f64` evaluation
that hits the composition's boundary case. The status table is the
credibility document for the crate; an overclaim there is far more
costly than a documented gap.

## Decision

1. **`cot`, `sec`, `csc`, `cbrt`, `hypot`, `rootn` are not added to
   pfloat as aliases over existing kernels.** The pfloat surface for
   the v1.0 tag does not gain them via composition. The Phase 1
   correctness sweep's frozen unary surface lists them as
   `absent, deferred to the libm phase`, not as
   `correctly-rounded via composition`. The surface gap is honest;
   the alias would be cosmetic.

2. **Each lands in the libm phase (Phase 2 in the Track B roadmap)
   as a direct primary kernel.** Each derives from a cited primary
   source (DLMF §4.14 for the trig reciprocals' range reduction
   shape; IEEE 754-2019 §9.2.1 for `hypot`; the relevant Newton or
   table-driven approach for `cbrt` and `rootn`, with the initial
   approximation derived not recalled). Each pins worked examples
   from the source before kernel code. Each is verified against the
   Phase 1 oracle harness as the regression net.

3. **The Phase 1 status table carries a `kernel_kind` column** with
   values `primary | derived_alias`. The `correctly-rounded` verdict
   is structurally unavailable for `derived_alias` rows; the highest
   verdict an alias can earn is `faithful`. This is a belt-and-
   braces guard: a future contributor who reaches for the alias
   shortcut anyway will produce a row that cannot claim
   `correctly-rounded`, regardless of how the sweep returns. The
   column is documentation in the type system.

## Consequences

- The frozen Phase 1 surface has explicit gaps for these six
  functions. The status table makes the gaps visible rather than
  papering them over with composition-aliased rows. A user
  comparing pfloat to libm or MPFR reads "not yet shipped" instead
  of "shipped, claimed correct, actually faithful." That is the
  disclosure-standards-consistent surface (per the global
  CLAUDE.md disclosure rules: no overclaiming on rounding rigor;
  the per-function table is the load-bearing artifact and is only
  as trustworthy as its least-honest row).

- The libm phase has six more kernels to write directly. Each
  derives from a cited source and verifies against the harness.
  This is more work than a one-line
  `pub fn cot(&self, mode) -> (Self, Status) { ... }` aliased over
  `tan` would be. The work earns the claim.

- **Future-reviewer reminder.** When the libm phase begins (Phase 2
  in the Track B roadmap), this ADR is the gate against the
  trivial-alias reflex. Anyone proposing a one-line wrapper such
  as `1 / self.tan(mode).0` for `cot`, or analogous compositions
  for `sec`, `csc`, `cbrt`, `hypot`, `rootn`, is steered to this
  ADR and the direct-kernel requirement. The libm phase's first
  commit should be a per-function kernel-list document that
  explicitly cites ADR-0032 against each of the six functions as
  `direct kernel required, not aliased; see ADR-0032 for the
  double-rounding rationale`.

- This is the same lesson as the slice 8a fork resolutions,
  generalized to the API surface: a recalled mathematical identity
  is not a correct implementation, and the cheap composition is
  rarely the correctly rounded one for adversarial inputs. The
  `derive don't recall` discipline applies at the kernel surface,
  not only at the coefficient and recurrence level. ADRs 0027
  through 0031 each closed a recalled shortcut found within an
  implementation; this one closes the recalled shortcut at the API
  composition level before it ships.

## Related

- The Phase 1 correctness sweep work breakdown (to be filed under
  `docs/decisions/plans/` when work on Phase 1 begins); the
  pre-freeze "trivial aliases" step is the one this ADR rules out.
- The Track B numerics ecosystem roadmap, Phase 2 (libm spinoff):
  the home for the six direct kernels this ADR defers there.
- ADR-0022 (`pow` Ziv interval test). The closest cousin:
  the original `pow` Ziv "recompute and compare" composition
  false-converged on hard-to-round inputs and was replaced by the
  interval test. The same property the alias trap exhibits.
- Slice 8a forks: `beta_case4` reciprocal product as `O(m)`
  (replaced by the O(1) `lgamma` factorial form, ADR-0030
  Correction); parser cap as `i64::MAX / log2(10)` (replaced by
  the pow5 storage budget, ADR-0031). The same recalled-shortcut
  pattern at the algorithm level; this ADR is the API-surface
  analog.
- Primary references for the direct-kernel posture: DLMF §4.14
  (trig reciprocals); IEEE 754-2019 §9.2.1 (`hypot`), §9.2
  (`rootn`); CRlibm (INRIA) and Sun `fdlibm` as the
  state-of-the-art ships of `cot`, `sec`, `csc`, `cbrt` as direct
  primary kernels with the double-rounding reasoning documented.
