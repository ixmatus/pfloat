//! [`Mag`]: an unsigned, upward-rounded magnitude — the ball radius.
//!
//! `Mag` is a single-limb binary float `m · 2^(e − 63)` with a `u64`
//! mantissa `m` (top bit set, so `m ∈ [2^63, 2^64)`) and an `i64`
//! exponent `e`, plus the two special values `0` and `+∞`. The exponent
//! width matches pfloat's, so a radius exponent and a midpoint exponent
//! compose without a second overflow regime. ADR-0074.
//!
//! # Why the type carries the soundness
//!
//! A ball radius must never be negative and must never round *inward*
//! (toward zero), or the enclosure it certifies becomes a lie. `Mag`
//! makes both unrepresentable:
//!
//! - The type is `{ finite ≥ 0, +∞ }` — no sign field and no NaN
//!   variant — so a negative or not-a-number radius cannot be written.
//! - Every arithmetic operation ([`add`](Mag::add), [`mul`](Mag::mul))
//!   and every conversion into `Mag` ([`from_bigfloat_ceil`](Mag::from_bigfloat_ceil))
//!   rounds the result *up* (toward `+∞`) to the single-limb mantissa,
//!   so an inward-rounded radius cannot be written either.
//!
//! # Mantissa width and the resolution cap
//!
//! A single 64-bit limb dominates Arb's 30-bit `mag_t` on both axes: it
//! gives markedly tighter radii, and being `Copy`, alloc-free, and
//! `Vec`-free is exactly what lets the round-up invariants discharge
//! under Kani. The one documented consequence: at midpoint precision
//! above 64 bits the radius has a relative resolution floor near
//! `2^-64`. This is sound because the radius is only ever an
//! upward-rounded *upper* bound, never an equality — the midpoint
//! carries the precision, the radius carries the certified slack, so the
//! radius is not the accuracy bottleneck in practice.

#[cfg(feature = "big")]
use pfloat::{BigFloat, Parts, RoundingMode, Sign};

/// An unsigned, upward-rounded magnitude: `0`, a finite `m · 2^(e − 63)`,
/// or `+∞`.
///
/// `Mag` is the radius half of a [`Ball`](crate::Ball). It is totally
/// ordered (`0 < finite < +∞`, finites by value), `Copy`, and never
/// negative or NaN by construction. Every operation rounds toward `+∞`.
///
/// The `Finite` fields are ordered `exponent` then `mantissa` so the
/// derived [`Ord`] is value order on the canonical (top-bit-set) form.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Mag {
    /// The exact magnitude zero (a radius of `0` is an *exact* ball).
    Zero,
    /// A finite positive magnitude `mantissa · 2^(exponent − 63)`.
    ///
    /// Invariant: `mantissa` has its top bit set (`mantissa ≥ 2^63`), so
    /// the representation is unique and `exponent` is `floor(log2 value)`.
    Finite {
        /// Binary exponent: `floor(log2(value))`.
        exponent: i64,
        /// Top-bit-set 64-bit significand.
        mantissa: u64,
    },
    /// An unbounded magnitude (`+∞`): a radius known only to be finite is
    /// not representable, so an over-the-ceiling or indeterminate result
    /// saturates here.
    Infinity,
}

impl Mag {
    /// The zero magnitude.
    pub const ZERO: Mag = Mag::Zero;

    /// The unbounded magnitude `+∞`.
    pub const INFINITY: Mag = Mag::Infinity;

    /// The exact magnitude `2^k`.
    #[must_use]
    pub const fn from_pow2(k: i64) -> Mag {
        Mag::Finite {
            exponent: k,
            mantissa: 1u64 << 63,
        }
    }

    /// `true` for the zero magnitude.
    #[must_use]
    pub const fn is_zero(self) -> bool {
        matches!(self, Mag::Zero)
    }

    /// `true` for `+∞`.
    #[must_use]
    pub const fn is_infinite(self) -> bool {
        matches!(self, Mag::Infinity)
    }

    /// `true` for `0` or a finite magnitude.
    #[must_use]
    pub const fn is_finite(self) -> bool {
        !matches!(self, Mag::Infinity)
    }

    /// Sum, rounded up: the smallest representable `Mag` greater than or
    /// equal to the true sum `self + rhs`.
    ///
    /// `+∞` absorbs; `0` is the identity. The soundness contract is
    /// `result ≥ self + rhs` exactly.
    // Deliberately not `core::ops::Add`: `+` implies exact addition, but
    // this rounds the result up. Naming it `add` keeps call sites in the
    // ball kernels readable (`a.rad.add(b.rad)`).
    #[allow(clippy::should_implement_trait)]
    #[must_use]
    pub fn add(self, rhs: Mag) -> Mag {
        match (self, rhs) {
            (Mag::Infinity, _) | (_, Mag::Infinity) => Mag::Infinity,
            (Mag::Zero, x) | (x, Mag::Zero) => x,
            (
                Mag::Finite {
                    exponent: ea,
                    mantissa: ma,
                },
                Mag::Finite {
                    exponent: eb,
                    mantissa: mb,
                },
            ) => {
                // Order so the larger exponent is `ea`.
                let (ea, ma, eb, mb) = if ea >= eb {
                    (ea, ma, eb, mb)
                } else {
                    (eb, mb, ea, ma)
                };

                // Frame: place `ma` at bit 63 so `acc = value · 2^(126 − ea)`
                // with 63 guard bits below the mantissa LSB. `<<63` (not
                // `<<64`) keeps `acc` below `2^128` even when both operands
                // have equal exponents and full mantissas.
                let acc_hi = (ma as u128) << 63;
                let shift = i128::from(ea) - i128::from(eb); // ≥ 0

                let (contribution, sticky) = if shift >= 127 {
                    // `mb` falls entirely below `acc`'s least bit.
                    (0u128, true)
                } else {
                    let shift = shift as u32; // 0..=126
                    if shift <= 63 {
                        ((mb as u128) << (63 - shift), false)
                    } else {
                        let down = shift - 63; // 1..=63
                        let lost = (mb & ((1u64 << down) - 1)) != 0;
                        ((mb as u128) >> down, lost)
                    }
                };

                let acc = acc_hi + contribution;
                pack_round_up(acc, sticky, i128::from(ea) - 126)
            }
        }
    }

    /// Product, rounded up: the smallest representable `Mag` greater than
    /// or equal to the true product `self · rhs`.
    ///
    /// `+∞` absorbs (including `0 · ∞`, which saturates to `+∞` rather
    /// than claiming an exact zero against an unbounded factor); `0`
    /// otherwise annihilates. The soundness contract is
    /// `result ≥ self · rhs` exactly.
    // Deliberately not `core::ops::Mul`: `*` implies exact multiplication,
    // but this rounds the result up.
    #[allow(clippy::should_implement_trait)]
    #[must_use]
    pub fn mul(self, rhs: Mag) -> Mag {
        match (self, rhs) {
            // `+∞` absorbs first, so `0 · ∞` saturates to `+∞` (the sound
            // over-estimate of an indeterminate radius product).
            (Mag::Infinity, _) | (_, Mag::Infinity) => Mag::Infinity,
            (Mag::Zero, _) | (_, Mag::Zero) => Mag::Zero,
            (
                Mag::Finite {
                    exponent: ea,
                    mantissa: ma,
                },
                Mag::Finite {
                    exponent: eb,
                    mantissa: mb,
                },
            ) => {
                let product = (ma as u128) * (mb as u128); // ∈ [2^126, 2^128)
                pack_round_up(product, false, i128::from(ea) + i128::from(eb) - 126)
            }
        }
    }

    /// The larger of two magnitudes (exact; no rounding).
    #[must_use]
    pub fn max(self, rhs: Mag) -> Mag {
        core::cmp::Ord::max(self, rhs)
    }

    /// The smallest `Mag` greater than or equal to `|x|`: the upward
    /// narrowing of a full-precision pfloat magnitude into a single-limb
    /// radius.
    ///
    /// `±0 → 0`. A non-finite `x` (`±∞` or NaN) maps to `+∞`, the sound
    /// over-estimate of an unbounded or indeterminate magnitude. A finite
    /// non-zero `x` keeps its top 64 mantissa bits and rounds up by one
    /// unit in the last place if any lower bit is set.
    #[cfg(feature = "big")]
    #[must_use]
    pub fn from_bigfloat_ceil(x: &BigFloat) -> Mag {
        match x.parts() {
            Parts::Zero { .. } => Mag::Zero,
            Parts::Infinity { .. } | Parts::Nan { .. } => Mag::Infinity,
            Parts::Normal {
                exponent, mantissa, ..
            } => {
                // Canonical storage is little-endian, top-bit-set in the
                // most-significant limb, low (limbs*64 − precision) bits
                // zero. The top 64 bits of |x| are the top limb; |x| =
                // top · 2^(exponent − 63) + (the lower limbs' tail).
                let limbs = mantissa.len();
                let top = mantissa[limbs - 1];
                let tail_nonzero = mantissa[..limbs - 1].iter().any(|&l| l != 0);
                if tail_nonzero {
                    // One ulp up bounds the tail: tail < 2^(exponent−63),
                    // exactly the weight of the mantissa LSB.
                    let (m1, carry) = top.overflowing_add(1);
                    if carry {
                        match exponent.checked_add(1) {
                            Some(e) => Mag::Finite {
                                exponent: e,
                                mantissa: 1u64 << 63,
                            },
                            None => Mag::Infinity,
                        }
                    } else {
                        Mag::Finite {
                            exponent,
                            mantissa: m1,
                        }
                    }
                } else {
                    // Exact: the top limb already is |x| (precision ≤ 64, or
                    // a wider mantissa whose tail is all zero).
                    Mag::Finite {
                        exponent,
                        mantissa: top,
                    }
                }
            }
        }
    }

    /// Exact conversion to a 64-bit [`BigFloat`].
    ///
    /// `0 → +0`, `+∞ → +∞`, and a finite `m · 2^(e − 63)` converts with
    /// no rounding (a `Mag` mantissa fits in 64 bits, so the value is
    /// exactly representable at precision 64). Used to subtract/add a
    /// radius against a midpoint for the exact `lower`/`upper` endpoints.
    #[cfg(feature = "big")]
    #[must_use]
    pub fn to_bigfloat(self) -> BigFloat {
        match self {
            Mag::Zero => BigFloat::try_new_zero(Sign::Positive, 64).expect("64 ≥ 1"),
            Mag::Infinity => BigFloat::try_new_infinity(Sign::Positive, 64).expect("64 ≥ 1"),
            Mag::Finite { exponent, mantissa } => {
                // Reassemble the 64-bit mantissa exactly from two halves
                // (each ≤ 32 bits, so each fits an i64 and precision 64),
                // then scale by 2^(exponent − 63).
                let hi = (mantissa >> 32) as i64;
                let lo = (mantissa & 0xFFFF_FFFF) as i64;
                let hi_bf = BigFloat::try_from_i64_exact(hi, 64).expect("≤ 32 bits fits p=64");
                let (hi_scaled, _) = hi_bf.scale_by_pow2(32);
                let lo_bf = BigFloat::try_from_i64_exact(lo, 64).expect("≤ 32 bits fits p=64");
                // hi·2^32 + lo is a 64-bit integer: exact at precision 64.
                let (m_bf, _) = hi_scaled.add(&lo_bf, RoundingMode::NearestEven);
                // exponent − 63 saturates only for radii within 63 of the
                // i64 exponent floor (unreachable at realistic scales).
                let (scaled, _) = m_bf.scale_by_pow2(exponent.saturating_sub(63));
                scaled
            }
        }
    }
}

/// Round a `u128` accumulator up to a single-limb `Mag`.
///
/// `acc ∈ [2^126, 2^128)` is the magnitude scaled so that its bit 0 has
/// weight `2^value_lsb_exp` (the true value is `acc · 2^value_lsb_exp`).
/// `sticky` carries any nonzero bits already dropped below bit 0. The
/// result is the smallest `Mag ≥` the true value: the top 64 bits of
/// `acc`, incremented by one ulp if any lower bit (or `sticky`) is set.
fn pack_round_up(acc: u128, sticky: bool, value_lsb_exp: i128) -> Mag {
    debug_assert!(
        acc >= (1u128 << 126),
        "accumulator must be normalized to the top two bits"
    );
    let top_bit = 127 - acc.leading_zeros(); // 126 or 127
    let drop = top_bit - 63;
    let mut mantissa = (acc >> drop) as u64; // top 64 bits, top-bit-set
    let remainder = acc & ((1u128 << drop) - 1);
    let mut exponent = i128::from(top_bit) + value_lsb_exp; // E = top_bit + value_lsb_exp

    if remainder != 0 || sticky {
        let (m1, carry) = mantissa.overflowing_add(1);
        if carry {
            // Mantissa was all ones: 1 ulp carries it to 2^64, renormalize
            // to the top bit at the next exponent.
            mantissa = 1u64 << 63;
            exponent += 1;
        } else {
            mantissa = m1;
        }
    }

    if exponent > i128::from(i64::MAX) {
        Mag::Infinity
    } else {
        // Underflow past the i64 floor saturates the exponent down, which
        // keeps the mantissa and so over-estimates (still sound). This is
        // unreachable for any radius at realistic exponent scales.
        let exponent = if exponent < i128::from(i64::MIN) {
            i64::MIN
        } else {
            exponent as i64
        };
        Mag::Finite { exponent, mantissa }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOP: u64 = 1u64 << 63;

    #[test]
    fn ordering_across_variants() {
        assert!(Mag::Zero < Mag::from_pow2(-100));
        assert!(Mag::from_pow2(100) < Mag::Infinity);
        assert!(Mag::Zero < Mag::Infinity);
        // Finites order by value: larger exponent wins, then larger mantissa.
        assert!(Mag::from_pow2(3) < Mag::from_pow2(4));
        assert!(
            Mag::Finite {
                exponent: 5,
                mantissa: TOP
            } < Mag::Finite {
                exponent: 5,
                mantissa: TOP + 1
            }
        );
    }

    #[test]
    fn from_pow2_is_a_power_of_two() {
        let m = Mag::from_pow2(7);
        assert_eq!(
            m,
            Mag::Finite {
                exponent: 7,
                mantissa: TOP
            }
        );
    }

    #[test]
    fn infinity_absorbs() {
        let one = Mag::from_pow2(0);
        assert_eq!(one.add(Mag::Infinity), Mag::Infinity);
        assert_eq!(Mag::Infinity.add(one), Mag::Infinity);
        assert_eq!(one.mul(Mag::Infinity), Mag::Infinity);
        assert_eq!(Mag::Infinity.mul(one), Mag::Infinity);
        // 0 · ∞ is the conservative +∞ (never a false exact zero).
        assert_eq!(Mag::Zero.mul(Mag::Infinity), Mag::Infinity);
        assert_eq!(Mag::Infinity.mul(Mag::Zero), Mag::Infinity);
    }

    #[test]
    fn zero_is_additive_identity_and_multiplicative_annihilator() {
        let x = Mag::from_pow2(5);
        assert_eq!(x.add(Mag::Zero), x);
        assert_eq!(Mag::Zero.add(x), x);
        assert_eq!(x.mul(Mag::Zero), Mag::Zero);
        assert_eq!(Mag::Zero.mul(x), Mag::Zero);
    }

    #[test]
    fn add_exact_power_sum() {
        // 2^0 + 2^0 = 2^1, exactly representable.
        assert_eq!(Mag::from_pow2(0).add(Mag::from_pow2(0)), Mag::from_pow2(1));
        // 2^10 + 2^2 = mantissa with two set bits, exact at 64-bit width.
        let s = Mag::from_pow2(10).add(Mag::from_pow2(2));
        match s {
            Mag::Finite { exponent, mantissa } => {
                assert_eq!(exponent, 10);
                // 2^10 + 2^2 = 1.0000000001b × 2^10; mantissa top bit + bit (63-8).
                assert_eq!(mantissa, TOP | (1u64 << (63 - 8)));
            }
            _ => panic!("expected finite"),
        }
    }

    #[test]
    fn mul_exact_powers() {
        // 2^3 · 2^4 = 2^7.
        assert_eq!(Mag::from_pow2(3).mul(Mag::from_pow2(4)), Mag::from_pow2(7));
        // 3 · 5 = 15 (mantissas with set low bits), exact at 64-bit width.
        let three = Mag::Finite {
            exponent: 1,
            mantissa: 3u64 << 62,
        };
        let five = Mag::Finite {
            exponent: 2,
            mantissa: 5u64 << 61,
        };
        let fifteen = three.mul(five);
        // 15 = 1.111b × 2^3.
        assert_eq!(
            fifteen,
            Mag::Finite {
                exponent: 3,
                mantissa: 15u64 << 60
            }
        );
    }

    #[test]
    fn add_rounds_up_when_inexact() {
        // (2^64 − 1) at exponent 0, plus a value far below: the tail is
        // entirely lost, so the result must round UP, never down.
        let big = Mag::Finite {
            exponent: 0,
            mantissa: u64::MAX,
        };
        let tiny = Mag::from_pow2(-200);
        let s = big.add(tiny);
        // The true sum is just above `big`; rounding up carries the
        // all-ones mantissa to 2^1.
        assert!(s > big, "rounding must go up");
        assert_eq!(s, Mag::from_pow2(1));
    }

    #[test]
    fn mul_rounds_up_when_inexact() {
        // Two mantissas whose product needs more than 64 bits.
        let a = Mag::Finite {
            exponent: 0,
            mantissa: u64::MAX,
        };
        let b = Mag::Finite {
            exponent: 0,
            mantissa: u64::MAX,
        };
        let p = a.mul(b);
        // (2 − 2^-63)^2 ≈ 4 − …; the exact product is below 4 but rounds
        // up. The result must be ≥ the true product.
        // True product = (2^64-1)^2 · 2^-126 = (2^128 - 2^65 + 1)·2^-126.
        // Rounded up to 64-bit mantissa.
        match p {
            Mag::Finite { exponent, .. } => assert!(exponent == 1),
            _ => panic!("finite"),
        }
        // Soundness checked exactly in the big-feature tests below.
    }
}

#[cfg(all(test, feature = "big"))]
mod soundness_tests {
    use super::*;

    // Exact high-precision oracle: convert a Mag to BigFloat and compare.
    const ORACLE_PREC: u32 = 300;

    fn mag_finite(exponent: i64, mantissa: u64) -> Mag {
        assert!(
            mantissa & (1u64 << 63) != 0,
            "test mantissa must be normalized"
        );
        Mag::Finite { exponent, mantissa }
    }

    fn ge(result: &BigFloat, truth: &BigFloat) -> bool {
        // result >= truth
        !matches!(result.partial_cmp(truth).0, Some(core::cmp::Ordering::Less))
    }

    // The exact sum/product of two 64-bit-mantissa values with a bounded
    // exponent gap fits in ORACLE_PREC bits, so the BigFloat oracle is
    // exact and the >= check is a true soundness assertion.
    fn exact_add(a: Mag, b: Mag) -> BigFloat {
        let (s, st) = a
            .to_bigfloat()
            .add_round(&b.to_bigfloat(), ORACLE_PREC, RoundingMode::NearestEven)
            .unwrap();
        assert!(
            !st.inexact(),
            "oracle add must be exact for the test inputs"
        );
        s
    }

    fn exact_mul(a: Mag, b: Mag) -> BigFloat {
        let (p, st) = a
            .to_bigfloat()
            .mul_round(&b.to_bigfloat(), ORACLE_PREC, RoundingMode::NearestEven)
            .unwrap();
        assert!(st.is_ok(), "oracle mul must be exact for the test inputs");
        p
    }

    #[test]
    fn add_never_under_estimates() {
        // A spread of mantissas and bounded exponent gaps. The exact sum
        // of two 64-bit values with gap ≤ ~120 fits in 300 bits.
        let mantissas = [
            1u64 << 63,
            (1u64 << 63) | 1,
            (1u64 << 63) | (1u64 << 40),
            u64::MAX,
            0xC000_0000_0000_0001,
            0xFFFF_FFFF_0000_0001,
        ];
        for &ma in &mantissas {
            for &mb in &mantissas {
                for &ea in &[-50i64, -1, 0, 7, 64, 130] {
                    for &eb in &[-50i64, -1, 0, 7, 64, 130] {
                        if (ea - eb).abs() > 120 {
                            continue;
                        }
                        let a = mag_finite(ea, ma);
                        let b = mag_finite(eb, mb);
                        let got = a.add(b);
                        let truth = exact_add(a, b);
                        assert!(
                            ge(&got.to_bigfloat(), &truth),
                            "add under-estimated: a={a:?} b={b:?} got={got:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn mul_never_under_estimates() {
        let mantissas = [
            1u64 << 63,
            (1u64 << 63) | 1,
            (1u64 << 63) | (1u64 << 40),
            u64::MAX,
            0xC000_0000_0000_0001,
            0xFFFF_FFFF_0000_0001,
        ];
        for &ma in &mantissas {
            for &mb in &mantissas {
                for &ea in &[-50i64, -1, 0, 7, 64, 130] {
                    for &eb in &[-30i64, -1, 0, 11, 90] {
                        let a = mag_finite(ea, ma);
                        let b = mag_finite(eb, mb);
                        let got = a.mul(b);
                        let truth = exact_mul(a, b);
                        assert!(
                            ge(&got.to_bigfloat(), &truth),
                            "mul under-estimated: a={a:?} b={b:?} got={got:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn add_is_minimal_over_estimate_when_inexact() {
        // When the result is finite and inexact, the next-smaller Mag must
        // be strictly below the true sum (tightness, not just soundness).
        let a = mag_finite(0, u64::MAX);
        let b = Mag::from_pow2(-10);
        let got = a.add(b);
        let truth = exact_add(a, b);
        assert!(ge(&got.to_bigfloat(), &truth));
        // got is the rounded-up neighbour; predecessor (got with mantissa-1)
        // must be < truth.
        if let Mag::Finite { exponent, mantissa } = got {
            let pred = if mantissa == (1u64 << 63) {
                Mag::Finite {
                    exponent: exponent - 1,
                    mantissa: u64::MAX,
                }
            } else {
                Mag::Finite {
                    exponent,
                    mantissa: mantissa - 1,
                }
            };
            assert!(
                pred.to_bigfloat().partial_cmp(&truth).0 == Some(core::cmp::Ordering::Less),
                "over-estimate is not minimal"
            );
        }
    }

    #[test]
    fn from_bigfloat_ceil_is_sound_and_exact_when_it_can_be() {
        // p ≤ 64: exact. Build 3 at precision 53.
        let three = BigFloat::try_from_i64_exact(3, 53).unwrap();
        let m = Mag::from_bigfloat_ceil(&three);
        assert!(ge(&m.to_bigfloat(), &three));
        assert_eq!(
            m.to_bigfloat().partial_cmp(&three).0,
            Some(core::cmp::Ordering::Equal),
            "p<=64 conversion must be exact"
        );

        // p > 64 with a nonzero tail: must round up (>= and within one ulp).
        // Build a 200-bit value with low bits set: 2^199 + 1.
        let big = BigFloat::try_from_i64_exact(1, 200).unwrap();
        let (two199, _) = big.scale_by_pow2(199);
        let (val, _) = two199.add(
            &BigFloat::try_from_i64_exact(1, 200).unwrap(),
            RoundingMode::NearestEven,
        );
        let m = Mag::from_bigfloat_ceil(&val);
        assert!(ge(&m.to_bigfloat(), &val), "ceil must over-cover");
    }

    #[test]
    fn to_bigfloat_round_trips_for_exact_values() {
        for &(e, m) in &[(0i64, 1u64 << 63), (5, (1u64 << 63) | 7), (-3, u64::MAX)] {
            let mag = mag_finite(e, m);
            let bf = mag.to_bigfloat();
            // Narrowing an exact 64-bit value back gives the same Mag.
            assert_eq!(Mag::from_bigfloat_ceil(&bf), mag);
        }
    }
}
