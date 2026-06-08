//! Enumerated C99/C11 Annex G special-value tables for the complex elementary
//! core, exercised through the PUBLIC `Complex` API (ADR-0092, the C5
//! verification pass for ADR-0091).
//!
//! This is the primary branch-cut guard. Where the kernels' own `#[cfg(test)]`
//! modules check a representative subset against the internal `*_big`
//! functions, this lane enumerates every Annex G special-value row through the
//! public `Complex::sqrt`/`exp`/`log`/`mul`/`div` bridge, across precisions and
//! across all five rounding modes. The mode sweep is the regression guard for
//! the C4 signed-zero defect class (ADR-0091's `resolve` warning): a result
//! zero whose sign is fixed by an input sign was wrong in four of five modes
//! when routed through the rounding-mode-stamped `resolve` path; the
//! `copysign(0, input)` fix must hold under every mode, so every signed-zero
//! row is asserted under all of them.
//!
//! Rows the standard leaves a sign genuinely unspecified are marked `[REP]` in
//! comments: the assertion pins magnitude and NaN-ness and the crate's chosen
//! representative sign, with the note that the standard permits either.

#![cfg(feature = "trig")]

mod common;

use common::{
    bf, cbf, expect_cls, expect_int, ninf, nz, pinf, pz, qnan, snan, ResultCls::*, ALL_MODES,
};
use core::cmp::Ordering;
use pfloat::RoundingMode;
use pfloat_complex::Complex;

const PRECS: [u32; 2] = [53, 113];

// ============================ csqrt (§G.6.4.2) ============================

#[test]
fn csqrt_infinite_imaginary_dominates_every_real() {
    // The y = ±inf rows dominate even x = ±inf and x = NaN:
    // csqrt(any + inf*i) = +inf + inf*i; csqrt(any - inf*i) = +inf - inf*i.
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            for x in [bf(3, p), bf(-3, p), pinf(p), ninf(p), qnan(p), pz(p), nz(p)] {
                let (up, _) = Complex::new(x.clone(), pinf(p)).sqrt(mode);
                expect_cls(&up.re, PosInf, "csqrt(any + inf*i).re");
                expect_cls(&up.im, PosInf, "csqrt(any + inf*i).im");
                let (lo, _) = Complex::new(x, ninf(p)).sqrt(mode);
                expect_cls(&lo.re, PosInf, "csqrt(any - inf*i).re");
                expect_cls(&lo.im, NegInf, "csqrt(any - inf*i).im");
            }
        }
    }
}

#[test]
fn csqrt_positive_infinity_real() {
    // csqrt(+inf + y) for finite y: +inf + copysign(0, y); for NaN y: +inf + NaN.
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            let (a, _) = Complex::new(pinf(p), bf(3, p)).sqrt(mode);
            expect_cls(&a.re, PosInf, "csqrt(+inf + 3i).re");
            expect_cls(&a.im, PosZero, "csqrt(+inf + 3i).im");
            let (b, _) = Complex::new(pinf(p), bf(-3, p)).sqrt(mode);
            expect_cls(&b.im, NegZero, "csqrt(+inf - 3i).im");
            let (c, _) = Complex::new(pinf(p), qnan(p)).sqrt(mode);
            expect_cls(&c.re, PosInf, "csqrt(+inf + NaN*i).re");
            expect_cls(&c.im, Nan, "csqrt(+inf + NaN*i).im");
        }
    }
}

#[test]
fn csqrt_negative_infinity_real() {
    // csqrt(-inf + y) for finite y: +0 + copysign(inf, y); for NaN y: NaN + inf*i
    // ([REP] imaginary sign +inf).
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            let (a, _) = Complex::new(ninf(p), bf(3, p)).sqrt(mode);
            expect_cls(&a.re, PosZero, "csqrt(-inf + 3i).re");
            expect_cls(&a.im, PosInf, "csqrt(-inf + 3i).im");
            let (b, _) = Complex::new(ninf(p), bf(-3, p)).sqrt(mode);
            expect_cls(&b.re, PosZero, "csqrt(-inf - 3i).re");
            expect_cls(&b.im, NegInf, "csqrt(-inf - 3i).im");
            // [REP]: csqrt(-inf + NaN*i) = NaN + (+inf)i; standard permits
            // either imaginary sign, the crate picks +inf.
            let (c, _) = Complex::new(ninf(p), qnan(p)).sqrt(mode);
            expect_cls(&c.re, Nan, "csqrt(-inf + NaN*i).re");
            expect_cls(&c.im, PosInf, "csqrt(-inf + NaN*i).im [REP +inf]");
        }
    }
}

#[test]
fn csqrt_nan_without_infinity() {
    // A quiet NaN with no infinity propagates without INVALID; a signaling NaN
    // raises INVALID.
    for &p in &PRECS {
        let (a, sa) = Complex::new(qnan(p), bf(2, p)).sqrt(RoundingMode::NearestEven);
        expect_cls(&a.re, Nan, "csqrt(NaN + 2i).re");
        expect_cls(&a.im, Nan, "csqrt(NaN + 2i).im");
        assert!(!sa.invalid(), "quiet NaN does not signal");
        let (b, sb) = Complex::new(bf(2, p), qnan(p)).sqrt(RoundingMode::NearestEven);
        expect_cls(&b.re, Nan, "csqrt(2 + NaN*i).re");
        assert!(!sb.invalid());
        let (_, sc) = Complex::new(snan(p), bf(2, p)).sqrt(RoundingMode::NearestEven);
        assert!(sc.invalid(), "signaling NaN raises INVALID");
    }
}

#[test]
fn csqrt_real_axis_signed_zero_under_all_modes() {
    // The C4 regression guard. The imaginary-zero sign on the real axis follows
    // the input y, NOT the rounding mode, in ALL five modes (the defect: routing
    // it through `resolve` gave +0 in four of five modes).
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            // x > 0: csqrt(4 + 0i) = 2 + 0i, csqrt(4 - 0i) = 2 - 0i.
            let (up, _) = Complex::new(bf(4, p), pz(p)).sqrt(mode);
            expect_int(&up.re, 2, "csqrt(4 + 0i).re");
            expect_cls(&up.im, PosZero, "csqrt(4 + 0i).im");
            let (dn, _) = Complex::new(bf(4, p), nz(p)).sqrt(mode);
            expect_int(&dn.re, 2, "csqrt(4 - 0i).re");
            expect_cls(
                &dn.im,
                NegZero,
                "csqrt(4 - 0i).im (must be -0 under every mode)",
            );

            // x < 0 (the branch cut): csqrt(-4 + 0i) = +0 + 2i,
            // csqrt(-4 - 0i) = +0 - 2i.
            let (cu, _) = Complex::new(bf(-4, p), pz(p)).sqrt(mode);
            expect_cls(&cu.re, PosZero, "csqrt(-4 + 0i).re");
            expect_int(&cu.im, 2, "csqrt(-4 + 0i).im");
            let (cd, _) = Complex::new(bf(-4, p), nz(p)).sqrt(mode);
            expect_cls(&cd.re, PosZero, "csqrt(-4 - 0i).re");
            expect_int(&cd.im, -2, "csqrt(-4 - 0i).im");

            // Origin: csqrt(±0 ± 0i) = +0 + copysign(0, y).
            let (o1, _) = Complex::new(pz(p), pz(p)).sqrt(mode);
            expect_cls(&o1.re, PosZero, "csqrt(+0 + 0i).re");
            expect_cls(&o1.im, PosZero, "csqrt(+0 + 0i).im");
            let (o2, _) = Complex::new(nz(p), nz(p)).sqrt(mode);
            expect_cls(&o2.re, PosZero, "csqrt(-0 - 0i).re");
            expect_cls(&o2.im, NegZero, "csqrt(-0 - 0i).im");
        }
    }
}

#[test]
fn csqrt_gaussian_integer_roots_are_exact() {
    // Exact algebraic outputs report OK (no forced INEXACT), one per Kahan sign
    // branch: (3+4i)->(2+i); (5-12i)->(3-2i); (-5-12i)->(2-3i); (-5+12i)->(2+3i);
    // (-7+24i)->(3+4i).
    for &p in &PRECS {
        for &(zr, zi, wr, wi) in &[
            (3i64, 4i64, 2i64, 1i64),
            (5, -12, 3, -2),
            (-5, -12, 2, -3),
            (-5, 12, 2, 3),
            (-7, 24, 3, 4),
        ] {
            let (r, s) = cbf(zr, zi, p).sqrt(RoundingMode::NearestEven);
            expect_int(&r.re, wr, "csqrt gaussian .re");
            expect_int(&r.im, wi, "csqrt gaussian .im");
            assert!(!s.inexact(), "csqrt({zr}+{zi}i) is an exact Gaussian root");
        }
    }
}

// ============================ cexp (§G.6.3.1) ============================

#[test]
fn cexp_real_axis_signed_zero_under_all_modes() {
    // cexp(x ± 0i) = (e^x, copysign(0, y)); the imaginary-zero sign follows y in
    // every mode (stamped, not formed by e^x * 0).
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            let (a, sa) = Complex::new(pz(p), pz(p)).exp(mode);
            expect_int(&a.re, 1, "cexp(+0 + 0i).re = 1");
            expect_cls(&a.im, PosZero, "cexp(+0 + 0i).im");
            assert!(!sa.inexact(), "cexp(0) is exact");
            let (b, _) = Complex::new(bf(2, p), nz(p)).exp(mode);
            expect_cls(
                &b.im,
                NegZero,
                "cexp(2 - 0i).im (must be -0 under every mode)",
            );
        }
    }
}

#[test]
fn cexp_pos_infinity_real() {
    for &p in &PRECS {
        // y = ±0: cexp(+inf + 0i) = +inf + copysign(0, y).
        let (a, _) = Complex::new(pinf(p), pz(p)).exp(RoundingMode::NearestEven);
        expect_cls(&a.re, PosInf, "cexp(+inf + 0i).re");
        expect_cls(&a.im, PosZero, "cexp(+inf + 0i).im");
        let (b, _) = Complex::new(pinf(p), nz(p)).exp(RoundingMode::NearestEven);
        expect_cls(&b.im, NegZero, "cexp(+inf - 0i).im");
        // y = ±inf: cexp(+inf + inf*i) = +inf + NaN*i + INVALID [REP +inf].
        let (c, sc) = Complex::new(pinf(p), pinf(p)).exp(RoundingMode::NearestEven);
        expect_cls(&c.re, PosInf, "cexp(+inf + inf*i).re [REP +inf]");
        expect_cls(&c.im, Nan, "cexp(+inf + inf*i).im");
        assert!(sc.invalid(), "cexp(+inf + inf*i) raises INVALID");
        // y = NaN: cexp(+inf + NaN*i) = +inf + NaN*i, no INVALID (quiet).
        let (d, sd) = Complex::new(pinf(p), qnan(p)).exp(RoundingMode::NearestEven);
        expect_cls(&d.re, PosInf, "cexp(+inf + NaN*i).re [REP +inf]");
        expect_cls(&d.im, Nan, "cexp(+inf + NaN*i).im");
        assert!(!sd.invalid());
        // y finite nonzero: signs from sign(cos y), sign(sin y).
        // cexp(+inf + 1i) = (+inf, +inf) (cos 1 > 0, sin 1 > 0).
        let (e, _) = Complex::new(pinf(p), bf(1, p)).exp(RoundingMode::NearestEven);
        expect_cls(&e.re, PosInf, "cexp(+inf + 1i).re");
        expect_cls(&e.im, PosInf, "cexp(+inf + 1i).im");
        // cexp(+inf + 3i) = (-inf, +inf) (cos 3 < 0, sin 3 > 0).
        let (f, _) = Complex::new(pinf(p), bf(3, p)).exp(RoundingMode::NearestEven);
        expect_cls(&f.re, NegInf, "cexp(+inf + 3i).re");
        expect_cls(&f.im, PosInf, "cexp(+inf + 3i).im");
    }
}

#[test]
fn cexp_neg_infinity_real_dominates_without_invalid() {
    // The load-bearing asymmetry: e^{-inf} = +0 dominates the indeterminate
    // angle, so the x = -inf rows are (+0, +0) with NO INVALID.
    for &p in &PRECS {
        let (a, _) = Complex::new(ninf(p), pz(p)).exp(RoundingMode::NearestEven);
        expect_cls(&a.re, PosZero, "cexp(-inf + 0i).re");
        expect_cls(&a.im, PosZero, "cexp(-inf + 0i).im");
        let (b, _) = Complex::new(ninf(p), nz(p)).exp(RoundingMode::NearestEven);
        expect_cls(&b.im, NegZero, "cexp(-inf - 0i).im");
        // [REP +0]: x = -inf, y = ±inf or NaN -> (+0, +0), NO INVALID.
        for y in [pinf(p), ninf(p), qnan(p)] {
            let (c, sc) = Complex::new(ninf(p), y.clone()).exp(RoundingMode::NearestEven);
            expect_cls(&c.re, PosZero, "cexp(-inf + {inf|NaN}*i).re [REP +0]");
            expect_cls(&c.im, PosZero, "cexp(-inf + {inf|NaN}*i).im [REP +0]");
            assert!(!sc.invalid(), "-inf real part suppresses INVALID");
        }
        // y finite nonzero: (sign(cos y)*0, sign(sin y)*0).
        // cexp(-inf + 3i): cos 3 < 0, sin 3 > 0 -> (-0, +0).
        let (d, _) = Complex::new(ninf(p), bf(3, p)).exp(RoundingMode::NearestEven);
        expect_cls(&d.re, NegZero, "cexp(-inf + 3i).re");
        expect_cls(&d.im, PosZero, "cexp(-inf + 3i).im");
    }
}

#[test]
fn cexp_finite_real_with_special_imaginary() {
    for &p in &PRECS {
        // y = ±inf: cexp(1 + inf*i) = NaN + NaN*i + INVALID.
        for y in [pinf(p), ninf(p)] {
            let (a, sa) = Complex::new(bf(1, p), y).exp(RoundingMode::NearestEven);
            expect_cls(&a.re, Nan, "cexp(1 + inf*i).re");
            expect_cls(&a.im, Nan, "cexp(1 + inf*i).im");
            assert!(sa.invalid(), "cos/sin of inf is INVALID");
        }
        // y = NaN: cexp(2 + NaN*i) = NaN + NaN*i, no INVALID (quiet).
        let (b, sb) = Complex::new(bf(2, p), qnan(p)).exp(RoundingMode::NearestEven);
        expect_cls(&b.re, Nan, "cexp(2 + NaN*i).re");
        assert!(!sb.invalid());
    }
}

#[test]
fn cexp_nan_real() {
    // cexp(NaN + 0i) = NaN + copysign(0, y); cexp(NaN - 0i) = NaN - 0i (the
    // conjugation symmetry pins the imaginary-zero sign); cexp(NaN + 2i) = NaN+NaN.
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            let (a, _) = Complex::new(qnan(p), pz(p)).exp(mode);
            expect_cls(&a.re, Nan, "cexp(NaN + 0i).re");
            expect_cls(&a.im, PosZero, "cexp(NaN + 0i).im");
            let (b, _) = Complex::new(qnan(p), nz(p)).exp(mode);
            expect_cls(&b.im, NegZero, "cexp(NaN - 0i).im (conjugation-pinned)");
        }
        let (c, sc) = Complex::new(qnan(p), bf(2, p)).exp(RoundingMode::NearestEven);
        expect_cls(&c.re, Nan, "cexp(NaN + 2i).re");
        expect_cls(&c.im, Nan, "cexp(NaN + 2i).im");
        assert!(!sc.invalid());
        let (_, sd) = Complex::new(snan(p), pz(p)).exp(RoundingMode::NearestEven);
        assert!(sd.invalid(), "signaling NaN raises INVALID");
    }
}

// ============================ clog (§G.6.3.2) ============================

#[test]
fn clog_four_poles_under_all_modes() {
    // clog(±0 ± 0i) = -inf + i*atan2(±0, ±0), all four carrying DIV_BY_ZERO.
    // atan2(+0,+0)=+0; atan2(+0,-0)=+pi; atan2(-0,+0)=-0; atan2(-0,-0)=-pi.
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            let cases = [
                (pz(p), pz(p), PosZero, "clog(+0 + 0i).im = +0"),
                (nz(p), pz(p), PosFin, "clog(-0 + 0i).im = +pi"),
                (pz(p), nz(p), NegZero, "clog(+0 - 0i).im = -0"),
                (nz(p), nz(p), NegFin, "clog(-0 - 0i).im = -pi"),
            ];
            for (x, y, im_cls, who) in cases {
                let (r, s) = Complex::new(x, y).log(mode);
                expect_cls(&r.re, NegInf, "clog pole .re = -inf");
                expect_cls(&r.im, im_cls, who);
                assert!(s.div_by_zero(), "clog pole raises DIV_BY_ZERO");
            }
        }
    }
}

#[test]
fn clog_one_is_exact_zero() {
    // clog(1 + 0i) = +0 + 0i, exact (hypot=1, ln(1)=+0; atan2(0,1)=+0).
    for &p in &PRECS {
        let (r, s) = Complex::new(bf(1, p), pz(p)).log(RoundingMode::NearestEven);
        expect_cls(&r.re, PosZero, "clog(1 + 0i).re = +0");
        expect_cls(&r.im, PosZero, "clog(1 + 0i).im = +0");
        assert!(!s.inexact(), "clog(1) is exact");
    }
}

#[test]
fn clog_infinite_part_makes_real_part_pos_inf() {
    for &p in &PRECS {
        // The subtle row: clog(NaN + inf*i) = +inf + NaN*i (hypot(NaN,inf)=+inf
        // dominates, so the REAL part is +inf, not NaN).
        let (a, _) = Complex::new(qnan(p), pinf(p)).log(RoundingMode::NearestEven);
        expect_cls(&a.re, PosInf, "clog(NaN + inf*i).re = +inf (not NaN)");
        expect_cls(&a.im, Nan, "clog(NaN + inf*i).im = NaN");
        // clog(-inf + 3i): re = +inf, im = atan2(3, -inf) = +pi.
        let (b, _) = Complex::new(ninf(p), bf(3, p)).log(RoundingMode::NearestEven);
        expect_cls(&b.re, PosInf, "clog(-inf + 3i).re");
        expect_cls(&b.im, PosFin, "clog(-inf + 3i).im = +pi");
        // clog(+inf + 3i): re = +inf, im = atan2(3, +inf) = +0.
        let (c, _) = Complex::new(pinf(p), bf(3, p)).log(RoundingMode::NearestEven);
        expect_cls(&c.re, PosInf, "clog(+inf + 3i).re");
        expect_cls(&c.im, PosZero, "clog(+inf + 3i).im = +0");
        // clog(2 + inf*i): re = +inf, im = atan2(inf, 2) = +pi/2.
        let (d, _) = Complex::new(bf(2, p), pinf(p)).log(RoundingMode::NearestEven);
        expect_cls(&d.re, PosInf, "clog(2 + inf*i).re");
        expect_cls(&d.im, PosFin, "clog(2 + inf*i).im = +pi/2");
    }
}

#[test]
fn clog_nan_without_infinity() {
    for &p in &PRECS {
        let (a, sa) = Complex::new(qnan(p), bf(2, p)).log(RoundingMode::NearestEven);
        expect_cls(&a.re, Nan, "clog(NaN + 2i).re");
        expect_cls(&a.im, Nan, "clog(NaN + 2i).im");
        assert!(!sa.invalid(), "quiet NaN does not signal");
        let (_, sb) = Complex::new(snan(p), bf(2, p)).log(RoundingMode::NearestEven);
        assert!(sb.invalid(), "signaling NaN raises INVALID");
    }
}

#[test]
fn clog_negative_real_axis_branch_cut_under_all_modes() {
    // clog(-2 + 0i) = ln 2 + i*pi, clog(-2 - 0i) = ln 2 - i*pi. Real parts equal
    // and positive (ln 2 > 0); imaginary parts opposite signs, under every mode.
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            let (u, _) = Complex::new(bf(-2, p), pz(p)).log(mode);
            let (l, _) = Complex::new(bf(-2, p), nz(p)).log(mode);
            expect_cls(&u.re, PosFin, "clog(-2 + 0i).re = ln 2 > 0");
            expect_cls(&l.re, PosFin, "clog(-2 - 0i).re = ln 2 > 0");
            assert_eq!(
                u.re.partial_cmp(&l.re).0,
                Some(Ordering::Equal),
                "both real parts equal ln 2"
            );
            expect_cls(&u.im, PosFin, "clog(-2 + 0i).im = +pi");
            expect_cls(&l.im, NegFin, "clog(-2 - 0i).im = -pi");
        }
    }
}

// ===================== §G.5.1 complex-infinity div =====================

#[test]
fn div_finite_over_complex_zero_is_directed_infinity() {
    // D1: a finite nonzero dividend over a complex-zero divisor is a directed
    // complex infinity; the direction is from the divisor real part c ONLY.
    for &p in &PRECS {
        // c = +0 -> +inf direction.
        let (a, _) = cbf(1, 1, p).div(&Complex::new(pz(p), pz(p)), RoundingMode::NearestEven);
        expect_cls(&a.re, PosInf, "(1+1i)/(+0+0i).re");
        expect_cls(&a.im, PosInf, "(1+1i)/(+0+0i).im");
        // c = -0 -> -inf direction (d is never consulted).
        let (b, _) = cbf(1, 1, p).div(&Complex::new(nz(p), pz(p)), RoundingMode::NearestEven);
        expect_cls(&b.re, NegInf, "(1+1i)/(-0+0i).re (direction from c=-0)");
        expect_cls(&b.im, NegInf, "(1+1i)/(-0+0i).im");
    }
}

#[test]
fn div_zero_over_zero_is_nan_invalid() {
    for &p in &PRECS {
        let (r, s) =
            Complex::new(pz(p), pz(p)).div(&Complex::new(pz(p), pz(p)), RoundingMode::NearestEven);
        expect_cls(&r.re, Nan, "0/0 .re");
        expect_cls(&r.im, Nan, "0/0 .im");
        assert!(s.invalid(), "0/0 raises INVALID");
    }
}

#[test]
fn div_partial_zero_dividend_over_zero() {
    // D1: (0 + 1i)/(0 + 0i): re = inf*0 = NaN, im = inf*1 = +inf; still a
    // complex infinity (§G.3) via the imaginary +inf; INVALID from the inf*0.
    for &p in &PRECS {
        let (r, s) = Complex::new(pz(p), bf(1, p))
            .div(&Complex::new(pz(p), pz(p)), RoundingMode::NearestEven);
        expect_cls(&r.re, Nan, "(0+1i)/(0+0i).re = inf*0 = NaN");
        expect_cls(&r.im, PosInf, "(0+1i)/(0+0i).im = +inf");
        assert!(s.invalid(), "inf*0 raises INVALID");
    }
}

#[test]
fn div_infinite_over_finite_and_finite_over_infinite() {
    for &p in &PRECS {
        // D2: (inf + 0i)/(2 + 0i) = +inf + NaN*i.
        let (a, _) = Complex::new(pinf(p), pz(p)).div(&cbf(2, 0, p), RoundingMode::NearestEven);
        expect_cls(&a.re, PosInf, "(inf+0i)/(2+0i).re");
        expect_cls(&a.im, Nan, "(inf+0i)/(2+0i).im");
        // D3: (2 + 0i)/(inf + 0i) = +0 + 0i (signed zeros).
        let (b, _) = cbf(2, 0, p).div(&Complex::new(pinf(p), pz(p)), RoundingMode::NearestEven);
        expect_cls(&b.re, PosZero, "(2+0i)/(inf+0i).re");
        expect_cls(&b.im, PosZero, "(2+0i)/(inf+0i).im");
        // D3 fires on an (inf + NaN*i) divisor too (the NaN part boxes to zero).
        let (c, _) = cbf(2, 0, p).div(&Complex::new(pinf(p), qnan(p)), RoundingMode::NearestEven);
        expect_cls(&c.re, PosZero, "(2+0i)/(inf+NaN*i).re");
        expect_cls(&c.im, PosZero, "(2+0i)/(inf+NaN*i).im");
        // inf/inf is NOT recovered: (inf+inf*i)/(inf+inf*i) = NaN + NaN*i.
        let (d, _) = Complex::new(pinf(p), pinf(p))
            .div(&Complex::new(pinf(p), pinf(p)), RoundingMode::NearestEven);
        expect_cls(&d.re, Nan, "inf/inf .re");
        expect_cls(&d.im, Nan, "inf/inf .im");
    }
}

#[test]
fn div_finite_nonzero_divisor_falls_through_to_ziv() {
    // The case the C3 directed-pair Ziv divide owns (no §G.5.1 dispatch):
    // (1 + 0i)/(3 + 0i) real part = 1/3 correctly rounded, INEXACT, equal to the
    // scalar 1/3 bit-for-bit under every mode.
    for &p in &PRECS {
        for &mode in &ALL_MODES {
            let (r, s) = cbf(1, 0, p).div(&cbf(3, 0, p), mode);
            let scalar = bf(1, p).div(&bf(3, p), mode).0;
            assert_eq!(
                r.re.partial_cmp(&scalar).0,
                Some(Ordering::Equal),
                "(1+0i)/(3+0i).re diverged from scalar 1/3 under {mode:?}"
            );
            // The imaginary part is an exact CANCELLATION zero (bc - ad = 0), so
            // its sign is mode-determined per the IEEE 754 same-sign-difference
            // rule (-0 under TowardNegative, +0 otherwise) -- NOT input-fixed.
            // Only the zero-ness is asserted; the sign is the C3 `resolve`
            // mode-sign behavior (ADR-0090), correct, not a branch-cut zero.
            assert!(r.im.is_zero(), "(1+0i)/(3+0i).im is a (cancellation) zero");
            assert!(s.inexact(), "1/3 is INEXACT");
        }
    }
}

// ===================== §G.5.1 complex-infinity mul =====================

#[test]
fn mul_recovers_complex_infinity() {
    for &p in &PRECS {
        // M1: (1 + 0i)*(inf + inf*i) = inf + inf*i (naive cross products give
        // (NaN, NaN); the recovery restores the infinity).
        let (a, _) = cbf(1, 0, p).mul(&Complex::new(pinf(p), pinf(p)), RoundingMode::NearestEven);
        expect_cls(&a.re, PosInf, "(1+0i)*(inf+inf*i).re");
        expect_cls(&a.im, PosInf, "(1+0i)*(inf+inf*i).im");
        // Sign follows the boxed parts: (1 + 0i)*(-inf + 0i): boxed (-1, +0);
        // n_re = -1 -> -inf; n_im = 0 -> inf*0 = NaN.
        let (b, _) = cbf(1, 0, p).mul(&Complex::new(ninf(p), pz(p)), RoundingMode::NearestEven);
        expect_cls(&b.re, NegInf, "(1+0i)*(-inf+0i).re");
        expect_cls(&b.im, Nan, "(1+0i)*(-inf+0i).im");
    }
}

#[test]
fn mul_infinity_times_zero_stays_nan() {
    for &p in &PRECS {
        // (inf + inf*i)*(0 + 0i) is a genuine inf*0: stays NaN + NaN*i + INVALID.
        let (r, s) = Complex::new(pinf(p), pinf(p))
            .mul(&Complex::new(pz(p), pz(p)), RoundingMode::NearestEven);
        expect_cls(&r.re, Nan, "(inf+inf*i)*(0+0i).re");
        expect_cls(&r.im, Nan, "(inf+inf*i)*(0+0i).im");
        assert!(s.invalid(), "inf*0 raises INVALID");
    }
}

#[test]
fn mul_finite_is_exact_with_cancellation() {
    // The finite path is unchanged (Annex G correct for finite inputs):
    // (3 + 4i)*(3 - 4i) = 25 + 0i, the imaginary part exactly +0 (no spurious
    // INEXACT on the cancelling component).
    for &p in &PRECS {
        let z = cbf(3, 4, p);
        let (r, s) = z.mul(&z.conj(), RoundingMode::NearestEven);
        expect_int(&r.re, 25, "(3+4i)*(3-4i).re = 25");
        expect_cls(&r.im, PosZero, "(3+4i)*(3-4i).im = +0");
        assert!(!s.inexact(), "exact product, no spurious INEXACT");
    }
}

// ============== Gap-closers from the C5 independent re-derivation ==============
// The annex-g-rederive-refute workflow (primary-source re-derivation +
// refutation, all rows SOUND) flagged a handful of rows the table had not
// pinned explicitly. Each was re-derived here before adding (verify-the-verdict;
// one of the workflow's six suggestions -- a `(2+3i)*(inf+inf*i)` "boxing
// asymmetry" test -- was itself wrong and is corrected below).

#[test]
fn clog_negative_infinity_with_nan_imaginary() {
    // The companion to clog(NaN + inf*i): BOTH signs of a real infinity make the
    // real part +inf (hypot(-inf, NaN) = +inf dominates), so
    // clog(-inf + NaN*i) = +inf + NaN*i. im = atan2(NaN, -inf) = NaN.
    for &p in &PRECS {
        let (r, _) = Complex::new(ninf(p), qnan(p)).log(RoundingMode::NearestEven);
        expect_cls(&r.re, PosInf, "clog(-inf + NaN*i).re = +inf");
        expect_cls(&r.im, Nan, "clog(-inf + NaN*i).im = NaN");
    }
}

#[test]
fn div_zero_divisor_direction_ignores_imaginary_part() {
    // §G.5.1 D1: the directed infinity's sign comes from the divisor REAL part
    // c ONLY; the divisor imaginary part d must NOT flip it. (1 + 1i)/(+0 - 0i)
    // is +inf + inf*i (c = +0 wins; d = -0 is never consulted), identical to
    // (1 + 1i)/(+0 + 0i).
    for &p in &PRECS {
        let (a, _) = cbf(1, 1, p).div(&Complex::new(pz(p), nz(p)), RoundingMode::NearestEven);
        expect_cls(&a.re, PosInf, "(1+1i)/(+0-0i).re (d=-0 must not flip)");
        expect_cls(&a.im, PosInf, "(1+1i)/(+0-0i).im");
        // And c = -0 with d = +0 still goes -inf (c wins, not d).
        let (b, _) = cbf(1, 1, p).div(&Complex::new(nz(p), pz(p)), RoundingMode::NearestEven);
        expect_cls(&b.re, NegInf, "(1+1i)/(-0+0i).re");
    }
}

#[test]
fn mul_recovery_keeps_finite_part_then_partial_nan_is_complex_infinity() {
    // Two related rows the table had not pinned:
    //
    // (a) (2 + 0i)*(inf + inf*i): naive = (NaN, NaN) (the 0 part makes both cross
    //     products carry 0*inf), so recovery fires. Boxing keeps the finite 2
    //     (NOT over-zeroed): n_re = 2*1 - 0*1 = 2 -> +inf, n_im = 2*1 + 0*1 = 2
    //     -> +inf. Result (+inf, +inf). Were the 2 zeroed, n_re would be NaN.
    //
    // (b) (2 + 3i)*(inf + inf*i): naive re = 2*inf - 3*inf = NaN but naive
    //     im = 2*inf + 3*inf = +inf, so the (NaN, NaN) recovery guard is NOT
    //     met -- pfloat (like compiler-rt __muldc3, which recovers only when
    //     BOTH parts are NaN) returns the PARTIAL complex infinity (NaN, +inf),
    //     itself a valid complex infinity by §G.3. (The workflow suggested this
    //     as a "boxing asymmetry" test; re-derivation shows it never reaches the
    //     boxing path, so it pins the recovery-firing CONDITION instead.)
    for &p in &PRECS {
        let (a, _) = cbf(2, 0, p).mul(&Complex::new(pinf(p), pinf(p)), RoundingMode::NearestEven);
        expect_cls(&a.re, PosInf, "(2+0i)*(inf+inf*i).re (finite 2 kept)");
        expect_cls(&a.im, PosInf, "(2+0i)*(inf+inf*i).im");

        let (b, _) = cbf(2, 3, p).mul(&Complex::new(pinf(p), pinf(p)), RoundingMode::NearestEven);
        expect_cls(
            &b.re,
            Nan,
            "(2+3i)*(inf+inf*i).re = NaN (recovery not fired)",
        );
        expect_cls(
            &b.im,
            PosInf,
            "(2+3i)*(inf+inf*i).im = +inf (partial complex inf)",
        );
    }
}
