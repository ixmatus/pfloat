//! Differential for `BigFloat::i0`/`i1`/`in_` and `k0`/`k1`/`kn`
//! (modified Bessel functions, integer order). rug 1.30 / MPFR
//! expose **no** modified-Bessel `I`/`K` primitive (only `j*`/`y*`),
//! so — exactly as for `Bi`/`Ai′`/`Bi′` (slice 6n) — there is no
//! bit-exact MPFR oracle. The oracle is the Airy-style **tiered**
//! one ([`feedback_differential_lane_cost`]):
//!
//! 1. a checked-in authoritative reference table (`mpmath`
//!    `besseli`/`besselk` at 86 digits, the slice-6m/6n
//!    primary-source recipe; treated as a fact) parsed by the
//!    bit-exact decimal parser;
//! 2. the **DLMF 10.28.2 I/K cross-tie**
//!    `I_ν(x)·K_{ν+1}(x) + I_{ν+1}(x)·K_ν(x) = 1/x` over the table —
//!    the strong binding identity (a **plus** and `1/x`, the
//!    modified-Bessel analog of 6p's J/Y Wronskian, and π-free, so
//!    no transcendental constant string is needed);
//! 3. a high-precision self-consistency sweep.
//!
//! Tolerance `p − 2` bits (the `differential_si`/`differential_bi`
//! posture; pfloat has no Ziv correct-rounding yet). Precisions
//! capped at **256**: `K` composes the `bessel_i` kernel, `γ`/`ln`,
//! and (for `n ≥ 2`) an upward recurrence, so the full table at
//! `p = 1024` is impractical here — the `p = 1024` path is pinned
//! cheaply by the in-module `bessel_i`/`bessel_k` unit tests
//! (`high_precision_pin` and `ik_wronskian_10_28_2`). All sweep
//! arguments are strictly positive (`K`'s domain is `x > 0`; `I`'s
//! negative-argument parity is covered by `tests/property_ik.rs`).
//! Validate a fast subset locally; the full lane is the CI slow
//! tier.

#![cfg(all(unix, feature = "differential-mpfr"))]

mod differential;

use differential::{next_i64_in, sweep_size, TRANSCENDENTAL_PRECISIONS};
use pfloat::{BigFloat, RoundingMode};

fn close_within(a: &BigFloat, b: &BigFloat, bits: u32) -> bool {
    let (diff, _) = a.sub(b, RoundingMode::NearestEven);
    let abs = diff.abs();
    if abs.is_zero() {
        return true;
    }
    let p = a.precision().max(b.precision());
    let two = BigFloat::try_from_i64_exact(2, p).unwrap();
    let mut bound = b.abs();
    if bound.is_zero() {
        bound = BigFloat::try_from_i64_exact(1, p).unwrap();
    }
    for _ in 0..bits {
        bound = bound.div(&two, RoundingMode::NearestEven).0;
    }
    matches!(
        abs.partial_cmp(&bound).0,
        Some(core::cmp::Ordering::Less | core::cmp::Ordering::Equal)
    )
}

fn at(num: i64, den: i64, p: u32) -> BigFloat {
    let n = BigFloat::try_from_i64_exact(num, p).unwrap();
    if den == 1 {
        n
    } else {
        n.div(
            &BigFloat::try_from_i64_exact(den, p).unwrap(),
            RoundingMode::NearestEven,
        )
        .0
    }
}

/// `(num, den, I0, I1, I2, K0, K1, K2)`. Source: `mpmath`
/// `besseli(k, S(num)/den)` / `besselk(k, …)` for `k = 0,1,2` at
/// 86 digits (mpmath backend); treated as a mathematical fact. All
/// arguments strictly positive (the `K` domain), spanning the tiny
/// (`< 1`), Miller/series moderate, the recurrence-deep regime, and
/// the large-`|x|` asymptotic regime past
/// `bessel_j_threshold(256) = 9` (`|x| ≥ 512`). The `(1024, 1)`
/// entry is the sub-slice 2b.2.a boundary-input gate at `e_x = 10`,
/// firmly in the shared-Hankel-asymptotic regime that
/// `src/math/bessel_j.rs:455-464` gates; values cross `e+442` /
/// `e-447`, parsed via the `parse_str` scientific-notation path.
#[allow(clippy::type_complexity)]
const IK_TABLE: &[(i64, i64, &str, &str, &str, &str, &str, &str)] = &[
    (1,2,"1.0634833707413235192631844154453565293295231748211049891695720746879267185056918544346","0.25789430539089631636247965952320963418774314964079457273094519087056586338943968672536","0.031906149177738253813265777352517992578550576257926698245791311205663264947933107533115","0.92441907122766586178192416753021698953876831195352968481501974063291996009501604867818","1.6564411200033008936964454031740915115341007594640774460554278145261965895145190494663","7.5501835512408694365677057802265830356751713498098394690367309987377063181530922465434"),
    (1,4,"1.01568614122360792335473449128722628948859494723791000292242694549907460720596030242","0.12597910894546792590770717563144950440655320941484934241078967980523634693606980424941","0.007853269659864516093077086235630254236169271919115263636109507057183831717401868424694","1.5415067512483028161695161405755016073986235364162416059412414417118922119016847235276","3.7470259744407116380341332708455740527700797690543374780677528606756512207755890693591","31.5177145467739959204425823073400940295592616888509414304832643271171019781063972784"),
    (3,4,"1.1456467780440013276475708003815183884962797819343354555070703582733562067723482962809","0.4019924615809222052521049056996070505153984698592972704812894555246799471561429689757","0.073666880494875446975291051849232920455217195642876067556965143540876347689300379012368","0.6105824221164641193509450911332914482419101029149004148029243131962773203955609697592","0.94958046696214023217778111231532620122029232808382516471914648956829188399592197509067","3.142797000682171405158361390640827984829356311138434187387314952045055677718019570001"),
    (1,1,"1.2660658777520083355982446252147175376076703113549622068081353312135750161227754703948","0.56515910399248502720769602760986330732889962162109200948029448947925564096437113409266","0.13574766976703828118285256999499092294987106811277818784754635225506373419403320220949","0.42102443824070833333562737921260903613621974822666047229896955145521267813810183909213","0.60190723019723457473754000153561733926158688996810645601776795916855358294623784016886","1.6248388986351774828107073822838437146593935281628733843345054697923198440305775194299"),
    (5,4,"1.430468717721829681772097928201953563980792068528956610750170148450999517181357861264","0.75528141834074719156429791740864825034310074139339329629585564898259063931735934848301","0.22201844837663417526922126034811636343183088229952733667680111007885449427358290369123","0.29760308908410591656278136339460034058956456879290612746443898664985555452794952603367","0.40212407978419541259110293163983248397758555135712657056499675288601561719695371601447","0.94100161673881857670854605401833231495370145096430864036843379126748054204307547165682"),
    (5,2,"3.2898391440501230357059082299060560261118015753483941612552870534405381281069497332407","2.5167162452886984415281917481223776723889473033969492189945395536862123621998489805479","1.2764661478191642824833548314081538882006437326308347860596554104915682383470705488024","0.062347553200366186029169529476013925996005578743445303838599172078372612679992899458805","0.073890816347747063648993540591217582101975744210962005306723010118640815423264789079999","0.1214602062785638369483643619489879916775861741122149080839775801732852650186047307228"),
    (7,2,"7.3782034322254796603438378461630199141446755167300286852469178705210901192317856330658","6.2058349222583654736226092526488319834977439218096555304605353036851866558756860777232","3.8320120480778422468452039875065444950031075614102255249837548398438406015885364457954","0.019598897170368489107572958926796913651978892162537233013413853914027774831153866757959","0.022239392925923833738563671982262212515964595467046296883362451005581469265985484609583","0.032307121699467822672466485773803892232530089572277974089620968774360042983145572249149"),
    (7,1,"168.59390851028969885732662718750084037652267923453171419319405566855416412467567826523","156.03909286995545346239058066071115563003105204154494317062487461487033616940013352004","124.01131054744528358235788985586908162508523579409030185872980577859121093341849725951","0.00042479574186923180685159865280657229397931752463243151357157944598225095002003175131102","0.00045418248688489697123995940271024363126799741919564224934062310511312502497023217346402","0.00055456216669348808434872991072378476005588821583118644195461461887171524286866951515789"),
    (9,2,"17.481171855609276043133119086994202726544377972749492371368218396062680521664252597626","15.389222753735923892693846937799842641505917109182904499784434445773227093923134251074","10.641517298393309868602520447972050441430637035334868149241803086830135146587304041593","0.0063998572432339750456378847251674384712947706746218510098038247661037271875769704281262","0.0070780949089680896929249404627387344537687526864660329763418439751577634660183310916723","0.0095456772027753482424934138197179871174142163130511989992890887550627331724740064688695"),
    (13,1,"49444.489582217572999254597698562867074686667334837161976037001700405441569237798726476","47502.987358995860929359580262337492005309074796158425017948043759592553456050602346928","42136.337680833594394737739196664791381562194289274327357891148814314279499076167596179","0.00000077845438614204963208218355523693605491657583401293492878803535575270982340413516015053","0.00000080785884122023473318192958112402262010635144523861761679180896720044658054003530438107","0.00000090274036171439343718709579848678568877909144097272225444831365839893237425644828390146"),
    (1024,1,"6.5067271678114070755797890594060035801280561916534415277733928563923658932384234534700e+442","6.5035492785155352684873291632180149619812567949280665965851640035476673838390974443838e+442","6.4940249231268064207585247446340933946554365494758476477019374579479368553793627162739e+442","7.5042536073015053508129957935075339188727241263790318594900600327676517501121119963267e-447","7.5079168999287341801191223296006867495055759964316901178012281786107305822650684334375e-447","7.5189175074966786597585409543075352601803522044970625042513905565540008332805984581107e-447"),
];

#[test]
fn ik_family_matches_authoritative_table() {
    for &p in TRANSCENDENTAL_PRECISIONS.iter().filter(|&&p| p <= 256) {
        for &(num, den, i0s, i1s, i2s, k0s, k1s, k2s) in IK_TABLE {
            let x = at(num, den, p);
            let (i0, _) = x.i0(RoundingMode::NearestEven);
            let (i1, _) = x.i1(RoundingMode::NearestEven);
            let (i2, _) = x.in_(2, RoundingMode::NearestEven);
            let (k0, _) = x.k0(RoundingMode::NearestEven);
            let (k1, _) = x.k1(RoundingMode::NearestEven);
            let (k2, _) = x.kn(2, RoundingMode::NearestEven);
            for (got, refstr, name) in [
                (&i0, i0s, "I0"),
                (&i1, i1s, "I1"),
                (&i2, i2s, "I2"),
                (&k0, k0s, "K0"),
                (&k1, k1s, "K1"),
                (&k2, k2s, "K2"),
            ] {
                let want = BigFloat::parse_str(refstr, p, RoundingMode::NearestEven)
                    .unwrap()
                    .0;
                assert!(
                    close_within(got, &want, p - 2),
                    "{name}({num}/{den}) at p={p}: got {got}, want {want}"
                );
            }
        }
    }
}

#[test]
fn ik_cross_tie_holds() {
    // DLMF 10.28.2: I_ν(x)·K_{ν+1}(x) + I_{ν+1}(x)·K_ν(x) = 1/x for
    // every x > 0 and every order. The primary cross-tie binding the
    // I and K kernels, with no independent oracle. RHS is the exact
    // rational 1/x (π-free), so unlike the Airy Wronskian no
    // transcendental constant string is needed.
    for &p in TRANSCENDENTAL_PRECISIONS.iter().filter(|&&p| p <= 256) {
        for &(num, den, ..) in IK_TABLE {
            let x = at(num, den, p);
            let (inv_x, _) = BigFloat::try_from_i64_exact(1, p)
                .unwrap()
                .div(&x, RoundingMode::NearestEven);
            for nu in [0i32, 1, 2, 3] {
                let (i_nu, _) = x.in_(nu, RoundingMode::NearestEven);
                let (i_nu1, _) = x.in_(nu + 1, RoundingMode::NearestEven);
                let (k_nu, _) = x.kn(nu, RoundingMode::NearestEven);
                let (k_nu1, _) = x.kn(nu + 1, RoundingMode::NearestEven);
                let (a, _) = i_nu.mul(&k_nu1, RoundingMode::NearestEven);
                let (b, _) = i_nu1.mul(&k_nu, RoundingMode::NearestEven);
                let (w, _) = a.add(&b, RoundingMode::NearestEven);
                assert!(
                    close_within(&w, &inv_x, p - 8),
                    "10.28.2 nu={nu} at {num}/{den}, p={p}: got {w}"
                );
            }
        }
    }
}

#[test]
fn i_self_consistent() {
    // I at p must agree with I at p+64 (rounded to p) to p−8 bits.
    // Dyadic arguments only (denominator a power of two), so the
    // constructed argument is bit-identical at p and p+64 — the
    // pf-ok9 lesson (a non-dyadic argument reconstructs differently
    // between the two precisions and spuriously fails).
    let mut state: u64 = u64::from_le_bytes(*b"pf6qisc_");
    for &p in &[53u32, 113] {
        for _ in 0..sweep_size().min(20) {
            let a = next_i64_in(&mut state, 1, 32);
            let n = next_i64_in(&mut state, 0, 4) as i32;
            let lo = at(a, 4, p).in_(n, RoundingMode::NearestEven).0;
            let hi = {
                let r = at(a, 4, p + 64).in_(n, RoundingMode::NearestEven).0;
                r.round_to_precision(p, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            };
            assert!(
                close_within(&lo, &hi, p - 8),
                "I_{n}({a}/4) self-consistency at p={p}: lo={lo} hi={hi}"
            );
        }
    }
}

#[test]
fn k_self_consistent() {
    let mut state: u64 = u64::from_le_bytes(*b"pf6qksc_");
    for &p in &[53u32, 113] {
        for _ in 0..sweep_size().min(20) {
            let a = next_i64_in(&mut state, 1, 32);
            let n = next_i64_in(&mut state, 0, 4) as i32;
            let lo = at(a, 4, p).kn(n, RoundingMode::NearestEven).0;
            let hi = {
                let r = at(a, 4, p + 64).kn(n, RoundingMode::NearestEven).0;
                r.round_to_precision(p, RoundingMode::NearestEven)
                    .unwrap()
                    .0
            };
            assert!(
                close_within(&lo, &hi, p - 8),
                "K_{n}({a}/4) self-consistency at p={p}: lo={lo} hi={hi}"
            );
        }
    }
}
