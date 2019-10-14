#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;

use criterion::{Benchmark, Criterion};

use zerocaf::backend::u64::{field, scalar};
use zerocaf::edwards::EdwardsPoint;
use zerocaf::ristretto::RistrettoPoint;
use zerocaf::traits::ops::*;
use zerocaf::constants::*;
#[allow(unused_imports)]
use zerocaf::traits::Identity;

mod field_benches {

    use super::*;
    use subtle::Choice;
    use zerocaf::field::FieldElement;
    use zerocaf::scalar::Scalar;

    /// `B = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static B: field::FieldElement = field::FieldElement([
        2766226127823335,
        4237835465749098,
        4503599626623787,
        4503599627370493,
        2199023255551,
    ]);

    /// `A = 182687704666362864775460604089535377456991567872`
    pub static A: FieldElement = FieldElement([0, 0, 0, 2, 0]);

    pub fn bench_field_element_ops(c: &mut Criterion) {
        c.bench(
            "Field Element",
            Benchmark::new("Addition", |b| b.iter(|| &B + &A)),
        );

        c.bench(
            "Field Element",
            Benchmark::new("Subtraction", |b| b.iter(|| &B - &A)),
        );

        c.bench(
            "Field Element",
            Benchmark::new("Mul", |b| b.iter(|| &B * &A)),
        );

        c.bench(
            "Field Element",
            Benchmark::new("Squaring", |b| b.iter(|| B.square())),
        );

        c.bench(
            "Field Element",
            Benchmark::new("Half", |b| b.iter(|| A.half())),
        );

        c.bench(
            "Field Element",
            Benchmark::new("Two Pow k (2^k)", |b| {
                b.iter(|| FieldElement::two_pow_k(213u64))
            }),
        );

        c.bench(
            "Field Element",
            Benchmark::new("Pow (a^b (mod l))", |b| b.iter(|| A.pow(&B))),
        );

        c.bench(
            "Field Element",
            Benchmark::new("Modular Sqrt", |b| b.iter(|| A.mod_sqrt(Choice::from(1u8)))),
        );
    }

    pub fn bench_modular_inverse(c: &mut Criterion) {
        c.bench(
            "Modular Inverse",
            Benchmark::new("Savas & Koç Modular Inverse algorithm", |b| {
                b.iter(|| field::FieldElement::inverse(&B))
            }),
        );
    }
}

mod scalar_benches {
    use super::*;
    use zerocaf::field::FieldElement;
    use zerocaf::scalar::Scalar;

    /// `C = 182687704666362864775460604089535377456991567872`.
    pub static C: Scalar = Scalar([0, 0, 0, 2, 0]);

    /// `D = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static D: Scalar = Scalar([
        2766226127823335,
        4237835465749098,
        4503599626623787,
        4503599627370493,
        2199023255551,
    ]);

    pub fn bench_scalar_element_ops(c: &mut Criterion) {
        c.bench("Scalar", Benchmark::new("Addition", |b| b.iter(|| &C + &D)));

        c.bench(
            "Scalar",
            Benchmark::new("Subtraction", |b| b.iter(|| &C - &D)),
        );

        c.bench("Scalar", Benchmark::new("Mul", |b| b.iter(|| &C * &D)));

        c.bench(
            "Scalar",
            Benchmark::new("Squaring", |b| b.iter(|| C.square())),
        );

        c.bench("Scalar", Benchmark::new("Half", |b| b.iter(|| C.half())));
    }
}

mod edwards_benches {

    use super::*;
    use subtle::Choice;
    use zerocaf::edwards::{CompressedEdwardsY, EdwardsPoint, ProjectivePoint};
    use zerocaf::field::FieldElement;
    use zerocaf::scalar::Scalar;

    pub static P1_EXTENDED: EdwardsPoint = EdwardsPoint {
        X: FieldElement([13, 0, 0, 0, 0]),
        Y: FieldElement([
            606320128494542,
            1597163540666577,
            1835599237877421,
            1667478411389512,
            3232679738299,
        ]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([
            2034732376387996,
            3922598123714460,
            1344791952818393,
            3662820838581677,
            6840464509059,
        ]),
    };

    pub static P2_EXTENDED: EdwardsPoint = EdwardsPoint {
        X: FieldElement([67, 0, 0, 0, 0]),
        Y: FieldElement([
            2369245568431362,
            2665603790611352,
            3317390952748653,
            1908583331312524,
            8011773354506,
        ]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([
            3474019263728064,
            2548729061993416,
            1588812051971430,
            1774293631565269,
            9023233419450,
        ]),
    };

    pub static P1_PROJECTIVE: ProjectivePoint = ProjectivePoint {
        X: FieldElement([13, 0, 0, 0, 0]),
        Y: FieldElement([
            606320128494542,
            1597163540666577,
            1835599237877421,
            1667478411389512,
            3232679738299,
        ]),
        Z: FieldElement([1, 0, 0, 0, 0]),
    };

    pub static P2_PROJECTIVE: ProjectivePoint = ProjectivePoint {
        X: FieldElement([67, 0, 0, 0, 0]),
        Y: FieldElement([
            2369245568431362,
            2665603790611352,
            3317390952748653,
            1908583331312524,
            8011773354506,
        ]),
        Z: FieldElement([1, 0, 0, 0, 0]),
    };

    /// `A = 182687704666362864775460604089535377456991567872`.
    pub static A: Scalar = Scalar([0, 0, 0, 2, 0]);

    /// `P1_EXTENDED on `CompressedEdwardsY` format.
    pub(self) static P1_COMPRESSED: CompressedEdwardsY = CompressedEdwardsY([
        206, 11, 225, 231, 113, 39, 18, 141, 213, 215, 201, 201, 90, 173, 14, 134, 192, 119, 133,
        134, 164, 26, 38, 1, 201, 94, 187, 59, 186, 170, 240, 2,
    ]);

    pub fn bench_extended_point_ops(c: &mut Criterion) {
        c.bench(
            "Extended Coordinates Point Addition",
            Benchmark::new("2008 Hisil–Wong–Carter–Dawson, Section 3.1.", |b| {
                b.iter(|| &P1_EXTENDED + &P2_EXTENDED)
            }),
        );

        c.bench(
            "Extended Coordinates Point Subtraction",
            Benchmark::new("2008 Hisil–Wong–Carter–Dawson, Section 3.1.", |b| {
                b.iter(|| &P1_EXTENDED - &P2_EXTENDED)
            }),
        );

        c.bench(
            "Extended Coordinates Point Doubling",
            Benchmark::new("2008 Hisil–Wong–Carter–Dawson, Section 3.1.", |b| {
                b.iter(|| P1_EXTENDED.double())
            }),
        );

        c.bench(
            "Extended Coordinates Scalar Mul",
            Benchmark::new("Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004) - Guide to Elliptic Curve Cryptography. ",
            |b| b.iter(|| &P1_EXTENDED * &A))
        );
    }

    pub fn bench_projective_point_ops(c: &mut Criterion) {
        c.bench(
            "Projective Coordinates Point Addition",
            Benchmark::new("D. J. Bernstein, P. Birkner, M. Joye, T. Lange, C. Peters. AFRICACRYPT 2008 - Section 6.", |b| b.iter(|| &P1_PROJECTIVE + &P2_PROJECTIVE))
        );

        c.bench(
            "Projective Coordinates Point Subtraction",
            Benchmark::new("D. J. Bernstein, P. Birkner, M. Joye, T. Lange, C. Peters. AFRICACRYPT 2008 - Section 6.", |b| b.iter(|| &P1_PROJECTIVE - &P2_PROJECTIVE))
        );

        c.bench(
            "Projective Coordinates Point Doubling",
            Benchmark::new("D. J. Bernstein, P. Birkner, M. Joye, T. Lange, C. Peters. AFRICACRYPT 2008 - Section 6.", |b| b.iter(|| P1_PROJECTIVE.double()))
        );

        c.bench(
            "Projective Coordinates Scalar Mul",
            Benchmark::new("Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004) - Guide to Elliptic Curve Cryptography. ",
            |b| b.iter(|| &P1_PROJECTIVE * &A))
        );

        c.bench(
            "Projective Coordinates Point Generation.",
            Benchmark::new("From y coordinate.", |b| {
                b.iter(|| {
                    ProjectivePoint::new_from_y_coord(
                        &FieldElement([
                            2369245568431362,
                            2665603790611352,
                            3317390952748653,
                            1908583331312524,
                            8011773354506,
                        ]),
                        Choice::from(1u8),
                    )
                })
            }),
        );
    }

    pub fn bench_point_compression_decompression(c: &mut Criterion) {
        c.bench(
            "Compress/Decompress",
            Benchmark::new("Point compression.", |b| b.iter(|| P1_EXTENDED.compress())),
        );

        c.bench(
            "Compress/Decompress",
            Benchmark::new("Point decompression.", |b| {
                b.iter(|| P1_COMPRESSED.decompress().unwrap())
            }),
        );
    }
}

mod ecdh_benches {
    use super::*;
    use rand::rngs::OsRng;
    use scalar::Scalar;

    fn generate_kp() -> (Scalar, RistrettoPoint) {
        let sk = Scalar::random(&mut OsRng);
        let pk = RISTRETTO_BASEPOINT * sk;

        (sk, pk)
    } 

    fn ecdh() -> () {
        let alice_kp = generate_kp();
        let bob_kp = generate_kp();

        let alice_computes_S = bob_kp.1 * alice_kp.0;
        let bob_computes_S = alice_kp.1 * bob_kp.0;
    }

    pub fn bench_ecdh(c: &mut Criterion) {
        c.bench(
                "ECDH key exchange",
                Benchmark::new("Key Exchange.", |b| b.iter(|| ecdh())),
            );
    }

    criterion_group! {
        name = ecdh_benches;
        config = Criterion::default();
        targets =
        bench_ecdh
    }
}

criterion_group!(
    benchmarks,
    field_benches::bench_field_element_ops,
    field_benches::bench_modular_inverse,
    scalar_benches::bench_scalar_element_ops,
    edwards_benches::bench_extended_point_ops,
    edwards_benches::bench_projective_point_ops,
    edwards_benches::bench_point_compression_decompression,
);
criterion_main!(
    ecdh_benches::ecdh_benches,
    benchmarks,
);
