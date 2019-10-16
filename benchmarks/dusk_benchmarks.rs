#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;

use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

use zerocaf::edwards::*;
use zerocaf::ristretto::*;
use zerocaf::traits::ops::*;
use zerocaf::constants::*;
use zerocaf::field::*;
use zerocaf::scalar::*;
#[allow(unused_imports)]
use zerocaf::traits::Identity;

use rand::rngs::OsRng;

mod field_benches {

    use super::*;
    use subtle::Choice;

    /// `B = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static B: FieldElement = FieldElement([
        2766226127823335,
        4237835465749098,
        4503599626623787,
        4503599627370493,
        2199023255551,
    ]);

    /// `A = 182687704666362864775460604089535377456991567872`
    pub static A: FieldElement = FieldElement([0, 0, 0, 2, 0]);

    pub fn bench_field_element_ops(c: &mut Criterion) {
        let inp = (A, B);
        let inp2 = (A, Choice::from(1u8));
        let pow = 249u64;

        c.bench_with_input(
            BenchmarkId::new("Addition", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.1 + inp.0);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Subtraction", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.0 - inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Multiplication", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.0 * inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Division", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.0/inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Square", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.1.square());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Half", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.0.half());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Two pow K", "Fixed constant"), &pow , |b, &pow| {
                b.iter(|| FieldElement::two_pow_k(pow))
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Power", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.0.pow(&inp.1));
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Legendre Symbol", "Fixed FieldElement"), &inp , |b, &inp| {
                b.iter(|| inp.0.legendre_symbol());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Modular Sqrt", "Fixed FieldElements"), &inp2 , |b, &inp| {
                b.iter(|| inp2.0.mod_sqrt(inp2.1));
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Modular inverse", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.0.inverse());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Sqrt ratio i", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.0.sqrt_ratio_i(&inp.1));
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Inverse sqrt", "Fixed FieldElements"), &inp , |b, &inp| {
                b.iter(|| inp.0.inv_sqrt());
            }
        );

        c.bench_function("Random FieldElement generation", |b| b.iter(|| FieldElement::random(&mut OsRng)));
    }
}

mod scalar_benches {
    use super::*;
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
        let inp = (C, D, 4u8);
        let pow = 215u64;

        c.bench_with_input(
            BenchmarkId::new("Addition", "Fixed Scalars"), &inp , |b, &inp| {
                b.iter(|| inp.0 + inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Subtraction", "Fixed Scalars"), &inp , |b, &inp| {
                b.iter(|| inp.0 - inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Mul", "Fixed Scalars"), &inp , |b, &inp| {
                b.iter(|| inp.0 * inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Squaring", "Fixed Scalars"), &inp , |b, &inp| {
                b.iter(|| inp.1.square());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Half", "Fixed Scalars"), &inp , |b, &inp| {
                b.iter(|| inp.0.half());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Two pow K", "Fixed Scalar"), &pow , |b, &pow| {
                b.iter(|| Scalar::two_pow_k(pow));
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Power", "Fixed Scalars"), &inp , |b, &inp| {
                b.iter(|| inp.0.pow(&inp.1));
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Into-bits conversion", "Fixed Scalar"), &inp , |b, &inp| {
                b.iter(|| inp.1.into_bits());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Shr", "Fixed Scalar"), &inp , |b, &inp| {
                b.iter(|| inp.1 >> 135);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Scalar mod 4", "Fixed Scalar"), &inp , |b, &inp| {
                b.iter(|| inp.1.mod_2_pow_k(inp.2));
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Compute NAF", "Fixed Scalar"), &inp , |b, &inp| {
                b.iter(|| inp.1.compute_NAF());
            }
        );

        c.bench_function("Random Scalar generation", |b| b.iter(|| Scalar::random(&mut OsRng)));
    }
}

mod edwards_benches {

    use super::*;
    use subtle::Choice;

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

    /// `D = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static D: Scalar = Scalar([
        2766226127823335,
        4237835465749098,
        4503599626623787,
        4503599627370493,
        2199023255551,
    ]);

    /// `P1_EXTENDED on `CompressedEdwardsY` format.
    pub(self) static P1_COMPRESSED: CompressedEdwardsY = CompressedEdwardsY([
        206, 11, 225, 231, 113, 39, 18, 141, 213, 215, 201, 201, 90, 173, 14, 134, 192, 119, 133,
        134, 164, 26, 38, 1, 201, 94, 187, 59, 186, 170, 240, 2,
    ]);

    pub fn bench_extended_point_ops(c: &mut Criterion) {

        let extend_inp = (P1_EXTENDED, P2_EXTENDED, D);
        let y_gen = (FieldElement([
                            2369245568431362,
                            2665603790611352,
                            3317390952748653,
                            1908583331312524,
                            8011773354506,
                        ]), Choice::from(1u8));
        
        c.bench_with_input(
            BenchmarkId::new("Extended Coordinates Point Addition", "Fixed Points"), &extend_inp , |b, &extend_inp| {
                b.iter(|| extend_inp.0 + extend_inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Extended Coordinates Point Subtraction", "Fixed Points"), &extend_inp , |b, &extend_inp| {
                b.iter(|| extend_inp.0 - extend_inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Extended Coordinates Point Doubling", "Fixed Points"), &extend_inp , |b, &extend_inp| {
                b.iter(|| extend_inp.0.double());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Extended Coordinates Point Multiplication", "Fixed Points"), &extend_inp , |b, &extend_inp| {
                b.iter(|| extend_inp.0 * extend_inp.2);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Extended Coordinates Point Generation", "Fixed y-coordinate"), &y_gen , |b, &y_gen| {
                b.iter(|| EdwardsPoint::new_from_y_coord(&y_gen.0, y_gen.1));
            }
        );

        c.bench_function("Random EdwardsPoint generation", |b| b.iter(|| EdwardsPoint::new_random_point(&mut OsRng)));
    }

    pub fn bench_projective_point_ops(c: &mut Criterion) {

        let proj_inp = (P1_PROJECTIVE, P2_PROJECTIVE, D);
        let y_gen = (FieldElement([
                            2369245568431362,
                            2665603790611352,
                            3317390952748653,
                            1908583331312524,
                            8011773354506,
                        ]), Choice::from(1u8));
        
        c.bench_with_input(
            BenchmarkId::new("Projective Coordinates Point Addition", "Fixed Points"), &proj_inp , |b, &proj_inp| {
                b.iter(|| proj_inp.0 + proj_inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Projective Coordinates Point Subtraction", "Fixed Points"), &proj_inp , |b, &proj_inp| {
                b.iter(|| proj_inp.0 - proj_inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Projective Coordinates Point Doubling", "Fixed Point"), &proj_inp , |b, &proj_inp| {
                b.iter(|| proj_inp.0.double());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Projective Coordinates Point Multiplication", "Fixed Points"), &proj_inp , |b, &proj_inp| {
                b.iter(|| proj_inp.0 * proj_inp.2);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Projective Coordinates Point Generation", "Fixed y-coordinate"), &y_gen , |b, &y_gen| {
                b.iter(|| ProjectivePoint::new_from_y_coord(&y_gen.0, y_gen.1));
            }
        );

        c.bench_function("Random ProjectivePoint generation", |b| b.iter(|| ProjectivePoint::new_random_point(&mut OsRng)));
    }

    pub fn bench_point_compression_decompression(c: &mut Criterion) {
        let cd_inp = (P1_COMPRESSED, P1_EXTENDED);
        
        c.bench_with_input(
            BenchmarkId::new("Point Compression", "Fixed Extended Point"), &cd_inp , |b, &cd_inp| {
                b.iter(|| cd_inp.1.compress());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Point Decompression", "Fixed Compressed-Point"), &cd_inp , |b, &cd_inp| {
                b.iter(|| cd_inp.0.decompress().unwrap());
            }
        );
    }
}

mod ristretto_benches {
    use super::*;

    /// `D = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static D: Scalar = Scalar([
        2766226127823335,
        4237835465749098,
        4503599626623787,
        4503599627370493,
        2199023255551,
    ]);

    pub fn bench_ristretto_point_ops(c: &mut Criterion) {

        let proj_inp = (RISTRETTO_BASEPOINT, D);
        
        c.bench_with_input(
            BenchmarkId::new("Ristretto Random Point Addition", "Fixed Points"), &proj_inp , |b, &proj_inp| {
                b.iter(|| proj_inp.0 + proj_inp.0);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Ristretto Random Point Subtraction", "Fixed Points"), &proj_inp , |b, &proj_inp| {
                b.iter(|| proj_inp.0 + -proj_inp.0);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Ristretto Random Point Doubling", "Fixed Point"), &proj_inp , |b, &proj_inp| {
                b.iter(|| proj_inp.0.double());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Ristretto Random Point Multiplication", "Fixed Points"), &proj_inp , |b, &proj_inp| {
                b.iter(|| proj_inp.0 * proj_inp.1);
            }
        );

        c.bench_function("Ristretto random Point generation", |b| b.iter(|| RistrettoPoint::new_random_point(&mut OsRng)));
    }

    pub fn bench_ristretto_protocol_impl(c: &mut Criterion) {

        let inputs = (RISTRETTO_BASEPOINT, 
                        FieldElement([
                            2369245568431362,
                            2665603790611352,
                            3317390952748653,
                            1908583331312524,
                            8011773354506,
                        ]),
                        RISTRETTO_BASEPOINT_COMPRESSED);

        c.bench_with_input(
            BenchmarkId::new("Ristretto Encoding", "Ristretto Basepoint"), &inputs , |b, &inputs| {
                b.iter(|| inputs.0.compress());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Ristretto Decoding", "Ristretto Basepoint"), &inputs , |b, &inputs| {
                b.iter(|| inputs.2.decompress().unwrap());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Ristretto Elligator", "Fixed FieldElement"), &inputs , |b, &inputs| {
                b.iter(|| RistrettoPoint::elligator_ristretto_flavor(&inputs.1));
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Ristretto Equalty", "Ristretto Basepoint"), &inputs , |b, &inputs| {
                b.iter(|| inputs.0 == inputs.0);
            }
        );
    }    
}

mod comparaisons {
    use super::*;
    use rand::rngs::OsRng;

    static P1_EXTENDED: EdwardsPoint = EdwardsPoint {
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

    /// `D = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static D: Scalar = Scalar([
        2766226127823335,
        4237835465749098,
        4503599626623787,
        4503599627370493,
        2199023255551,
    ]);

    pub fn bench_point_ops_impl(c: &mut Criterion) {
        let i = P1_EXTENDED;
        let mul = (P1_EXTENDED, D);

        // Equalty
        let mut group = c.benchmark_group("Equalty");

        group.bench_with_input(BenchmarkId::new("Compressing", "Fixed Point"), &i, 
            |b, &i| b.iter(|| i.compress() == i.compress()));
        group.bench_with_input(BenchmarkId::new("To Affine", "Fixed Point"), &i, 
            |b, &i| b.iter(|| i == i));
        
        group.finish();

        // Point Mul
        let mut group = c.benchmark_group("Point Multiplication");

        group.bench_with_input(BenchmarkId::new("Double And Add", "Fixed inputs"), &mul, 
            |b, &mul| b.iter(|| double_and_add(&mul.0, &mul.1)));
        group.bench_with_input(BenchmarkId::new("Left to right binary method", "Fixed inputs"), &mul, 
            |b, &mul| b.iter(|| ltr_bin_mul(&mul.0, &mul.1)));
        group.bench_with_input(BenchmarkId::new("NAF binary", "Fixed inputs"), &mul, 
            |b, &mul| b.iter(|| binary_naf_mul(&mul.0, &mul.1)));
        
        group.finish();
    }   

    // Helpers for benchmarking the whole ECDH process //

    // Left to right shift method.
    fn generate_kp_ltrs() -> (Scalar, RistrettoPoint) {
        let sk = Scalar::random(&mut OsRng);
        let pk = ltr_bin_mul(&RISTRETTO_BASEPOINT, &sk);

        (sk, pk)
    } 

    fn ecdh_ltrs() -> () {
        let alice_kp = generate_kp_ltrs();
        let bob_kp = generate_kp_ltrs();

        let alice_computes_S = ltr_bin_mul(&bob_kp.1, &alice_kp.0);
        let bob_computes_S = ltr_bin_mul(&alice_kp.1, &bob_kp.0);
    }

    // Binary NAF
    fn generate_kp_binary_naf() -> (Scalar, EdwardsPoint) {
        let sk = Scalar::random(&mut OsRng);
        let pk = binary_naf_mul(&RISTRETTO_BASEPOINT.0, &sk);

        (sk, pk)
    } 

    fn ecdh_binary_naf() -> () {
        let alice_kp = generate_kp_binary_naf();
        let bob_kp = generate_kp_binary_naf();

        let alice_computes_S = binary_naf_mul(&bob_kp.1, &alice_kp.0);
        let bob_computes_S = binary_naf_mul(&alice_kp.1, &bob_kp.0);
    }

    // Double and Add
    fn generate_kp_double_add() -> (Scalar, RistrettoPoint) {
        let sk = Scalar::random(&mut OsRng);
        let pk = double_and_add(&RISTRETTO_BASEPOINT, &sk);

        (sk, pk)
    } 

    fn ecdh_double_add() -> () {
        let alice_kp = generate_kp_double_add();
        let bob_kp = generate_kp_double_add();

        let alice_computes_S = double_and_add(&bob_kp.1, &alice_kp.0);
        let bob_computes_S = double_and_add(&alice_kp.1, &bob_kp.0);
    }

    

    /// ECDH bench function.
    pub fn bench_ecdh(c: &mut Criterion) {
        let mut group = c.benchmark_group("ECDH");

        group.bench_function("Double and Add", |b| b.iter(|| ecdh_double_add()));
        group.bench_function("Left to right binary method", |b| b.iter(|| ecdh_ltrs()));
        group.bench_function("Binary NAF method", |b| b.iter(|| ecdh_binary_naf()));
        
        group.finish();
    } 

}

criterion_group!(
    benchmarks,
    comparaisons::bench_point_ops_impl,
    comparaisons::bench_ecdh,
    field_benches::bench_field_element_ops,
    scalar_benches::bench_scalar_element_ops,
    edwards_benches::bench_extended_point_ops, 
    edwards_benches::bench_projective_point_ops, 
    edwards_benches::bench_point_compression_decompression,
    ristretto_benches::bench_ristretto_point_ops,
    ristretto_benches::bench_ristretto_protocol_impl,
);
criterion_main!(
    benchmarks,
);
