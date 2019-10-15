#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;

use criterion::{Benchmark, criterion_group, criterion_main, Criterion, BenchmarkId};

use zerocaf::edwards::*;
use zerocaf::ristretto::*;
use zerocaf::traits::ops::*;
use zerocaf::constants::*;
use zerocaf::field::*;
use zerocaf::scalar::*;
#[allow(unused_imports)]
use zerocaf::traits::Identity;

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
        let inp = (C, D);
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
    }
}

criterion_group!(
    benchmarks,
    field_benches::bench_field_element_ops,
    scalar_benches::bench_scalar_element_ops,
);
criterion_main!(
    benchmarks,
);
