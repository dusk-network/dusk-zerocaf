#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;

use criterion::{Benchmark, criterion_group, criterion_main, Criterion, BenchmarkId};

use zerocaf::backend::u64::{field, scalar};
use zerocaf::edwards::*;
use zerocaf::ristretto::RistrettoPoint;
use zerocaf::traits::ops::*;
use zerocaf::constants::*;
use zerocaf::field::*;
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
        let inp = (A, B);
        let inp2 = (A, Choice::from(1u8));
        let pow = 249u64;

        c.bench_with_input(
            BenchmarkId::new("Addition", "Fixed constants"), &inp , |b, &inp| {
                b.iter(|| inp.1 + inp.0);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Subtraction", "Fixed constants"), &inp , |b, &inp| {
                b.iter(|| inp.0 - inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Multiplication", "Fixed constants"), &inp , |b, &inp| {
                b.iter(|| inp.0 * inp.1);
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Square", "Fixed constants"), &inp , |b, &inp| {
                b.iter(|| inp.1.square());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Half", "Fixed constants"), &inp , |b, &inp| {
                b.iter(|| inp.0.half());
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Two pow K", "Fixed constant"), &pow , |b, &pow| {
                b.iter(|| FieldElement::two_pow_k(pow))
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Power", "Fixed constants"), &inp , |b, &inp| {
                b.iter(|| inp.0.pow(&inp.1));
            }
        );

        c.bench_with_input(
            BenchmarkId::new("Modular Sqrt", "Fixed constants"), &inp2 , |b, &inp| {
                b.iter(|| inp2.0.mod_sqrt(inp2.1));
            }
        );
        c.bench_with_input(
            BenchmarkId::new("Modular inverse", "Fixed constants"), &inp , |b, &inp| {
                b.iter(|| inp.0.inverse());
            }
        );
    }
}

criterion_group!(
    benchmarks,
    field_benches::bench_field_element_ops,
);
criterion_main!(
    benchmarks,
);
