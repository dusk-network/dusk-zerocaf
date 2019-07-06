#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;


use criterion::{Criterion, Benchmark};

use zerocaf::backend::u64::{scalar, field};
use zerocaf::edwards::EdwardsPoint;
#[allow(unused_imports)]
use zerocaf::traits::{Identity, Square};



mod field_benches {

    use super::*;
    use zerocaf::field::FieldElement;
    use zerocaf::scalar::Scalar;

    /// `B = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static B: field::FieldElement = field::FieldElement([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    /// `A = 182687704666362864775460604089535377456991567872`
    pub static A: FieldElement = FieldElement([0, 0, 0, 2, 0]);


    pub fn bench_field_element_ops(c: &mut Criterion) {
        c.bench(
            "Field Element",
            Benchmark::new("Addition", |b| b.iter(|| &B + &A))
        );

        c.bench(
            "Field Element",
            Benchmark::new("Subtraction", |b| b.iter(|| &B - &A))
        );

        c.bench(
            "Field Element",
            Benchmark::new("Mul", |b| b.iter(|| &B * &A))
        );

        c.bench(
            "Field Element",
            Benchmark::new("Squaring", |b| b.iter(|| B.square()))
        );

        c.bench(
            "Field Element",
            Benchmark::new("Half", |b| b.iter(|| A.half()))
        );

        c.bench(
            "Field Element",
            Benchmark::new("Two Pow k (2^k)", |b| b.iter(|| FieldElement::two_pow_k(&213u64)))
        );
    }

    pub fn bench_modular_inverse(c: &mut Criterion) {
        c.bench(
            "Modular Inverse",
            Benchmark::new("Savas & Koç Modular Inverse algorithm", |b| b.iter(|| field::FieldElement::inverse(&B)))
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
    pub static D: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    pub fn bench_scalar_element_ops(c: &mut Criterion) {
        c.bench(
            "Scalar",
            Benchmark::new("Addition", |b| b.iter(|| &C + &D))
        );

        c.bench(
            "Scalar",
            Benchmark::new("Subtraction", |b| b.iter(|| &C - &D))
        );

        c.bench(
            "Scalar",
            Benchmark::new("Mul", |b| b.iter(|| &C * &D))
        );

        c.bench(
            "Scalar",
            Benchmark::new("Squaring", |b| b.iter(|| C.square()))
        );

        c.bench(
            "Scalar",
            Benchmark::new("Half", |b| b.iter(|| C.half()))
        );
    }
}

mod edwards_benches {

    use super::*;

    pub fn bench_point_addition(c: &mut Criterion) {
    c.bench(
        "Extended Coordinates Point Addition",
        Benchmark::new("2008 Hisil–Wong–Carter–Dawson, Section 3.1.", |b| b.iter(|| &EdwardsPoint::identity() + &EdwardsPoint::identity()))
    );
    }
}

criterion_group!(benchmarks, 
                field_benches::bench_field_element_ops,
                field_benches::bench_modular_inverse,
                scalar_benches::bench_scalar_element_ops,
                edwards_benches::bench_point_addition);
//criterion_group!(benchmarks, field_benches::bench_modular_inverse);
criterion_main!(benchmarks);