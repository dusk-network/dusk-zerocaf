#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;


use criterion::{Criterion, Benchmark};

use zerocaf::backend::u64::{scalar, field};
use zerocaf::edwards::EdwardsPoint;
#[allow(unused_imports)]
use zerocaf::traits::Identity;



mod field_benches {

    use super::*;

    /// `B = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static B: field::FieldElement = field::FieldElement([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    pub fn bench_modular_inverse(c: &mut Criterion) {
    c.bench(
        "Modular Inverse",
        Benchmark::new("Kalinski inverse", |b| b.iter(|| field::FieldElement::inverse(&B)))
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

criterion_group!(benchmarks, field_benches::bench_modular_inverse, edwards_benches::bench_point_addition);
//criterion_group!(benchmarks, field_benches::bench_modular_inverse);
criterion_main!(benchmarks);