#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;


use criterion::{Criterion, Benchmark};

use zerocaf::backend::u64::{scalar, field};
use zerocaf::edwards::EdwardsPoint;
#[allow(unused_imports)]
use zerocaf::traits::Identity;



mod scalar_benches {   
    use super::*; 
    // B = 904625697166532776746648320197686575422163851717637391703244652875051672039
    pub static B: scalar::Scalar = scalar::Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    // BA = B - A = 904625697166532776746648320014998870755800986942176787613709275418060104167
    pub static BA: scalar::Scalar = scalar::Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370491, 2199023255551]);
    
    // Test if some implementation performs much better than the other even `inline`
    // forces to replace the code too. So both implementations should perform similarly.
    pub fn bench_mul_internal(c: &mut Criterion) {
    c.bench(
        "mul_internal",
        Benchmark::new("Function", |b| b.iter(|| scalar::Scalar::mul_internal(&B, &BA)))
            .with_function("Macro", |b| b.iter(|| scalar::Scalar::mul_internal_macros(&B, &BA))),
    );
    }
}

mod field_benches {

    use super::*;

    /// `B = 904625697166532776746648320197686575422163851717637391703244652875051672039`
    pub static B: field::FieldElement = field::FieldElement([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    pub fn bench_modular_inverse(c: &mut Criterion) {
    c.bench(
        "Modular Inverse",
        Benchmark::new("Kalinski inverse", |b| b.iter(|| field::FieldElement::kalinski_inverse(&B))).
            with_function("Savas & Koç inverse", |b| b.iter(|| field::FieldElement::savas_koc_inverse(&B))),
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

criterion_group!(benchmarks, scalar_benches::bench_mul_internal, field_benches::bench_modular_inverse, edwards_benches::bench_point_addition);
//criterion_group!(benchmarks, field_benches::bench_modular_inverse);
criterion_main!(benchmarks);