#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;


use criterion::{Criterion, Benchmark};

use corretto::backend::u64::scalar::Scalar;



mod scalar_benches {   
    use super::*; 
    // B = 904625697166532776746648320197686575422163851717637391703244652875051672039
    pub static B: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    // BA = B - A = 904625697166532776746648320014998870755800986942176787613709275418060104167
    pub static BA: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370491, 2199023255551]);
    
    // Test if some implementation performs much better than the other even `inline`
    // forces to replace the code too. So both implementations should perform similarly.
    pub fn bench_mul_internal(c: &mut Criterion) {
    c.bench(
        "mul_internal",
        Benchmark::new("Function", |b| b.iter(|| Scalar::mul_internal(&B, &BA)))
            .with_function("Macro", |b| b.iter(|| Scalar::mul_internal_macros(&B, &BA))),
    );
}

}
criterion_group!(benchmarks, scalar_benches::bench_mul_internal);
criterion_main!(benchmarks);