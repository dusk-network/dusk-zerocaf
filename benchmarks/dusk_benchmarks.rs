#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;
extern crate doppio_curve;

use criterion::{Criterion, Benchmark};

use doppio_curve::backend::u64::constants;
use doppio_curve::backend::u64::scalar::Scalar;



mod scalar_benches {   
    use super::*; 
    // B = 904625697166532776746648320197686575422163851717637391703244652875051672039
    pub static B: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370493, 2199023255551]);

    // BA = B - A = 904625697166532776746648320014998870755800986942176787613709275418060104167
    pub static BA: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370491, 2199023255551]);
    /*
    // benchmarking mul_internal function with function calls
    pub fn mul_int_func(c: &mut Criterion) {
        // B * BA computed in Sage limb by limb. (Since we don't have any other way to verify it.
        c.bench_function("mul_internal_function", move |b| b.iter(|| Scalar::mul_internal(&B, &BA)));
    }  
    // benchmarking mul_internal function with macro calls
    pub fn mul_int_macros(c: &mut Criterion) {
        // B * BA computed in Sage limb by limb. (Since we don't have any other way to verify it.
        c.bench_function("mul_internal_macros", move |b| b.iter(|| Scalar::mul_internal_macros(&B, &BA)));
    }
    */
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