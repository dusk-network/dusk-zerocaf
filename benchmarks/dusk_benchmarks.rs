#![allow(non_snake_case)]

#[macro_use]
extern crate criterion;


use criterion::{Criterion, Benchmark};

use zerocaf::backend::u64::{scalar, field};
use zerocaf::edwards::EdwardsPoint;

use zerocaf::traits::ops::*;

#[allow(unused_imports)]
use zerocaf::traits::Identity;



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
    use zerocaf::edwards::{EdwardsPoint, ProjectivePoint};
    use zerocaf::scalar::Scalar;
    use zerocaf::field::FieldElement;

    pub static P1_EXTENDED: EdwardsPoint = EdwardsPoint {
        X: FieldElement([23, 0, 0, 0, 0]),
        Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([4351986304670635, 4020128726404030, 674192131526433, 1158854437106827, 6468984742885])
    };

    pub static P2_EXTENDED: EdwardsPoint = EdwardsPoint {
        X: FieldElement([68, 0, 0, 0, 0]),
        Y: FieldElement([1799957170131195, 4493955741554471, 4409493758224495, 3389415867291423, 16342693473584]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([3505259403500377, 292342788271022, 2608000066641474, 796697979921534, 2995435405555])
    };

    pub static P1_PROJECTIVE: ProjectivePoint = ProjectivePoint {
        X: FieldElement([23, 0, 0, 0, 0]),
        Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
        Z: FieldElement([1, 0, 0, 0, 0])
    };

    pub static P2_PROJECTIVE: ProjectivePoint = ProjectivePoint {
        X: FieldElement([68, 0, 0, 0, 0]),
        Y: FieldElement([1799957170131195, 4493955741554471, 4409493758224495, 3389415867291423, 16342693473584]),
        Z: FieldElement([1, 0, 0, 0, 0])
    };



    /// `A = 182687704666362864775460604089535377456991567872`.
    pub static A: Scalar = Scalar([0, 0, 0, 2, 0]);

    pub fn bench_extended_point_ops(c: &mut Criterion) {
        c.bench(
            "Extended Coordinates Point Addition",
            Benchmark::new("2008 Hisil–Wong–Carter–Dawson, Section 3.1.", |b| b.iter(|| &P1_EXTENDED + &P2_EXTENDED))
        );

        c.bench(
            "Extended Coordinates Point Subtraction",
            Benchmark::new("2008 Hisil–Wong–Carter–Dawson, Section 3.1.", |b| b.iter(|| &P1_EXTENDED - &P2_EXTENDED))
        );

        c.bench(
            "Extended Coordinates Point Doubling",
            Benchmark::new("2008 Hisil–Wong–Carter–Dawson, Section 3.1.", |b| b.iter(|| P1_EXTENDED.double()))
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
    }
}

criterion_group!(benchmarks, 
                field_benches::bench_field_element_ops,
                field_benches::bench_modular_inverse,
                scalar_benches::bench_scalar_element_ops,
                edwards_benches::bench_extended_point_ops,
                edwards_benches::bench_projective_point_ops);
//criterion_group!(benchmarks, field_benches::bench_modular_inverse);
criterion_main!(benchmarks);