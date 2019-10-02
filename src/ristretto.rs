#![allow(dead_code)]
#![allow(non_snake_case)]
//! Implementation of the Ristretto Protocol over the
//! Doppio curve.
//!
//! Notes extracted from: https://ristretto.group/ristretto.html.
//! Go there for the full lecture or check the paper here:
//! https://tools.ietf.org/pdf/draft-hdevalence-cfrg-ristretto-00.pdf
//!
//! The code wa originaly created by Isis Agora Lovecruft and
//! Henry de Valence [here](https://github.com/dalek-cryptography/curve25519-dalek/blob/master/src/ristretto.rs)
//!
//! # What's Ristretto?
//! Ristretto is a construction of a prime-order group using a non-prime-order Edwards curve.
//! The Decaf paper suggests using a non-prime-order curve E\mathcal EE to implement a prime-order
//! group by constructing a quotient group. Ristretto uses the same idea, but with different formulas,
//! in order to allow the use of cofactor-888 curves such as Curve25519.
//!
//! Internally, a Ristretto point is represented by an Edwards point.
//! Two Edwards points `P, Q` may represent the same Ristretto point, in the same way that
//! different projective (X,Y,Z) coordinates may represent the same Edwards point.
//!
//! Group operations on Ristretto points are carried out with no overhead by performing the
//! operations on the representative Edwards points.

use crate::constants;
use crate::edwards::{double_and_add, EdwardsPoint};
use crate::field::FieldElement;
use crate::scalar::Scalar;
use crate::traits::ops::*;
use crate::traits::{Identity, ValidityCheck};

use core::ops::{Add, Index, Mul, Neg};

use std::fmt::Debug;

use rand::{CryptoRng, Rng};
use subtle::{Choice, ConditionallyNegatable, ConditionallySelectable, ConstantTimeEq};

/// Ristretto Point expressed in wire format.
/// Since the Ristretto bytes encoding is canonical,
/// two points are equal if their encodin form is equal.
#[derive(Debug, Clone, Copy)]
pub struct CompressedRistretto(pub [u8; 32]);

impl Index<usize> for CompressedRistretto {
    type Output = u8;
    fn index(&self, _index: usize) -> &u8 {
        &(self.0[_index])
    }
}

impl ConstantTimeEq for CompressedRistretto {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.as_bytes().ct_eq(&other.as_bytes())
    }
}

impl PartialEq for CompressedRistretto {
    fn eq(&self, other: &CompressedRistretto) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}

impl Eq for CompressedRistretto {}

impl Identity for CompressedRistretto {
    /// Returns the Identity point on `CompressedRistretto`
    /// format.
    fn identity() -> CompressedRistretto {
        CompressedRistretto([0u8; 32])
    }
}

impl CompressedRistretto {
    /// Get the bytes of the `CompressedRistretto` point.
    pub fn as_bytes(&self) -> [u8; 32] {
        self.0
    }

    pub fn copy_from_slice(bytes: &[u8]) -> CompressedRistretto {
        let mut inp = [0u8; 32];
        inp.copy_from_slice(bytes);
        CompressedRistretto(inp)
    }

    #[allow(non_snake_case)]
    /// Attempt to decompress a `CompressedRistretto` point.
    /// This proces is done following the formulas derived from the
    /// isogenies that suit for our curve selection.
    ///
    /// # Returns
    /// - If the decompression/decoding succeeds -> `Some(RistrettoPoint)`.
    /// - If the decompression/decoding fails -> `None`.
    pub fn decompress(&self) -> Option<RistrettoPoint> {
        // Step 1: Check that the byte-string is a valid FieldElement.

        // As Ristretto paper says: "If the implementation's field element
        // encoding function produces canonical outputs, one way to check
        // that s_bytes is a canonical encoding (in step 1) is to decode
        // s_bytes into sss, then re-encode sss into s_bytes_check, and ensure
        // that s_bytes == s_bytes_check.
        let s: FieldElement = FieldElement::from_bytes(&self.as_bytes());
        let s_check = s.to_bytes();
        let s_correct_enc = s_check.ct_eq(&self.as_bytes());
        let s_is_positive = s.is_positive();

        // If the byte-encoding was incorrect or the representation is
        // a negative `FieldElement` (according to the definition of
        // positive found on Decaf paper), return `None`.
        if s_is_positive.unwrap_u8() == 0u8 || s_correct_enc.unwrap_u8() == 0u8 {
            return None;
        };
        // Step 2: Attempt to decompress the CompressedRistretto.
        let one = FieldElement::one();

        // u1 = 1 + as² with a = -1.
        let u1 = one - s.square();
        // u2 = 1 - as² with a = -1.
        let u2 = one + s.square();
        let u2_sq = u2.square();

        // v = a*d*u1² - u2²
        let v = -(constants::EDWARDS_D * u1.square()) - u2_sq;
        // I = 1/sqrt(v*u2²), returns `None` if the sqrt does not exist.
        let (ok, I) = (v * u2_sq).inv_sqrt();
        if ok.unwrap_u8() == 0 {
            return None;
        };

        // Compute the Extended Point Coordinates Y & T
        let Dx = I * u2;
        let Dy = I * Dx * v;

        // Compute ABS(2*s*Dx) and negate if it is negative.
        let mut x = (s + s) * Dx;
        let x_is_pos = x.is_positive();
        x.conditional_negate(!x_is_pos);
        // Compute Y and T coordinates.
        let y = u1 * Dy;
        let t = x * y;

        if t.is_positive().unwrap_u8() == 0u8 || y == FieldElement::zero() {
            return None;
        };

        Some(RistrettoPoint(EdwardsPoint {
            X: x,
            Y: y,
            Z: one,
            T: t,
        }))
    }
}

#[derive(Clone, Copy)]
pub struct RistrettoPoint(pub(crate) EdwardsPoint);

impl Debug for RistrettoPoint {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "{:?}", &self.coset4())
    }
}

impl ConstantTimeEq for RistrettoPoint {
    /// As specified on the Ristretto protocol docs:
    /// https://ristretto.group/formulas/equality.html
    /// and we are on the twisted case, we compare
    /// `X1*Y2 == Y1*X2 | X1*X2 == Y1*Y2`.
    fn ct_eq(&self, other: &RistrettoPoint) -> Choice {
        let a = (self.0.X * other.0.Y).ct_eq(&(self.0.Y * other.0.X));
        let b = (self.0.X * other.0.X).ct_eq(&(self.0.Y * other.0.Y));
        a | b
    }
}

impl PartialEq for RistrettoPoint {
    fn eq(&self, other: &RistrettoPoint) -> bool {
        self.ct_eq(&other).unwrap_u8() == 1u8
    }
}

impl Eq for RistrettoPoint {}

impl Identity for RistrettoPoint {
    /// Gives back the Identity point for the Extended Edwards Coordinates
    /// which is endoded as a `RistrettoPoint` with coordinates:
    /// `(X, Y, Z, T)` = `(0, 1, 1, 0)`.
    fn identity() -> RistrettoPoint {
        RistrettoPoint(EdwardsPoint::identity())
    }
}

impl Default for RistrettoPoint {
    /// Gives back the Identity point for the Extended Edwards Coordinates
    /// which is endoded as a `RistrettoPoint` with coordinates:
    /// `(X, Y, Z, T)` = `(0, 1, 1, 0)`.
    fn default() -> RistrettoPoint {
        RistrettoPoint::identity()
    }
}

//TODO: Review RistrettoPoint original implementation correctness: #83.
impl ValidityCheck for RistrettoPoint {
    /// A valid `RistrettoPoint` should have exactly
    /// order `L` (Scalar Field Order) and also
    /// verify the curve equation.
    ///
    /// This trait is mostly implemented for debugging purposes.
    ///
    /// # Returns
    /// - `Choice(1) if the point has order L (not 2L, 4L or 8L) &
    /// satisfies the curve equation.
    /// - `Choice(0) if the point does not satisfy one of the conditions
    /// mentioned avobe.
    fn is_valid(&self) -> Choice {
        // Verify that the point has order `L` (Sub group order).
        let has_order_l = (self.0 * constants::L).ct_eq(&EdwardsPoint::identity());
        has_order_l & self.0.is_valid()
    }
}

impl<'a> Neg for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    /// Negates a `RistrettoPoint` giving it's negated representation
    /// as a result.
    ///
    /// Since the negative of a point is (-X:Y:Z:-T), it
    /// gives as a result: `(-X:Y:Z:-T)`.
    fn neg(self) -> RistrettoPoint {
        RistrettoPoint(-self.0)
    }
}

impl Neg for RistrettoPoint {
    type Output = RistrettoPoint;
    /// Negates a `RistrettoPoint` giving it's negated representation
    /// as a result.
    ///
    /// Since the negative of a point is (-X:Y:Z:-T), it
    /// gives as a result: `(-X:Y:Z:-T)`.
    fn neg(self) -> RistrettoPoint {
        RistrettoPoint(-&self.0)
    }
}

impl<'a, 'b> Add<&'a RistrettoPoint> for &'b RistrettoPoint {
    type Output = RistrettoPoint;
    /// Performs the addition of two RistrettoPoints following the
    /// Twisted Edwards Extended Coordinates formulae but using as
    /// `d` parameter the `RISTRETTO_D` instead of the `EDWARDS_D`.
    ///
    /// This implementation is specific for curves with `a = -1` as
    /// the isomorphic twist is for Doppio.
    ///
    /// [Source: 2008 Hisil–Wong–Carter–Dawson],
    /// (http://eprint.iacr.org/2008/522), Section 3.1.
    fn add(self, other: &'a RistrettoPoint) -> RistrettoPoint {
        RistrettoPoint(&self.0 + &other.0)
    }
}

impl Add<RistrettoPoint> for RistrettoPoint {
    type Output = RistrettoPoint;
    /// Performs the addition of two RistrettoPoints following the
    /// Twisted Edwards Extended Coordinates formulae but using as
    /// `d` parameter the `RISTRETTO_D` instead of the `EDWARDS_D`.
    ///
    /// This implementation is specific for curves with `a = -1` as
    /// the isomorphic twist is for Doppio.
    ///
    /// [Source: 2008 Hisil–Wong–Carter–Dawson],
    /// (http://eprint.iacr.org/2008/522), Section 3.1.
    fn add(self, other: RistrettoPoint) -> RistrettoPoint {
        &self + &other
    }
}

impl<'a> Double for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    /// Performs the point doubling operation
    /// ie. `2*P` over the Twisted Edwards Extended
    /// Coordinates.
    ///
    /// This implementation is specific for curves with `a = -1` as
    /// the isomorphic twist is.
    /// Source: 2008 Hisil–Wong–Carter–Dawson,
    /// http://eprint.iacr.org/2008/522, Section 3.1.
    /// Cost: 4M+ 4S+ 1D
    fn double(self) -> RistrettoPoint {
        RistrettoPoint(self.0.double())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a RistrettoPoint {
    type Output = RistrettoPoint;
    /// Scalar multiplication: compute `self * Scalar`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    ///
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004).
    /// Guide to Elliptic Curve Cryptography.
    /// Springer Professional Computing. New York: Springer-Verlag.
    fn mul(self, scalar: &'b Scalar) -> RistrettoPoint {
        double_and_add(self, scalar)
    }
}

impl Mul<Scalar> for RistrettoPoint {
    type Output = RistrettoPoint;
    /// Scalar multiplication: compute `self * Scalar`.
    /// This implementation uses the algorithm:
    /// `add_and_doubling` which is the standard one for
    /// this operations and also adds less constraints on
    /// R1CS.
    ///
    /// Hankerson, Darrel; Vanstone, Scott; Menezes, Alfred (2004).
    /// Guide to Elliptic Curve Cryptography.
    /// Springer Professional Computing. New York: Springer-Verlag.
    fn mul(self, scalar: Scalar) -> RistrettoPoint {
        &self * &scalar
    }
}

impl RistrettoPoint {
    /// Encode a Ristretto point represented by the point `(X:Y:Z:T)`
    /// in extended coordinates.
    #[allow(non_snake_case)]
    pub fn compress(&self) -> CompressedRistretto {
        let u1 = (self.0.Z + self.0.Y) * (self.0.Z - self.0.Y);
        let u2 = self.0.X * self.0.Y;
        let (_, I) = (u1 * u2.square()).inv_sqrt();
        let D1 = u1 * I;
        let D2 = u2 * I;
        let Zinv = D1 * D2 * self.0.T;
        let mut xy;
        let D;
        if (self.0.T * Zinv).is_positive().unwrap_u8() == 0u8 {
            xy = (
                constants::SQRT_MINUS_ONE * self.0.Y,
                constants::SQRT_MINUS_ONE * self.0.X,
            );
            D = D1 * constants::INV_SQRT_A_MINUS_D;
        } else {
            xy = (self.0.X, self.0.Y);
            D = D2;
        };

        xy.1.conditional_negate(!(xy.0 * Zinv).is_positive());
        // We are on the Twisted case, so a = -1.
        // Then s = ABS((Z-Y) * D)
        let mut s = (self.0.Z - xy.1) * D;
        s.conditional_negate(!s.is_positive());

        CompressedRistretto(s.to_bytes())
    }

    /// Computes the Ristretto Elligator map.
    /// This gets a `RistrettoPoint` from a given
    /// `FieldElement´.
    pub(crate) fn elligator_ristretto_flavor(r_0: &FieldElement) -> RistrettoPoint {
        let d = constants::EDWARDS_D;
        let one = FieldElement::one();
        let mut c = -one;
        // 1 - d^2
        let one_minus_d_sq = one - d.square();

        // r = i*r0^2
        let r = constants::SQRT_MINUS_ONE * r_0.square();
        // Ns = a(r+1)*(a+d)*(a-d)
        let N_s = (r + one) * one_minus_d_sq;
        // D = (d*r -a)*(a*r -d)
        let D = (c - (d * r)) * (r + d);
        // s = sqrt(Ns/D)
        let (Ns_D_is_sq, mut s) = N_s.sqrt_ratio_i(&D);

        //s' = -ABS(s*r0)
        let mut s_prim = &s * r_0;
        s_prim.conditional_negate(s_prim.is_positive());

        s.conditional_assign(&s_prim, !Ns_D_is_sq);
        c.conditional_assign(&r, !Ns_D_is_sq);
        // Nt = c(r-1)*(d-1)^2 - D
        let N_t = ((c * (r - one)) * (d - one).square()) - D;
        let s_square = s.square();

        // Get the `CompletePoint` coordinates.
        let W0 = (s + s) * D;
        let W1 = N_t * constants::SQRT_AD_MINUS_ONE;
        let W2 = one - s_square;
        let W3 = one + s_square;

        // Get the `EdwardsPoint` that comes from the
        // `CompletePoint` obtained by the original
        // algorithm.
        RistrettoPoint(EdwardsPoint {
            X: W0 * W3,
            Y: W2 * W1,
            Z: W1 * W3,
            T: W0 * W2,
        })
    }

    /// Debugging function used to get the 4coset where a point
    /// lives.
    pub(self) fn coset4(&self) -> [EdwardsPoint; 4] {
        self.0.coset4()
    }

    /// Construct a `RistrettoPoint` from 64 bytes of data.
    ///
    /// If the input bytes are uniformly distributed, the resulting
    /// point will be uniformly distributed over the group, and its
    /// discrete log with respect to other points should be unknown.
    ///
    /// # Implementation
    ///
    /// This function splits the input array into two 32-byte halves,
    /// takes the low 255 bits of each half mod p, applies the
    /// Ristretto-flavored Elligator map to each, and adds the results.
    ///
    /// This function is taken from the Ristretto255 implementation found
    /// in [curve25519-dalek](https://github.com/dalek-cryptography/curve25519-dalek/blob/cf03d39f0fc3e1c625b9f1e9be0473758b324526/src/ristretto.rs#L713)
    pub fn from_uniform_bytes(bytes: &[u8; 64]) -> RistrettoPoint {
        let mut r_1_bytes = [0u8; 32];
        r_1_bytes.copy_from_slice(&bytes[0..32]);
        let r_1 = FieldElement::from_bytes(&r_1_bytes);
        let R_1 = RistrettoPoint::elligator_ristretto_flavor(&r_1);

        let mut r_2_bytes = [0u8; 32];
        r_2_bytes.copy_from_slice(&bytes[32..64]);
        let r_2 = FieldElement::from_bytes(&r_2_bytes);
        let R_2 = RistrettoPoint::elligator_ristretto_flavor(&r_2);

        // Applying Elligator twice and adding the results ensures a
        // uniform distribution.
        R_1 + R_2
    }

    /// Generate a random `RistrettoPoint` from a 64-byte array generated
    /// with user-provided rng.
    ///
    /// The provided `rng` has to implement: `Rng` + `CryptoRng`.
    ///
    /// This function uses the elligator hash map twice, once for [0..31] &
    /// another for [32..64] giving a uniformly distributed random value.
    ///
    /// This implementation follows the idea pointed on the
    /// random point generation used in [curve25519-dalek](https://github.com/dalek-cryptography/curve25519-dalek).
    pub fn new_random_point<T: Rng + CryptoRng>(rand: &mut T) -> RistrettoPoint {
        let mut bytes = [0u8; 64];
        rand.try_fill(&mut bytes).unwrap();
        RistrettoPoint::from_uniform_bytes(&bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "rand")]
    use rand::rngs::OsRng;

    #[test]
    fn basepoint_compr_decompr() {
        let compress = RistrettoPoint(constants::BASEPOINT).compress();

        let decompress = compress.decompress().unwrap();
        assert!(decompress == RistrettoPoint(constants::BASEPOINT));
    }

    #[test]
    fn valid_encoding_test_vectors() {
        // The following are the byte encodings of small multiples
        //     [0]B, [1]B, ..., [15]B
        // of the basepoint, represented as hex strings.
        let encodings_of_small_multiples = [
            // This is the identity point
            "0000000000000000000000000000000000000000000000000000000000000000",
            // This is the basepoint
            "0200000000000000000000000000000000000000000000000000000000000000",
            // These are small multiples of the basepoint
            "abe4ea98eaaeda5a9c63879cb3c4d9b4a01ed31ac383acefd7ed49861e1a8002",
            "1064fe35b16525f90f1d2f7d3dc448ba31a118f136c53eed88c2e951f1832907",
            "a826cf66461dea21e51187dddd8753299b726a7d4217cb75758aefbf5a2d4f01",
            "4d2e0705a9b47d122f98bd74808d386cf1691bc5407af703dd0c4808038b7f07",
            "f3a3592fde5fa05a881b80b4e732b37c32c7f684a5be33cdb8b7bdaf53db6f04",
            "51626c7960da63010efc5e064e62962f158f59928914fc108257ec2653745e01",
            "d5f8144c1b04954291785be578633a79131752e82afb990bd4a25b41cbd49001",
            "1372ed81add54633970746cd4b38ceb8a3e538b916288ac3d7c0dfbd54a42b06",
            "a83d7a262a80926724a0beb75a5f26e9a622205e6a64730e14ce64c4b2acf704",
            "a6b2712a6e586ab552f7bcf438168304b8b8a3f3b2852a06ae183e6303406503",
            "7876266b939b889c1da827a76da5c220eb1ff934472d35de60c9e4c3528fcc06",
            "11a0f75ab351572b572c38bf073b076aa964cdff70d53ad7588174dae2729306",
            "64f2fb80b45fbf73793e9e8e509f98848ecdb452c98c83c55c5c31fb233d9907",
            "1de5afbe9fd279f1651306d8ac0f68f0cb2689609ccfe8db1636f9481a33e205",
        ];

        // Test the encodings of small multiples

        let B = constants::RISTRETTO_BASEPOINT;
        let mut P = RistrettoPoint::identity();
        for i in 0..16 {
            assert_eq!(
                hex::encode(P.compress().as_bytes()),
                encodings_of_small_multiples[i],
            );
            P = P + B;
        }
    }

    #[test]
    fn decompress_id() {
        use crate::edwards::CompressedEdwardsY;

        let compressed_id = CompressedRistretto::identity();
        let id = compressed_id.decompress().unwrap();
        let mut identity_in_coset = false;
        for P in &id.0.coset4() {
            if P.compress() == CompressedEdwardsY::identity() {
                identity_in_coset = true;
            }
        }
        assert!(identity_in_coset);
    }

    #[test]
    fn four_torsion_diff() {
        use crate::edwards::{mul_by_pow_2, CompressedEdwardsY};

        let bp_compr_decompr = constants::RISTRETTO_BASEPOINT
            .compress()
            .decompress()
            .unwrap()
            .0;

        // Check that bp_compr_decompr differs from the
        // original RistrettoBasepoint by a point of order 4.
        let point_order_4 = &constants::RISTRETTO_BASEPOINT.0 - &bp_compr_decompr;

        let verif = mul_by_pow_2(&point_order_4, 2);
        assert_eq!(verif.compress(), CompressedEdwardsY::identity());
    }

    #[cfg(feature = "rand")]
    #[test]
    fn four_torsion_diff_random() {
        let mut rng = OsRng::new().unwrap();
        let B = &constants::RISTRETTO_BASEPOINT;
        let P = B * &Scalar::random(&mut rng);
        let P_coset = P.coset4();
        for i in 0..4 {
            assert_eq!(P, RistrettoPoint(P_coset[i]));
        }

        // Check that P_compr_decompr differs from the
        // original P by a point of order 4
        let P_compr_decompr = P.compress().decompress().unwrap();
        let point_order_4 = &P + &-P_compr_decompr;
        assert!(mul_by_pow_2(&point_order_4, &2) == RistrettoPoint::identity())
    }

    #[test]
    fn four_coset_eq_basepoint() {
        let basepoint = constants::RISTRETTO_BASEPOINT;
        let basep_coset = basepoint.coset4();

        for coset_point in &basep_coset {
            assert!(RistrettoPoint(*coset_point) == basepoint);
        }
    }

    #[test]
    fn validity_check() {
        // RISTRETTO_BASEPOINT should be valid.
        assert!(constants::RISTRETTO_BASEPOINT.is_valid().unwrap_u8() == 1u8);
        // The identity and multiples of the basepoint should also be valid.
        let mut P = RistrettoPoint::identity();
        let basep = constants::RISTRETTO_BASEPOINT;
        for _i in 0..16 {
            assert!(P.is_valid().unwrap_u8() == 1u8);
            P = P + basep;
        }

        // This point has order `8L` is a valid `EdwardsPoint`
        // but it's not a valid `RistrettoPoint`.
        let y_coord_bytes_8L = FieldElement::from_bytes(&[
            177, 118, 250, 81, 30, 181, 58, 122, 224, 214, 112, 52, 50, 60, 95, 199, 213, 167, 143,
            108, 154, 218, 242, 27, 175, 111, 152, 152, 213, 211, 157, 15,
        ]);
        let point_8L =
            EdwardsPoint::new_from_y_coord(&y_coord_bytes_8L, Choice::from(0u8)).unwrap();
        assert!(point_8L.is_valid().unwrap_u8() == 1u8);
        assert!(RistrettoPoint(point_8L).is_valid().unwrap_u8() == 0u8);
    }

    #[cfg(feature = "rand")]
    #[test]
    fn random_point_validity() {
        let mut rng = OsRng::new().unwrap();
        for i in 0..100 {
            let P = RistrettoPoint::new_random_point(&mut rng);
            // Check that the resulting `EdwardsPoint` relies on the curve.
            assert!(P.0.is_valid().unwrap_u8() == 1u8);
            P.compress().decompress();
        }
    }

    #[test]
    fn elligator_vs_ristretto_sage() {
        // This test uses the Sage script `ristretto.sage` located in the
        // `curve25519-dalek` repository in order to get test vectors of the
        // ristretto_elligator algorithm.
        let expected_point = RistrettoPoint(EdwardsPoint {
            X: FieldElement([
                520984263488427,
                2866053035698784,
                356812350072736,
                1177086814167286,
                17585355348321,
            ]),
            Y: FieldElement([
                2224110940152212,
                767723869121786,
                2519083920383090,
                3478258567033985,
                6072297619626,
            ]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([
                3761248848988017,
                3474827148739807,
                3137090891116602,
                1521420215868592,
                8052069914602,
            ]),
        });
        assert!(expected_point.is_valid().unwrap_u8() == 1u8);

        let raw_bytes =
            hex::decode("2e2d7c6f887c81c1593f32e2fa31a7b65d4fbbf38f8ab3045ead22fc45743219")
                .unwrap();
        let mut bytes = [0u8; 32];
        bytes.copy_from_slice(&raw_bytes);
        let point_from_ellig =
            RistrettoPoint::elligator_ristretto_flavor(&FieldElement::from_bytes(&bytes));

        assert!(point_from_ellig.0.is_valid().unwrap_u8() == 1u8);
        assert!(point_from_ellig == expected_point);
        assert!(point_from_ellig.compress() == expected_point.compress())
    }
}
