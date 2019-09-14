//! Implementation of the Ristretto Protocol over the
//! Doppio curve.
//! 
//! Notes extracted from: https://ristretto.group/ristretto.html. 
//! Go there for the full lecture or check the paper here: 
//! https://tools.ietf.org/pdf/draft-hdevalence-cfrg-ristretto-00.pdf
//! 
//! # What's Ristretto?
//! Ristretto is a construction of a prime-order group using a non-prime-order Edwards curve.
//! The Decaf paper suggests using a non-prime-order curve E\mathcal EE to implement a prime-order
//! group by constructing a quotient group. Ristretto uses the same idea, but with different formulas,
//! in order to allow the use of cofactor-888 curves such as Curve25519.
//! 
//! Internally, a Ristretto point is represented by an Edwards point. 
//! Two Edwards points P,QP, QP,Q may represent the same Ristretto point, in the same way that 
//! different projective (X,Y,Z) coordinates may represent the same Edwards point. 
//! 
//! Group operations on Ristretto points are carried out with no overhead by performing the 
//! operations on the representative Edwards points.



use crate::constants;
use crate::edwards::{EdwardsPoint, double_and_add};
use crate::field::FieldElement;
use crate::scalar::Scalar;
use crate::traits::ops::*;
use crate::traits::Identity;

use core::ops::{Add, Mul, Neg, Index};

use std::fmt::Debug;

use subtle::{Choice, ConstantTimeEq, ConditionallyNegatable};


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
            return None
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
        let (ok, I) = (v*u2_sq).inv_sqrt();
        if ok.unwrap_u8() == 0 {return None};
        
        // Compute the Extended Point Coordinates Y & T
        let Dx = I*u2;
        let Dy = I*Dx*v;

        // Compute ABS(2*s*Dx) and negate if it is negative.
        let mut x = (&s + &s) * Dx;
        let x_is_pos = x.is_positive();
        x.conditional_negate(!x_is_pos);
        // Compute Y and T coordinates. 
        let y = u1 * Dy;
        let t = x * y;

        if t.is_positive().unwrap_u8() == 0u8 || y == FieldElement::zero() {
            return None
        };

        Some(RistrettoPoint(EdwardsPoint {
            X: x,
            Y: y,
            Z: one,
            T: t
        }))
    }
}

#[derive(Clone, Copy)]
pub struct RistrettoPoint (pub(crate) EdwardsPoint);

impl Debug for RistrettoPoint {
    fn fmt(&self, f: &mut ::core::fmt::Formatter) -> ::core::fmt::Result {
        write!(f, "{:?}", &self.0)
    }
}

impl ConstantTimeEq for RistrettoPoint {
    /// As specified on the Ristretto protocol docs: 
    /// https://ristretto.group/formulas/equality.html
    /// and we are on the twisted case, we compare 
    /// X1*Y2 == Y1*X2. 
    fn ct_eq(&self, other: &RistrettoPoint) -> Choice {
        (self.0.X * other.0.Y).to_bytes().ct_eq(
            &(self.0.Y * other.0.X).to_bytes())
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
            xy = 
                (constants::SQRT_MINUS_ONE * self.0.Y,
                constants::SQRT_MINUS_ONE * self.0.X);
            D = D1 * constants::INV_SQRT_A_MINUS_D;
        } else {
            xy = (self.0.X, self.0.Y);
            D = D2;
        };

        xy.1.conditional_negate(!(xy.0*Zinv).is_positive());
        // We are on the Twisted case, so a = -1. 
        // Then s = ABS((Z-Y) * D)
        let mut s = (self.0.Z - xy.1) * D;
        s.conditional_negate(!s.is_positive());

        CompressedRistretto(s.to_bytes()) 
    }

    /// Computes the Ristretto Elligator map.
    ///
    /// # Note
    ///
    /// This method is not public because it's just used for hashing
    /// to a point -- proper elligator support is deferred for now.
    pub(crate) fn elligator_ristretto_flavor(r_0: &FieldElement) -> RistrettoPoint {
        unimplemented!()
    }

    /// Debugging function used to get the 4coset where a point
    /// lives.
    pub(self) fn coset4(&self) -> [EdwardsPoint; 4] {
        self.0.coset4()
    }
}

mod tests {
    use hex;
    use super::*;

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
}

