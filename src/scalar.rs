use crate::backend;

use subtle::Choice;
use subtle::ConstantTimeEq;



/// A `Scalar` represents an element of the field GF(l), optimized for speed.
///
/// This is a type alias for one of the Scalar types in the `backend`
/// module.
#[cfg(feature = "u64_backend")]
pub type Scalar = backend::u64::scalar::Scalar;

impl PartialEq for Scalar {
    fn eq(&self, other: &Scalar) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}

impl ConstantTimeEq for Scalar {
    /// Test equality between two `Scalar`s.  Since the
    /// internal representation is not canonical, the field elements
    /// are normalized to wire format before comparison.
    fn ct_eq(&self, other: &Scalar) -> Choice {
        self.to_bytes().ct_eq(&other.to_bytes())
    }
}

impl Eq for Scalar {}

/// This is a type alias for the Scalar type in the `curve25519-dalek` lib.
pub type Ristretto255Scalar = curve25519_dalek::scalar::Scalar;
