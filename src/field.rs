//! A `FieldElement` represents an element of the finite field 
//! modulo `2^252 + 27742317777372353535851937790883648493`.
//! 
//! The `FieldElement` type is an alias for one of the platform-specific
//! implementations. 
//! 
//! The most commonly used on all architectures will be the `u64` backend
//! which is the one that has been implemented at this moment.

use core::cmp::PartialEq;

use subtle::Choice;
use subtle::ConstantTimeEq;

use crate::backend;


#[cfg(feature = "u64_backend")]
pub use backend::u64::field::*;
/// A `FieldElement` represents an element of the field 
/// `2^252 + 27742317777372353535851937790883648493`
/// 
/// The `FieldElement` type is an alias for one of the platform-specific
/// implementations.
#[cfg(feature = "u64_backend")]
pub type FieldElement = backend::u64::field::FieldElement;


impl PartialEq for FieldElement {
    fn eq(&self, other: &FieldElement) -> bool {
        self.ct_eq(other).unwrap_u8() == 1u8
    }
}

impl ConstantTimeEq for FieldElement {
    /// Test equality between two `FieldElement`s.  Since the
    /// internal representation is not canonical, the field elements
    /// are normalized to wire format before comparison.
    fn ct_eq(&self, other: &FieldElement) -> Choice {
        self.to_bytes().ct_eq(&other.to_bytes())
    }
}

