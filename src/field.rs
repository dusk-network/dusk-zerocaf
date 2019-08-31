//! A `FieldElement` represents an element of the finite field 
//! modulo `2^252 + 27742317777372353535851937790883648493`.
//! 
//! The `FieldElement` type is an alias for one of the backend
//! implementations. 
//! 
//! `ConstantTimeEq` and `PartialEq` traits have been implemented 
//! here since they will be the samme across all of the different backends.
//!  
//! # Examples
//! ```rust
//! use zerocaf::field::FieldElement;
//! use zerocaf::traits::ops::*;
//! use zerocaf::constants::EDWARDS_D;
//! use subtle::Choice;
//! 
//! // You can create a FieldElement from a byte-array as follows:
//! let a = FieldElement::from_bytes(&[0u8;32]); 
//! 
//! // You ca also create a FieldElement from an uint type as follows:
//! let b = FieldElement::from(&86649u128);
//! let c = FieldElement::from(&86650u64);
//! 
//! // You can create random FieldElements by calling: 
//! let rand = FieldElement::generate_random();
//! 
//! // The last way of creating a FieldElement it by calling the
//! // constructor. THIS IS NOT RECOMMENDED since NO checks about
//! // the correctness of the input will be done at all. 
//! // It can be done as follows: 
//! let d: FieldElement = FieldElement([0, 1, 0, 0, 0]); // d = 2^52.
//! assert!(d == FieldElement::two_pow_k(&52u64));
//! 
//! // All of the basuc modular operations are implemented 
//! // for FieldElement type:  
//! let mut res = &a + &b; // Performs a + b (mod l).
//! res = a - b; // Performs a - b (mod l).
//! res = a * b; // Performs a * b (mod l).
//! res = a.square(); // Performs a^2 (mod l).
//! res = -&a; // Performs Negation over the modulo l.
//! res = a.pow(&b); // Performs Modular exponentiation.
//! res = a.mod_sqrt(Choice::from(1u8)).unwrap(); //Performs
//! // modular sqrt.
//! // Returs `None` if the input is not a QR on the field.
//! // Returns Some(result) if everything is correct.
//! 
//! // Division has been also implemented. Remember that when we write
//! // a/b (mod l), we are indeed performing a * inverse_mod(b, l) (mod l).
//! assert!(-&(b / c) == EDWARDS_D);
//! 
//! // Dividing by two even FieldElements is recommended through the `Half`
//! // trait implmementation since it's much faster.
//! if a.is_even() {
//!     let half_a = &a.half(); // This will panic if a isn't even.
//! };
//! 
//! // We can finally perform inversion modulo l for a FieldElement:
//! let inv_a = &c.inverse(); // Performs a^-1 (mod l).
//! 
//! // You can export your `FieldElement` as an slice of 32 bytes in Little
//! // Endian encoding by:
//! let c_bytes: [u8; 32] = c.to_bytes();
//! ```
//! 
//! `PartialOrd`, `Ord`, `PartialEq` and `Eq` are also implemented for
//! `FieldElement` type. 
//! 
//! All `std::core::ops traits -> (Add, Sub, Mul, Div)` are implemented
//! for both, `&FieldElement` and `FieldElement`.

use core::cmp::PartialEq;

use subtle::{Choice, ConstantTimeEq, ConditionallySelectable}; 

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

impl ConditionallySelectable for FieldElement {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        FieldElement([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
            u64::conditional_select(&a.0[2], &b.0[2], choice),
            u64::conditional_select(&a.0[3], &b.0[3], choice),
            u64::conditional_select(&a.0[4], &b.0[4], choice)
        ])
    }
}

