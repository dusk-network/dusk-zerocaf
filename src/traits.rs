//! Module for Public Trait implementations.

use subtle::Choice;

/// Gives the Identity element for the
/// type which it has been implemented on.
///
/// This trait is implemented following the rules that
/// mandate over the Type that is being implemented.
pub trait Identity {
    #[must_use]
    /// Returns the identity element for the implemented
    /// type, following it's mathematical rules.
    fn identity() -> Self;
}

/// This trait pretends to be a verification in ct_time
/// about a point correctness.
///
/// This is done through checking that the (X, Y) coordinates
/// of the point are valid and satisfy the curve equation.
pub trait ValidityCheck {
    #[must_use]
    /// Checks the point coordinates agains the curve equation
    /// to validate that the point relies on the curve and its
    /// valid.
    ///
    /// # Returns
    /// - `Choice(0)` if the point isn't a valid point.
    /// - `Choice(1)` if the point is valid.
    fn is_valid(&self) -> Choice;
}

pub mod ops {
    use super::*;

    /// Trait that represents the `^2` operation for any
    /// kind of element on the library.
    ///
    /// This trait is implemented following the rules that
    /// mandate over the Type that is being implemented.
    pub trait Square {
        type Output;

        #[must_use]
        /// Returns the square of the input: `x^2`.
        fn square(self) -> Self::Output;
    }

    /// Trait that represents the Point doubling operation
    /// for any type of Point that is used on the lib.
    ///
    /// This trait is implemented following the rules that
    /// mandate over the Type that is being implemented.
    pub trait Double {
        type Output;

        #[must_use]
        /// Performs the point-doubling operation over the
        /// coordinates which this trait has been implemented
        /// for.
        fn double(self) -> Self::Output;
    }

    /// Trait that represents the `/2` operation for any
    /// kind of element on the library.
    ///
    /// This is a more performant way of performing the
    /// division by 2 that dividing by 2 with the `Div`
    /// trait implementation.
    ///
    /// This trait is implemented following the rules that
    /// mandate over the Type that is being implemented.
    pub trait Half {
        type Output;

        #[must_use]
        /// Returns the half of the input: `x/2`.
        fn half(self) -> Self::Output;
    }

    /// Trait that represents the modular exponentiation
    /// operation, ie.`a ^ b (mod l)`, for any
    /// kind of element on the library (except points).
    ///
    /// This trait is implemented following the rules that
    /// mandate over the Type that is being implemented.
    pub trait Pow<T> {
        type Output;

        #[must_use]
        /// Returns  `a^b (mod l)`.
        fn pow(self, exp: T) -> Self::Output;
    }

    pub trait ModSqrt {
        type Output;

        #[must_use]
        /// Performs the modular Square Root operation over a finite
        /// field ie. `sqrt(x) (mod l)`.
        ///
        /// With the given `Choice`, the impl is able to provide the
        /// result that corresponds to the positive or negative sign choosen.
        ///
        /// # Returns
        ///
        /// `Some(symb_choosen_result)` if the input is a QR for the prime modulo.
        /// Otherways it returns `None`
        fn mod_sqrt(self, choice: Choice) -> Self::Output;
    }

    pub trait InvSqrt {
        type Output;

        #[must_use]
        /// Performs the Inverse Square root of a given value.
        ///
        /// This operation returns always the positive result of the
        /// modular sqrt, understanding positive as the definition that
        /// appears on the Decaf paper: 0 < result < (P - 1)/2.  
        fn inv_sqrt(self) -> Self::Output;
    }

    pub trait SqrtRatioI<T> {
        type Output;

        #[must_use]
        /// Using the same trick as in ed25519 decoding, we merge the
        /// inversion, the square root, and the square test.AsMut
        ///
        /// The first part of the return value signals whether u/v was square,
        /// and the second part contains a square root.
        /// Specifically, it returns:
        ///
        ///- (true, +sqrt(u/v)) if v is nonzero and u/v is square;
        ///- (true, zero) if u is zero;
        ///- (false, zero) if v is zero and uuu is nonzero;
        ///- (false, +sqrt(i*u/v)) if u/v is nonsquare (so iu/v is square).
        fn sqrt_ratio_i(&self, v: T) -> Self::Output;
    }
}
