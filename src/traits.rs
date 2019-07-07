//! Module for Public Trait implementations.


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

pub mod ops {
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

    pub trait Double {
    type Output;

    #[must_use]
    /// Performs the point-doubling operation over the
    /// coordinates which this trait has been implemented
    /// for.
    fn double(self) -> Self::Output;
}
}