//! This module contains backend-specific constant values, such as the 64-bit limbs of curve constants.

use crate::backend::u64::field::FieldElement;
use crate::backend::u64::scalar::Scalar;

/// `R` is the order of base point, i.e. 2^249 - 15145038707218910765482344729778085401
pub(crate) const R: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370495, 2199023255551]);