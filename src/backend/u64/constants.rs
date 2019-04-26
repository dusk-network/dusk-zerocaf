//! This module contains backend-specific constant values, such as the 64-bit limbs of curve constants.

use crate::backend::u64::field::FieldElement;
use crate::backend::u64::scalar::Scalar;

/// `R` is the order of base point, i.e. 2^249 - 15145038707218910765482344729778085401
pub(crate) const R: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370495, 2199023255551]);

/// `L` * `LFACTOR` = -1 (mod 2^52)
pub(crate) const LFACTOR: u64 = 0x51da312547e1b;

/// `R^2` = (2^260)^2
pub(crate) const RR: Scalar = Scalar([1682248870925813, 4078880264703668, 2289123149127681, 4169238435752846, 2104335664921]);