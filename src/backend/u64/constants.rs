//! This module contains backend-specific constant values, such as the 64-bit limbs of curve constants.

use crate::backend::u64::field::FieldElement;
use crate::backend::u64::scalar::Scalar;

/// `L` is the order of base point for Doppio, i.e. 2^249 - 15145038707218910765482344729778085401
pub(crate) const L: Scalar = Scalar([2766226127823335, 4237835465749098, 4503599626623787, 4503599627370495, 2199023255551]);

/// `L` * `LFACTOR` = -1 (mod 2^52)
pub(crate) const LFACTOR: u64 = 547593343082025;

/// `R^2` = (2^260)^2 % R
pub(crate) const RR: Scalar = Scalar([1682248870925813, 4078880264703668, 2289123149127681, 4169238435752846, 2104335664921]);

/// `FIELD_L` is the order of the Prime field for Doppio, i.e. 2^252 + 27742317777372353535851937790883648493`
pub(crate) const FIELD_L: FieldElement = FieldElement([671914833335277, 3916664325105025, 1367801, 0, 17592186044416]);

/// `L_FIELD` * `LFACTOR_FIELD` = -1 (mod 2^52)
pub(crate) const LFACTOR_FIELD: u64 = 1439961107955227;