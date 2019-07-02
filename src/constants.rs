//! This module contains curve-specific constant values, 
//! needed to either testing or implementing algorithms.

use crate::scalar::Scalar;
use crate::field::FieldElement;

/// Edwards `a` variable value = `-1 (mod l)` equals:
/// `7237005577332262213973186563042994240857116359379907606001950938285454250988`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub (crate) static EDWARDS_A: FieldElement = FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]);

/// Edwards `d` variable value = `-1 (mod l)` equals:
/// `1201935917638644956968126114584555454358623906841733991436515590915937358637`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub (crate) static EDWARDS_D: FieldElement = FieldElement([939392471225133, 587442007554368, 4497154776428662, 4184267646867733, 2921744366591]);

/// Edwards helper value mentioned on: 2008 Hisil–Wong–Carter–Dawson, 
/// http://eprint.iacr.org/2008/522, Section 3.1.
/// EDWARDS_2_D_PRIME = `2d' = 2*(-a/d) = 9694442075178136385882072623406230049103477363454219155949387174793353019892`.
pub (crate) static EDWARDS_2_D_PRIME: FieldElement = FieldElement([4470874449035764, 1422240272386778, 4036862279045373, 3331298855343795, 23565883259442]);