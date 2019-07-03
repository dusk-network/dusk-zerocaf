//! Contains the curve-constants needed by the different algorithm implementations.

use crate::scalar::Scalar;
use crate::field::FieldElement;

/// Edwards `a` variable value = `-1 (mod l)` equals:
/// `7237005577332262213973186563042994240857116359379907606001950938285454250988`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_A: FieldElement = FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]);

/// Edwards `d` variable value = `86649/86650 (mod l)` equals:
/// `6035069659693617257005060448458438786498492452538173614565435347369516892352`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_D: FieldElement = FieldElement([4236121989480640, 3329222317550656, 6444852309635, 319331980502762, 14670441677824]);

/// Edwards helper value mentioned on: [2008 Hisil–Wong–Carter–Dawson, ]
/// (http://eprint.iacr.org/2008/522), Section 3.1.
/// EDWARDS_2_D_PRIME = `2d' = 2*(-a/d) = 4833133742054972300036934333873883332139868545696439623128919756453579533715`.
pub static EDWARDS_2_D_PRIME: FieldElement = FieldElement([3296729518255507, 2741780309996288, 12889703251469, 638663961005524, 11748697311232]);