//! This module contains backend-specific constant values as the 64-bit limbs of curve constants.

use crate::backend::u64::field::FieldElement;
use crate::backend::u64::scalar::Scalar;
use crate::edwards::*;
use crate::ristretto::{RistrettoPoint, CompressedRistretto};


/// `L` is the order of base point for Sonny, in this case it is equivalent to 2^249 + 14490550575682688738086195780655237219
pub const L: Scalar = Scalar([1129677152307299, 1363544697812651, 714439, 0, 2199023255552]);

/// `(L - 1) / 2` used to check positiveness of a `FieldElement` on the Decaf paper. 
pub(crate) const POS_RANGE: FieldElement = FieldElement([2587757230352886, 4210131976237760, 683900, 0, 8796093022208]);

/// Scalar-LFACTOR is the value that satisfies the equation: `L * LFACTOR = -1 (mod 2^52)`
/// In this case, `LFACTOR` is the one used for the Montgomery Reduction algorithm,
/// implemented on Scalar Arithmetics module.
pub const LFACTOR: u64 = 1331240223835829;

/// Montgomery modulus defined for Scalar arithmetics, `R^2 = (2^260)^2 % L`
pub const RR: Scalar = Scalar([137682194168839, 3209056245311277, 1480926248458276, 2533620989757837, 1314911199310]);

/// `FIELD_L` is the order of the Prime field for Sonny, in this case it is equivalent to 2^252 + 27742317777372353535851937790883648493`
pub const FIELD_L: FieldElement = FieldElement([671914833335277, 3916664325105025, 1367801, 0, 17592186044416]);

/// Montgomery modulus defined for FieldElement arithmetics, `R^2 = (2^260)^2 % FIELD_L`
pub const RR_FIELD: FieldElement = FieldElement([2764609938444603, 3768881411696287, 1616719297148420, 1087343033131391, 10175238647962]);

/// FieldElement-LFACTOR is the value that satisfies the equation: `L * LFACTOR = -1 (mod 2^52)`
/// In this case, `LFACTOR` is the one used for the Montgomery Reduction algorithm,
/// implemented on FieldElement Arithmetics module.
pub const LFACTOR_FIELD: u64 = 1439961107955227;

/// Montgomery modulus defined for FieldElements on `inverse()` functions scope. 
/// It is used for the Montgomery Mul operation that takes place on the `Inversion
// operation`. It's defined as: `R^2 = (2^253)^2 % L`
pub const INV_RR: FieldElement = FieldElement([2210115751650724, 3809421927348411, 2357176729341513, 3420097284349172, 7483527818736]);

/// Edwards `a` variable value = `-1 (mod l)` equals:
/// `7237005577332262213973186563042994240857116359379907606001950938285454250988`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_A: FieldElement = FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]);

/// Edwards `d` variable value = `-126296/126297 (mod l)` equals:
/// `951605751702391019481481818669129158712512026257330939079110344917983315091`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_D: FieldElement = FieldElement([3304133203739795, 2446467598308289, 1534112949566882, 2032729967918914, 2313225441931]);

/// Holds the value of one of both `sqrt(-1 (mod p)) values. 
/// `SQRT_MINUS_ONE = 3034649101460298094273452163494570791663566989388331537498831373842135895065`. 
pub const SQRT_MINUS_ONE: FieldElement = FieldElement([3075585030474777, 2451921961843096, 1194333869305507, 2218299809671669, 7376823328646]); 

/// `(+)1/SQRT(a) (mod l)` equals: `4202356475871964119699734399548423449193549369991576068503119564443318355924`.
pub static INV_SQRT_A: FieldElement = FieldElement([2099929430230996, 1464742363261928, 3309265759432790, 2285299817698826, 10215362715769]);

/// `(-)SQRT(a) (mod l)` equals: `4202356475871964119699734399548423449193549369991576068503119564443318355924`. 
pub static MINUS_SQRT_A: FieldElement = FieldElement([2099929430230996, 1464742363261928, 3309265759432790, 2285299817698826, 10215362715769]);

/// `INV_SQRT_A_MINUS_D = 482283834104289360917429750399313974390948281833312135312952165682596457149`.
pub const INV_SQRT_A_MINUS_D: FieldElement = FieldElement([550050132044477, 3953042081665262, 2971403105229349, 212915494370164, 1172367057772]);

/// `SQRT_AD_MINUS_ONE = `.
pub const SQRT_AD_MINUS_ONE : FieldElement = FieldElement([3601277882726560, 1817821323014817, 1726005090908779, 2111284621343800, 648674458156]);

/// 4Coset of a RistrettoPoint. 
pub(crate) const FOUR_COSET_GROUP: [EdwardsPoint; 4] = 
    [
        EdwardsPoint {
            X: FieldElement([1, 0, 0, 0, 0]),
            Y: FieldElement([0, 0, 0, 0, 0]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([0, 0, 0, 0, 0])
        }, 

        EdwardsPoint {
            X: FieldElement([2099929430230996, 1464742363261928, 3309265759432790, 2285299817698826, 10215362715769]),
            Y: FieldElement([0, 0, 0, 0, 0]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([0, 0, 0, 0, 0])
        }, 

        EdwardsPoint {
            X: FieldElement([0, 0, 0, 0, 0]),
            Y: FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([0, 0, 0, 0, 0])
        },
        
        EdwardsPoint {
            X: FieldElement([3075585030474777, 2451921961843096, 1194333869305507, 2218299809671669, 7376823328646]),
            Y: FieldElement([0, 0, 0, 0, 0]),
            Z: FieldElement([1, 0, 0, 0, 0]),
            T: FieldElement([0, 0, 0, 0, 0])
        },
    ];

/// Holds the value of the Curve basepoint, which has been constructed
/// from taking `y-coodrinate = 3/5 (mod l)`.
pub const BASEPOINT: EdwardsPoint = EdwardsPoint {
    X: FieldElement([276718085098056, 1646536057461434, 2704687245600312, 2630386667454967, 13476148227069]),
    Y: FieldElement([1303868825475266, 3250718520537114, 2702159777242978, 2702159776422297, 10555311626649]),
    Z: FieldElement([1, 0, 0, 0, 0]),
    T: FieldElement([3634527586288175, 2006028620404053, 3424252198034825, 2478951925947079, 4567251727358])
};

/// Ristretto Basepoint. 
pub const RISTRETTO_BASEPOINT: RistrettoPoint = RistrettoPoint(BASEPOINT);
