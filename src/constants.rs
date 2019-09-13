//! Contains the curve-constants needed by different algorithm implementations.

use crate::field::FieldElement;
use crate::edwards::{EdwardsPoint, CompressedEdwardsY};
use crate::ristretto::{CompressedRistretto, RistrettoPoint};

#[cfg(feature = "u64_backend")]
pub use crate::backend::u64::constants::*;

/// Edwards `a` variable value = `-1 (mod l)` equals:
/// `7237005577332262213973186563042994240857116359379907606001950938285454250988`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_A: FieldElement = FieldElement([671914833335276, 3916664325105025, 1367801, 0, 17592186044416]);

/// Edwards `d` variable value = `-126296/126297 (mod l)` equals:
/// `951605751702391019481481818669129158712512026257330939079110344917983315091`
/// where `l = Prime of the field = 2^252 + 27742317777372353535851937790883648493`
pub static EDWARDS_D: FieldElement = FieldElement([3304133203739795, 2446467598308289, 1534112949566882, 2032729967918914, 2313225441931]);

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