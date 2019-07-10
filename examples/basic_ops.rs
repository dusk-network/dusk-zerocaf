#![allow(non_snake_case)]
/// The purpose of this implementation is to provide support for one of the most 
/// commonly used operations over EC which is Random Scalar Mul. 
/// 
/// This is one of the first steps to take on the implementation of algoritnms
/// like Diffie-Hellman Key Exchange.
/// 
/// Let G be a point over Doppio, let k = Scalar random value over the Doppio Sub-group. 
/// Let P = G*k;

extern crate zerocaf;
extern crate rand;

use zerocaf::field::FieldElement;
use zerocaf::scalar::Scalar;
use zerocaf::edwards::EdwardsPoint;

use rand::{Rng, thread_rng};

fn main() -> () {

    // Let G be an `EdwardsPoint` which is a point over the Twisted Eds Extended Coordinates. 
    let G: EdwardsPoint = EdwardsPoint {
        X: FieldElement([23, 0, 0, 0, 0]),
        Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
        Z: FieldElement([1, 0, 0, 0, 0]),
        T: FieldElement([4351986304670635, 4020128726404030, 674192131526433, 1158854437106827, 6468984742885])
    };

    let scalar: Scalar = rand_scalar_generation();
    println!("{:?}", scalar);

    // Perform G*k, Point mul uses the `add_and_double` standard algorithm.
    let P = &G  * &scalar;
    println!("{:?}", P);
}

/// Generate a random `Scalar` defined over the sub-group field
/// modulo: `2^249 - 15145038707218910765482344729778085401`
pub fn rand_scalar_generation() -> Scalar {
    // Gen random 32-byte array. 
    let mut bytes = [0u8;32];

    // Fill the bytes varible with random bytes. We can use the 32 bytes co give
    // total randomness but then we will need to be aware because we can generate
    // values greater than `L = 2^252 + 27742317777372353535851937790883648493` and
    // the program will panic if we don't catch the error correctly on the 
    // `from_bytes()` Scalar method call.
    thread_rng().try_fill(&mut bytes[..31]).expect("Error getting the random bytes");

    Scalar::from_bytes(&bytes)
}

