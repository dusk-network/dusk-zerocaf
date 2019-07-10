#![doc(html_logo_url = "https://lh3.googleusercontent.com/SmwswGxtgIANTbDrCOn5EKcRBnVdHjmYsHYxLq2HZNXWCQ9-fZyaea-bNgdX9eR0XGSqiMFi=w128-h128-e365")]
#![doc(html_favicon_url = "https://dusk.network/lib/img/favicon-16x16.png")]
//! <a href="https://codecov.io/gh/dusk-network/Dusk-Zerocaf">
//! <img src="https://codecov.io/gh/dusk-network/Dusk-Zerocaf/branch/master/graph/badge.svg" />
//! </a>
//! <a href="https://travis-ci.com/dusk-network/dusk-zerocaf">
//! <img src="https://travis-ci.com/dusk-network/dusk-zerocaf.svg?branch=master" />
//! </a>
//! <a href="https://github.com/dusk-network/dusk-zerocaf">
//! <img alt="GitHub closed issues" src="https://img.shields.io/github/issues-closed/dusk-network/dusk-zerocaf.svg?color=4C4CFF">
//! </a>
//! <a href="https://crates.io/crates/zerocaf">
//! <img alt="Crates.io" src="https://img.shields.io/crates/v/zerocaf.svg?color=f07900">
//! </a>
//! </hr>
//! 
//! <div>
//! <img src="https://camo.githubusercontent.com/db129d98b9686d0db27a9fd27c8e54086b14a6a7/68747470733a2f2f692e696d6775722e636f6d2f496a61645a50592e6a7067" width="100%" />
//! </div>
//! <br>
//! 
//! # What is Zerocaf?
//! Zerocaf is a pure Rust cryptographic library constructed to define operations for an elliptic curve embedded
//! into the Ristretto scalar field, which allows the construction of prime order groups
//! from an otherwise non prime order curve. <br>
//! 
//! The ultimate purpose of defining operations is for set inclusion proofs - where it is shown, in zero-knowledge, 
//! that a private key exists in a set of many public keys. <br>
//! 
//! Additionally, the zero-knowledge proofs use Bulletproofs as the argument for arithmetic circuits
//! that are used to form arbitrary constraint systems.
//! 
//! # What can it be used for?
//! The main goal of the library, as said before, is to be the base of operations over set inclusion proofs
//! and other Zero-Knowledge protocols.<br>
//! 
//! But since Zerocaf is build upon the Doppio Curve using the Ristretto protocol, it allows other devs to build
//! cryptographic protocols over it without needing to take care about the co-factor of the curve. <br>
//! 
//! This, brings to developers, a good mathematical backend library which can be used as a `mid-level` API
//! for building all kinds of cryptographic protocols over it such as key agreement, signatures, 
//! anonymous credentials, rangeproofs... <br>
//! 
//! # Usage
//! To import the library as a dependency of your project, just add on your `Cargo.toml`s project file:
//! ```toml
//! zerocaf = "0.1.1"
//! ```
//! 
//! Then import the crate as:
//! ```rust
//! extern crate zerocaf;
//! ```
//! 
//! # Backends.
//! Zerocaf has been built following the [Curve25519-dalek](https://docs.rs/curve25519-dalek/1.2.1/curve25519_dalek/) library structure, which allows for multiple
//! backend implementations. All of the works are built to enable modularity.<br>
//! 
//! Currently, `Zerocaf` has implemented the u64 backend.
//! By default, the `u64` backend is the one which is used to perform all of
//! the operations.
//! Additionly, for future works, we would like to implement a `u32` backend aswell. <br>
//! 
//! To select a backend type, the following method can be used:
//! ```sh
//! // For unoptimized builds:
//! cargo build --features "u64_backend"
//! 
//! // For optimized/release builds:
//! cargo build --release --features "u64_backend"
//! ``` 
//! <br>
//! 
//! NOTE: If no backend is selected, the compilation will fail!<br>
//! 
//! # Security and features of Zerocaf
//! 
//! As is previously mentioned, zerocaf is designed to host the fastest possible curve operations whilst
//! simultaneously avoiding all of the drawbacks associated with having a cofactor such that h > 1.<br>
//! 
//! To achieve this we make use of Ristretto, which is a technique to construct prime order elliptic curve groups.
//! The Ristretto protocol compresses the cofactor by adding a thin abstraction layer to allow small changes
//! in code to ultimately omit the cofactor issues. <br>
//! 
//! This is achieved by having defining the twisted edwards curve over the ristretto scalar field, 
//! which means to perform every operation on the curve in modulo L, 
//! where L is the order of the ristretto scalar field. <br>
//! 
//! `L =  2^252 + 27742317777372353535851937790883648493`. <br> 
//! 
//! By expounding the operations in this manner, we can benefit from the speed of a non-prime order twisted 
//! edwards curve whilst not suffering the pitfalls of a cofactor greater than one.
//! 
//! 
//! # Performance & Benchmarks
//! Benchmarks have been implemented using [Criterion.rs](https://docs.rs/criterion/0.2.11/criterion/).
//! To run them just execute `cargo bench` on the repository root.<br>
//! 
//! All of the operatons have been implemented using bit-shifting techniques to allow better performance
//! and a significant reduction in execution time.  
//! 
//! # Examples
//! We are planning to add some examples about tha basics of the `Zerocaf` library usage.<br>
//! They will be uploaded to the [examples](https://github.com/dusk-network/dusk-zerocaf/tree/master/docs) folder.<br>
//! 
//! This is a very basic usage example of the Zerocaf lib:
//! ```rust
//! extern crate zerocaf;
//! extern crate rand;
//! 
//! use zerocaf::field::FieldElement;
//! use zerocaf::scalar::Scalar;
//! use zerocaf::edwards::EdwardsPoint;
//! 
//! use rand::{Rng, thread_rng};
//! 
//! fn main() -> () {
//! 
//!     // Let G be an `EdwardsPoint` which is a point over the Twisted Eds Extended Coordinates. 
//!     let G: EdwardsPoint = EdwardsPoint {
//!         X: FieldElement([23, 0, 0, 0, 0]),
//!         Y: FieldElement([1664892896009688, 132583819244870, 812547420185263, 637811013879057, 13284180325998]),
//!         Z: FieldElement([1, 0, 0, 0, 0]),
//!         T: FieldElement([4351986304670635, 4020128726404030, 674192131526433, 1158854437106827, 6468984742885])
//!     };
//! 
//!     let scalar: Scalar = rand_scalar_generation();
//!     println!("{:?}", scalar);
//! 
//!     // Perform G*k, Point mul uses the `add_and_double` standard algorithm.
//!     let P = &G  * &scalar;
//!     println!("{:?}", P);
//! }
//! 
//! /// Generate a random `Scalar` defined over the sub-group field
//! /// modulo: `2^249 - 15145038707218910765482344729778085401`
//! pub fn rand_scalar_generation() -> Scalar {
//!     // Gen random 32-byte array. 
//!     let mut bytes = [0u8;32];
//! 
//!     // Fill the bytes varible with random bytes. We can use the 32 bytes co give
//!     // total randomness but then we will need to be aware because we can generate
//!     // values greater than `L = 2^252 + 27742317777372353535851937790883648493` and
//!     // the program will panic if we don't catch the error correctly on the 
//!     // `from_bytes()` Scalar method call.
//!     thread_rng().try_fill(&mut bytes[..31]).expect("Error getting the random bytes");
//! 
//!     Scalar::from_bytes(&bytes)
//! }
//! ```
//! <br>
//! We will also publish some videos talking about how is the library built and
//! the maths that are happening behind the scenes.<br>
//! Videos can also include programming examples using `Zerocaf` as a dependency.<br>
//! You can check them on the [Dusk Network Youtube Channel](https://www.youtube.com/channel/UCAfY3VcuaxAelPp44B253Rw).
//! 


// Used for traits related to constant-time code.
extern crate subtle;
// Used for Ristretto255Scalar trait.
extern crate curve25519_dalek;
extern crate num;


pub mod backend;
pub mod constants;
pub mod edwards;
pub mod field;
pub mod montgomery;
pub mod scalar;
pub mod traits;
