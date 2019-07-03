//! <a href="https://codecov.io/gh/dusk-network/Dusk-Zerocaf">
//! <img src="https://codecov.io/gh/dusk-network/Dusk-Zerocaf/branch/master/graph/badge.svg" />
//! </a>
//! <a href="https://travis-ci.com/dusk-network/dusk-zerocaf">
//! <img src="https://travis-ci.com/dusk-network/dusk-zerocaf.svg?branch=master" />
//! </a>
//! </hr>
//! 
//! <div>
//! <img src="https://camo.githubusercontent.com/db129d98b9686d0db27a9fd27c8e54086b14a6a7/68747470733a2f2f692e696d6775722e636f6d2f496a61645a50592e6a7067" width="100%" />
//! </div>
//! <hr/>


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
