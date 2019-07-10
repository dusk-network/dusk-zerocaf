# Dusk-Zerocaf [![Build Status](https://travis-ci.com/dusk-network/Dusk-Zerocaf.svg?branch=master)](https://travis-ci.com/dusk-network/Dusk-Zerocaf)  [![codecov](https://codecov.io/gh/dusk-network/dusk-corretto/branch/master/graph/badge.svg)](https://codecov.io/gh/dusk-network/dusk-corretto) ![GitHub closed issues](https://img.shields.io/github/issues-closed/dusk-network/dusk-zerocaf.svg?color=4C4CFF) ![Crates.io](https://img.shields.io/crates/v/zerocaf.svg)

<img
 width="100%"
 src="https://camo.githubusercontent.com/db129d98b9686d0db27a9fd27c8e54086b14a6a7/68747470733a2f2f692e696d6775722e636f6d2f496a61645a50592e6a7067">

# WARNING: WIP Repo.

## Fast, efficient and bulletproof-friendly cryptographic operations.

This repository contains an implementation of the Doppio curve over the `Ristretto Scalar field`: a pure Rust implementation designed by [Dusk](https://dusk.network) team.

### Ristretto curve 

Ristretto is a technique for constructing prime order elliptic curve groups with non-malleable encodings. The [Ristretto protocol](https://ristretto.group/ristretto.html) arose as an extension of [Mike Hamburg's Decaf](https://www.shiftleft.org/papers/decaf/decaf.pdf) approach to cofactor elimination, which is applicable to curves of
cofactor 4, whereas the Ristretto is designed for non-prime-order curves of cofactor 8 or 4.

### Ristretto Scalar Field And Bulletproofs.

Originally designed to abstract _non-prime-order curves into prime-order scalar fields_, the `Ristretto` abstraction would have been far too inefficient to implement for Bulletproofs zero-knowledge proof. Therefore the `Ristretto scalar field` is used to **solve all negative impacts of using cofactors equalling 8 on the Ristretto curve.**. The strategy is to use a _Ristretto embedded curve_ (also called `Doppio Curve`), as the initial operations within `zerocaf` are performed therein. `zerocaf` opens up new opportunities for the use cases of **zero-knowledge proofs** inside the Dusk Network protocol as well as making a _Bulletproof-integrated ring signature substitute possible_, with orders of magnitude performance improvements compared to the fastest ringsig implementation.

Within this library, the implementation of the Ristretto to construct the curve with desired properties is made possible by 
defining the curve over the scalar field, using only a thin abstraction layer, which in turn allows for systems that use signatures to be safely extended with zero-knowledge protocols. These zero-knowledge protocols are utilised with no additional cryptographic assumptions and minimal changes in the code. The Ristretto scalar field is Bulletproof friendly, which makes it possible to use both cryptographic protocols in tandem with one another, as they are centric to contemporary applications of elliptic curve operations.

Special thanks to [@ebfull](https://github.com/ebfull) who triggered this work with the following tweet:

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">Here&#39;s an &quot;embedded&quot; curve over ristretto255&#39;s scalar field<br><br>-x^2 + y^2 = 1 - (86649/86650)x^2y^2<br><br>which is Ristretto-ready and birationally equivalent to<br><br>y^2 = x^3 + 346598x^2 + x (and it&#39;s twist secure)<br><br>Any other suggestions?</p>&mdash; Sean Bowe (@ebfull) <a href="https://twitter.com/ebfull/status/1087571257057406976?ref_src=twsrc%5Etfw">January 22, 2019</a></blockquote>

## Details

### Curve parameters:

| Variable | Value | Explanation |
|--|--|--|
| Equation | Edwards -x²+y²=1-$`\frac{86649}{86650}`$x²y² | -|
| a | -1 | - |
| d | $`\frac{86649}{86650}`$ | - |
| B | $`\frac{8}{9}`$ | Edwards Basepoint Y-coordinate With X > 0 | 

<br/>

| Montgomery | y²=x³+346598*x²+x |
|--|--|
| u(P) | 17 | `u` coordinate of the Montgomery Basepoint, X-coordinate | \
| A | 346598 | |

<br/>

| Weierstrass | y²=x³+ax+b |
|--|--|
| a | 2412335192444087404657728854347664746952372119793302535333983646055108025796 | |
| b | 1340186218024493002587627141304258192751317844329612519629993998710484804961 | |
| x | 2412335192444087404657728854347664746952372119793302535333983646095151532546 | |
| y | 6222320563903764848551041754877036140234555813488015858364752483591799173948 | |

| Variable | Value | Explanation |
|--|--|--|
| G | 2²⁵² - 121160309657751286123858757838224683208 | Curve order |
| p | 2²⁵² + 27742317777372353535851937790883648493 | Prime of the field |
| r | 2²⁴⁹ - 15145038707218910765482344729778085401 | Prime of the Sub-Group |\

<br/>

### Encoding / Decoding tools
In order to work with our points along the curve, or any non trivial computuaions, for example those with tough notations - there has been a set of tools and examples which have been created to make facilitate the Encoding/Decoding processes. These can be found at: `tools/src/main.rs` 

### Examples

```rust
num_from_bytes_le(&[76, 250, 187, 243, 105, 92, 117, 70, 234, 124, 126, 180, 87, 149, 62, 249, 16, 149, 138, 56, 26, 87, 14, 76, 251, 39, 168, 74, 176, 202, 26, 84]);
// Prints: 38041616210253564751207933125345413214423929536328854382158537130491690875468
    
let res = to_field_elem_51(&"1201935917638644956968126114584555454358623906841733991436515590915937358637");
println!("{:?}", res);
// Gives us: [939392471225133, 1174884015108736, 2226020409917912, 1948943783348399, 46747909865470]

hex_bytes_le("120193591763864495696812611458455545435862390684173399143651559091593735863735685683568356835683");
// Prints: Encoding result -> [63, 41, b7, c, b, 79, 94, 7b, 21, d2, fe, 7b, c8, 89, c9, 7f, 76, c8, 9b, a3, 58, 18, 39, a, f2, d2, 7c, 17, ed, 7f, 6, c4, 9d, 44, f3, 7c, 85, c2, 67, e]
// Put the 0x by yourseleves and if there's any value alone like `c` padd it with a 0 on the left like: `0x0c`

from_radix_to_radix_10("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", 16u32);
// Prints: 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787

```

> When performing operations with large values, such as: `2²⁵² - 121160309657751286123858757838224683208`, it is recomended to compute them through `SageMath`, as the user interface adheres to these types of functions. From `SageMath`, they can be converted in a consistent format and easily compiled into Rust.

### Roadmap:

Note: the refactoring relations are expressed as indentations


- [x] Build Scalar Arithmetics and Scalar Struct definition.
    - [x] Find the proper radix value for FieldElement.
    - [x] Add the required constants for computation.
      - [x] Implement Addition.
      - [x] implement Subtraction.
      - [x] Implement Byte-encoding/decoding.
      - [x] Implement From for uint native types.
      - [x] Implement Ord, PartialOrd, Eq & PartialEq.
      - [x] Implement Multiplication on u64-backend with u128 usage.
      - [x] Implement Squaring.
      - [x] Implement Half for even numbers.
      - [x] Implement Modular Negation.
      - [x] Implement Montgomery_reduction.
      - [x] Define Montgomery_reduction algorithm.
- [x] Create FieldElement Struct and implement the basic operations we need on a u64 backend.
  - [x] Find the proper radix value for FieldElement.
  - [x] Add basic and needed constants.
  - [x] Implement Reduce function to make the FieldElements fit on a 5 u64-bit limbs.
    - [x] Implement Addition.
    - [x] Implement Subtraction.
    - [x] Implement Byte-encoding/decoding.
    - [x] Implement From for uint native types.
    - [x] Implement Ord, PartialOrd, Eq & PartialEq.
    - [x] Implement Multiplication on u64-backend with u128 usage.
    - [x] Implement Squaring.
    - [x] Implement Half for even numbers
    - [x] Implement Modular Negation.
    - [x] Implement Montgomery_reduction.
    - [x] Define Montgomery_reduction algorithm.
    - [x] Implement Modular inversion.
    - [x] Research about addition chains inversion methods.
  - [x] Add proper tests for every function.
- [ ] Implement Edwards Points
    - [ ] Implement Twisted Edwards Extended Coordiantes.
       - [x] Implement Point Addition.
       - [x] Implement Point Subtraction.
       - [x] Implement Point Doubling.
       - [x] Implement Scalar Mul.
       - [ ] Implement from_bytes conversions.
       - [ ] Implement to byte conversions.
       - [ ] Implement compressed Edwards point Y-coordinate.
    - [ ] Implement Twisted Edwards Projective Coordiates.
       - [x] Implement Point Addition.
       - [x] Implement Point Subtraction.
       - [x] Implement Point Doubling.
       - [x] Implement Scalar Mul.
       - [ ] Implement from_bytes conversions.
       - [ ] Implement to byte conversions.
       - [ ] Implement compressed Edwards point Y-coordinate.
    - [ ] Represent Edwards points as Ristretto points using wrapping function (research).
    - [x] Cargo doc testing and improvement.
    - [ ] Decide the best use cases of the various Edwards coordinate types (compressed, standard, extended, projective).
    - [ ] Benchmark different implementations and algorithms.
    - [ ] Research About Niels and ProjectiveNiels coordinates usage.
    - [ ] Implement Ristretto Mapping.
    - [ ] Build and test torsion points.
    - [ ] Test all point operations for Edwards Points.
 - [ ] Implement Montgomery and Edwards operations & functions.
 
 
> Operations with large numbers are recommended to be done in `SageMath`, from which they can be converted in a continuous format with rust and compiled easily after each computation. 
