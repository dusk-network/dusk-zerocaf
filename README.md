# Coretto: our cryptographic protocol for set inclusion using elliptic curve operations

**This repository contains the implementation of the `Doppio Curve` over the `Ristretto Scalar field`. This is a pure Rust implementation designed by the Dusk-Team.**

**WIP**


# Curve parameters:

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

# TODO:

### The refactoring relations are expressed as indentations.
- [ ] Create FieldElement Struct and implement the basic operations we need on a u64 backend.
  - [x] Find the proper radix value for FieldElement.
  - [ ] Add basic and needed constants.
  - [ ] Implement Reduce function to make the FieldElements fit on a 5 u64-bit limbs.
    - [x] Implement Addition. (Testing needed)
    - [ ] Implement Subtraction.
    - [ ] Implement Byte-encoding/decoding. (Encoding done)
    - [ ] Implement Multiplication on u64-backend with u128 usage.
  - [ ] Add proper tests for every function.
- [ ] Build Scalar Arithmetics and Scalar Struct definition.
    - [x] Find the proper radix value for FieldElement.
    - [ ] Add the required constants for computation.
      - [x] Implement Addition.
      - [x] implement Subtraction.
      - [x] Implement Inner_Multiplication.
      - [x] Implement Inner_Squaring.
      - [ ] Implement Montgomery_reduction. (BigUint implementation done on tools/src/main.rs to test it and get an idea.)
        - [ ] Implement Montgomery_Muliplication.
        - [ ] Implement Montgomery_Squaring.
- [ ] Create Conversions from Montgomery points to Weierstrass ones. (Not clear if necessary yet.)

# Encoding / Decoding tools and examples.
In order to work with our points along the curve, or any non trivial computuaions, for example those with tough notations - there has been a set of tools and examples has been created to make facilitiate the Encoding/Decoding processes. Thye can be found at: `tools/src/main.rs` 

**Examples**
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
> When performing operations with large values, such as: `2²⁵² - 121160309657751286123858757838224683208` it is recomended to compute through `SageMath`, as the user interface is adheres to these types of fucntions.

<br/>

> Operations with large numbers are recommended to be done in `SageMath`, where they can be converted in a continuouse format into rust and easily compiled each time. 
