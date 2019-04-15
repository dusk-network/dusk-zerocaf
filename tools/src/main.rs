extern crate num_bigint;
extern crate num;

use num_bigint::BigUint;
use num::{Zero, One, Integer, Num};
use std::str::FromStr;
use num_traits::cast::ToPrimitive;


fn main() {
    /*Examples of function calls
    num_from_bytes_le(&[76, 250, 187, 243, 105, 92, 117, 70, 234, 124, 126, 180, 87, 149, 62, 249, 16, 149, 138, 56, 26, 87, 14, 76, 251, 39, 168, 74, 176, 202, 26, 84]);
    
    let res = to_field_elem_51(&"1201935917638644956968126114584555454358623906841733991436515590915937358637");
    println!("{:?}", res);

    let scalar_res = to_scalar_base_52(&"1201935917638644956968126114584555454358623906841733991436515590915937358637");
    println!("{:?}", res);

    hex_bytes_le("120193591763864495696812611458455545435862390684173399143651559091593735863735685683568356835683");
    from_radix_to_radix_10("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab", 16u32);

    let m = BigUint::from_str("750791094644726559640638407699").unwrap();
    let x1 = BigUint::from_str("540019781128412936473322405310").unwrap();
    let x2 = BigUint::from_str("515692107665463680305819378593").unwrap();

    reduction(m, x1, x2).unwrap();
    */
    let a = to_scalar_base_52("182687704666362864775460604089535377456991567872");
    let b = to_scalar_base_52("904625697166532776746648320197686575422163851717637391703244652875051672039");
    let b_minus_a = to_scalar_base_52("904625697166532776746648320014998870755800986942176787613709275418060104167");
    let a_minus_b = to_scalar_base_52("365375409332725729550921208179070754913983135744");

    println!("A: {:?}\nB: {:?}\nA - B: {:?}\nB - A: {:?}\n",a, b, a_minus_b, b_minus_a);

    let r = to_scalar_base_52("904625697166532776746648320380374280088526716493097995792780030332043239911");
    println!("R: {:?}", r);

}

/// The num has to be positive! Otherways it will fail
pub fn hex_bytes_le(num: &str) {
    let num: BigUint = BigUint::from_str(num).unwrap();
    println!("Encoding result -> {:x?}", num.to_bytes_le());
}

/// Prints a number in radix 10 given a reference of a slice of Little Endian bytes of it.
pub fn num_from_bytes_le(num: &[u8]){
    println!("{}", BigUint::from_bytes_le(num));
}

/// Changes a number from radix X to radix 10
pub fn from_radix_to_radix_10(num: &str, radix: u32) {
    println!("{}",BigUint::from_str_radix(num, radix).unwrap());
}

/// Gets a number as a &str and encodes it like a FieldElement51 representation of corretto repo.
pub fn to_field_elem_51(num: &str) -> [u64; 5] {
    let mut resp_as_array = [0u64;5];
    let mut response = vec!();
    let two_pow_51 = BigUint::from_str("2251799813685248").unwrap();
    let mut num = BigUint::from_str(num).unwrap();

    for _i in 0..5 {
        response.push(&num % &two_pow_51);
        num = &num / &two_pow_51;
    }
    
    let mut res2: Vec<u64> = response.iter().map(|x| u64::from_str(&x.to_str_radix(10u32)).unwrap()).collect();
    resp_as_array.swap_with_slice(&mut res2);
    resp_as_array
}

/// Gets a number as a &str and encodes it like a Scalar representation of corretto repo.
pub fn to_scalar_base_52(num: &str) -> [u64; 5] {
    let mut resp_as_array = [0u64;5];
    let mut response = vec!();
    let two_pow_52 = BigUint::from_str("4503599627370496").unwrap();
    let mut num = BigUint::from_str(num).unwrap();

    for _i in 0..5 {
        response.push(&num % &two_pow_52);
        num = &num / &two_pow_52;
    }
    
    let mut res2: Vec<u64> = response.iter().map(|x| u64::from_str(&x.to_str_radix(10u32)).unwrap()).collect();
    resp_as_array.swap_with_slice(&mut res2);
    resp_as_array
}

/// Montgomery struct
#[derive(Debug)]
pub struct Montgomery {
    pub BASE: BigUint,
    pub m: BigUint,
    pub rrm: BigUint,
    pub n: i32
}

impl Montgomery {
    pub fn new(m: BigUint) -> Result<Self, &'static str> {
        if m < BigUint::zero() || m.is_even() {
            Err("Invalid input!")
        } else {
            let mont = Montgomery {
                BASE: BigUint::from(2u32),
                m: m.clone(),
                n: m.bits() as i32,
                rrm: (BigUint::one() << (m.bits() * 2)) % m
            }; 
            return Ok(mont)
        }
    }

    pub fn reduce(&self, prop: BigUint) -> BigUint {
        let mut a: BigUint = prop.clone();
        for i in 0..self.n {
            if !a.is_even() {a+=self.m.clone()}
            a = a.clone() >> 1;
        }
        if a >= self.m {a -= self.m.clone()}
        a
    }
}

pub fn reduction(m: BigUint, x1: BigUint, x2: BigUint) -> Result<(), &'static str> {
    let mont = Montgomery::new(m.clone()).unwrap();
    let t1 = x1.clone() * mont.rrm.clone();
    let t2 = x2.clone() * mont.rrm.clone();

    let r1 = mont.reduce(t1.clone());
    let r2 = mont.reduce(t2.clone());

    let r = BigUint::one() << mont.n.clone() as usize;
    println!("b = {}\n", mont.BASE);
    println!("n = {}\n", mont.n);
    println!("r = {}\n", r);
    println!("m = {}\n", m);
    println!("t1 = {}\n", t1);
    println!("t2 = {}\n", t2);
    println!("r1 = {}\n", r1);
    println!("r2 = {}\n", r2);
    println!("Original x1: {}\n", x1.clone());
    println!("Recovered from r1: {}\n", mont.reduce(r1.clone()));
    println!("Original x2: {}\n", x2.clone());
    println!("Recovered from r2: {}\n", mont.reduce(r2.clone()));
    println!("Montgomery computation of x1 ^ x2 mod m :");
    let prod = x1.modpow(&x2, &mont.m);

    println!("{}\n", prod);

    Ok(())
}