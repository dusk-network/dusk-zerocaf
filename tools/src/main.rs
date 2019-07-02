extern crate num_bigint;
extern crate num;

use num_bigint::BigUint;
use num::{Zero, One, Integer, Num};
use std::str::FromStr;
use num_traits::cast::ToPrimitive;
use num_traits::pow;



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
    

    let x = to_scalar_base_52("1809251394333065553493296640760748560207343510400633813116524750123642650623");
    println!("R: {:?}", x);
    let y = to_scalar_base_52("6145104759870991071742105800796537629880401874866217824609283457819451087098");
    println!("R: {:?}", y);
    let x_y_mont = to_scalar_base_52("532695731852302855127194384594826090862673352940952884676179882511233227548");
    println!("R: {:?}", x_y_mont);
    let x_y = to_scalar_base_52("890263784947025690345271110799906008759402458672628420828189878638015362081");
    println!("R: {:?}", x_y);
    let mont_x = to_scalar_base_52("658448296334113745583381664921721413881518248721417041768778176391714104386");
    println!("{:?}", mont_x);
    */
   
    //let _ = from_scalar_base_52(&[1682248870925813, 4078880264703668, 2289123149127681, 4169238435752846, 2104335664921]);
    //let mut x = from_scalar_base_52(&[0, 1, 0, 0, 0]);
    //let x = from_scalar_base_52(&[3986048181299919, 3661573232122943, 4503599627233715, 4503599627370495, 17592186044415]);
    //let x = phase1(&BigUint::from_str("182687704666362864775460604089535377456991567872").unwrap());
    let x = to_scalar_base_52("3289281484154527543902855200365966432305982921603192322330710365097124688423");
    println!("{:?}", x);

    //let res = phase1(&BigUint::from_str("2009874587549").unwrap());
    //println!("res: {} as FieldELEM: {:?}", res, to_scalar_base_52(&res.to_str_radix(10)));

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

/// Gets a number as a &str and encodes it like a FieldElement51 representation of zerocaf repo.
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

/// Gets a number as a &str and encodes it like a Scalar representation of zerocaf repo.
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

/// Gets a Scalar in base 52 and transforms it into a normal representation
pub fn from_scalar_base_52(limbs: &[u64; 5]) -> () {
    let mut res: BigUint = BigUint::zero();
    let two_pow_52 = BigUint::from_str("4503599627370496").unwrap();
    for i in 0..5 {
        res = res + (pow(two_pow_52.clone(), i) * limbs[i]);
    }
    println!("{}", res);
    ()
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

pub fn phase1(input: &BigUint) -> BigUint {
    let p = BigUint::from_str("7237005577332262213973186563042994240832148273380272487819684194273658967345").unwrap();
    let mut u = p.clone();
    let mut v = input.clone();
    let mut r = BigUint::zero();
    let mut s = BigUint::one();
    let two = &s + &s;
    let mut k: u64 = 0;
    while v > BigUint::zero() {
       
        match(u.is_even(), v.is_even(), u > v, v >= u) {
            // u is even
            (true, _, _, _) => {

                u = &u / &two;
                s = (&s * &two);
            },
            // u isn't even but v is even
            (false, true, _, _) => {

                v = &v / &two;
                r = (&r * &two);
            },
            // u and v aren't even and u > v
            (false, false, true, _) => {

                u = (&u - &v);
                u = &u / &two;
                r = (&r + &s);
                s = (&s * &two);
            },
            // u and v aren't even and v > u
            (false, false, false, true) => {

                v = (&v - &u);
                v = &v / &two;
                s = (&r + &s);
                r = (&r * &two);
            },
            (false, false, false, false) => panic!("InverseMod does not exist"),
        }
         k+=1;
        println!("Values on iteration: {}: \nr = {:?}\ns = {:?}\nv = {:?}\nu = {:?}", k, to_scalar_base_52(&r.to_str_radix(10)), to_scalar_base_52(&s.to_str_radix(10)),to_scalar_base_52(&v.to_str_radix(10)),to_scalar_base_52(&u.to_str_radix(10)));
    }
    if r >= p {
        println!("Inside if: {:?}", to_scalar_base_52(&(&r - &p).to_str_radix(10)));
        r = &r - &p;
    }
    println!("Outside if: {}", (&p - &r));
    &p - &r
 }