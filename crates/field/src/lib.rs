use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub trait Field:
    Copy
    + Eq
    + std::fmt::Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
{
    const ZERO: Self;
    const ONE: Self;
    const TWO_ADICITY: u32;

    fn from_u64(x: u64) -> Self;
    fn inverse(self) -> Option<Self>;
    fn root_of_unity(log_n: u32) -> Self;

    fn pow(self, mut e: u64) -> Self {
        let (mut base, mut acc) = (self, Self::ONE);
        while e > 0 {
            if e & 1 == 1 {
                acc *= base;
            }
            base *= base;
            e >>= 1;
        }
        acc
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Goldilocks(pub u64);

impl Goldilocks {
    pub const PRIME: u64 = 0xFFFFFFFF00000001; // 2^64 - 2^32 + 1
    const GENERATOR: u64 = 7;
}

impl Field for Goldilocks {
    const ZERO: Self = Self(0);
    const ONE: Self = Self(1);
    const TWO_ADICITY: u32 = 32;

    fn from_u64(x: u64) -> Self {
        Self(x % Self::PRIME)
    }

    fn inverse(self) -> Option<Self> {
        if self == Self::ZERO {
            None
        } else {
            // Fermat: a^(p-2) = a^(-1) mod p
            Some(self.pow(Self::PRIME - 2))
        }
    }

    fn root_of_unity(log_n: u32) -> Self {
        assert!(log_n <= Self::TWO_ADICITY, "log_n exceeds 2-adicity");
        let exponent = (Self::PRIME - 1) >> log_n;
        Self(Self::GENERATOR).pow(exponent)
    }
}

impl Add for Goldilocks {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let sum = (self.0 as u128 + rhs.0 as u128) % (Self::PRIME as u128);
        Self(sum as u64)
    }
}

impl Sub for Goldilocks {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let diff = ((self.0 as u128 + Self::PRIME as u128) - rhs.0 as u128) % (Self::PRIME as u128);
        Self(diff as u64)
    }
}

impl Mul for Goldilocks {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let prod = (self.0 as u128 * rhs.0 as u128) % (Self::PRIME as u128);
        Self(prod as u64)
    }
}

impl Neg for Goldilocks {
    type Output = Self;
    fn neg(self) -> Self::Output {
        if self.0 == 0 {
            Self(0)
        } else {
            Self(Self::PRIME - self.0)
        }
    }
}

impl AddAssign for Goldilocks {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Goldilocks {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for Goldilocks {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    prop_compose! {
        fn field_element()(val in 0..Goldilocks::PRIME) -> Goldilocks {
            Goldilocks(val)
        }
    }

    proptest! {
        #[test]
        fn test_associativity(a in field_element(), b in field_element(), c in field_element()) {
            assert_eq!((a+b) + c, a + (b+c));
            assert_eq!((a*b) * c, a * (b*c));
        }

        #[test]
        fn test_distributivity(a in field_element(), b in field_element(), c in field_element()) {
            assert_eq!(a * (b+c), (a*b) + (a*c));
        }

        #[test]
        fn test_inverse(a in field_element()) {
            if a != Goldilocks::ZERO {
                let inv = a.inverse().unwrap();
                assert_eq!(a * inv, Goldilocks::ONE);
            }
        }
    }

    #[test]
    fn test_generator() {
        let g = Goldilocks(7);
        let p_minus_1 = Goldilocks::PRIME - 1;
        // p-1 = 2^32 * (2^32-1) and 2^32 - 1 = 3 * 5 * 17 * 257 * 65537
        let factors = [2, 3, 5, 17, 257, 65537];

        for q in factors.iter() {
            let exponent = p_minus_1 / q;
            assert_ne!(
                g.pow(exponent),
                Goldilocks::ONE,
                "7 is not a generator, factor {}",
                q
            );
        }
    }

    #[test]
    fn test_roots_of_unity() {
        for k in 1..=Goldilocks::TWO_ADICITY {
            let root = Goldilocks::root_of_unity(k);
            assert_eq!(root.pow(1 << k), Goldilocks::ONE); // order divides 2^k
            assert_ne!(root.pow(1 << (k - 1)), Goldilocks::ONE); // order is exactly 2^k
        }
    }
}
