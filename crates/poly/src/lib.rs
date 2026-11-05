use field::{Field, Goldilocks};
use std::ops::{Mul, Sub};

#[derive(Debug, PartialEq)]
struct DivResult<F: Field> {
    quotient: Polynomial<F>,
    remainder: Polynomial<F>,
}

#[derive(Clone, Debug, PartialEq)]
struct Polynomial<F: Field> {
    coefficients: Vec<F>,
}

impl<F: Field> Polynomial<F> {
    fn new(mut coefficients: Vec<F>) -> Self {
        // Ignore trailing zeros
        let end = coefficients
            .iter()
            .rposition(|&c| c != F::ZERO)
            .map_or(0, |i| i + 1);
        coefficients.drain(end..);
        Self { coefficients }
    }
    fn degree(&self) -> Option<usize> {
        self.coefficients.len().checked_sub(1)
    }
    fn coefficient(&self, i: usize) -> Option<&F> {
        self.coefficients.get(i)
    }

    // Another impl is:
    // 1. store the last element (highest degree) in an accumulator
    // 2. highest degree coefficient gets multiplied by x
    // 3. next coefficient gets added to the value in the accumulator
    // 4. repeat until list of coefficients is complete
    // eg:
    //      x=5
    //      5x^3 + 2x^2 + 3x + 7 = 5*125 + 2*25 + 22 = 625 + 50 + 22 = 697
    //      [7,3,2,5]
    //      1. acc is 5 - we enter the loop
    //      2. i=2, acc is 5 * x (25) + 2 (from 2*x^2) = 27
    //      3. i=1, acc is 27 * x (135) + 3 (from 3*x) = 138
    //      4. i=0, acc is 138 * x (690) + 7 (constant) = 697
    fn eval(&self, x: F) -> F {
        let mut acc = F::ZERO;
        for i in (0..self.coefficients.len()).rev() {
            acc = acc.mul(x).add(*self.coefficients.get(i).unwrap());
        }
        acc
    }

    fn div(&self, g: Polynomial<F>) -> DivResult<F> {
        if self.coefficients.len() == 0 || g.degree() > self.degree() {
            DivResult {
                quotient: Polynomial::new(vec![]),
                remainder: self.clone(),
            }
        } else {
            let mut quotient = Polynomial {
                coefficients: vec![Field::ZERO; self.coefficients.len()],
            };
            let mut remainder_degree = self.degree().unwrap();
            let mut remainder = self.clone();
            let g_degree = g.degree().unwrap();

            while remainder_degree >= g_degree {
                let lead_f = remainder.coefficients.get(remainder_degree).unwrap();
                let lead_g = g.coefficients.get(g_degree);
                let div_result = lead_f.mul(lead_g.unwrap().inverse().unwrap());
                let deg_diff = remainder_degree.sub(g_degree);
                // div_result * X^deg_diff
                quotient.coefficients[deg_diff] = div_result;
                let mut temp_quotient = Polynomial {
                    coefficients: vec![Field::ZERO; quotient.coefficients.len()],
                };
                for i in 0..g.coefficients.len() {
                    temp_quotient.coefficients[i + deg_diff] = div_result.mul(g.coefficients[i]);
                }
                for i in 0..remainder.coefficients.len() {
                    remainder.coefficients[i] =
                        remainder.coefficients[i].sub(temp_quotient.coefficients[i]);
                }
                remainder = Polynomial::new(remainder.coefficients);
                remainder_degree = remainder.degree().unwrap();
            }
            DivResult {
                quotient: Polynomial::new(quotient.coefficients),
                remainder,
            }
        }
    }

    fn add(self, g: Polynomial<F>) -> Polynomial<F> {
        let lhs_degree = match self.degree() {
            Some(v) => v,
            None => return g,
        };

        let rhs_degree = match g.degree() {
            Some(v) => v,
            None => return self,
        };

        let (shorter, longer) = if rhs_degree < lhs_degree {
            (&g, &self)
        } else {
            (&self, &g)
        };

        let mut result: Polynomial<F> = Polynomial {
            coefficients: vec![F::ZERO; longer.coefficients.len()],
        };

        for i in 0..longer.coefficients.len() {
            if i < shorter.coefficients.len() {
                result.coefficients[i] = shorter.coefficients[i] + longer.coefficients[i];
            } else {
                result.coefficients[i] = longer.coefficients[i];
            }
        }

        Polynomial::new(result.coefficients)
    }

    fn neg(mut self) -> Self {
        for c in self.coefficients.iter_mut() {
            *c = c.neg()
        }

        self
    }

    fn sub(self, g: Polynomial<F>) -> Polynomial<F> {
        self.add(g.neg())
    }

    fn scalar_mul(mut self, s: F) -> Self {
        for c in self.coefficients.iter_mut() {
            *c = c.mul(s);
        }

        Polynomial::new(self.coefficients)
    }

    fn mul(self, g: Polynomial<F>) -> Polynomial<F> {
        let f_len = self.coefficients.len();
        let g_len = g.coefficients.len();

        let mut result = Polynomial {
            coefficients: vec![F::ZERO; f_len + g_len]
        };

        for i in 0..f_len {
            for j in 0..g_len {
                result.coefficients[i+j] += self.coefficients[i].mul(g.coefficients[j]);
            }
        }
        Polynomial::new(result.coefficients)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use field::Goldilocks;
    use std::ops::Neg;
    #[test]
    fn test_degree() {
        // 5x^3 + 2x^2 + 3x
        let degree: usize = 3;
        let coefficients: Vec<Goldilocks> = vec![
            Goldilocks::ZERO,
            Goldilocks(3),
            Goldilocks(2),
            Goldilocks(5),
            Goldilocks(0),
            Goldilocks(0),
        ];
        let poly = Polynomial::new(coefficients);
        assert_eq!(poly.degree().unwrap(), degree);
    }
    #[test]
    fn test_empty() {
        let empty_coefficients: Vec<Goldilocks> = vec![];
        let poly = Polynomial::new(empty_coefficients);
        assert_eq!(poly.degree(), None)
    }
    #[test]
    fn test_zeros_only() {
        let zero_coefficients: Vec<Goldilocks> = vec![Goldilocks::ZERO, Goldilocks::ZERO];
        let poly = Polynomial::new(zero_coefficients);
        assert_eq!(poly.degree(), None)
    }

    #[test]
    fn test_coefficient() {
        // 5x^3 + 2x^2 + 3x
        let coefficients: Vec<Goldilocks> = vec![
            Goldilocks::ZERO,
            Goldilocks(3),
            Goldilocks(2),
            Goldilocks(5),
            Goldilocks(0),
            Goldilocks(0),
        ];
        let poly = Polynomial::new(coefficients);
        assert_eq!(poly.coefficient(3).unwrap(), &Goldilocks(5));
    }
    #[test]
    fn test_out_of_range_coefficient() {
        // 5x^3 + 2x^2 + 3x
        let coefficients: Vec<Goldilocks> = vec![
            Goldilocks::ZERO,
            Goldilocks(3),
            Goldilocks(2),
            Goldilocks(5),
            Goldilocks(0),
            Goldilocks(0),
        ];
        let poly = Polynomial::new(coefficients);
        assert_eq!(poly.coefficient(5), None);
    }

    #[test]
    fn test_eval() {
        // x=5
        // 5x^3 + 2x^2 + 3x + 7 = 5*125 + 2*25 + 22 = 625 + 50 + 22 = 697
        let coefficients: Vec<Goldilocks> = vec![
            Goldilocks(7),
            Goldilocks(3),
            Goldilocks(2),
            Goldilocks(5),
            Goldilocks(0),
            Goldilocks(0),
        ];
        let poly = Polynomial::new(coefficients);
        assert_eq!(poly.eval(Goldilocks(5)), Goldilocks(697));
    }
    #[test]
    fn test_eval_at_zero() {
        let coefficients: Vec<Goldilocks> = vec![
            Goldilocks(7),
            Goldilocks(3),
            Goldilocks(2),
            Goldilocks(5),
            Goldilocks(0),
            Goldilocks(0),
        ];
        let poly = Polynomial::new(coefficients);
        assert_eq!(poly.eval(Goldilocks(0)), Goldilocks(7));
    }

    #[test]
    fn test_eval_zero_poly() {
        let coefficients: Vec<Goldilocks> = vec![Goldilocks(0), Goldilocks(0), Goldilocks(0)];
        let poly = Polynomial::new(coefficients);
        assert_eq!(poly.eval(Goldilocks(5)), Goldilocks::ZERO);
    }

    #[test]
    fn test_monic_divisor_nonzero_remainder() {
        let f_coeffs: Vec<Goldilocks> =
            vec![Goldilocks(4), Goldilocks(3), Goldilocks(2), Goldilocks(1)];
        let g_coeffs: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(1)];
        let quotient_coeffs: Vec<Goldilocks> = vec![Goldilocks(2), Goldilocks(1), Goldilocks(1)];
        let remainder_coeffs: Vec<Goldilocks> = vec![Goldilocks(2)];

        let f = Polynomial::new(f_coeffs);
        let g = Polynomial::new(g_coeffs);
        let q = Polynomial::new(quotient_coeffs);
        let r = Polynomial::new(remainder_coeffs);

        assert_eq!(
            f.div(g),
            DivResult {
                quotient: q,
                remainder: r
            }
        );
    }

    #[test]
    fn test_add() {
        let f_coeffs: Vec<Goldilocks> =
            vec![Goldilocks(4), Goldilocks(3), Goldilocks(2), Goldilocks(1)];
        let g_coeffs: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(1)];

        let f = Polynomial::new(f_coeffs);
        let g = Polynomial::new(g_coeffs);

        let result: Vec<Goldilocks> =
            vec![Goldilocks(5), Goldilocks(4), Goldilocks(2), Goldilocks(1)];
        let r = Polynomial::new(result);

        assert_eq!(f.add(g), r);
    }

    #[test]
    fn test_add_zero_poly() {
        let f_coeffs: Vec<Goldilocks> =
            vec![Goldilocks(4), Goldilocks(3), Goldilocks(2), Goldilocks(1)];
        let f = Polynomial::new(f_coeffs);
        let result = f.clone();
        let zero_poly: Polynomial<Goldilocks> = Polynomial::new(vec![Field::ZERO]);

        assert_eq!(f.add(zero_poly), result);
    }

    #[test]
    fn test_add_to_zero_poly() {
        let f_coeffs: Vec<Goldilocks> =
            vec![Goldilocks(4), Goldilocks(3), Goldilocks(2), Goldilocks(1)];
        let f = Polynomial::new(f_coeffs);
        let result = f.clone();
        let zero_poly: Polynomial<Goldilocks> = Polynomial::new(vec![Field::ZERO]);

        assert_eq!(zero_poly.add(f), result);
    }

    #[test]
    fn test_cancellation() {
        let f_coeffs: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(0), Goldilocks(1)];
        let g_coeffs: Vec<Goldilocks> = vec![Goldilocks(0), Goldilocks(0), Goldilocks::ONE.neg()];

        let one_poly = Polynomial::new(vec![Goldilocks::ONE]);

        let f = Polynomial::new(f_coeffs);
        let g = Polynomial::new(g_coeffs);

        assert_eq!(f.add(g), one_poly);
    }

    #[test]
    fn test_complete_cancellation() {
        let f_coeffs: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(0), Goldilocks(1)];
        let empty_poly = Polynomial::new(vec![]);

        let f = Polynomial::new(f_coeffs);
        let f_clone = f.clone();

        assert_eq!(f.add(f_clone.neg()), empty_poly);
    }

    #[test]
    fn test_subtraction() {
        let f_coeffs: Vec<Goldilocks> =
            vec![Goldilocks(4), Goldilocks(3), Goldilocks(2), Goldilocks(1)];
        let g_coeffs: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(1)];

        let f = Polynomial::new(f_coeffs);
        let g = Polynomial::new(g_coeffs);

        let result: Vec<Goldilocks> =
            vec![Goldilocks(3), Goldilocks(2), Goldilocks(2), Goldilocks(1)];
        let r = Polynomial::new(result);

        assert_eq!(f.sub(g), r);
    }

    #[test]
    fn test_poly_minus_itself() {
        let f_coeffs: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(0), Goldilocks(1)];
        let empty_poly = Polynomial::new(vec![]);
        let f = Polynomial::new(f_coeffs);
        let f_clone = f.clone();
        assert_eq!(f.sub(f_clone), empty_poly);
    }

    #[test]
    fn test_poly_minus_zero() {
        let f_coeffs: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(0), Goldilocks(1)];
        let f = Polynomial::new(f_coeffs);
        let empty_poly = Polynomial::new(vec![]);
        assert_eq!(f.clone().sub(empty_poly), f);
    }

    #[test]
    fn test_zero_minus_poly() {
        let f_coeffs: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(0), Goldilocks(1)];
        let f = Polynomial::new(f_coeffs);
        let empty_poly = Polynomial::new(vec![]);
        assert_eq!(empty_poly.sub(f.clone()), f.neg());
    }

    #[test]
    fn test_minus_then_add() {
        let f_coeffs: Vec<Goldilocks> = vec![Goldilocks(6), Goldilocks(1), Goldilocks(1)];
        let f = Polynomial::new(f_coeffs);
        let g = Polynomial::new(vec![Goldilocks(5)]);
        let result: Vec<Goldilocks> = vec![Goldilocks(1), Goldilocks(1), Goldilocks(1)];
        let r = Polynomial::new(result);
        assert_eq!(f.clone().sub(g.clone()), r);
        assert_eq!(r.add(g), f);
    }
}
