use nalgebra::*;
use num_traits::Float;

use crate::adaptive_strategy::AdaptiveStrategy;

pub struct IntegralController<F> {
    init: F,
    atol: F,
    rtol: F,
    order: u32,
}

impl<F> IntegralController<F> {
    pub fn new(init: F, atol: F, rtol: F, order: u32) -> Self {
        Self {
            init,
            atol,
            rtol,
            order,
        }
    }
}

impl<F: Float + PartialOrd> AdaptiveStrategy<F, F> for IntegralController<F> {
    fn init_dt(&self) -> F {
        self.init
    }

    fn try_accept(&self, cur_dt: F, error: F, y: DVectorView<F>) -> Result<F, F> {
        let tol =
            self.atol + self.rtol * y.iter().map(|x| x.abs()).reduce(|x, y| x.max(y)).unwrap();
        let new_dt = cur_dt * (tol / error).powf(F::from(self.order + 1).unwrap().recip());
        Ok(new_dt)
    }
}
