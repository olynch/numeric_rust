use nalgebra::*;
use num_traits::Float;

use crate::adaptive_strategy::AdaptiveStrategy;

pub struct ProportionalIntegralController<F> {
    init: F,
    atol: F,
    rtol: F,
    order: u32,
    beta1: F,
    beta2: F,
}

impl<F> ProportionalIntegralController<F> {
    pub fn new(init: F, atol: F, rtol: F, order: u32) -> Self {
        todo!()
        // Self {
        //     init,
        //     atol,
        //     rtol,
        //     order,
        // }
    }
}

impl<F: Float + PartialOrd> AdaptiveStrategy<F, F> for ProportionalIntegralController<F> {
    fn init_dt(&self) -> F {
        self.init
    }

    // TODO: implement PI controller, taking code from: https://github.com/SciML/OrdinaryDiffEq.jl/blob/8d830c9cc6326904d54896971194010c421ee890/lib/OrdinaryDiffEqCore/src/integrators/controllers.jl#L73
    // It seems like it might make sense to factor out init, atol, rtol, etc. to be
    // common to PI controller and P controller
    fn try_accept(&self, cur_dt: F, error: F, y: DVectorView<F>) -> Result<F, F> {
        let tol =
            self.atol + self.rtol * y.iter().map(|x| x.abs()).reduce(|x, y| x.max(y)).unwrap();
        let new_dt = cur_dt * (tol / error).powf(F::from(self.order + 1).unwrap().recip());
        Ok(new_dt)
    }
}
