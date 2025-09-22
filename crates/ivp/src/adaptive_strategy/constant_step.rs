use num_traits::Float;

use super::*;

pub struct ConstantStep<F>(pub F);

impl<F: Float, E> AdaptiveStrategy<F, E> for ConstantStep<F> {
    fn init_dt(&self) -> F {
        self.0
    }

    fn try_accept(&self, _cur_dt: F, _error: E, _y: DVectorView<F>) -> Result<F, F> {
        Ok(self.0)
    }
}
