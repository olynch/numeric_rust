use nalgebra::*;

pub mod adaptive_step;
pub mod constant_step;

pub use adaptive_step::*;
pub use constant_step::*;

pub trait AdaptiveStrategy<F, E> {
    fn init_dt(&self) -> F;

    /// Return Ok(new_dt) if the step succeeds, and Err(new_dt) if the step
    /// should be retried with a smaller dt
    fn try_accept(&self, cur_dt: F, error: E, y: DVectorView<F>) -> Result<F, F>;
}
