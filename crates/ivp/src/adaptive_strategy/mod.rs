use nalgebra::*;

pub mod constant_step;
pub mod integral_controller;
pub mod proportional_integral_controller;

pub use constant_step::*;
pub use integral_controller::*;
pub use proportional_integral_controller::*;

pub trait AdaptiveStrategy<F, E> {
    fn init_dt(&self) -> F;

    /// Return Ok(new_dt) if the step succeeds, and Err(new_dt) if the step
    /// should be retried with a smaller dt
    fn try_accept(&self, cur_dt: F, error: E, y: DVectorView<F>) -> Result<F, F>;
}
