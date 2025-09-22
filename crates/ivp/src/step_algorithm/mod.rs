use crate::system::*;
use nalgebra::*;
use num_traits::Float;

pub mod euler;
pub mod rosenbrock23;
pub mod tsit5;

pub use euler::*;
pub use rosenbrock23::*;
pub use tsit5::*;

pub trait StepAlgorithm<F: Scalar + Float> {
    type Cache;
    type Interpolant;
    type ErrorEstimate;

    fn init_cache<S: OdeSystem<F>>(&self, sys: &S) -> Self::Cache;

    fn step<S: OdeSystem<F>>(
        &self,
        cache: &mut Self::Cache,
        system: &S,
        y1: DVectorViewMut<F>,
        y0: DVectorView<F>,
        t: F,
        dt: F,
    ) -> (Self::Interpolant, Self::ErrorEstimate);

    fn interpolate(
        &self,
        y0: DVectorView<F>,
        y1: DVectorView<F>,
        interpolant: &Self::Interpolant,
        dt: F,
        s: F,
    ) -> DVector<F>;
}
