use super::*;
use num_traits::Float;

pub struct Euler;

impl<F: Float + Scalar + ComplexField<RealField = F>> StepAlgorithm<F> for Euler {
    type Cache = DVector<F>;
    type Interpolant = ();
    type ErrorEstimate = ();

    fn init_cache<S: OdeSystem<F>>(&self, sys: &S) -> Self::Cache {
        DVector::zeros(sys.dimension())
    }

    fn step<S: OdeSystem<F>>(
        &self,
        cache: &mut Self::Cache,
        system: &S,
        mut y1: DVectorViewMut<F>,
        y0: DVectorView<F>,
        t: F,
        dt: F,
    ) -> ((), ()) {
        system.vfield(cache.as_view_mut(), y0, t);
        y1.copy_from(&y0);
        y1.axpy(dt, &cache, F::one());
        ((), ())
    }

    fn interpolate(
        &self,
        y0: DVectorView<F>,
        y1: DVectorView<F>,
        _interpolant: &(),
        _dt: F,
        t: F,
    ) -> DVector<F> {
        y0.scale(F::from(1.).unwrap() - t) + y1.scale(t)
    }
}
