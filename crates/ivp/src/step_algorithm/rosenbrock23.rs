use crate::OdeSystem;

use super::StepAlgorithm;
use nalgebra::*;
use num_traits::Float;

pub struct Rosenbrock23;

pub struct Rosenbrock23Cache<F> {
    e32: F,
    d: F,
    mass_matrix: DMatrix<F>,
    f0: DVector<F>,
    f1: DVector<F>,
    f2: DVector<F>,
}

impl<F: Float + ComplexField<RealField = F>> StepAlgorithm<F> for Rosenbrock23 {
    type Cache = Rosenbrock23Cache<F>;
    type Interpolant = [DVector<F>; 3];
    type ErrorEstimate = F;

    fn init_cache<S: OdeSystem<F>>(&self, sys: &S) -> Self::Cache {
        let n = sys.dimension();
        Rosenbrock23Cache {
            e32: F::from(6. + f64::sqrt(2.0)).unwrap(),
            d: F::from(1. / (2. + f64::sqrt(2.0))).unwrap(),
            mass_matrix: sys.mass_matrix(),
            f0: DVector::zeros(n),
            f1: DVector::zeros(n),
            f2: DVector::zeros(n),
        }
    }

    fn interpolate(
        &self,
        y0: DVectorView<F>,
        y1: DVectorView<F>,
        _interpolant: &Self::Interpolant,
        _dt: F,
        s: F,
    ) -> DVector<F> {
        y0.scale(F::one() - s) + y1.scale(s)
    }

    fn step<S: OdeSystem<F>>(
        &self,
        cache: &mut Self::Cache,
        system: &S,
        mut y1: DVectorViewMut<F>,
        y0: DVectorView<F>,
        t: F,
        dt: F,
    ) -> (Self::Interpolant, Self::ErrorEstimate) {
        let n = system.dimension();
        // let d = cache.d;
        let dto2 = dt / F::from(2.).unwrap();
        let dtd = dt * cache.d;
        let neginvdtd = -Float::recip(dtd);
        let mut w = DMatrix::zeros(n, n);
        system.jacobian(w.as_view_mut(), y0, t);
        w += cache.mass_matrix.scale(neginvdtd);
        let wlu = w.lu();
        system.vfield(cache.f0.as_view_mut(), y0, t);
        // TODO: for non-autonomous system, we need f0 + dtd * dT
        let mut k1 = wlu.solve(&cache.f0).unwrap();
        k1.scale_mut(neginvdtd);
        y1.copy_from(&y0);
        y1.axpy(dto2, &k1, F::one());
        system.vfield(cache.f1.as_view_mut(), y1.as_view(), t + dto2);
        let mut k2 = &cache.f1 - &cache.mass_matrix * &k1;
        wlu.solve_mut(&mut k2);
        k2.axpy(F::one(), &k1, neginvdtd);
        y1.copy_from(&y0);
        y1.axpy(dt, &k2, F::one());
        system.vfield(cache.f2.as_view_mut(), y1.as_view(), t + dt);
        let mut k3 = DVector::zeros(n);
        for i in 0..n {
            k3[i] = cache.f2[i]
                - cache.e32 * (k2[i] - cache.f1[i])
                - F::from(2.).unwrap() * (k1[i] - cache.f0[i]);
        }
        wlu.solve_mut(&mut k3);
        let error = EuclideanNorm
            .norm(&(&k1 - k2.scale(F::from(2.).unwrap()) + &k3).scale(dt / F::from(6.).unwrap()));
        ([k1, k2, k3], error)
    }
}
